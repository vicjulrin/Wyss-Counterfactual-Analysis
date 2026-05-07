# Functions used by Reproducible_MC2.
# This file keeps only functions required by script.R or by another function in this file.
# Numbered development versions were renamed to final names.

# data_driven_qres ####
# What it does: Computes randomized quantile residuals for binomial or Poisson GLM objects.
# Inputs:
# - model: fitted glm object with binomial or Poisson family.
# Output:
# - Numeric vector of standardized quantile residuals, one value per model observation.
data_driven_qres <- function(model) {
  y <- model$y
  mu <- model$fitted.values
  eta <- model$linear.predictors
  fam <- family(model)$family
  n <- length(y)
  
  if (fam == "poisson") {
    p_lower <- ppois(y - 1, lambda = mu)
    p_upper <- ppois(y, lambda = mu)
  } else if (fam == "binomial") {
    size <- if (is.null(model$prior.weights)) rep(1, n) else model$prior.weights
    p_lower <- pbinom(y - 1, size = size, prob = mu)
    p_upper <- pbinom(y, size = size, prob = mu)
  } else {
    stop("data_driven_qres() only supports Poisson or Binomial GLMs.")
  }
  
  u_intra_data <- abs(eta) %% 1
  u_final <- p_lower + u_intra_data * (p_upper - p_lower)
  z <- qnorm(pmin(pmax(u_final, 1e-16), 1 - 1e-16))
  
  as.numeric(scale(z, center = TRUE, scale = TRUE))
}


# multicol_analysis ####
# What it does: Removes aliased variables, clusters correlated covariates, and selects one covariate per correlation group.
# Inputs:
# - data_x: analysis data.frame or data.table.
# - treatment_col_x: treatment column name.
# - covars_x: candidate covariate names.
# - covars_x_key: covariates protected or prioritized during selection.
# - cor_threshold_x: maximum correlation threshold used to define groups.
# - family_distribution: GLM family used for the preliminary model; default is binomial.
# Output:
# - List with correlation matrix, grouped VIF table, dendrogram plot, and selected covariates.
multicol_analysis <- function (data_x, treatment_col_x, covars_x, covars_x_key, cor_threshold_x, family_distribution = "binomial") 
{
    covars_x2 <- covars_x
    covars_x_key2 <- if (is.null(covars_x_key)) {
        NULL
    }
    else {
        intersect(covars_x_key, covars_x)
    }
    repeat {
        formula_test_multicor <- as.formula(paste0(treatment_col_x, "~", paste0(covars_x2, collapse = "+")))
        test_multicor <- glm(formula_test_multicor, data = data_x, family = family_distribution)
        aliased <- alias(test_multicor)$Complete
        vars_to_remove <- {
            if (!is.null(aliased)) {
                names(which(apply(aliased, 1, any)))
            }
            else {
                NULL
            }
        }
        if (length(vars_to_remove) == 0) {
            break
        }
        vars_to_remove <- c(intersect(covars_x_key2, vars_to_remove), setdiff(vars_to_remove, covars_x_key2))
        covars_x2 <- setdiff(covars_x2, tail(vars_to_remove, 1))
    }
    vif_data <- car::vif(test_multicor) %>% as.data.frame() %>% {
        data.frame(Var = rownames(.), VIF = .[, 1])
    } %>% arrange(VIF)
    cordataR <- summary(test_multicor, correlation = T)[["correlation"]] %>% as.data.frame.matrix()
    cordataR[, "(Intercept)"] <- NULL
    cordataR <- cordataR[2:nrow(cordataR), ]
    NACol <- names(which(rowSums(is.na(cordataR)) > (ncol(cordataR)/2)))
    corMatrix <- cordataR %>% {
        .[!names(.) %in% NACol, ]
    } %>% {
        .[, !colnames(.) %in% NACol]
    }
    corMatrix[is.na(corMatrix)] <- 0
    corhclust <- hclust(as.dist(1 - abs(corMatrix)))
    cordend <- as.dendrogram(corhclust)
    cordend_data <- dendro_data(cordend)
    group_VIF_covars <- cutree(corhclust, h = 1 - cor_threshold_x) %>% as.data.frame %>% rownames_to_column("Var") %>% dplyr::rename(group = ".") %>% dplyr::filter(!Var %in% "(Intercept)") %>% list(vif_data) %>% join_all() %>% arrange(group, !Var %in% covars_x_key2, match(Var, covars_x_key2, nomatch = Inf), VIF)
    covars_no_multicol <- dplyr::filter(group_VIF_covars, !duplicated(group))$Var
    cordend_data$labels$label <- ifelse(cordend_data$labels$label %in% covars_no_multicol, paste0("* ", cordend_data$labels$label), cordend_data$labels$label)
    var_table <- with(cordend_data$labels, data.frame(y_center = x, y_min = x - 0.5, y_max = x + 0.5, Variable = as.character(label), height = 1))
    is.odd <- function(x) {
        x%%2 == 0
    }
    var_table$col <- rep_len(c("#EBEBEB", "white"), length.out = length(var_table$Variable)) %>% {
        if (is.odd(length(.))) {
            rev(.)
        }
        else {
            .
        }
    }
    segment_data <- with(segment(cordend_data), data.frame(x = y, y = x, xend = yend, yend = xend, cor = 1 - yend))
    multicol_dendro_plot <- ggplot() + annotate("rect", xmin = -0.05, xmax = 1.04, fill = var_table$col, ymin = var_table$y_min, ymax = var_table$y_max, alpha = 0.75) + geom_segment(data = segment_data, aes(x = 1 - x, y = y, xend = 1 - xend, yend = yend, label = cor), size = 0.3) + scale_y_continuous(breaks = cordend_data$labels$x, labels = cordend_data$labels$label) + coord_cartesian(expand = F) + labs(x = "Correlation", y = "Variables") + geom_vline(xintercept = cor_threshold_x, linetype = "dashed", 
        col = "red") + theme(panel.grid.major = element_line(color = "gray"))
    return(list(corMatrix = corMatrix, group_VIF_covars = group_VIF_covars, multicol_dendro_plot = multicol_dendro_plot, covars_no_multicol = covars_no_multicol))
}


# bestmodels_test ####
# What it does: Searches subsets of selected covariates and ranks candidate GLMs by AIC/BIC-style criteria.
# Inputs:
# - data_x: analysis data.frame or data.table.
# - treatment_col_x: treatment column name.
# - covars_x: candidate covariate names.
# - covars_x_key: covariates forced into the subset search when present.
# - cor_threshold_x: retained for compatibility; not used inside this function.
# - criteria_x: ranking criterion column, usually AIC or BIC.
# - family_distribution: GLM family used to score candidate models; default is binomial.
# Output:
# - List with model-variable table, ranked models, split model list, and model-selection plot.
bestmodels_test <- function (data_x, treatment_col_x, covars_x, covars_x_key, cor_threshold_x, criteria_x, family_distribution = "binomial") 
{
    covars_x_key2 <- if (is.null(covars_x_key)) {
        NULL
    }
    else {
        intersect(covars_x_key, covars_x)
    }
    covars_x2 <- c(intersect(covars_x_key2, covars_x), setdiff(covars_x, covars_x_key2))
    pre_formula_glm <- as.formula(paste0(treatment_col_x, " ~ ", paste0(covars_x2, collapse = "+")))
    model <- regsubsets(pre_formula_glm, data = data_x, nvmax = length(covars_x2), method = "seqrep", force.in = seq_along(covars_x_key2))
    summ_model <- summary(model)[c("rsq", "rss", "adjr2", "cp", "bic")] %>% as.data.frame() %>% dplyr::mutate(model = seq(nrow(.)))
    AIC_models <- lapply(seq(nrow(summ_model)), function(m) {
        coefs <- coef(model, id = m)
        vars <- names(coefs)[-1]
        form_test <- as.formula(paste0(treatment_col_x, "~", paste0(vars, collapse = "+")))
        glm_test <- glm(form_test, data = data_x, family = family_distribution)
        data_AIC <- data.frame(model = m, AIC = extractAIC(glm_test)[2])
        data_vars <- data.frame(model = m, vars = vars)
        data_form <- data.frame(model = m, formula = gsub("\\s+", " ", paste(deparse(form_test), collapse = "")))
        list(data_AIC = data_AIC, data_vars = data_vars, data_form = data_form)
    })
    forms_models <- rbind.fill(purrr::map(AIC_models, "data_form"))
    rank_models <- rbind.fill(purrr::map(AIC_models, "data_AIC")) %>% list(summ_model) %>% join_all() %>% dplyr::arrange(bic) %>% dplyr::mutate(rank_BIC = seq(nrow(.))) %>% dplyr::arrange(AIC) %>% dplyr::mutate(rank_AIC = seq(nrow(.)))
    data_vars <- rbind.fill(purrr::map(AIC_models, "data_vars")) %>% list(rank_models) %>% join_all() %>% dplyr::group_by(vars) %>% dplyr::mutate(freq_var = n()) %>% dplyr::arrange(freq_var) %>% dplyr::mutate(vars = factor(vars, levels = unique(.$vars)))
    vars_models <- data_vars %>% dplyr::arrange(!!sym(criteria_x)) %>% dplyr::mutate(model = factor(model, levels = unique(.$model))) %>% as.data.frame()
    list_models <- vars_models %>% split(.$model)
    better_models_plot <- ggplot() + geom_tile(data = vars_models, aes(x = model, y = vars, fill = !!sym(criteria_x)), color = "black", alpha = 0.5, size = 0.2) + scale_fill_gradientn(criteria_x, colors = brewer.pal(11, "Spectral"))
    return(list(vars_models = vars_models, rank_models = rank_models, list_models = list_models, better_models_plot = better_models_plot))
}


# moran_by_distance ####
# What it does: Estimates global Moran's I for treatment-model residuals at one distance band.
# Inputs:
# - data_x: data with PUID, x, y, treatment column, and covariates.
# - distance_x: neighborhood radius in map units.
# - treatment_col_x: treatment column name.
# - covars_x: covariates used in the treatment GLM.
# - resamp: number of permutations used for the p-value confidence interval.
# Output:
# - Data.frame with Moran.I.statistic, permutation expectation, variance, p-value, and distance.
moran_by_distance <- function(data_x, distance_x, treatment_col_x, covars_x, resamp= 999){
  
  #parameters
  ## organize data
  covars_x2<- intersect(covars_x, names(data_x))
  data_x<- data_x %>%   dplyr::mutate(treatment = !!sym(treatment_col_x)) %>%  tidyr::drop_na(all_of(covars_x2)) %>% tibble::remove_rownames()
  
  ## Initial values
  data_x_filter <- as.data.table(data_x)[, from := .I]
  block_size <- 10000
  coords <- data_x_filter[, c("x","y")]
  band_guide <- as.data.table(data_x_filter)[, `:=`(from = .I, block = ((.I - 1L) %/% block_size) + 1L)]
  band_data <- band_guide[, {
    fr  <- dbscan::frNN(x = coords, query = coords[from, , drop = FALSE], eps = distance_x, sort = FALSE)
    data.table(from     = rep.int(from, lengths(fr$id)), to = unlist(fr$id,   use.names = FALSE), distance = unlist(fr$dist, use.names = FALSE))
  }, by = block][from != to] 
  
  ## Loop
  current_from_indices <- data_x_filter$from
  band_data_current <- band_data[from %in% current_from_indices & to %in% current_from_indices]
  data_x_band <- band_guide %>% dplyr::filter(from %in% band_data_current[, unique(from)]) %>% dplyr::mutate(nb_from = seq_len(n()))
  
  map_indices <- data_x_band[, .(from, nb_from, treatment)]
  band_data2 <- band_data_current[map_indices, on = "from", nomatch = 0L][map_indices, on = c(to = "from"), nomatch = 0L] %>% 
    dplyr::rename(nb_to = i.nb_from) %>% dplyr::select(-to, -from, -i.treatment)
  
  pairs_test<- table(data_x_band$treatment) %>% {.[.>0]}; 
  
  form_sp_autocor<- as.formula( paste0(treatment_col_x, "~", paste0(c(covars_x2), collapse = "+")) )
  glm_initial <- glm(form_sp_autocor, family = binomial(), data = data_x_band, na.action = na.pass)
  data_x_band[,"qresiduals"]<- data_driven_qres(glm_initial)
  
  h <- max(band_data2$distance)
  band_data2[, dist_kernel_gauss := exp(-(distance / h)^2)]
  
  W_sparse <- sparseMatrix(i = band_data2$nb_from, j = band_data2$nb_to, x = band_data2$dist_kernel_gauss, dims =  rep(nrow(data_x_band), 2) )
  nodes_W <- nrow(W_sparse); row_mean_W <- Matrix::rowMeans(W_sparse); col_mean_W <- Matrix::colMeans(W_sparse)
  global_mean_W <- mean(W_sparse@x); ones_vec <- rep(1, nodes_W)
  
  centered_operator <- function(x, args = NULL) {
    sum_x <- sum(x); weighted_sum <- sum(col_mean_W * x); Wx <- W_sparse %*% x
    as.numeric(Wx - row_mean_W * sum_x - ones_vec * weighted_sum + global_mean_W * sum_x)
  }
  
  eig_res <- RSpectra::eigs_sym(A = centered_operator, k = min(20, nodes_W-1), n = nodes_W, which = "LM")
  sdj <- apply(eig_res$vectors, 2, sd); 
  keep_col<- (sdj >= 1e-10 & eig_res$values>0)
  W_meig <- list(sf = scale(eig_res$vectors[,keep_col], center=TRUE, scale=FALSE), ev = eig_res$values[keep_col],other = list(interact = FALSE, fast = FALSE, coords_z = NULL))
  
  mem_data<- {if(ncol(W_meig$sf)>1){
    fit_fast <- lm.fit( y= data_x_band$qresiduals, x=eig_res$vectors)
    as.numeric(fit_fast$fitted.values) 
  }else {eig_res$vectors}} %>% as.data.frame() %>% setNames(paste0("MEM_",  distance_x , "_", ncol(.)))
  
  data_x_band<- cbind(dplyr::select(data_x_band, c(names(data_x), "nb_from")), mem_data)  %>% as.data.frame()  
  names(data_x_band) <- make.unique(names(data_x_band), sep = "_")
  mem_cols <- grep("^MEM", names(data_x_band), value = TRUE)
  
  form_sp_autocor2<- as.formula( paste0(treatment_col_x, "~", paste0(c(covars_x2, mem_cols), collapse = "+")) )
  
  glm_MEM <- glm(form_sp_autocor2, family = binomial(), data = data_x_band, na.action = na.pass)
  data_x_band[,"qresiduals"]<-  data_driven_qres(glm_MEM)
  
  n <- nrow(data_x_band)
  z <- data_x_band$qresiduals 
  
  m <- mean(z, na.rm = TRUE)
  s <- sd(z, na.rm = TRUE)
  zscal <- ((z - m) / s) / sqrt((n - 1) / n)
  
  rowsum_vec <- band_data2[, sum(dist_kernel_gauss), by = nb_from][order(nb_from)]$V1
  band_data2[, w_W := dist_kernel_gauss / rowsum_vec[nb_from]]
  Wz <- band_data2[, .(Wz = sum(w_W * zscal[nb_to])), by = nb_from]
  
  Ww <- Matrix::sparseMatrix(i = band_data2$nb_from, j = band_data2$nb_to, x = band_data2$w_W, dims = c(n, n))
  w<- sum(band_data2$w_W)
  
  I_point <- mean(zscal[Wz$nb_from] * Wz$Wz)
  
  perm <- replicate(resamp, { zi <- zscal[sample.int(n)]
  drop(crossprod(zi, Ww %*% zi)) / w
  })
  
  p_greater <- (sum(perm >= I_point) + 1L) / (resamp + 1L)
  p_greater_IC <- binom.test(sum(perm >= I_point) + 1L, resamp + 1L)$conf.int
  
  moran_data <- data.frame(
    Moran.I.statistic = I_point,
    Expectation = mean(perm, na.rm = TRUE),
    Variance = var(perm, na.rm = TRUE),
    pval = min(p_greater_IC),
    distance = distance_x
  )
  
  moran_data
}


# moran_plot ####
# What it does: Compares initial and adjusted global Moran's I values across distance bands.
# Inputs:
# - list_autocor: output from reduce_autocor_by_distance().
# - treatment_col_x: treatment column name.
# - covars_x: covariates used in the treatment GLM.
# - distances: vector of distance bands to evaluate.
# Output:
# - ggplot object with initial and adjusted Moran's I by distance.
moran_plot <- function(list_autocor, treatment_col_x, covars_x, distances = c(600, 1000, 5000)){
  
  list_iterations<- list(list_autocor$PUID_initial, list_autocor$PUID_noautocor) %>% 
    setNames(c("1", list_autocor$listw_v))
  
  test_moran_list<- lapply(names(list_iterations), function(i) {
    
    pblapply(distances, function(d) {
      moran_by_distance(data_x = list_iterations[[i]], treatment_col_x = treatment_col_x, covars_x = covars_x,
                          distance_x = d) 
    }) %>% plyr::rbind.fill() %>% dplyr::mutate(iteration= i)
    
  })
  
  test_moran<- plyr::rbind.fill(test_moran_list) %>% 
    arrange(iteration, distance)  %>% dplyr::mutate(sign= ifelse(pval <= 0.05 |is.na(pval), "Sign.", "No sign."), iteration= as.numeric(iteration)) %>% 
    dplyr::mutate(name= if_else(iteration==1, "Initial Autocor", "Adjusted Autocor")) %>% 
    dplyr::mutate(name= factor(name, levels= unique(.$name)))
  
  
  moran_plot<-    ggplot(data= test_moran, aes(x = distance/1000, y = Moran.I.statistic, color= iteration, group= iteration)) +
    geom_line(method="loess", alpha=0.5)+
    geom_point(alpha= 0.5, size = 3, aes(shape= sign))+
    scale_shape_manual(values = c("Sign." = 17, "No sign." = 16), "Moran\np value") +
    scale_color_gradient(low = "red", high = "darkblue", breaks= unique(test_moran$iteration), "Pruning\niteration")+
    scale_x_log10(breaks = function(lims) {b <- scales::breaks_pretty(n = 5)(lims); b <- b[b != 0]   ;sort(  unique(c(b,0.5, 1, 5))    )   })+
    labs(x = "Distance (km)", y = "Global Moran's I Statistic")+
    theme(
      axis.text.x = element_text(angle= 45,hjust = 1, vjust = 1),
      axis.line.x = element_line(color = "black", linewidth=0.5),
      axis.line.y = element_line(color = "black", linewidth=0.5))+
    facet_wrap(~ name)+
    coord_cartesian(ylim = c(0,1))
  
  moran_plot
  
}


# reduce_autocor_by_distance ####
# What it does: Iteratively removes locally autocorrelated pixels and adds a MEM covariate for residual spatial structure.
# Inputs:
# - data_x: data with PUID, x, y, treatment column, and covariates.
# - distance_x: neighborhood radius used to build spatial weights.
# - treatment_col_x: treatment column name.
# - covars_x: covariates used in the treatment GLM.
# Output:
# - List with initial data, autocorrelation-reduced data, and final pruning iteration.
reduce_autocor_by_distance <- function(data_x, distance_x, treatment_col_x, covars_x){
  
  #parameters
  resamp <- 999
  list_prev<- list()
  
  ## organize data
  covars_x2<- intersect(covars_x, names(data_x))
  data_x<- data_x %>%   dplyr::mutate(treatment = !!sym(treatment_col_x)) %>%  tidyr::drop_na(all_of(covars_x2)) %>% tibble::remove_rownames()
  
  ## Initial values
  data_x_filter <- as.data.table(data_x)[, from := .I]
  block_size <- 10000
  coords <- data_x_filter[, c("x","y")]
  band_guide <- as.data.table(data_x_filter)[, `:=`(from = .I, block = ((.I - 1L) %/% block_size) + 1L)]
  band_data <- band_guide[, {
    fr  <- dbscan::frNN(x = coords, query = coords[from, , drop = FALSE], eps = distance_x, sort = FALSE)
    data.table(from     = rep.int(from, lengths(fr$id)), to = unlist(fr$id,   use.names = FALSE), distance = unlist(fr$dist, use.names = FALSE))
  }, by = block][from != to] 
  
  
  ## Loop
  for(i in 1:nrow(data_x_filter)){
    
    current_from_indices <- data_x_filter$from
    band_data_current <- band_data[from %in% current_from_indices & to %in% current_from_indices]
    data_x_band <- band_guide %>% dplyr::filter(from %in% band_data_current[, unique(from)]) %>% dplyr::mutate(nb_from = seq_len(n()))
    
    map_indices <- data_x_band[, .(from, nb_from, treatment)]
    band_data2 <- band_data_current[map_indices, on = "from", nomatch = 0L][map_indices, on = c(to = "from"), nomatch = 0L] %>% 
      dplyr::rename(nb_to = i.nb_from) %>% dplyr::select(-to, -from, -i.treatment)
    
    pairs_test<- table(data_x_band$treatment) %>% {.[.>0]}; print("pairs_test"); print(pairs_test)
    if(length(pairs_test)<2){break}
    
    form_sp_autocor<- as.formula( paste0(treatment_col_x, "~", paste0(c(covars_x2), collapse = "+")) )
    glm_initial <- glm(form_sp_autocor, family = binomial(), data = data_x_band, na.action = na.pass)
    data_x_band[,"qresiduals"]<- data_driven_qres(glm_initial)
    
    
    
    h <- max(band_data2$distance)
    band_data2[, dist_kernel_gauss := exp(-(distance / h)^2)]
    
    W_sparse <- sparseMatrix(i = band_data2$nb_from, j = band_data2$nb_to, x = band_data2$dist_kernel_gauss, dims =  rep(nrow(data_x_band), 2) )
    nodes_W <- nrow(W_sparse); row_mean_W <- Matrix::rowMeans(W_sparse); col_mean_W <- Matrix::colMeans(W_sparse)
    global_mean_W <- mean(W_sparse@x); ones_vec <- rep(1, nodes_W)
    
    centered_operator <- function(x, args = NULL) {
      sum_x <- sum(x); weighted_sum <- sum(col_mean_W * x); Wx <- W_sparse %*% x
      as.numeric(Wx - row_mean_W * sum_x - ones_vec * weighted_sum + global_mean_W * sum_x)
    }
    
    eig_res <- RSpectra::eigs_sym(A = centered_operator, k = min(20, nodes_W-1), n = nodes_W, which = "LM")
    sdj <- apply(eig_res$vectors, 2, sd); 
    keep_col<- (sdj >= 1e-10 & eig_res$values>0)
    W_meig <- list(sf = scale(eig_res$vectors[,keep_col], center=TRUE, scale=FALSE), ev = eig_res$values[keep_col],other = list(interact = FALSE, fast = FALSE, coords_z = NULL))
    
    mem_data<- {if(ncol(W_meig$sf)>1){
      fit_fast <- lm.fit( y= data_x_band$qresiduals, x=eig_res$vectors)
      as.numeric(fit_fast$fitted.values) 
    }else {eig_res$vectors}} %>% as.data.frame() %>% setNames(paste0("MEM", ncol(.)))
    
    
    
    mem_cols <- grep("^MEM", names(data_x_band), value = TRUE)
    
    form_sp_autocor2<- as.formula( paste0(treatment_col_x, "~", paste0(c(covars_x2, mem_cols), collapse = "+")) )
    
    
    
    data_x_band<- cbind(dplyr::select(data_x_band, c(names(data_x), "nb_from")), mem_data)  %>% as.data.frame()  
    
    glm_MEM <- glm(form_sp_autocor2, family = binomial(), data = data_x_band, na.action = na.pass)
    data_x_band[,"qresiduals"]<-  data_driven_qres(glm_MEM)
    
    
    n <- nrow(data_x_band)
    z <- data_x_band$qresiduals 
    
    m <- mean(z, na.rm = TRUE)
    s <- sd(z, na.rm = TRUE)
    zscal <- ((z - m) / s) / sqrt((n - 1) / n)
    
    
    rowsum_vec <- band_data2[, sum(dist_kernel_gauss), by = nb_from][order(nb_from)]$V1
    band_data2[, w_W := dist_kernel_gauss / rowsum_vec[nb_from]]
    Wz <- band_data2[, .(Wz = sum(w_W * zscal[nb_to])), by = nb_from]
    
    Ww <- Matrix::sparseMatrix(i = band_data2$nb_from, j = band_data2$nb_to, x = band_data2$w_W, dims = c(n, n))
    w<- sum(band_data2$w_W)
    
    print(w)    
    
    
    
    
    
    I_point <- mean(zscal[Wz$nb_from] * Wz$Wz)
    print("moran"); print(I_point)
    
    list_prev[[i]]<- list()
    list_prev[[i]]$moran<- data.frame(moran= I_point, pval=NA, it= i)
    list_prev[[i]]$data<- data_x_filter
    
    if(I_point<= 0.01){
      
      perm <- replicate(resamp, { zi <- zscal[sample.int(n)]
      drop(crossprod(zi, Ww %*% zi)) / w
      })
      
      p_greater <- sum(perm >= I_point)/(resamp + 1)
      p_greater_IC <- binom.test(sum(perm >= I_point) + 1L, resamp + 1L)$conf.int
      print("moran_p"); print(p_greater_IC)
      
      list_prev[[i]]$moran$pval<- min(p_greater_IC)
      
      if(min(p_greater_IC)>0.05){break}
    }
    
    
    
    
    ## moran local
    local_I <- band_data2[, .(lisa_I = zscal[nb_from] * sum(w_W * zscal[nb_to])), by = .(nb_from, treatment)]
    local_I_perm <- replicate(resamp, { zi <- zscal[sample.int(n)]; zi * as.numeric(Ww %*% zi) })
    
    ge <- sweep(local_I_perm, 1, local_I$lisa_I, `>=`)
    local_I[, lisa_p := (rowSums(ge) + 1) / (resamp + 1)]
    local_I[, lisa_q := p.adjust(lisa_p, method = "BH")]
    
    
    pixels_autocor<- filter(local_I, !is.na(lisa_I)  & lisa_q <= 0.05)
    
    if(nrow(pixels_autocor)<1){
      print("thr_I")
      
      thr_I <- quantile(abs(local_I$lisa_I), 0.99, na.rm=TRUE)
      S1 <- local_I[!is.na(lisa_I) & abs(lisa_I) >= thr_I]
      # thr_q <- quantile(S1$lisa_q, 0.05, na.rm = TRUE)
      # pixels_autocor <- dplyr::filter(S1, lisa_q <= thr_q)    
      pixels_autocor<- filter(S1, lisa_p <= 0.05)
    }
    if(nrow(pixels_autocor)<1){break}
    
    
    ### Pixels maximum autocor by neigborhood    
    nb_autocor<- band_data2[pixels_autocor, on = .(nb_from), nomatch = 0L]
    
    edges <- unique(nb_autocor[, .(from = nb_from, to = nb_to)])
    g <- graph_from_data_frame(edges, directed = FALSE)
    comp <- components(g)
    groups_dt <- data.table(
      nb_from  = as.integer(names(comp$membership)),
      group_id = comp$membership
    )
    nb_autocor2 <- as.data.table(nb_autocor)[groups_dt, on = "nb_from", nomatch = 0L]
    
    node_ids <- as.integer(V(g)$name)
    deg_vec  <- degree(g)
    
    deg_dt <- data.table(
      nb_from = node_ids,
      deg     = deg_vec
    )
    
    nb_autocor3 <- deg_dt[nb_autocor2, on = "nb_from"]
    nb_autocor3[, score := abs(lisa_I) * deg]
    
    
    if(nrow(nb_autocor)<1){break}
    
    ids_autocor_node<- nb_autocor3  %>% dplyr::arrange(treatment, dplyr::desc(score), desc(abs(lisa_I)), lisa_q, lisa_p) %>% dplyr::group_by(nb_to) %>% slice(1)
    pixels_autocor2 <- data_x_band  %>% dplyr::filter(!duplicated(PUID))  %>%  dplyr::filter(nb_from %in% ids_autocor_node$nb_from) 
    
    
    
    
    tabx<- table(pixels_autocor2[["treatment"]])
    print("autocor"); print(tabx)
    
    
    data_x_filter<- data_x_filter %>% dplyr::filter(!PUID %in% pixels_autocor2$PUID)
    print("no_autocor");
    print(table(data_x_filter[[treatment_col_x]]))
    
  }
  
  
  PUID_noautocor<- data_x_filter %>% dplyr::select(names(data_x)) %>% dplyr::left_join(dplyr::select(data_x_band, c("PUID", "MEM1")), by = "PUID")
  PUID_noautocor$MEM1[is.na(PUID_noautocor$MEM1)] <- 0
  
  
  
  
  # Results
  return(list(PUID_initial= data_x, PUID_noautocor=PUID_noautocor,  listw_v= i ))  
}


# counterfactual_function ####
# What it does: Estimates propensity scores, runs nearest-neighbor matching, and prunes matched pairs with high imbalance.
# Inputs:
# - data_x: analysis data after autocorrelation adjustment.
# - treatment_col_x: treatment column name.
# - covars_x: covariates used for propensity-score estimation.
# - exact_factor_vars: optional factor covariates for exact matching.
# - caliper_test: maximum propensity-score caliper.
# - replace: whether controls can be reused in matching.
# Output:
# - List with matched data, discarded matched rows, balance table, balance plot, propensity-score plot, and iteration summary.
counterfactual_function <- function(data_x, treatment_col_x, covars_x, exact_factor_vars= NULL, caliper_test=0.1, replace= F) {
  
  data_x2<- data.table(data_x) 
  
  #### Set treatment = 1 as the minority class 
  tab_treatment <- table(data_x2[[treatment_col_x]])
  if(names(which.max(tab_treatment))=="1"){
    data_x3<- data_x2
    data_x3[,treatment_col_x] <- ifelse(data_x2[[treatment_col_x]] == 1, 0, 1)
  } else {data_x3<-data_x2}
  
  formula_ps<- as.formula( paste0(treatment_col_x, "~", paste0(covars_x, collapse = "+")) )
  data_x_test<- data_x3 
  
  # Estimate propensity scores
  ps <- glm(formula_ps, data =  data_x_test , family = binomial())
  
  data_x_test2<- data_x_test %>% 
    dplyr::filter(complete.cases(model.frame(formula_ps, data = ., na.action = na.pass))) %>% 
    dplyr::mutate(ps = predict(ps, type = "response") )
  
  exact_vars<- {if(is.null(exact_factor_vars)){NULL} else {
    intersect(exact_factor_vars, names(data_x_test2)[sapply(data_x_test2, is.factor)]) %>%  paste0(collapse = "+")    }}
  if (length(exact_vars) == 0 || all(exact_vars == "")) {exact_vars <- NULL}
  
  
  
  
  
  # Matching
  matching_test <- matchit(formula_ps, data= data_x_test2, method= "nearest", ratio = 1, distance= data_x_test2$ps, 
                           exact= exact_vars,
                           replace = replace,
                           caliper = caliper_test,  std.caliper = F,
                           m.order = "largest", discard = "both")
  
  matched_df <- match.data(matching_test)
  
  if(names(which.max(tab_treatment))=="1"){matched_df[[treatment_col_x]] <- ifelse(matched_df[[treatment_col_x]] == 1, 0, 1)}
  
  prev_mean_imb<- Inf
  matched_df2<- matched_df
  
  balance_improvement_list<- list()
  
  for(j in seq(nrow(matched_df2))){
    
    # Balance improvement
    list_versions <- list(PreMatch=data_x_test2, PosMatch= matched_df2)
    
    balance_improvement <- imap(list_versions, function(v, name) {
      
      covars_imb<- intersect(names(v), covars_x)
      levels_count<- sapply(v[, ..covars_imb], function(x) {if (is.factor(x)) nlevels(droplevels(x)) else NA}) %>% {.[!is.na(.) & . < 2]} %>% names()
      
      if(length(levels_count)>0){ covars_imb<- covars_imb[!covars_imb %in% levels_count]}
      
      balance_data<- bal.tab(x = v[, ..covars_imb], treat = v[[treatment_col_x]], un = TRUE)[["Balance"]] %>% as.data.frame() %>% 
        mutate(!!name := abs(Diff.Un)*100) %>% rownames_to_column("variable") %>% dplyr::select(variable, !!name)
      
      
      if(length(levels_count)>0){ 
        balance_data<- plyr::rbind.fill(c(list(balance_data), lapply(levels_count, function(j) data.frame(variable= paste0(j, "_", unique(v[[j]]))))))
      } else {
        balance_data
      }
      
    }) %>% plyr::join_all() %>% mutate(across(everything(), ~replace_na(., 0))) %>% 
      dplyr::mutate(improve_abs = PosMatch-PreMatch)
    
    balance_improvement<- balance_improvement %>% dplyr::mutate(PreMatch= if_else(PreMatch>100,100,PreMatch))
    
    
    caliper_used <- matched_df2 %>% dplyr::group_by(subclass) %>% 
      dplyr::summarise(caliper_pair = abs(ps[treatment == 1] - ps[treatment == 0])) %>%
      dplyr::summarise( mean_calliper= mean(caliper_pair), max_calliper= max(caliper_pair))
    
    
    balance_improvement_list[[j]]<- balance_improvement %>% dplyr::mutate(n_pairs= n_distinct(matched_df2$subclass)) %>% cbind(caliper_used)
    
    current_imbalance<- dplyr::filter(balance_improvement, PosMatch>25)
    
    if(nrow(current_imbalance)>0){
      
      print(nrow(matched_df2)); print(current_imbalance)
      current_mean_imb <- mean(current_imbalance$PosMatch, na.rm=TRUE)
      
      vars_fix <- intersect(current_imbalance$variable, names(matched_df2))
      
      # 2. Calcular "Puntaje de Culpa" (Badness Score) para cada fila
      # LÃ³gica: Si Media_Tratados > Media_Control, penaliza a Tratados altos y Controles bajos (y viceversa).
      badness_matrix <- sapply(vars_fix, function(v) {
        x <- matched_df2[[v]]
        # DirecciÃ³n del desbalance: 1 si Tratados > Control, -1 si Control > Tratados
        dir <- sign(mean(x[matched_df2[[treatment_col_x]]==1]) - mean(x[matched_df2[[treatment_col_x]]==0]))
        # Pone puntaje alto a quien contribuya a esa direcciÃ³n
        scale(x) * dir * ifelse(matched_df2[[treatment_col_x]] == 1, 1, -1)
      })
      
      # Sumar culpas (rowSums maneja si hay 1 o varias variables)
      total_badness <- if(length(vars_fix) > 1) rowSums(badness_matrix) else badness_matrix
      
      # 3. Identificar y eliminar al peor 1% de los subclasses (o mÃ­nimo 1)
      n_bad <- max(1, ceiling(nrow(matched_df2) * 0.01))
      worst_subs <- unique(matched_df2$subclass[order(total_badness, decreasing = TRUE)[1:n_bad]])
      
      matched_df2 <- matched_df2[!matched_df2$subclass %in% worst_subs, ]
      
      
    } else {break}
    
    
  }
  
  
  balance_long <- balance_improvement %>%
    pivot_longer(cols = c(PreMatch, PosMatch), names_to = "type", values_to = "value") %>% 
    dplyr::mutate(type= factor(type, levels = unique(.$type)), 
                  improve_abs= if_else(type== "PreMatch", NA, improve_abs))
  
  imbalance_colors<- balance_long %>% dplyr::filter(type == "PosMatch") %>%
    dplyr::mutate(direction_label = ifelse(improve_abs <= 0, "â†“ Decreased ", "â†‘ Increased"),
                  direction_color = ifelse(improve_abs <=0, "forestgreen", "red")) %>% 
    dplyr::mutate(value= if_else(value>100,100,0))
  
  imbalance_colors_pallete<- dplyr::select(imbalance_colors, direction_color, direction_label) %>% dplyr::distinct()
  
  
  # Imbalance plot
  vars_imbalance_plot<-   ggplot(data= balance_long)+
    geom_point(aes(x=value, y= variable, color= type), size= 1)+
    scale_color_manual("Match", values = c(PreMatch = "#4F79B6", PosMatch = "#D04F4F"), labels= c(PreMatch = "Unmatched", PosMatch = "Matched") )+
    geom_vline(aes(xintercept= 25),  size= 0.5, linetype="dashed", color = "black")+
    labs(x= "Index of covariate imbalance", y= "Variables")+
    ggnewscale::new_scale_color()+
    geom_label(data = imbalance_colors,
               aes(y= variable, x= Inf, label = paste0(round(improve_abs, 2), "%"),
                   color= direction_color),
               alpha= 0.5,label.size = 0,
               hjust= 1, vjust= -0.1, inherit.aes = F)+
    scale_color_manual("Imbalance\nchange",
                       values = setNames(imbalance_colors_pallete$direction_color, imbalance_colors_pallete$direction_color),
                       labels = setNames(imbalance_colors_pallete$direction_label, imbalance_colors_pallete$direction_color)
    ) +
    coord_cartesian(xlim = c(0,100)) +
    theme(panel.background = element_rect(fill = NA), panel.grid.major = element_line(color = "gray"),
          axis.line = element_line(size = 0.5, colour = "black") )
  
  # Propensity scores plot
  join_versions <- imap(list_versions, function(v, name) {
    dplyr::mutate(v, type= name)}) %>% plyr::rbind.fill() %>% dplyr::mutate(type= factor(type, levels = unique(.$type)))
  
  join_subsets<- join_versions %>% split(join_versions[[treatment_col_x]])
  
  ps_plot<- ggplot(join_versions, aes(x = ps)) +
    geom_histogram(
      data = join_subsets[["1"]],
      aes(y = after_stat(count), fill = "1", color = "1"),
      bins = 30,
      alpha = 0.6,
      position = "identity"
    ) +
    geom_histogram(
      data = join_subsets[["0"]],
      aes(y = -after_stat(count), fill = "0", color = "0"),
      bins = 30,
      alpha = 0.6,
      position = "identity"
    ) +
    facet_wrap(~ type, scales = "free_y") +
    scale_fill_manual(name = "Window Pixels",
                      values = c("1" = "goldenrod", "0" = "royalblue"),
                      labels = c("1" = "In", "0" = "Out")
    )+
    scale_color_manual(name = "Window Pixels",
                       values = c("1" = "goldenrod", "0" = "royalblue"),
                       labels = c("1" = "In", "0" = "Out")
    )+
    scale_y_continuous(labels = abs) + 
    labs(x= "Propensity Score", y= "Number of Units")+
    guides(size = "none", fill = guide_legend(title.position="top", title.hjust = 0.5))+
    theme(panel.background = element_rect(fill = NA), 
          panel.grid.major = element_line(color = "gray"),
          legend.position =  "bottom", text = element_text(size = 10),
          panel.border = element_rect(color = NA,fill = NA),
          axis.line = element_line(colour = "black", size = 0.05, linetype = "solid"))
  
  
  matched_discard<-matched_df %>% dplyr::filter(!subclass %in% matched_df2$subclass)
  
  
  max_distance<- matched_df2 %>% dplyr::group_by(subclass) %>% 
    dplyr::summarise(diff= abs(diff(ps)))
  
  matching_interation_adjust <- balance_improvement_list %>% 
    imap(~ {
      
      imbalance_covars<- .x %>% dplyr::filter(PosMatch>25) %>% dplyr::mutate(label= paste0(variable, " (", round(PosMatch, 2), ")"))
      
      data.frame(iteration= .y,
                 n_pairs= unique(.x$n_pairs), mean_calliper= unique(.x$mean_calliper),  
                 max_calliper=unique(.x$max_calliper), 
                 n_imbalance_covars = nrow(imbalance_covars),
                 imbalance_vars= paste0(imbalance_covars$label, collapse = "; ")
      ) %>% cbind((dplyr::summarise(.x, across(c("PreMatch","PosMatch", "improve_abs"), mean, na.rm = TRUE)) %>% setNames(paste0("mean_", names(.)))))
      
    }) %>% plyr::rbind.fill()
  
  # Results
  return(list(matched_df=matched_df2, 
              matched_discard= matched_discard, 
              balance_improvement=balance_improvement, 
              vars_imbalance_plot=vars_imbalance_plot,
              ps_plot=ps_plot,
              matching_interation_adjust= matching_interation_adjust
  ))
  
}


# effectivenes_analysis_glm ####
# What it does: Summarizes pre/post outcome change by treatment and tests treatment differences with a GLM.
# Inputs:
# - data_x: matched analysis data.
# - treatment_col_x: treatment column name.
# - col_change_t1_x: baseline outcome column.
# - col_change_t2_x: final outcome column.
# - col_change: response variable used in the GLM significance test.
# - family_distribution: GLM family, for example binomial or gaussian.
# - col_site: optional site column retained for compatibility.
# - transform_value: optional multiplier applied to numeric effect summaries.
# Output:
# - List with data_change by treatment and summ_sign_change with GLM coefficient and effectiveness metrics.
effectivenes_analysis_glm <- function(data_x, treatment_col_x, col_change_t1_x, col_change_t2_x,
                                        col_change = "deforest",
                                        family_distribution = "binomial",
                                        col_site = NULL,
                                        transform_value = NULL) {
  
  data_x <- as.data.frame(data_x)
  
  list_treatments <- {
    if(treatment_col_x %in% names(data_x)) {
      data_x %>% split(.[, treatment_col_x])
    } else {
      setNames(list(data_x), treatment_col_x)
    }
  }
  
  list_treatments_change <- imap(list_treatments, function(x, name) {
    
    change_data <- x %>% 
      dplyr::summarise(
        mean_t1 = mean(!!sym(col_change_t1_x), na.rm = TRUE),
        mean_t2 = mean(!!sym(col_change_t2_x), na.rm = TRUE)
      ) %>% 
      dplyr::mutate(change_t1_t2 = mean_t1 - mean_t2)
    
    CI_test <- if(nrow(x) > 2) {
      DescTools::MeanDiffCI(
        x = x[[col_change_t1_x]],
        y = x[[col_change_t2_x]],
        conf.level = 0.95,
        sides = "two.sided",
        paired = TRUE
      ) %>% t()
    } else {
      t(data.frame(rep(NA, 3)))
    } %>% 
      as.data.frame() %>% 
      dplyr::mutate(across(everything(), ~ ifelse(is.nan(.), 0, .)))
    
    lims <- range(x[[col_change_t1_x]] - x[[col_change_t2_x]], na.rm = TRUE)
    
    CI_test <- CI_test %>% 
      as.data.frame() %>% 
      setNames(c("change", "change_lwr", "change_upr")) %>% 
      dplyr::mutate(
        across(everything(), ~ ifelse(. < lims[1], lims[1], .)),
        across(everything(), ~ ifelse(. > lims[2], lims[2], .))
      )
    
    cbind(change_data, CI_test) %>% 
      dplyr::mutate(
        change = change_t1_t2,
        ndata = nrow(x),
        !!sym(treatment_col_x) := name
      )
  })
  
  data_change <- plyr::rbind.fill(list_treatments_change) %>% 
    dplyr::mutate(treatment = !!sym(treatment_col_x)) %>% 
    dplyr::relocate(tail(names(.), 2))
  
  if(all(c("1", "0") %in% names(list_treatments_change))) {
    
    z <- qnorm(0.975)
    
    n_vals <- data_x %>% 
      dplyr::group_by(!!sym(treatment_col_x)) %>% 
      dplyr::summarise(n = n_distinct(col_change)) %>% 
      {sum(.$n, na.rm = TRUE)}
    
    form_sign_change <- as.formula(
      paste0(
        if(col_change == "deforest2") {
          paste0("cbind(", col_change_t1_x, ", ", col_change_t2_x, ")")
        } else {
          col_change
        },
        "~",
        treatment_col_x
      )
    )
    
    model_change <- glm(
      form_sign_change,
      data = data_x,
      family = family_distribution
    )
    
    model_change_sum <- if(n_vals > 2) {
      broom.mixed::tidy(model_change, conf.int = TRUE)
    } else {
      broom.mixed::tidy(model_change, conf.int = FALSE) %>% 
        dplyr::mutate(
          df = df.residual(model_change),
          conf.low = estimate - z * std.error,
          conf.high = estimate + z * std.error
        )
    } %>% 
      dplyr::mutate(
        term = sub("[0-9]+$", "", term),
        p.value = if_else(df <= 0, NA_real_, p.value)
      ) %>% 
      dplyr::filter(term %in% treatment_col_x)
    
    models_test <- list(Fixed = model_change_sum)
    
    model_sign <- models_test %>%
      imap(function(x, name) {
        as.data.frame(x) %>% 
          setNames(gsub("\\.+", "", make.names(names(.), unique = TRUE))) %>% 
          dplyr::mutate(
            pvalue_round = ifelse(
              pvalue > 0.001,
              formatC(pvalue, format = "f", digits = 3),
              formatC(pvalue, format = "e", digits = 0)
            ),
            sign = dplyr::case_when(
              is.na(pvalue) ~ "",
              pvalue < 0.001 ~ "***",
              pvalue < 0.05 ~ "**",
              pvalue < 0.1 ~ "*",
              TRUE ~ "ns"
            )
          ) %>% 
          setNames(ifelse(names(.) == "term", "term", paste0("mod", name, "_", names(.))))
      }) %>% 
      plyr::join_all() %>% 
      setNames(gsub("\\.+", "", make.names(names(.), unique = TRUE)))
    
    delta_change_treatment <- list_treatments_change[["1"]][, "change_t1_t2"] -
      list_treatments_change[["0"]][, "change_t1_t2"]
    
    delta_change_control <- list_treatments_change[["0"]][, "change_t1_t2"] -
      list_treatments_change[["1"]][, "change_t1_t2"]
    
    change_treatment <- list_treatments_change[["1"]][, "change_t1_t2"]
    change_control <- list_treatments_change[["0"]][, "change_t1_t2"]
    
    AD_treatment <- change_treatment - change_control
    AD_control <- change_control - change_treatment
    
    RAD_treatment <- case_when(
      change_control > 0 ~ AD_treatment / change_control,
      change_control == 0 & change_treatment == 0 ~ 0,
      change_control == 0 & change_treatment > 0 ~ -change_treatment,
      TRUE ~ NA_real_
    )
    
    RAD_control <- case_when(
      change_treatment > 0 ~ AD_control / change_treatment,
      change_treatment == 0 & change_control == 0 ~ 0,
      change_treatment == 0 & change_control > 0 ~ -change_control,
      TRUE ~ NA_real_
    )
    
    max_den <- pmax(abs(change_treatment), abs(change_control))
    max_den <- if_else(max_den == 0, 100, max_den)
    
    RAD_adjusted_treatment <- (AD_treatment / max_den) * 100
    RAD_adjusted_control <- (AD_control / max_den) * 100
    
    data_change_sign <- dplyr::mutate(
      model_sign,
      delta_change_treatment = delta_change_treatment,
      change_treatment = change_treatment,
      AD_treatment = AD_treatment,
      RAD_treatment = RAD_treatment,
      RAD_adjusted_treatment = RAD_adjusted_treatment,
      delta_change_control = delta_change_control,
      change_control = change_control,
      AD_control = AD_control,
      RAD_control = RAD_control,
      RAD_adjusted_control = RAD_adjusted_control
    )
    
    if(!is.null(transform_value)) {
      data_change_sign <- data_change_sign %>% 
        dplyr::mutate(
          dplyr::across(
            where(is.numeric) & !matches("_stderror|pvalue", ignore.case = TRUE),
            ~ .x * transform_value
          )
        )
    }
    
  } else {
    data_change_sign <- NULL
  }
  
  return(list(data_change = data_change, summ_sign_change = data_change_sign))
}


# make_change_plot ####
# What it does: Plots treatment/control change bars with confidence intervals and a treatment-effect label.
# Inputs:
# - data_change: treatment-level change data with term, treatment, change_prop, change_prop_lwr, and change_prop_upr.
# - data_change_sign: model/effect summary with term and the metric columns used for labels/colors.
# - metric_change_plot: column used as the plotted effect label.
# - color_metric_change_plot: column whose sign controls label color.
# - model_plot: model prefix retained for compatibility with previous plotting code.
# - aes_x_plot: data.frame mapping term to x-axis labels.
# - aes_fill_plot: data.frame mapping treatment values to legend labels and colors.
# - geom_vline, xlab_title, ylab_title, long_sec_axis: retained plotting arguments.
# Output:
# - ggplot object.
make_change_plot <- function(data_change, data_change_sign,
                             metric_change_plot,
                             color_metric_change_plot,
                             model_plot = "modFixed",
                             aes_x_plot,
                             aes_fill_plot,
                             geom_vline = NULL,
                             xlab_title = "",
                             ylab_title = "",
                             long_sec_axis = 1.5) {
  
  data_change <- data_change %>% 
    dplyr::mutate(
      term = as.character(term),
      treatment = as.character(treatment)
    ) %>% 
    dplyr::filter(treatment %in% aes_fill_plot$treatment) %>% 
    dplyr::mutate(
      term = factor(term, levels = aes_x_plot$term),
      treatment = factor(treatment, levels = aes_fill_plot$treatment)
    )
  
  data_change_sign <- data_change_sign %>% 
    dplyr::mutate(
      term = as.character(term),
      term = factor(term, levels = aes_x_plot$term),
      label_change = paste0(round(!!rlang::sym(metric_change_plot), 2), "%"),
      label_change_color = dplyr::if_else(
        !!rlang::sym(color_metric_change_plot) < 0,
        "red",
        "darkgreen"
      )
    )
  
  dodge_x <- ggplot2::position_dodge(width = 0.8)
  
  y_top <- max(data_change$change_prop_upr, data_change$change_prop, na.rm = TRUE)
  y_bottom <- min(data_change$change_prop_lwr, data_change$change_prop, 0, na.rm = TRUE)
  y_pad <- (y_top - y_bottom) * 0.18
  
  change_plot <- ggplot2::ggplot(
    data_change,
    ggplot2::aes(x = term, y = change_prop, fill = treatment)
  ) +
    ggplot2::geom_col(
      position = dodge_x,
      width = 0.72,
      color = "gray70",
      linewidth = 0.25
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = change_prop_lwr, ymax = change_prop_upr),
      position = dodge_x,
      width = 0.18,
      linewidth = 0.45,
      color = "black"
    ) +
    ggplot2::geom_text(
      data = data_change_sign,
      ggplot2::aes(x = term, y = y_top + y_pad, label = label_change),
      inherit.aes = FALSE,
      color = data_change_sign$label_change_color,
      size = 3.5
    ) +
    ggplot2::scale_fill_manual(
      "",
      values = stats::setNames(aes_fill_plot$color_fill, aes_fill_plot$treatment),
      labels = stats::setNames(aes_fill_plot$label_fill, aes_fill_plot$treatment),
      na.translate = FALSE
    ) +
    ggplot2::scale_x_discrete(
      labels = stats::setNames(aes_x_plot$label_x, aes_x_plot$term)
    ) +
    ggplot2::coord_cartesian(
      ylim = c(y_bottom, y_top + y_pad * 1.8),
      clip = "off"
    ) +
    ggplot2::labs(x = xlab_title, y = ylab_title) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "bottom",
      text = ggplot2::element_text(size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text.y = ggplot2::element_text(angle = 90, vjust = 0.5),
      axis.line.x = ggplot2::element_line(color = "black"),
      axis.line.y = ggplot2::element_line(color = "black"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  change_plot
}


