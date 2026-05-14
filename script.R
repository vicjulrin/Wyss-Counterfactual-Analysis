# Set environment ####
packages_list <- c(
  "magrittr", "data.table", "dplyr", "tidyr", "tibble", "purrr", "plyr",
  "ggplot2", "ggdendro", "ggnewscale", "ggbeeswarm", "ggpubr", "scales",
  "RColorBrewer", "MatchIt", "cobalt", "igraph", "leaps", "pbapply",
  "rlang", "car", "DescTools", "broom.mixed", "dbscan", "Matrix",
  "RSpectra", "this.path"
)


## Check - install libraries ####
packagesPrev<- .packages(all.available = TRUE)
lapply(packages_list, function(x) {   if ( ! x %in% packagesPrev ) { install.packages(x, force=T)}    })

## Load libraries ####
lapply(packages_list, library, character.only = TRUE)

## Set directory work ####
dir_work<- this.path::this.path() %>% dirname()
if(F){
  path_enviroment<- file.path(dir_work, paste0(basename(dir_work), ".RData"))
  load(file = path_enviroment)
}


## Set inputs ####
setwd(dir_work)
path_data<- file.path(dir_work, "inputs", "data_MC.rds")
treatment_col<- "MC"
col_forest_t1<- "Forest2005"; col_forest_t2<- "Forest2021"
col_carbon_t1<- "Carbon2005"; col_carbon_t2<- "Carbon2021"
cor_threshold<- 0.65 # Correlation threshold

col_site <- "id_site"
col_area_site <- "area_site"



covars_list<- list(numeric= c("Anual_Prec","D17Set10","D17Set1000","D17Set10000","D17Set5000","D7Set10",
              "D7Set1000","D7Set10000","D7Set5000","defor_2000_2005_10km2","defor_2000_2005_1km2",
              "Departme_R","Dis_Def","Dis_Rivers","District_R",
              "Elevation","National_R","Pop2000","Pop2020","Prec_Seas","Slope",
              "Tra_Time00","Tra_Time15"),
              factor= c("Department", "Ecoregions"))
  

## Set functions ####
path_matching_functions<- file.path(dir_work, "functions", "functions.R")
source(path_matching_functions, encoding = "UTF-8")


# Script ####
## Load data ####
data_PUID<-readRDS(file= path_data)
covars_test<- unlist(covars_list)


## Variables selection ####
### Compress categorical covariates ####
covars_factor<- covars_list$factor

if(!is.null(covars_factor) | length(covars_factor)<1){
  
  data_PUID2<-data_PUID
  
  for(i in covars_factor){
    data_PUID2<- data_PUID2 %>%  mutate(!!sym(i) := factor(!!sym(i), exclude = NULL))
    dummy <- model.matrix(as.formula(paste0("~", i, "-1")), data=data_PUID2)
    pca_1 <- prcomp(dummy, center=TRUE, scale.=FALSE)$x[,1]
    data_PUID2<- data_PUID2 %>%  mutate(!!sym(i) := pca_1)
  }
} else {
  data_PUID2<- data_PUID
}

### Multicollinearity test ####
multicol_test<- multicol_analysis(data_x= data_PUID2, treatment_col_x= treatment_col,  cor_threshold_x= cor_threshold, 
                                  covars_x= covars_test, covars_x_key= c("defor_2000_2005_1km2", "defor_2000_2005_10km2"))

covars_no_multicol<- multicol_test$covars_no_multicol

### Model test ####
model_selection<- bestmodels_test(data_x= data_PUID2, treatment_col_x= treatment_col, 
                                  covars_x= covars_no_multicol, covars_x_key= c("defor_2000_2005_1km2", "defor_2000_2005_10km2"),
                                  criteria_x= "AIC")

covars_better_model<-as.character(model_selection[["list_models"]][[1]]$vars) 

### Spatial autocor ####
selected_vars<- c(col_forest_t1, covars_better_model)

data_adjust_autocor<- reduce_autocor_by_distance(data_x = data_PUID2, treatment_col_x = treatment_col,
                                                        covars_x = selected_vars, 
                                                        distance_x = 1000)

moran_plot_autocor<- moran_plot(list_autocor=data_adjust_autocor, treatment_col_x = treatment_col,
                                covars_x = selected_vars, distances= c(600,1000,5000))


###  Matching Analysis ####
matching_analysis <- counterfactual_function(data_x = data_adjust_autocor$PUID_noautocor,
                                              treatment_col_x = treatment_col,
                                              covars_x = selected_vars, caliper_test = 0.2)

###  Post-Matching Analysis ####

matched_data <- matching_analysis$matched_df %>% 
  dplyr::filter(!!sym(col_forest_t1) > 0) %>% 
  dplyr::mutate(deforest = ifelse(!!sym(col_forest_t1) == 1 & !!sym(col_forest_t2) == 0, 1, 0),
                forest_persistence = ifelse(!!sym(col_forest_t1) == 1 & !!sym(col_forest_t2) == 0, 0, 1),
                carbon_emmisions = !!sym(col_carbon_t1) - !!sym(col_carbon_t2)
  ) %>% as.data.frame() 



#### Pixel-Level Analysis ####

deforest_test <- effectivenes_analysis_glm(
  data_x = matched_data,
  treatment_col_x = treatment_col,
  col_change_t1_x = col_forest_t1,
  col_change_t2_x = col_forest_t2,
  col_change = "deforest",
  family_distribution = "binomial",
  transform_value = -1
)

deforest_change_plot_data <- deforest_test$data_change %>% 
  dplyr::mutate(term = treatment_col) %>% 
  dplyr::rename(change_prop = change,change_prop_lwr = change_lwr,change_prop_upr = change_upr)

deforest_change_sign_plot_data <- deforest_test$summ_sign_change %>%  dplyr::mutate(term = treatment_col)


aes_x_mc <- data.frame(
  term = treatment_col,
  label_x = "Mining concessions"
)


deforest_plot <- make_change_plot(
  data_change = deforest_change_plot_data,
  data_change_sign = deforest_change_sign_plot_data,
  metric_change_plot = "RAD_adjusted_treatment",
  color_metric_change_plot = "RAD_adjusted_treatment",
  model_plot = "modFixed",
  aes_x_plot = aes_x_mc,
  aes_fill_plot = data.frame(
    treatment = c("0", "1", "modFixed_estimate"),
    label_fill = c("Control post-matching", "Treatment", "Model coefficient"),
    color_fill = c("rosybrown1", "lightgoldenrodyellow", "#CCCCCC")
  ),
  xlab_title = "",
  ylab_title = "Forest change"
)

CO2emissions_test <- effectivenes_analysis_glm(
  data_x = matched_data,
  treatment_col_x = treatment_col,
  col_change_t1_x = col_carbon_t1,
  col_change_t2_x = col_carbon_t2,
  col_change = "carbon_emmisions",
  family_distribution = "gaussian",
  transform_value = -1
)

CO2emissions_change_plot_data <- CO2emissions_test$data_change %>% 
  dplyr::mutate(
    term = treatment_col,
    treatment = as.character(treatment),
    change = change_t1_t2
  ) %>% 
  dplyr::rename(
    change_prop = change,
    change_prop_lwr = change_lwr,
    change_prop_upr = change_upr
  )

CO2emissions_change_sign_plot_data <- CO2emissions_test$summ_sign_change %>% 
  dplyr::mutate(term = treatment_col)

CO2emissions_plot <- make_change_plot(
  data_change = CO2emissions_change_plot_data,
  data_change_sign = CO2emissions_change_sign_plot_data,
  metric_change_plot = "RAD_adjusted_treatment",
  color_metric_change_plot = "RAD_adjusted_treatment",
  model_plot = "modFixed",
  aes_x_plot = aes_x_mc,
  aes_fill_plot = data.frame(
    treatment = c("0", "1", "modFixed_estimate"),
    label_fill = c("Control post-matching", "Treatment", "Model coefficient"),
    color_fill = c("rosybrown1", "lightgoldenrodyellow", "#CCCCCC")
  ),
  xlab_title = "",
  ylab_title = "Carbon emissions"
)


#### Mean Pixel Analysis by Site ####
matched_data_site <- matched_data %>% 
  dplyr::filter(!is.na(!!sym(col_site)))

matched_subclasses <- matched_data %>% 
  dplyr::group_by(subclass) %>% dplyr::filter(dplyr::n_distinct(!!sym(treatment_col)) > 1) %>% 
  dplyr::ungroup() %>% dplyr::pull(subclass) %>% unique()

list_sites <- matched_data_site %>% 
  dplyr::filter(subclass %in% matched_subclasses,!!sym(treatment_col) == 1) %>% 
  split(.[[col_site]])

data_sites <- pbapply::pblapply(list_sites, function(site_data) {
  
  site_x <- dplyr::first(site_data[[col_site]])
  area_x <- dplyr::first(site_data[[col_area_site]])
  vals <- as.character(unique(site_data$subclass))
  
  matched_data_y <- matched_data %>% 
    dplyr::filter(subclass %in% vals)
  
  deforest_test_site <- tryCatch({
    effectivenes_analysis_glm(
      data_x = matched_data_y,treatment_col_x = treatment_col,
      col_change_t1_x = col_forest_t1,col_change_t2_x = col_forest_t2,col_change = "deforest",
      family_distribution = "binomial",transform_value = -1
    ) %>% 
      lapply(function(x) {
        dplyr::mutate(x,term = treatment_col,id_site = site_x,area_site = area_x)})}, error = function(e) NULL)
  
  CO2emissions_test_site <- tryCatch({
    effectivenes_analysis_glm(
      data_x = matched_data_y,
      treatment_col_x = treatment_col,
      col_change_t1_x = col_carbon_t1,
      col_change_t2_x = col_carbon_t2,
      col_change = "carbon_emmisions",
      family_distribution = "gaussian",
      transform_value = -1
    ) %>% 
      lapply(function(x) {
        dplyr::mutate(
          x,
          term = treatment_col,
          id_site = site_x,
          area_site = area_x
        )
      })
  }, error = function(e) NULL)
  
  list(
    deforest_test = deforest_test_site,
    CO2emissions_test = CO2emissions_test_site
  )
})

forest_sites_data_change <- purrr::map(data_sites, "deforest_test") %>% 
  purrr::compact() %>% purrr::map("data_change") %>% plyr::rbind.fill()

forest_sites_summ_sign_change <- purrr::map(data_sites, "deforest_test") %>% 
  purrr::compact() %>% purrr::map("summ_sign_change") %>% plyr::rbind.fill()

forest_sites_data_change <- forest_sites_data_change %>% 
  dplyr::mutate(treatment = as.character(treatment),change_prop = change,change_prop_lwr = change_lwr,change_prop_upr = change_upr)

mean_forest_sites <- forest_sites_data_change %>% 
  dplyr::filter(!is.na(change_prop)) %>% 
  dplyr::group_by(id_site, term, treatment) %>% 
  dplyr::summarise(change_prop = mean(change_prop, na.rm = TRUE),.groups = "drop") %>% 
  dplyr::group_by(term, treatment) %>% 
  dplyr::summarise(
    change_prop_lwr = ifelse(dplyr::n() > 1, t.test(change_prop)$conf.int[1], NA_real_),
    change_prop_upr = ifelse(dplyr::n() > 1, t.test(change_prop)$conf.int[2], NA_real_),
    change_prop = mean(change_prop, na.rm = TRUE),
    ndata = dplyr::n_distinct(id_site), .groups = "drop")

aes_x_sites <- mean_forest_sites %>% 
  dplyr::distinct(term, ndata) %>% 
  dplyr::mutate(label_x = paste0("Mining concessions\n(", ndata, ")"))

delta_mean_forest_sites <- mean_forest_sites %>% 
  dplyr::group_by(term) %>% 
  dplyr::summarise(
    change_mean_treatment = change_prop[treatment == "1"],
    change_mean_control = change_prop[treatment == "0"],
    RAD_adjusted_mean_treatment = (
      (change_mean_treatment - change_mean_control) /
        pmax(abs(change_mean_treatment), abs(change_mean_control), 1e-9)
    ) * -100,
    .groups = "drop"
  )

test_forest_sites_sign <- forest_sites_data_change %>% 
  split(.$term) %>% 
  purrr::imap(function(x, name) {
    effectivenes_analysis_glm(
      data_x = x,
      treatment_col_x = "treatment",
      col_change_t1_x = "mean_t1",
      col_change_t2_x = "mean_t2",
      col_change = "change_prop",
      family_distribution = "gaussian",
      transform_value = -1
    )[["summ_sign_change"]] %>% 
      dplyr::mutate(term = name)
  }) %>% 
  plyr::rbind.fill() %>% 
  list(delta_mean_forest_sites) %>% 
  plyr::join_all()

forest_mean_sites_plot <- make_change_plot(
  data_change = mean_forest_sites,
  data_change_sign = test_forest_sites_sign,
  metric_change_plot = "RAD_adjusted_mean_treatment",
  color_metric_change_plot = "RAD_adjusted_mean_treatment",
  model_plot = "modFixed",
  aes_x_plot = aes_x_sites,
  aes_fill_plot = data.frame(
    treatment = c("0", "1", "modFixed_estimate"),
    label_fill = c("Control post-matching", "Treatment", "Model coefficient"),
    color_fill = c("rosybrown1", "lightgoldenrodyellow", "#CCCCCC")
  ),
  xlab_title = "",
  ylab_title = "Forest loss (%)"
)


carbon_sites_data_change <- purrr::map(data_sites, "CO2emissions_test") %>% 
  purrr::compact() %>% 
  purrr::map("data_change") %>% 
  plyr::rbind.fill()

carbon_sites_data_change <- carbon_sites_data_change %>% 
  dplyr::mutate(
    treatment = as.character(treatment),
    change_prop = change_t1_t2,
    change_prop_lwr = change_lwr,
    change_prop_upr = change_upr
  )

mean_carbon_sites <- carbon_sites_data_change %>% 
  dplyr::filter(!is.na(change_prop)) %>% 
  dplyr::group_by(id_site, term, treatment) %>% 
  dplyr::summarise(change_prop = mean(change_prop, na.rm = TRUE), .groups = "drop") %>% 
  dplyr::group_by(term, treatment) %>% 
  dplyr::summarise(
    change_prop_lwr = ifelse(dplyr::n() > 1, t.test(change_prop)$conf.int[1], NA_real_),
    change_prop_upr = ifelse(dplyr::n() > 1, t.test(change_prop)$conf.int[2], NA_real_),
    change_prop = mean(change_prop, na.rm = TRUE),
    ndata = dplyr::n_distinct(id_site),
    .groups = "drop"
  )

delta_mean_carbon_sites <- mean_carbon_sites %>% 
  dplyr::group_by(term) %>% 
  dplyr::summarise(
    change_mean_treatment = change_prop[treatment == "1"],
    change_mean_control = change_prop[treatment == "0"],
    RAD_adjusted_mean_treatment = (
      (change_mean_control - change_mean_treatment) /
        pmax(abs(change_mean_treatment), abs(change_mean_control), 1e-9)
    ) * 100,
    .groups = "drop"
  )

test_carbon_sites_sign <- carbon_sites_data_change %>% 
  split(.$term) %>% 
  purrr::imap(function(x, name) {
    effectivenes_analysis_glm(
      data_x = x,
      treatment_col_x = "treatment",
      col_change_t1_x = "mean_t1",
      col_change_t2_x = "mean_t2",
      col_change = "change_prop",
      family_distribution = "gaussian",
      transform_value = -1
    )[["summ_sign_change"]] %>% 
      dplyr::mutate(term = name)
  }) %>% 
  plyr::rbind.fill() %>% 
  list(delta_mean_carbon_sites) %>% 
  plyr::join_all()

carbon_mean_sites_plot <- make_change_plot(
  data_change = mean_carbon_sites,
  data_change_sign = test_carbon_sites_sign,
  metric_change_plot = "RAD_adjusted_mean_treatment",
  color_metric_change_plot = "RAD_adjusted_mean_treatment",
  model_plot = "modFixed",
  aes_x_plot = aes_x_sites,
  aes_fill_plot = data.frame(
    treatment = c("0", "1", "modFixed_estimate"),
    label_fill = c("Control post-matching", "Treatment", "Model coefficient"),
    color_fill = c("rosybrown1", "lightgoldenrodyellow", "#CCCCCC")
  ),
  xlab_title = "",
  ylab_title = "Carbon emissions"
)

#### Histogram - sites ####
area_sites <- forest_sites_data_change %>% 
  dplyr::filter(treatment == "1") %>% 
  dplyr::distinct(term, id_site, area_site, ndata) %>% 
  dplyr::mutate(
    area_km2 = area_site / 1e6,
    term = as.character(term)
  )

hist_fill_values <- setNames("#8D8D8D", unique(area_sites$term))
hist_fill_labels <- setNames(
  as.character(aes_x_sites$label_x[match(unique(area_sites$term), aes_x_sites$term)]),
  unique(area_sites$term)
)

Hist_PixelsPerSiteArea_detail <- ggplot2::ggplot(
  area_sites,
  ggplot2::aes(x = area_km2, fill = term, weight = ndata)
) +
  ggplot2::geom_histogram(position = "identity", alpha = 1) +
  ggplot2::scale_fill_manual(
    "Governance\ntype",
    values = hist_fill_values,
    labels = hist_fill_labels,
    drop = FALSE
  ) +
  ggplot2::ylab("Pixel count") +
  ggplot2::xlab(expression(Area~(km^2))) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

Hist_SitesPerArea_detail <- ggplot2::ggplot(
  area_sites,
  ggplot2::aes(x = area_km2, fill = term)
) +
  ggplot2::geom_histogram(position = "identity", alpha = 1) +
  ggplot2::scale_fill_manual(
    "Governance\ntype",
    values = hist_fill_values,
    labels = hist_fill_labels,
    drop = FALSE
  ) +
  ggplot2::ylab("Sites count") +
  ggplot2::xlab(expression(Area~(km^2))) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

#### Boxplots Change by site ####

forest_sites_summ_sign_change <- forest_sites_summ_sign_change %>% 
  dplyr::mutate(
    term = as.character(term),
    area_km2 = area_site / 1e6
  )

site_size_breaks <- c(0, 10, 100, 1000, Inf)
site_size_labels <- c("< 10 km2", "10 - 100 km2", "100 - 1000 km2", ">= 1000 km2")

avoidDefor_sites <- forest_sites_summ_sign_change %>% 
  dplyr::transmute(
    term,
    id_site,
    area_km2,
    site_size = cut(
      area_km2,
      breaks = site_size_breaks,
      labels = site_size_labels,
      include.lowest = TRUE,
      right = FALSE
    ),
    avoid_forest_loss = AD_treatment * 100
  )

eff_avoidDefor_sites <- forest_sites_summ_sign_change %>% 
  dplyr::transmute(
    term,
    id_site,
    area_km2,
    site_size = cut(
      area_km2,
      breaks = site_size_breaks,
      labels = site_size_labels,
      include.lowest = TRUE,
      right = FALSE
    ),
    RAD_adjusted_treatment
  )

avoidDefor_sites_area <- avoidDefor_sites %>% 
  dplyr::filter(!is.na(site_size))

eff_avoidDefor_sites_area <- eff_avoidDefor_sites %>% 
  dplyr::filter(!is.na(site_size))

site_labels <- setNames(aes_x_sites$label_x, aes_x_sites$term)
site_fill <- setNames("#8D8D8D", unique(avoidDefor_sites$term))

avoidDefor_sites_detailPlot <- ggplot2::ggplot(
  avoidDefor_sites,
  ggplot2::aes(x = term, y = avoid_forest_loss, fill = term)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.2) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 1.5, width = 0.2, color = "gray40") +
  ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.4) +
  ggplot2::scale_fill_manual(values = site_fill, labels = site_labels, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = site_labels) +
  ggplot2::xlab("") +
  ggplot2::ylab("Avoided forest loss (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

avoidDefor_sites_generalPlot <- ggplot2::ggplot(
  avoidDefor_sites_area,
  ggplot2::aes(x = site_size, y = avoid_forest_loss)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.2) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 1.5, width = 0.2, color = "gray40") +
  ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.4, fill = "#8D8D8D") +
  ggplot2::xlab("Site area") +
  ggplot2::ylab("Avoided forest loss (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

eff_avoidDefor_sites_detailPlot <- ggplot2::ggplot(
  eff_avoidDefor_sites,
  ggplot2::aes(x = term, y = RAD_adjusted_treatment, fill = term)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.2) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 1.5, width = 0.2, color = "gray40") +
  ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.4) +
  ggplot2::scale_fill_manual(values = site_fill, labels = site_labels, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = site_labels) +
  ggplot2::xlab("") +
  ggplot2::ylab("Relative avoided forest loss (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
  )

eff_avoidDefor_sites_generalPlot <- ggplot2::ggplot(
  eff_avoidDefor_sites_area,
  ggplot2::aes(x = site_size, y = RAD_adjusted_treatment)
) +
  ggplot2::geom_hline(yintercept = 0, color = "gray50", linewidth = 0.2) +
  ggbeeswarm::geom_quasirandom(alpha = 0.5, size = 1.5, width = 0.2, color = "gray40") +
  ggplot2::geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.4, fill = "#8D8D8D") +
  ggplot2::xlab("Site area") +
  ggplot2::ylab("Relative avoided forest loss (%)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

avoidDefor_sites_plot <- ggpubr::ggarrange(
  avoidDefor_sites_detailPlot,
  avoidDefor_sites_generalPlot,
  ncol = 2
)

eff_avoidDefor_sites_plot <- ggpubr::ggarrange(
  eff_avoidDefor_sites_detailPlot,
  eff_avoidDefor_sites_generalPlot,
  ncol = 2
)



