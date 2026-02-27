library(tidyverse)
library(censusapi)
library(sf)
library(tidycensus)
library(augsynth)
library(ggcheck)
library(patchwork)
library(gt)

## Spending (Biweekly)

### BID

#### East Harlem

bid_data <- read_rds("data/bid_analysis_dat.rds")

harlem_bid_spend <- bid_data$h_spend %>%
  rename(name = f_all_bi_2)

harlem_bid_spend_u <- harlem_bid_spend %>%
  select(name, biweekly_marker, YEAR, prop_rest:treated)

bid_spend_dat <- read_rds("data/bid_spending_inf.rds") %>%
  rename(name = f_all_bi_2)

bid_census_dat <- read_rds("data/bid_census.rds") %>%
  rename(name = f_all_bi_2)

harlem_bid_spend_cl <- harlem_bid_spend_u %>%
  left_join(bid_spend_dat, by = c("name", "biweekly_marker")) %>%
  left_join(bid_census_dat, by = c("name"))

unadj_form_spend <- "median_total_spend ~ treated"

hspend_bid_ascm <- single_augsynth(form = as.formula(unadj_form_spend),
                                   unit = name,
                                   time = biweekly_marker,
                                   data = harlem_bid_spend_cl,
                                   scm = T, t_int = 0,
                                   progfunc = "ridge")

##### Adjusted Ridge aSCM

adj_form_spend_bid <- "median_total_spend ~ treated | med_age + prop_fam_pov + prop_unemploy + med_income + prop_hs_diploma + prop_rest + prop_groc + prop_gent + overall_mean_hosprate + arrests_area"

hspend_bid_ascm_adj <- single_augsynth(form = as.formula(adj_form_spend_bid),
                                       unit = name,
                                       time = biweekly_marker,
                                       data = harlem_bid_spend_cl,
                                       scm = T, t_int = 0,
                                       progfunc = "ridge",
                                       residualize = T)



#### Washington Heights

##### Ridge aSCM
  
wh_bid_spend <- bid_data$wh_spend %>%
  rename(name = f_all_bi_2)

wh_bid_spend_u <- wh_bid_spend %>%
  select(name, biweekly_marker, YEAR, prop_rest:treated)

wh_bid_spend_cl <- wh_bid_spend_u %>%
  left_join(bid_spend_dat, by = c("name", "biweekly_marker")) %>%
  left_join(bid_census_dat, by = c("name"))

whspend_bid_ascm <- single_augsynth(form = as.formula(unadj_form_spend),
                                    unit = name,
                                    time = biweekly_marker,
                                    data = wh_bid_spend_cl,
                                    scm = T, t_int = 0,
                                    progfunc = "ridge")

##### Adjusted Ridge aSCM

whspend_bid_ascm_adj <- single_augsynth(form = as.formula(adj_form_spend_bid),
                                        unit = name,
                                        time = biweekly_marker,
                                        data = wh_bid_spend_cl,
                                        scm = T, t_int = 0,
                                        progfunc = "ridge",
                                        residualize = T)


# Functions for permutation testing with spatial overlap checking

assign_treatment <- function(data, unit_name){
  out_df <- data %>%
    mutate(treated = 0) %>%
    mutate(treated = case_when(biweekly_marker >= 0 & name == unit_name ~ 1,
                               TRUE ~ 0))
  return(out_df)
}

create_permutation <- function(original_dat, unit_name, spatial_dat, buffer_size, custom_form,
                               progfunc_opt, residualize_opt, overlap = TRUE){
  new_dat <- assign_treatment(original_dat, unit_name)
  if(overlap == TRUE){
    new_overlaps <- check_for_overlaps(spatial_dat, unit_name, buffer_size)
  } else {new_overlaps <- NULL}
  permute_results <- run_augsynth(new_dat, new_overlaps,
                                  custom_form,
                                  progfunc_opt,
                                  residualize_opt)
  return(permute_results)
}

check_for_overlaps <- function(spatial_dat, unit_name, buffer_size){
  buffers <- spatial_dat %>%
    st_transform(crs = "EPSG:32118") %>%
    st_buffer(buffer_size)
  
  site_only <- buffers %>%
    filter(name == unit_name) %>%
    mutate(target_site = 1)
  
  overlaps <- st_join(buffers, site_only) %>%
    drop_na(target_site) %>%
    filter(name.x != unit_name) %>%
    mutate(name_to_remove = name.x)
  
  return(overlaps)
}

run_augsynth <- function(dataset, overlaps, formula, progfunc_opt, residualize_opt){
  df <- dataset %>%
    filter(name %in% overlaps$name_to_remove == FALSE)
  
  res <- single_augsynth(form = as.formula(formula),
                         unit = name,
                         time = biweekly_marker,
                         data = df,
                         scm = T, t_int = 0, 
                         progfunc = progfunc_opt,
                         residualize = residualize_opt)
  
  graph_data <- plot(res) %>%
    get_data()
  
  rmse_out <- graph_data %>%
    mutate(timepoint = case_when(Time >= 0 ~ "Post",
                                 Time < 0 ~ "Pre")) %>%
    group_by(timepoint) %>%
    summarize(RMSPE = sqrt(sum((Estimate)^2)/n())) %>%
    pivot_wider(names_from = timepoint,
                values_from = RMSPE) %>%
    mutate(ratio = Post/Pre) %>%
    mutate(graph_dat = list(graph_data))
  
  return(rmse_out)
} 

calculate_p <- function(df, obs_stat){
  p_val <- df %>%
    mutate(as_extreme = case_when(ratio >= obs_stat ~ 1,
                                  TRUE ~ 0)) %>%
    group_by(as_extreme) %>%
    summarize(n = n()) %>%
    mutate(per = n/sum(n)) %>%
    filter(as_extreme == 1) %>%
    pull(per)
  return(p_val) 
}

get_permutation_results <- function(fit_obj, rel_form, spatial_data, buffer_size,
                                    progfunc_opt, residualize_opt){
  
  treated_unit <- fit_obj$raw_data %>%
    filter(treated == 1) %>%
    pull(name) %>%
    unique()
  
  standard_res <- create_permutation(fit_obj$raw_data, treated_unit, spatial_data, buffer_size, rel_form,
                                     progfunc_opt, residualize_opt) %>%
    mutate(unit = treated_unit) %>%
    mutate(treatment = 1)
  
  permutation_subjects <- fit_obj$raw_data %>%
    filter(name != treated_unit) %>%
    pull(name) %>%
    unique()
  
  permutation_out <- map_dfr(permutation_subjects, function(a){
    create_permutation(fit_obj$raw_data, a, spatial_data, buffer_size, rel_form,
                       progfunc_opt, residualize_opt) %>%
      mutate(unit = a)
  }) %>%
    bind_rows(standard_res)
  
  perm_results <- permutation_out %>%
    mutate(p_val = map_dbl(ratio, function(a) calculate_p(permutation_out, a))) %>%
    ungroup() %>%
    mutate(sig_status = case_when(p_val < 0.05 ~ "Significant",
                                  p_val >= 0.05 ~ "Not Significant") %>%
             factor(levels = c("Significant", "Not Significant")))
  
}


### 10-Minutes

#### East Harlem

sites <- read_rds("data/site_cd_lookup.rds")

highres10min <- read_rds("data/highres10min.rds")

sites_to_remove <- sites %>% 
  filter(boro_cd %in% 110:112) %>% 
  filter(str_detect(name, "OnPoint") == FALSE)

site_data_10 <- read_rds("data/site_analysis_dat_10min.rds")

harlem_10_spend <- site_data_10$h_spend %>%
  filter(name %in% sites_to_remove$name == FALSE)

unadj_form_spend <- "median_total_spend ~ treated"

hspend_10_ascm <- single_augsynth(form = as.formula(unadj_form_spend), 
                                  unit = name,
                                  time = biweekly_marker,
                                  data = harlem_10_spend,
                                  scm = T, t_int = 0, 
                                  progfunc = "ridge")

rel_sites <- highres10min %>%
  filter(name %in% harlem_10_spend$name)

hspend_perm_10_ascm <- get_permutation_results(hspend_10_ascm, unadj_form_spend, rel_sites, 0,
                                     "ridge", TRUE)


##### Adjusted Ridge aSCM

adj_form_spend <- "median_total_spend ~ treated | med_age + prop_fam_pov + prop_unemploy + med_income + prop_hs_diploma + prop_rest + prop_groc + gentrifying + overall_mean_hosprate + n_arrests"

hspend_10_ascm_adj <- single_augsynth(form = as.formula(adj_form_spend), 
                                      unit = name,
                                      time = biweekly_marker,
                                      data = harlem_10_spend,
                                      scm = T, t_int = 0, 
                                      progfunc = "ridge",
                                      residualize = T)

hspend_perm_10_ascm_adj <- get_permutation_results(hspend_10_ascm_adj, adj_form_spend, rel_sites, 0, "ridge", TRUE)
 

#### Washington Heights

##### Ridge aSCM
  
wh_10_spend <- site_data_10$wh_spend %>% 
  filter(name %in% sites_to_remove$name == FALSE)

whspend_10_ascm <- single_augsynth(form = as.formula(unadj_form_spend),
                                   unit = name,
                                   time = biweekly_marker,
                                   data = wh_10_spend,
                                   scm = T, t_int = 0,
                                   progfunc = "ridge")

rel_sites <- highres10min %>%
  filter(name %in% wh_10_spend$name)

wh_perm_10_ascm <- get_permutation_results(whspend_10_ascm, unadj_form_spend, rel_sites, 0,
                                        "ridge", TRUE)

whspend_perm_10_ascm <- read_rds("data/whspend_perm_10_ascm.rds")


##### Adjusted Ridge aSCM

whspend_10_ascm_adj <- single_augsynth(form = as.formula(adj_form_spend),
                                       unit = name,
                                       time = biweekly_marker,
                                       data = wh_10_spend,
                                       scm = T, t_int = 0,
                                       progfunc = "ridge",
                                       residualize = T)

wh_spend_perm_10_ascm_adj <- get_permutation_results(whspend_10_ascm_adj, adj_form_spend, rel_sites, 0,
                                        "ridge", TRUE)


### 5-Minutes

#### East Harlem

sites <- read_rds("data/site_cd_lookup.rds")

highres5min <- read_rds("data/highres5min.rds")

sites_to_remove <- sites %>% 
  filter(boro_cd %in% 110:112) %>% 
  filter(str_detect(name, "OnPoint") == FALSE)

site_data_5 <- read_rds("data/site_analysis_dat_5min.rds")

harlem_5_spend <- site_data_5$h_spend %>% 
  filter(name %in% sites_to_remove$name == FALSE)

all_med_spending <- harlem_5_spend %>%
  group_by(name) %>%
  summarize(median_spend = median(median_total_spend)) %>%
  filter(median_spend <= 1840)

red_harlem_5_spend <- harlem_5_spend %>%
  filter(name %in% all_med_spending$name)

unadj_form_spend <- "median_total_spend ~ treated"

hspend_5_ascm <- single_augsynth(form = as.formula(unadj_form_spend), 
                                 unit = name,
                                 time = biweekly_marker,
                                 data = red_harlem_5_spend,
                                 scm = T, t_int = 0, 
                                 progfunc = "ridge")

rel_sites <- highres5min %>%
  filter(name %in% red_harlem_5_spend$name)

hspend_perm_5_ascm <- get_permutation_results(hspend_5_ascm, unadj_form_spend, rel_sites, 0,
                                              "ridge", TRUE)

##### Adjusted Ridge aSCM

adj_form_spend <- "median_total_spend ~ treated | med_age + prop_fam_pov + prop_unemploy + med_income + prop_hs_diploma + prop_rest + prop_groc + gentrifying + overall_mean_hosprate + n_arrests"


all_med_spending <- harlem_5_spend %>%
  group_by(name) %>%
  summarize(median_spend = median(median_total_spend)) %>%
  filter(median_spend <= 1840)

red_harlem_5_spend <- harlem_5_spend %>%
  filter(name %in% all_med_spending$name)

hspend_5_ascm_adj <- single_augsynth(form = as.formula(adj_form_spend), 
                                     unit = name,
                                     time = biweekly_marker,
                                     data = red_harlem_5_spend,
                                     scm = T, t_int = 0, 
                                     progfunc = "ridge",
                                     residualize = T)

rel_sites <- highres5min %>%
  filter(name %in% red_harlem_5_spend$name)

hspend_perm_5_ascm_adj <- get_permutation_results(hspend_5_ascm_adj, adj_form_spend, rel_sites, 0,
                                                  "ridge", TRUE)


#### Washington Heights


##### Ridge aSCM
  
wh_5_spend <- site_data_5$wh_spend %>%
  filter(name %in% sites_to_remove$name == FALSE)

all_med_spending_wh <- wh_5_spend %>%
  group_by(name) %>%
  summarize(median_spend = median(median_total_spend)) %>%
  filter(median_spend <= 816)

red_wh_5_spend <- wh_5_spend %>%
  filter(name %in% all_med_spending_wh$name)

whspend_5_ascm <- single_augsynth(form = as.formula(unadj_form_spend),
                                  unit = name,
                                  time = biweekly_marker,
                                  data = red_wh_5_spend,
                                  scm = T, t_int = 0,
                                  progfunc = "ridge")

rel_sites <- highres5min %>%
  filter(name %in% red_wh_5_spend$name)

wh_perm_5_ascm <- get_permutation_results(whspend_5_ascm, unadj_form_spend, rel_sites, 0,
                                          "ridge", TRUE)

##### Adjusted Ridge aSCM

whspend_5_ascm_adj <- single_augsynth(form = as.formula(adj_form_spend),
                                      unit = name,
                                      time = biweekly_marker,
                                      data = red_wh_5_spend,
                                      scm = T, t_int = 0,
                                      progfunc = "ridge",
                                      residualize = T)
rel_sites <- highres5min %>%
  filter(name %in% red_wh_5_spend$name)

wh_spend_perm_5_ascm_adj <- get_permutation_results(whspend_5_ascm_adj, adj_form_spend, rel_sites, 0,
                                                    "ridge", TRUE)


## Foot Traffic (Biweekly)

### BID

#### East Harlem

  ##### Ridge aSCM
  
harlem_bid_foot <- bid_data$h_foot %>%
  rename(name = f_all_bi_2)

unadj_form_foot <- "median_total_foot ~ treated"

hfoot_bid_ascm <- single_augsynth(form = as.formula(unadj_form_foot),
                                  unit = name,
                                  time = biweekly_marker,
                                  data = harlem_bid_foot,
                                  scm = T, t_int = 0,
                                  progfunc = "ridge")


##### Adjusted Ridge aSCM

adj_form_foot_bid <- "median_total_foot ~ treated | med_age + prop_fam_pov + prop_unemploy + med_income + prop_hs_diploma + prop_rest + prop_groc + prop_gent + overall_mean_hosprate + arrests_area"

harlem_bid_foot_cl <- harlem_bid_foot %>%
  left_join(bid_census_dat, by = c("name"))

hfoot_bid_ascm_adj <- single_augsynth(form = as.formula(adj_form_foot_bid),
                                      unit = name,
                                      time = biweekly_marker,
                                      data = harlem_bid_foot_cl,
                                      scm = T, t_int = 0,
                                      progfunc = "ridge",
                                      residualize = T)

#### Washington Heights
  
##### Ridge aSCM
  
wh_bid_foot <- bid_data$wh_foot %>%
  rename(name = f_all_bi_2)

whfoot_bid_ascm <- single_augsynth(form = as.formula(unadj_form_foot),
                                   unit = name,
                                   time = biweekly_marker,
                                   data = wh_bid_foot,
                                   scm = T, t_int = 0,
                                   progfunc = "ridge")


##### Adjusted Ridge aSCM

wh_bid_foot_cl <- wh_bid_foot %>%
  left_join(bid_census_dat, by = c("name"))

whfoot_bid_ascm_adj <- single_augsynth(form = as.formula(adj_form_foot_bid),
                                       unit = name,
                                       time = biweekly_marker,
                                       data = wh_bid_foot_cl,
                                       scm = T, t_int = 0,
                                       progfunc = "ridge",
                                       residualize = T)

### 10-Minutes
  
#### East Harlem

##### Ridge aSCM
  
unadj_form_foot <- "median_total_foot ~ treated"

harlem_10_foot <- site_data_10$h_foot %>% 
  filter(name %in% sites_to_remove$name == FALSE)

hfoot_10_ascm <- single_augsynth(form = as.formula(unadj_form_foot), 
                                 unit = name,
                                 time = biweekly_marker,
                                 data = harlem_10_foot,
                                 scm = T, t_int = 0, 
                                 progfunc = "ridge")

rel_sites <- highres10min %>%
  filter(name %in% harlem_10_foot$name)

hfoot_perm_10_ascm <- get_permutation_results(hfoot_10_ascm, unadj_form_foot, rel_sites, 0,
                                     "ridge", TRUE)

##### Adjusted Ridge aSCM

adj_form_foot <- "median_total_foot ~ treated | med_age + prop_fam_pov + prop_unemploy + med_income + prop_rest + prop_groc + gentrifying + prop_hs_diploma + overall_mean_hosprate + n_arrests"

hfoot_10_ascm_adj <- single_augsynth(form = as.formula(adj_form_foot), 
                                     unit = name,
                                     time = biweekly_marker,
                                     data = harlem_10_foot,
                                     scm = T, t_int = 0, 
                                     progfunc = "ridge",
                                     residualize = T)

hfoot_perm_10_ascm_adj <- get_permutation_results(hfoot_10_ascm_adj, adj_form_foot, rel_sites, 0,
                                     "ridge", TRUE)

##### Ridge aSCM
  
wh_10_foot <- site_data_10$wh_foot %>%
  filter(name %in% sites_to_remove$name == FALSE) %>%
  filter(name != "Long Island Jewish Medical Center OTP")


whfoot_10_ascm <- single_augsynth(form = as.formula(unadj_form_foot),
                                  unit = name,
                                  time = biweekly_marker,
                                  data = wh_10_foot,
                                  scm = T, t_int = 0,
                                  progfunc = "ridge")

rel_sites <- highres10min %>%
  filter(name %in% wh_10_foot$name)

wh_perm_10_ascm <- get_permutation_results(whfoot_10_ascm, unadj_form_foot, rel_sites, 0,
                                        "ridge", TRUE)

##### Adjusted Ridge aSCM

whfoot_10_ascm_adj <- single_augsynth(form = as.formula(adj_form_foot),
                                      unit = name,
                                      time = biweekly_marker,
                                      data = wh_10_foot,
                                      scm = T, t_int = 0,
                                      progfunc = "ridge",
                                      residualize = T)

wh_foot_perm_ascm_adj <- get_permutation_results(whfoot_10_ascm_adj, adj_form_foot, rel_sites, 0,
                                        "ridge", TRUE)

### 5-Minutes

#### East Harlem

##### Ridge aSCM
  
harlem_5_foot <- site_data_5$h_foot %>%
  filter(name %in% sites_to_remove$name == FALSE)

hfoot_5_ascm <- single_augsynth(form = as.formula(unadj_form_foot), 
                                unit = name,
                                time = biweekly_marker,
                                data = harlem_5_foot,
                                scm = T, t_int = 0, 
                                progfunc = "ridge")

rel_sites <- highres5min %>%
  filter(name %in% harlem_5_foot$name)

hfoot_perm_5_ascm <- get_permutation_results(hfoot_5_ascm, unadj_form_foot, rel_sites, 0,
                                     "ridge", TRUE)

##### Adjusted Ridge aSCM

hfoot_5_ascm_adj <- single_augsynth(form = as.formula(adj_form_foot), 
                                    unit = name,
                                    time = biweekly_marker,
                                    data = harlem_5_foot,
                                    scm = T, t_int = 0, 
                                    progfunc = "ridge",
                                    residualize = T)

hfoot_perm_5_ascm_adj <- get_permutation_results(hfoot_5_ascm, adj_form_foot, rel_sites, 0,
                                     "ridge", TRUE)

hfoot_perm_5_ascm_adj <- read_rds("data/hfoot_perm_5_ascm_adj.rds")


#### Washington Heights

##### Ridge aSCM
  
wh_5_foot <- site_data_5$wh_foot %>%
  filter(name %in% sites_to_remove$name == FALSE)

whfoot_5_ascm <- single_augsynth(form = as.formula(unadj_form_foot),
                                 unit = name,
                                 time = biweekly_marker,
                                 data = wh_5_foot,
                                 scm = T, t_int = 0,
                                 progfunc = "ridge")

wh_perm_5_ascm <- get_permutation_results(whfoot_5_ascm, unadj_form_foot, rel_sites, 0,
                                         "ridge", TRUE)

##### Adjusted Ridge aSCM

whfoot_5_ascm_adj <- single_augsynth(form = as.formula(adj_form_foot),
                                     unit = name,
                                     time = biweekly_marker,
                                     data = wh_5_foot,
                                     scm = T, t_int = 0,
                                     progfunc = "ridge",
                                     residualize = T)

wh_foot_perm_ascm_adj <- get_permutation_results(whfoot_5_ascm_adj, adj_form_foot, rel_sites, 0,
                                        "ridge", TRUE)

