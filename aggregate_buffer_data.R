# Create Spending Data 5-minute buffers

library(tidyverse)
library(arrow)
library(sf)
library(lubridate)
library(tictoc)

spend_data <- open_dataset("data/spend")

site_data <- read_rds("data/site_cd_lookup.rds")

nys_poi <- read_csv("raw_data/dewey_data/poi/ny_places/all_nyc_places.csv")

nys_poi_spatial <- st_as_sf(nys_poi, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
  st_transform(crs = "EPSG:32118")

nys_shared_poi <- read_rds("data/nys_shared_poi.rds")

excluded_places <- read_csv("data/excluded_placekeys.csv")

buffer_5min_s <- read_rds("data/highres5min.rds") %>%
  st_transform(crs = "EPSG:32118")

create_spend_data <- function(df_input, df_index){
  df_rel <- df_input %>%
    slice(df_index)
  
  within_site <- st_join(nys_poi_spatial, df_rel) %>%
    drop_na(name)
  
  site_lu <- within_site %>%
    st_drop_geometry() %>%
    select(PLACEKEY, name)
  
  nys_spend_site <- spend_data %>%
    filter(PLACEKEY %in% within_site$PLACEKEY) %>%
    left_join(nys_poi, by = "PLACEKEY") %>%
    left_join(site_lu, by = "PLACEKEY") %>%
    select(PLACEKEY, name, LOCATION_NAME, RAW_TOTAL_SPEND, RAW_NUM_TRANSACTIONS,
           NAICS_CODE, TOP_CATEGORY, SUB_CATEGORY,
           RAW_NUM_CUSTOMERS, SPEND_BY_DAY, ONLINE_SPEND,
           SPEND_DATE_RANGE_START, SPEND_DATE_RANGE_END) %>%
    collect()
  
  full_data_site <- nys_spend_site %>%
    mutate(SPEND_DATE_RANGE_START = ymd(SPEND_DATE_RANGE_START),
           SPEND_DATE_RANGE_END = ymd(SPEND_DATE_RANGE_END)) %>%
    group_by(PLACEKEY, name) %>%
    filter(min(SPEND_DATE_RANGE_START) <= ymd("2021-06-01")) %>%
    filter(max(SPEND_DATE_RANGE_END) >= ymd("2022-07-01")) %>%
    filter(all(!is.na(RAW_TOTAL_SPEND))) %>%
    ungroup() %>%
    filter(SPEND_DATE_RANGE_END <= ymd("2022-07-01")) %>%
    filter(SPEND_DATE_RANGE_START >= ymd("2021-06-01")) %>%
    filter(PLACEKEY %in% excluded_places$placekey == FALSE) %>%
    filter(PLACEKEY %in% nys_shared_poi$PLACEKEY == FALSE) %>%
    drop_na(RAW_TOTAL_SPEND) %>%
    group_by(PLACEKEY, name) %>%
    mutate(num_obs = n()) %>%
    filter(num_obs == 13)
  
  return(full_data_site)
  
}

custom_toc <- function(tic, toc, msg, info){
  str_c("Time elapsed: ", round((toc - tic)/60,2), " minutes!")
}


all_sites_spend_data_5 <- map_dfr(seq_along(buffer_5min_s$name),
                           function(a){
                             tic()
                             lo <- create_spend_data(buffer_5min_s, a)
                             print(str_c("Site ", a, "/", nrow(buffer_5min_s), "!!!"))
                             toc(func.toc = custom_toc)
                             return(lo)
                           })

write_rds(all_sites_spend_data_5, "data/all_sites_spend_5min.rds")

# Create 10 minute buffer spend data
buffer_10min <- read_rds("data/highres10min.rds") %>%
  st_transform(crs = "EPSG:32118")

all_sites_spend_data_10 <- map_dfr(seq_along(buffer_10min$name),
                                function(a){
                                  tic()
                                  lo <- create_spend_data(buffer_10min, a)
                                  print(str_c("Site ", a, "/", nrow(buffer_10min), "!!!"))
                                  toc(func.toc = custom_toc)
                                  return(lo)
                                })

write_rds(all_sites_spend_data, "data/all_sites_spend_10min.rds")

# Create 5-minute foot data

foot_data <- open_dataset("data/foot")

create_foot_data <- function(df_input, df_index){
  df_rel <- df_input %>%
    slice(df_index)
  
  within_site <- st_join(nys_poi_spatial, df_rel, st_within) %>%
    drop_na(name)
  
  site_lu <- within_site %>%
    st_drop_geometry() %>%
    select(PLACEKEY, name)
  
  nys_foot_site <- foot_data %>%
    filter(PLACEKEY %in% within_site$PLACEKEY) %>%
    select(PLACEKEY, LOCATION_NAME, RAW_VISIT_COUNTS,  NAICS_CODE, TOP_CATEGORY, SUB_CATEGORY,
           NORMALIZED_VISITS_BY_STATE_SCALING, VISITS_BY_DAY,
           DATE_RANGE_START, DATE_RANGE_END) %>%
    filter(PLACEKEY %in% within_site$PLACEKEY) %>%
    left_join(site_lu, by = "PLACEKEY") %>%
    collect()
  
  full_data_site <- nys_foot_site %>%
    mutate(DATE_RANGE_START = ymd(DATE_RANGE_START),
           DATE_RANGE_END = ymd(DATE_RANGE_END)) %>%
    group_by(PLACEKEY, LOCATION_NAME, name) %>%
    filter(min(DATE_RANGE_START) <= ymd("2021-06-01")) %>%
    filter(max(DATE_RANGE_END) >= ymd("2022-07-01")) %>%
    filter(all(!is.na(RAW_VISIT_COUNTS))) %>%
    ungroup() %>%
    filter(DATE_RANGE_END <= ymd("2022-07-01")) %>%
    filter(DATE_RANGE_START >= ymd("2021-06-01")) %>%
    filter(PLACEKEY %in% excluded_places$placekey == FALSE) %>%
    filter(PLACEKEY %in% nys_shared_poi$PLACEKEY == FALSE) %>%
    drop_na(RAW_VISIT_COUNTS) %>%
    group_by(PLACEKEY, LOCATION_NAME, name) %>%
    mutate(num_obs = n()) %>%
    filter(num_obs == 13)
  
  return(full_data_site)
  
}

custom_toc <- function(tic, toc, msg, info){
  str_c("Time elapsed: ", round((toc - tic)/60,2), " minutes!")
}

all_sites_foot_data_5 <- map_dfr(seq_len(nrow(buffer_5min_s)),
                               function(a){
                                 tic()
                                 lo <- create_foot_data(buffer_5min_s, a)
                                 print(str_c("Site ", a, "/", nrow(buffer_5min_s), "!!!"))
                                 toc(func.toc = custom_toc)
                                 return(lo)
                               })

write_rds(all_sites_foot_data_5, "data/all_sites_foot_5min.rds")

# Create 10-minute foot data

buffer_10min <- read_rds("data/highres10min.rds") %>%
  st_transform(crs = "EPSG:32118")

all_sites_foot_data_10 <- map_dfr(seq_along(buffer_10min$name),
                               function(a){
                                 tic()
                                 lo <- create_foot_data(buffer_10min, a)
                                 print(str_c("Site ", a, "/", nrow(buffer_10min), "!!!"))
                                 toc(func.toc = custom_toc)
                                 return(lo)
                               })

write_rds(all_sites_foot_data_10, "data/all_sites_foot_10min.rds")

