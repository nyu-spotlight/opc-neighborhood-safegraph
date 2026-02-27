###### Code to summarize outcomes and add covariates

library(tidyverse)
library(sf)
library(censusapi)
library(tidycensus)
library(tictoc)
library(timetk)
library(quantmod)

spending_data <- read_rds("data/all_sites_spend_10min.rds")

site_data <- read_rds("data/site_cd_lookup.rds")

cd_lookup <- read_rds("data/site_cd_lookup.rds") %>%
  select(name, boro_cd) %>%
  st_drop_geometry()

site_lookup <- site_data %>%
  st_drop_geometry() %>%
  select(name, site_type)

#### Spending Data Setup

spend_daily <- spending_data %>%
  mutate(pre_process = str_remove_all(SPEND_BY_DAY, "\\[") %>%
           str_remove_all("\\]")) %>%
  rowwise() %>%
  mutate(date_vec = str_c(seq(SPEND_DATE_RANGE_START, SPEND_DATE_RANGE_END-1, by = "day"), collapse = ",")) %>%
  separate_longer_delim(c(pre_process, date_vec), delim = ",")

### Account for inflation:

getSymbols("CPIAUCSL", src = "FRED")

cpi_tbl <- tk_tbl(CPIAUCSL) %>%
  mutate(year_cl = year(index),
         month_cl = month(index) %>% as.numeric())  %>%
  mutate(my_cl = str_c(month_cl, "_", year_cl)) %>%
  filter(year_cl %in% 2021:2022)

calculate_inflation <- function(my_input){
  june_2022 <- cpi_tbl %>%
    filter(my_cl == "6_2022") %>%
    pull(CPIAUCSL)

  input_inf <- cpi_tbl %>%
    filter(my_cl == my_input) %>%
    pull(CPIAUCSL)
  
  return(june_2022/input_inf)
}

inflation_vals <- tibble(my_cl = c(str_c(1:12, "_2021"), str_c(1:12, "_2022"))) %>%
  rowwise() %>%
  mutate(inf_ratio = calculate_inflation(my_cl))

spend_daily_inf <- spend_daily %>%
  mutate(year_cl = year(date_vec),
         month_cl = month(date_vec) %>% as.numeric())  %>%
  mutate(my_cl = str_c(month_cl, "_", year_cl)) %>%
  left_join(inflation_vals, by = "my_cl") %>%
  mutate(pre_process = as.numeric(pre_process)) %>%
  mutate(pre_process_inf = inf_ratio*pre_process)

spend_biweekly <- spend_daily_inf %>%
  mutate(date_vec = ymd(date_vec)) %>%
  mutate(weeks_since = days(date_vec - ymd("2021-11-30"))/ days(14)) %>%
  select(PLACEKEY, LOCATION_NAME, name, pre_process_inf, RAW_TOTAL_SPEND, TOP_CATEGORY, NAICS_CODE, ONLINE_SPEND, SPEND_DATE_RANGE_START, num_obs, date_vec, weeks_since) %>%
  mutate(biweekly_marker = case_when(weeks_since < 0 ~ floor(weeks_since),
                                     weeks_since >= 0 ~ ceiling(weeks_since + 1/14) - 1)) %>%
  group_by(PLACEKEY, LOCATION_NAME, name, biweekly_marker) %>%
  mutate(marker_obs = n()) %>%
  filter(biweekly_marker %in% -13:13) %>%
  filter(num_obs == 13)

biweekly_summary <- spend_biweekly %>%
  mutate(prop_spending = (RAW_TOTAL_SPEND - ONLINE_SPEND)/RAW_TOTAL_SPEND) %>%
  mutate(day_spend = as.numeric(pre_process_inf)) %>%
  mutate(true_spending = prop_spending*day_spend) %>%
  group_by(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY, NAICS_CODE, biweekly_marker) %>%
  summarize(total_spending = sum(true_spending))

median_name_cd <- biweekly_summary %>%
  group_by(name, biweekly_marker) %>%
  summarize(mean_total_spend = mean(total_spending),
            median_total_spend = median(total_spending)) %>%
  left_join(site_lookup, by = "name")

### Foot data

foot_data <- read_rds("data/all_sites_foot_10min.rds")
  
foot_daily <- foot_data %>%
  filter(DATE_RANGE_START != ymd("2021-11-01")) %>%
  mutate(pre_process = str_remove_all(VISITS_BY_DAY, "\\[") %>%
           str_remove_all("\\]")) %>%
  rowwise() %>%
  mutate(date_vec = str_c(seq(DATE_RANGE_START, DATE_RANGE_END-1, by = "day"), collapse = ",")) %>%
  separate_longer_delim(c(pre_process, date_vec), delim = ",")

nov_problems <- foot_data %>%
  filter(DATE_RANGE_START == ymd("2021-11-01")) %>%
  mutate(pre_process = str_remove_all(VISITS_BY_DAY, "\\[") %>%
           str_remove_all("\\]")) %>%
  rowwise() %>%
  mutate(date_vec = str_c(seq(DATE_RANGE_START, DATE_RANGE_END-2, by = "day"), collapse = ",")) %>%
  separate_longer_delim(c(pre_process, date_vec), delim = ",") %>%
  mutate(date_vec = ymd(date_vec)) %>%
  mutate(pre_process = as.numeric(pre_process)) %>%
  mutate(date_vec = case_when(date_vec >= ymd("2021-11-07") ~ date_vec + 1,
                              TRUE ~ date_vec))

real_dst_vals <- nov_problems %>%
  group_by(PLACEKEY, LOCATION_NAME, name, DATE_RANGE_START, RAW_VISIT_COUNTS) %>%
  summarize(total_foot = sum(pre_process)) %>%
  mutate(pre_process = RAW_VISIT_COUNTS - total_foot) %>%
  mutate(date_vec = ymd("2021-11-07")) %>%
  select(PLACEKEY, LOCATION_NAME, name, DATE_RANGE_START, RAW_VISIT_COUNTS, date_vec, pre_process)

nov_fixed <- nov_problems %>%
  bind_rows(real_dst_vals) 

foot_daily_complete <- foot_daily %>%
  mutate(pre_process = as.numeric(pre_process)) %>%
  mutate(date_vec = ymd(date_vec)) %>%
  bind_rows(nov_fixed)

# Set up for biweekly analysis

foot_biweekly <- foot_daily_complete %>%
  mutate(weeks_since = days(date_vec - ymd("2021-11-30"))/ days(14)) %>%
  select(PLACEKEY, LOCATION_NAME, name, pre_process, TOP_CATEGORY, NAICS_CODE, DATE_RANGE_START, num_obs, date_vec, weeks_since) %>%
  mutate(biweekly_marker = case_when(weeks_since < 0 ~ floor(weeks_since),
                                     weeks_since >= 0 ~ ceiling(weeks_since + 1/14) - 1)) %>%
  group_by(PLACEKEY, LOCATION_NAME, name, biweekly_marker) %>%
  mutate(marker_obs = n()) %>%
  filter(biweekly_marker %in% -13:13) %>%
  filter(num_obs == 13)

biweekly_summary_foot <- foot_biweekly %>%
  group_by(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY, NAICS_CODE, biweekly_marker) %>%
  summarize(total_foot = sum(pre_process))

median_bid_foot <- biweekly_summary_foot %>%
  group_by(name, biweekly_marker) %>%
  summarize(mean_total_foot = mean(total_foot),
            median_total_foot = median(total_foot)) %>%
  left_join(site_lookup, by = "name")


# Linking covariates

biweekly_lookup <- spend_biweekly %>%
  ungroup() %>%
  select(biweekly_marker, SPEND_DATE_RANGE_START) %>%
  distinct()

get_census_vars <- function(year){
  name_code = str_c(year, "/acs/acs5")
  target_vars <- c("STATE",
                   "NAME",
                   "B02009_001E",
                   "B01003_001E", 
                   "B03002_012E",
                   "B01001_011E", 
                   "B01001_012E",
                   "B01001_013E", 
                   "B01001_014E", 
                   "B01001_035E", 
                   "B01001_036E", 
                   "B01001_037E", 
                   "B01001_038E", 
                   "B17010_001E", 
                   "B17010_002E", 
                   "B23001_006E","B23001_013E","B23001_020E","B23001_027E","B23001_034E","B23001_041E","B23001_048E","B23001_055E", "B23001_062E", "B23001_069E", 
                   "B23001_092E","B23001_099E","B23001_106E","B23001_113E","B23001_120E","B23001_127E","B23001_134E","B23001_141E", "B23001_148E", "B23001_155E", 
                   "B23001_008E","B23001_015E","B23001_022E","B23001_029E","B23001_036E","B23001_043E","B23001_050E","B23001_057E", "B23001_064E", "B23001_071E", 
                   "B23001_094E","B23001_101E","B23001_108E","B23001_115E","B23001_122E","B23001_129E","B23001_136E","B23001_143E", "B23001_150E", "B23001_157E",
                   "B19013_001E", 
                   "B21004_001E",
                   "B01002_001E",
                   "B01001_001E",
                   "B01001_026E",  
                   "B19058_001E",
                   "B19058_002E" 
  ) 
  if(year > 2008){
    target_vars <- c(target_vars, 
                     "B15002_001E",
                     "B15002_003E", "B15002_004E", "B15002_005E", "B15002_006E", "B15002_007E", "B15002_008E", "B15002_009E", "B15002_010E",
                     "B15002_020E", "B15002_021E", "B15002_022E", "B15002_023E", "B15002_024E", "B15002_025E", "B15002_026E", "B15002_027E"
    )
  }
  
  census_dat <- getCensus(
    name= name_code,
    region = "tract:*",
    regionin="state:36",
    vars = target_vars
  )
  
  output_df <- census_dat %>%
    group_by(STATE, NAME, county, tract) %>%
    mutate(YEAR = year) %>%
    mutate(prop_black = B02009_001E / B01003_001E,
           prop_hispanic = B03002_012E / B01003_001E,
           prop_25_44 = (B01001_011E + B01001_012E + B01001_013E + B01001_014E + B01001_035E + B01001_036E + B01001_037E + B01001_038E)/B01003_001E,
           prop_fam_pov = B17010_002E/B17010_001E,
           prop_unemploy =  sum(c_across(B23001_008E:B23001_157E))/sum(c_across(B23001_013E:B23001_155E)),
           med_hh_income = B19013_001E,
           med_income = B21004_001E,
           prop_hs_diploma = 1 - sum(c_across(B15002_003E:B15002_027E))/B15002_001E,
           med_age = B01002_001E,
           prop_female = B01001_026E/B01001_001E,
           prop_pa_fs = B19058_002E/B19058_001E) %>%
    select(STATE, NAME, county, tract, prop_black, prop_hispanic, med_age, prop_fam_pov, prop_unemploy, med_income, prop_female, YEAR,
           prop_pa_fs, prop_hs_diploma) %>%
    mutate(GEOID = str_c(STATE, county, tract)) %>%
    select(STATE, YEAR, everything())
  
  return(output_df)
  
}

Sys.setenv(CENSUS_KEY="XXX")

census_2021 <- get_census_vars(2021) %>%
  mutate(med_age = if_else(med_age <= 0, NA_real_, med_age)) %>%
  mutate(med_income = if_else(med_income <= 0, NA_real_, med_income))

highres10min <- read_rds("data/highres10min.rds")

shape21 <- st_read("shapefiles/tl_2021_36_tract/tl_2021_36_tract.shp", quiet = T)

shape21_t <- st_transform(shape21, crs = st_crs(highres10min)) %>%
  left_join(census_2021, by = "GEOID")

shape21_inter <- st_join(shape21_t, highres10min)

only_inter21 <- shape21_inter %>%
  drop_na(name)

covars <- c("prop_black",
            "prop_hispanic",
            "med_age",
            "prop_fam_pov",
            "prop_unemploy",
            "med_income",
            "prop_female",
            "prop_pa_fs",
            "prop_hs_diploma")

custom_toc <- function(tic, toc, msg, info){
  str_c("Time elapsed: ", round((toc - tic)/60,2), " minutes!")
}

interp_df_2021 <- map_dfc(covars, function(x){
  print(str_c("Starting areal interpolation of: ", x, "!"))
  tic()
  df_out <- st_interpolate_aw(shape21_t %>% select(all_of(x)), highres10min, extensive = FALSE,
                    keep_NA = TRUE, na.rm = TRUE) %>%
    st_drop_geometry() %>%
    select(x) 
  print(str_c("Finished with: ", x, "!"))
  toc(func.toc = custom_toc)
  return(df_out)
}) %>%
  bind_cols(highres10min %>% st_drop_geometry() %>% select(name))


tog21 <- only_inter21 %>%
  filter(med_age > 0) %>%
  mutate(med_income = if_else(med_income <= 0, NA_real_, med_income)) %>%
  st_drop_geometry() %>%
  group_by(name, YEAR) %>%
  summarize(across(any_of(covars), function(a) mean(a, na.rm = T), .names = "{.col}_mean"))

# POI type variables

prop_types <- spend_biweekly %>%
  ungroup() %>%
  select(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY) %>%
  distinct() %>%
  group_by(name, TOP_CATEGORY) %>%
  summarize(n_cat = n()) %>%
  group_by(name) %>%
  mutate(proportion = n_cat/sum(n_cat))

prop_restaurants <- prop_types %>%
  filter(TOP_CATEGORY == "Restaurants and Other Eating Places") %>%
  select(name, prop_rest = proportion)

prop_grocery <- prop_types %>%
  filter(TOP_CATEGORY == "Grocery Stores") %>%
  select(name, prop_groc = proportion)

prop_types_f <- foot_biweekly %>%
  ungroup() %>%
  select(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY) %>%
  distinct() %>%
  group_by(name, TOP_CATEGORY) %>%
  summarize(n_cat = n()) %>%
  group_by(name) %>%
  mutate(proportion = n_cat/sum(n_cat))

prop_restaurants_f <- prop_types_f %>%
  filter(TOP_CATEGORY == "Restaurants and Other Eating Places") %>%
  select(name, prop_rest = proportion)

prop_grocery_f <- prop_types_f %>%
  filter(TOP_CATEGORY == "Grocery Stores") %>%
  select(name, prop_groc = proportion)

# Gentrification variable setup

cbsa_2011 <- get_acs(geography = "cbsa",
                     variables = "B19326_001",
                     year = 2011) %>%
  filter(GEOID == 35620) %>%
  pull(estimate)

inc_2011 <- get_acs(geography = "tract",
                    state = "36",
                    year = 2011,
                    variables = "B19326_001")

tract_11 <- st_read("shapefiles/tl_2011_36_tract/tl_2011_36_tract.shp")

highres10min <- read_rds("data/highres10min.rds")

tract_11_t <- st_transform(tract_11, crs = st_crs(highres10min))

buffer_tract_inter_11 <- st_join(highres10min, tract_11_t)

rel_inc_2011 <- inc_2011 %>%
  filter(GEOID %in% buffer_tract_inter_11$GEOID) %>%
  drop_na() %>%
  select(GEOID, income_2011 = estimate)

cbsa_2021 <- get_acs(geography = "cbsa",
                     variables = "B19326_001",
                     year = 2021) %>%
  filter(GEOID == 35620) %>%
  pull(estimate)

inc_2021<- get_acs(geography = "tract",
                   state = "36",
                   year = 2021, 
                   variables = "B19326_001")

rel_inc_2021 <- inc_2021 %>%
  filter(GEOID %in% buffer_tract_inter_11$GEOID) %>%
  drop_na() %>%
  select(GEOID, income_2021 = estimate)

gent_var <- buffer_tract_inter_11 %>%
  left_join(rel_inc_2011, by = "GEOID") %>%
  left_join(rel_inc_2021, by = "GEOID") %>%
  mutate(ratio_2011 = income_2011/cbsa_2011,
         ratio_2021 = income_2021/cbsa_2021) %>%
  drop_na(ratio_2011, ratio_2021) %>%
  st_drop_geometry() %>%
  select(GEOID, name, income_2011, income_2021, ratio_2011, ratio_2021) %>%
  mutate(difference = ratio_2021 - ratio_2011) %>%
  mutate(gentrifying = as.numeric(difference >= 0.1))

gent_shapefile <- shape21_t %>%
  left_join(gent_var, by = "GEOID")

gent_interp <- st_interpolate_aw(gent_shapefile %>% select(gentrifying), highres10min, extensive = FALSE,
                                 keep_NA = TRUE, na.rm = TRUE) %>%
  st_drop_geometry() %>%
  select(gentrifying) %>%
  bind_cols(highres10min %>% st_drop_geometry() %>% select(name))

gent_summary <- gent_var %>%
  group_by(name) %>%
  summarize(prop_gent = mean(gentrifying))

# COVID-19 Hospitalization Data

covid_hosprate <- read_csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/refs/heads/master/trends/hosprate-by-modzcta.csv",
                           col_types = "c") %>%
  filter(str_detect(date, "2021"))

covid_long <- covid_hosprate %>%
  dplyr::select(date, HOSPRATE_10001:HOSPRATE_11697) %>%
  pivot_longer(cols = HOSPRATE_10001:HOSPRATE_11697, names_to = "zip", names_transform = function(a) str_remove(a, "HOSPRATE_"),
               values_to = "hosprate") %>%
  mutate(date_link = my(date)) 

covid_zips <- st_read("shapefiles/COVID ZCTA/MODZCTA_2010.shp")

covid_zips_t <- st_transform(covid_zips, crs = st_crs(highres10min)) %>%
  st_make_valid()

walking_zips <- st_join(highres10min, covid_zips_t) %>%
  drop_na(MODZCTA)

walking_hosprate <- walking_zips %>%
  st_drop_geometry() %>%
  left_join(covid_long, by = c("MODZCTA" = "zip")) %>%
  group_by(name, date, date_link) %>%
  summarize(mean_hosprate = mean(hosprate, na.rm = T)) %>%
  ungroup() %>%
  select(name, date_link, mean_hosprate)

hosprate_interp <- walking_zips %>%
  left_join(covid_long, by = c("MODZCTA" = "zip")) %>%
  group_by(date_link) %>%
  nest()

covid_to_interp <- covid_zips_t %>%
  left_join(covid_long, by = c("MODZCTA" = "zip")) %>%
  group_by(date_link) %>%
  nest() %>%
  drop_na(date_link) %>%
  mutate(site_list = map(data, function(x) st_interpolate_aw(x %>% select(hosprate), highres10min, extensive = FALSE,
                                                             keep_NA = TRUE, na.rm = TRUE) %>%
                           bind_cols(highres10min %>% st_drop_geometry() %>% select(name))))

covid_interp <- covid_to_interp %>%
  select(-data) %>%
  unnest(cols = site_list)
  
hosprate_walk10_biweekly <- covid_interp %>%
  left_join(biweekly_lookup, by = c("date_link" = "SPEND_DATE_RANGE_START")) %>%
  drop_na(biweekly_marker) %>%
  group_by(name, biweekly_marker) %>%
  summarize(overall_mean_hosprate = mean(hosprate, na.rm = T))

# NYPD Arrest Data

raw_arrest_data <- read_csv("raw_data/arrests/NYPD_Arrests_Data__Historic__20241112.csv")

cleaned_arrest_data <- raw_arrest_data %>%
  mutate(clean_date = mdy(ARREST_DATE)) %>%
  filter(clean_date >= ymd("2021-06-01")) %>%
  filter(clean_date <= ymd("2022-07-01"))

arrest_spatial <- st_as_sf(cleaned_arrest_data, coords = c("Longitude", "Latitude"),
                           crs = 4326)

arrest_t <- st_transform(arrest_spatial, crs = st_crs(highres10min))

arrest_j <- st_join(arrest_spatial, highres10min) %>%
  drop_na(name)

biweekly_arrest_data <- arrest_j %>%
  mutate(weeks_since = days(clean_date - ymd("2021-11-30"))/ days(14)) %>%
  mutate(biweekly_marker = case_when(weeks_since < 0 ~ floor(weeks_since),
                                     weeks_since >= 0 ~ ceiling(weeks_since + 1/14) - 1)) %>%
  st_drop_geometry() %>%
  group_by(name, biweekly_marker) %>%
  summarize(n_arrests = n()) %>%
  select(name, biweekly_marker, n_arrests)

mean_arrests <- biweekly_arrest_data %>%
  filter(biweekly_marker < 0) %>%
  group_by(name) %>%
  summarize(mean_arrests = mean(n_arrests, na.rm = T)) 

biweekly_arrest_data <- biweekly_arrest_data %>%
  left_join(mean_arrests, by = "name") %>%
  mutate(n_arrests = case_when(is.na(n_arrests) ~ mean_arrests,
                                  TRUE ~ n_arrests))

# Combine all data sources

# East Harlem

harlem_cd <- cd_lookup %>%
  filter(str_detect(name, "East Harlem"))

wh_cd <- cd_lookup %>%
  filter(str_detect(name, "OnPoint Washington Heights"))

harlem_comps <- cd_lookup %>%
  group_by(name) %>%
  filter(all(boro_cd %in% harlem_cd$boro_cd == FALSE)) %>%
  filter(all(boro_cd %in% wh_cd$boro_cd == FALSE))

harlem_spend <- median_name_cd %>%
  left_join(tog21, by = "name") %>%
  left_join(prop_restaurants, by = "name") %>%
  left_join(prop_grocery, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_summary, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint Washington Heights") %>%
  filter(name == "OnPoint East Harlem" | name %in% harlem_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint East Harlem" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 


# Washington Heights

wh_comps <- cd_lookup %>%
  group_by(name) %>%
  filter(all(boro_cd %in% wh_cd$boro_cd == FALSE)) %>%
  filter(all(boro_cd %in% harlem_cd$boro_cd == FALSE))

wh_spend <- median_name_cd %>%
  left_join(tog21, by = "name") %>%
  left_join(prop_restaurants, by = "name") %>%
  left_join(prop_grocery, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_summary, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint East Harlem") %>%
  filter(name == "OnPoint Washington Heights" | name %in% wh_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint Washington Heights" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 


# Foot Traffic Setup

harlem_foot <- median_bid_foot %>%
  left_join(tog21, by = "name") %>%
  left_join(prop_restaurants, by = "name") %>%
  left_join(prop_grocery, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_summary, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  filter(name != "Greenwich House East OTP",
         name != "Lafayette Medical Approach, LLC OTP",
         name != "Long Island Jewish Medical Center OTP") %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint Washington Heights") %>%
  filter(name == "OnPoint East Harlem" | name %in% harlem_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint East Harlem" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 

# Washington Heights

wh_foot <- median_bid_foot %>%
  left_join(tog21, by = "name") %>%
  left_join(prop_restaurants, by = "name") %>%
  left_join(prop_grocery, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_summary, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  filter(name != "Greenwich House East OTP",
         name != "Lafayette Medical Approach, LLC OTP",
         name != "Long Island Jewish Medical Center OTP") %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint East Harlem") %>%
  filter(name == "OnPoint Washington Heights" | name %in% wh_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint Washington Heights" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 

site_list <- list(h_spend = harlem_spend,
                  h_foot = harlem_foot,
                  wh_spend = wh_spend,
                  wh_foot = wh_foot)

write_rds(site_list, "data/site_analysis_dat_10min.rds")

# Same process for 5-minute buffers:

spending_data <- read_rds("data/all_sites_spend_5min.rds")

cd_lookup <- read_rds("data/site_cd_lookup.rds") %>%
  select(name, boro_cd) %>%
  st_drop_geometry()

spend_daily <- spending_data %>%
  mutate(pre_process = str_remove_all(SPEND_BY_DAY, "\\[") %>%
           str_remove_all("\\]")) %>%
  rowwise() %>%
  mutate(date_vec = str_c(seq(SPEND_DATE_RANGE_START, SPEND_DATE_RANGE_END-1, by = "day"), collapse = ",")) %>%
  separate_longer_delim(c(pre_process, date_vec), delim = ",")

### Account for inflation:

inflation_vals <- tibble(my_cl = c(str_c(1:12, "_2021"), str_c(1:12, "_2022"))) %>%
  rowwise() %>%
  mutate(inf_ratio = calculate_inflation(my_cl))

spend_daily_inf <- spend_daily %>%
  mutate(year_cl = year(date_vec),
         month_cl = month(date_vec) %>% as.numeric())  %>%
  mutate(my_cl = str_c(month_cl, "_", year_cl)) %>%
  left_join(inflation_vals, by = "my_cl") %>%
  mutate(pre_process = as.numeric(pre_process)) %>%
  mutate(pre_process_inf = inf_ratio*pre_process)

spend_biweekly <- spend_daily_inf %>%
  mutate(date_vec = ymd(date_vec)) %>%
  mutate(weeks_since = days(date_vec - ymd("2021-11-30"))/ days(14)) %>%
  select(PLACEKEY, LOCATION_NAME, name, pre_process_inf, RAW_TOTAL_SPEND, TOP_CATEGORY, NAICS_CODE, ONLINE_SPEND, SPEND_DATE_RANGE_START, num_obs, date_vec, weeks_since) %>%
  mutate(biweekly_marker = case_when(weeks_since < 0 ~ floor(weeks_since),
                                     weeks_since >= 0 ~ ceiling(weeks_since + 1/14) - 1)) %>%
  group_by(PLACEKEY, LOCATION_NAME, name, biweekly_marker) %>%
  mutate(marker_obs = n()) %>%
  filter(biweekly_marker %in% -13:13) %>%
  filter(num_obs == 13)

biweekly_summary <- spend_biweekly %>%
  mutate(prop_spending = (RAW_TOTAL_SPEND - ONLINE_SPEND)/RAW_TOTAL_SPEND) %>%
  mutate(day_spend = as.numeric(pre_process_inf)) %>%
  mutate(true_spending = prop_spending*day_spend) %>%
  group_by(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY, NAICS_CODE, biweekly_marker) %>%
  summarize(total_spending = sum(true_spending))

median_name_cd <- biweekly_summary %>%
  group_by(name, biweekly_marker) %>%
  summarize(mean_total_spend = mean(total_spending),
            median_total_spend = median(total_spending)) %>%
  left_join(site_lookup, by = "name")

# Foot Traffic Data

foot_data <- read_rds("data/all_sites_foot_5min.rds")

foot_daily <- foot_data %>%
  filter(DATE_RANGE_START != ymd("2021-11-01")) %>%
  mutate(pre_process = str_remove_all(VISITS_BY_DAY, "\\[") %>%
           str_remove_all("\\]")) %>%
  rowwise() %>%
  mutate(date_vec = str_c(seq(DATE_RANGE_START, DATE_RANGE_END-1, by = "day"), collapse = ",")) %>%
  separate_longer_delim(c(pre_process, date_vec), delim = ",")

nov_problems <- foot_data %>%
  filter(DATE_RANGE_START == ymd("2021-11-01")) %>%
  mutate(pre_process = str_remove_all(VISITS_BY_DAY, "\\[") %>%
           str_remove_all("\\]")) %>%
  rowwise() %>%
  mutate(date_vec = str_c(seq(DATE_RANGE_START, DATE_RANGE_END-2, by = "day"), collapse = ",")) %>%
  separate_longer_delim(c(pre_process, date_vec), delim = ",") %>%
  mutate(date_vec = ymd(date_vec)) %>%
  mutate(pre_process = as.numeric(pre_process)) %>%
  mutate(date_vec = case_when(date_vec >= ymd("2021-11-07") ~ date_vec + 1,
                              TRUE ~ date_vec))

real_dst_vals <- nov_problems %>%
  group_by(PLACEKEY, LOCATION_NAME, name, DATE_RANGE_START, RAW_VISIT_COUNTS) %>%
  summarize(total_foot = sum(pre_process)) %>%
  mutate(pre_process = RAW_VISIT_COUNTS - total_foot) %>%
  mutate(date_vec = ymd("2021-11-07")) %>%
  select(PLACEKEY, LOCATION_NAME, name, DATE_RANGE_START, RAW_VISIT_COUNTS, date_vec, pre_process)

nov_fixed <- nov_problems %>%
  bind_rows(real_dst_vals) 

foot_daily_complete <- foot_daily %>%
  mutate(pre_process = as.numeric(pre_process)) %>%
  mutate(date_vec = ymd(date_vec)) %>%
  bind_rows(nov_fixed)

# Set up for analysis

foot_biweekly <- foot_daily_complete %>%
  mutate(weeks_since = days(date_vec - ymd("2021-11-30"))/ days(14)) %>%
  select(PLACEKEY, LOCATION_NAME, name, pre_process, TOP_CATEGORY, NAICS_CODE, DATE_RANGE_START, num_obs, date_vec, weeks_since) %>%
  mutate(biweekly_marker = case_when(weeks_since < 0 ~ floor(weeks_since),
                                     weeks_since >= 0 ~ ceiling(weeks_since + 1/14) - 1)) %>%
  group_by(PLACEKEY, LOCATION_NAME, name, biweekly_marker) %>%
  mutate(marker_obs = n()) %>%
  filter(biweekly_marker %in% -13:13) %>%
  filter(num_obs == 13)

biweekly_summary_foot <- foot_biweekly %>%
  group_by(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY, NAICS_CODE, biweekly_marker) %>%
  summarize(total_foot = sum(pre_process))

median_bid_foot <- biweekly_summary_foot %>%
  group_by(name, biweekly_marker) %>%
  summarize(mean_total_foot = mean(total_foot),
            median_total_foot = median(total_foot)) %>%
  left_join(site_lookup, by = "name")

biweekly_lookup <- spend_biweekly %>%
  ungroup() %>%
  select(biweekly_marker, SPEND_DATE_RANGE_START) %>%
  distinct()

highres5min <- read_rds("data/highres5min.rds")

shape21 <- st_read("shapefiles/tl_2021_36_tract/tl_2021_36_tract.shp", quiet = T)

shape21_t <- st_transform(shape21, crs = st_crs(highres5min)) %>%
  left_join(census_2021, by = "GEOID")

interp_df_2021 <- map_dfc(covars, function(x){
  print(str_c("Starting areal interpolation of: ", x, "!"))
  tic()
  df_out <- st_interpolate_aw(shape21_t %>% select(all_of(x)), highres5min, extensive = FALSE,
                              keep_NA = TRUE, na.rm = TRUE) %>%
    st_drop_geometry() %>%
    select(x) 
  print(str_c("Finished with: ", x, "!"))
  toc(func.toc = custom_toc)
  return(df_out)
}) %>%
  bind_cols(highres5min %>% st_drop_geometry() %>% select(name))

prop_types <- spend_biweekly %>%
  ungroup() %>%
  select(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY) %>%
  distinct() %>%
  group_by(name, TOP_CATEGORY) %>%
  summarize(n_cat = n()) %>%
  group_by(name) %>%
  mutate(proportion = n_cat/sum(n_cat))

prop_restaurants <- prop_types %>%
  filter(TOP_CATEGORY == "Restaurants and Other Eating Places") %>%
  select(name, prop_rest = proportion)

prop_grocery <- prop_types %>%
  filter(TOP_CATEGORY == "Grocery Stores") %>%
  select(name, prop_groc = proportion)

prop_types_f <- foot_biweekly %>%
  ungroup() %>%
  select(PLACEKEY, LOCATION_NAME, name, TOP_CATEGORY) %>%
  distinct() %>%
  group_by(name, TOP_CATEGORY) %>%
  summarize(n_cat = n()) %>%
  group_by(name) %>%
  mutate(proportion = n_cat/sum(n_cat))

prop_restaurants_f <- prop_types_f %>%
  filter(TOP_CATEGORY == "Restaurants and Other Eating Places") %>%
  select(name, prop_rest = proportion)

prop_grocery_f <- prop_types_f %>%
  filter(TOP_CATEGORY == "Grocery Stores") %>%
  select(name, prop_groc = proportion)

cbsa_2011 <- get_acs(geography = "cbsa",
                     variables = "B19326_001",
                     year = 2011) %>%
  filter(GEOID == 35620) %>%
  pull(estimate)

inc_2011 <- get_acs(geography = "tract",
                    state = "36",
                    year = 2011,
                    variables = "B19326_001")

tract_11 <- st_read("shapefiles/tl_2011_36_tract/tl_2011_36_tract.shp")

tract_11_t <- st_transform(tract_11, crs = st_crs(highres5min))

buffer_tract_inter_11 <- st_join(highres5min, tract_11_t)

rel_inc_2011 <- inc_2011 %>%
  filter(GEOID %in% buffer_tract_inter_11$GEOID) %>%
  drop_na() %>%
  select(GEOID, income_2011 = estimate)

cbsa_2021 <- get_acs(geography = "cbsa",
                     variables = "B19326_001",
                     year = 2021) %>%
  filter(GEOID == 35620) %>%
  pull(estimate)

inc_2021<- get_acs(geography = "tract",
                   state = "36",
                   year = 2021, 
                   variables = "B19326_001")

rel_inc_2021 <- inc_2021 %>%
  filter(GEOID %in% buffer_tract_inter_11$GEOID) %>%
  drop_na() %>%
  select(GEOID, income_2021 = estimate)

gent_var <- buffer_tract_inter_11 %>%
  left_join(rel_inc_2011, by = "GEOID") %>%
  left_join(rel_inc_2021, by = "GEOID") %>%
  mutate(ratio_2011 = income_2011/cbsa_2011,
         ratio_2021 = income_2021/cbsa_2021) %>%
  drop_na(ratio_2011, ratio_2021) %>%
  st_drop_geometry() %>%
  select(GEOID, name, income_2011, income_2021, ratio_2011, ratio_2021) %>%
  mutate(difference = ratio_2021 - ratio_2011) %>%
  mutate(gentrifying = as.numeric(difference >= 0.1))

gent_shapefile <- shape21_t %>%
  left_join(gent_var, by = "GEOID")

gent_interp <- st_interpolate_aw(gent_shapefile %>% select(gentrifying), highres5min, extensive = FALSE,
                                 keep_NA = TRUE, na.rm = TRUE) %>%
  st_drop_geometry() %>%
  select(gentrifying) %>%
  bind_cols(highres5min %>% st_drop_geometry() %>% select(name))

covid_hosprate <- read_csv("https://raw.githubusercontent.com/nychealth/coronavirus-data/refs/heads/master/trends/hosprate-by-modzcta.csv",
                           col_types = "c") %>%
  filter(str_detect(date, "2021"))

covid_long <- covid_hosprate %>%
  dplyr::select(date, HOSPRATE_10001:HOSPRATE_11697) %>%
  pivot_longer(cols = HOSPRATE_10001:HOSPRATE_11697, names_to = "zip", names_transform = function(a) str_remove(a, "HOSPRATE_"),
               values_to = "hosprate") %>%
  mutate(date_link = my(date)) 

covid_zips <- st_read("shapefiles/COVID ZCTA/MODZCTA_2010.shp")

covid_zips_t <- st_transform(covid_zips, crs = st_crs(highres5min)) %>%
  st_make_valid()

walking_zips <- st_join(highres5min, covid_zips_t) %>%
  drop_na(MODZCTA)

walking_hosprate <- walking_zips %>%
  st_drop_geometry() %>%
  left_join(covid_long, by = c("MODZCTA" = "zip")) %>%
  group_by(name, date, date_link) %>%
  summarize(mean_hosprate = mean(hosprate, na.rm = T)) %>%
  ungroup() %>%
  select(name, date_link, mean_hosprate)

hosprate_interp <- walking_zips %>%
  left_join(covid_long, by = c("MODZCTA" = "zip")) %>%
  group_by(date_link) %>%
  nest()

covid_to_interp <- covid_zips_t %>%
  left_join(covid_long, by = c("MODZCTA" = "zip")) %>%
  group_by(date_link) %>%
  nest() %>%
  drop_na(date_link) %>%
  mutate(site_list = map(data, function(x) st_interpolate_aw(x %>% select(hosprate), highres5min, extensive = FALSE,
                                                             keep_NA = TRUE, na.rm = TRUE) %>%
                           bind_cols(highres5min %>% st_drop_geometry() %>% select(name))))

covid_interp <- covid_to_interp %>%
  select(-data) %>%
  unnest(cols = site_list)

buffer_hosprate_biweekly <- covid_interp %>%
  left_join(biweekly_lookup, by = c("date_link" = "SPEND_DATE_RANGE_START")) %>%
  drop_na(biweekly_marker) %>%
  group_by(name, biweekly_marker) %>%
  summarize(overall_mean_hosprate = mean(hosprate, na.rm = T))

mean_hosprates <- buffer_hosprate_biweekly %>%
  filter(biweekly_marker < 0) %>%
  group_by(name) %>%
  summarize(mean_hosprate_imp = mean(overall_mean_hosprate, na.rm = T))

buffer_hosprate_biweekly <- buffer_hosprate_biweekly %>%
  left_join(mean_hosprates, by = "name") %>%
  mutate(overall_mean_hosprate = case_when(is.na(overall_mean_hosprate) ~ mean_hosprate_imp,
                                           TRUE ~ overall_mean_hosprate))

harlem_cd <- cd_lookup %>%
  filter(str_detect(name, "East Harlem"))

wh_cd <- cd_lookup %>%
  filter(str_detect(name, "OnPoint Washington Heights"))

harlem_comps <- cd_lookup %>%
  group_by(name) %>%
  filter(all(boro_cd %in% harlem_cd$boro_cd == FALSE)) %>%
  filter(all(boro_cd %in% wh_cd$boro_cd == FALSE))

harlem_spend <- median_name_cd %>%
  left_join(interp_df_2021, by = "name") %>%
  left_join(prop_restaurants, by = "name") %>%
  left_join(prop_grocery, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_interp, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint Washington Heights") %>%
  filter(name == "OnPoint East Harlem" | name %in% harlem_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint East Harlem" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 

wh_comps <- cd_lookup %>%
  group_by(name) %>%
  filter(all(boro_cd %in% wh_cd$boro_cd == FALSE)) %>%
  filter(all(boro_cd %in% harlem_cd$boro_cd == FALSE))

wh_spend <- median_name_cd %>%
  left_join(interp_df_2021, by = "name") %>%
  left_join(prop_restaurants, by = "name") %>%
  left_join(prop_grocery, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_interp, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint East Harlem") %>%
  filter(name == "OnPoint Washington Heights" | name %in% wh_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint Washington Heights" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 

harlem_foot <- median_bid_foot %>%
  left_join(interp_df_2021, by = "name") %>%
  left_join(prop_restaurants_f, by = "name") %>%
  left_join(prop_grocery_f, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_interp, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint Washington Heights") %>%
  filter(name == "OnPoint East Harlem" | name %in% harlem_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint East Harlem" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 

wh_foot <- median_bid_foot %>%
  left_join(interp_df_2021, by = "name") %>%
  left_join(prop_restaurants_f, by = "name") %>%
  left_join(prop_grocery_f, by = "name") %>%
  mutate(prop_rest = if_else(is.na(prop_rest), 0, prop_rest),
         prop_groc = if_else(is.na(prop_groc), 0, prop_groc)) %>%
  left_join(gent_interp, by = "name") %>%
  left_join(buffer_hosprate_biweekly, by = c("name", "biweekly_marker")) %>%
  left_join(biweekly_arrest_data, by = c("name", "biweekly_marker")) %>%
  mutate(n_arrests = if_else(is.na(n_arrests), 0, n_arrests)) %>%
  filter(name != "OnPoint East Harlem") %>%
  filter(name == "OnPoint Washington Heights" | name %in% wh_comps$name) %>%
  mutate(treated = case_when(name == "OnPoint Washington Heights" & biweekly_marker >= 0 ~ 1,
                             TRUE ~ 0)) %>%
  ungroup() 

site_list <- list(h_spend = harlem_spend,
                  h_foot = harlem_foot,
                  wh_spend = wh_spend,
                  wh_foot = wh_foot)

write_rds(site_list, "data/site_analysis_dat_5min.rds")
