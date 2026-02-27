# Code to produce walking isochrones

library(osrm)
library(tidyverse)
library(sf)
library(tictoc)

site_data <- read_rds("data/site_cd_lookup.rds")

custom_toc <- function(tic, toc, msg, info){
  str_c("Time elapsed: ", round((toc - tic)/60,2), " minutes!")
}



all_site_walking <- map_dfr(seq_along(site_data$name), 
                         function(a){
                           tic()
                           
                           df <- osrmIsochrone(site_data[a,], breaks = seq(from = 0, to = 10, by = 5), res = 200, osrm.profile = "foot") %>%
                             mutate(name = site_data$name[a])
                            toc(func.toc = custom_toc)   
                            return(df)
                         }
                         )
  
  
highres_10min <- all_site_walking %>%
  filter(isomax == 10) %>%
  st_make_valid() %>%
  geos_concave_hull(allow_holes = F, ratio = 0.2) %>%
  st_as_sf() %>%
  bind_cols(all_site_walking %>% st_drop_geometry())

write_rds(highres_10min, "data/highres10min.rds")

highres_5min <- all_site_walking %>%
  filter(isomax == 5) %>%
  st_make_valid() %>%
  geos_concave_hull(allow_holes = F, ratio = 0.2) %>%
  st_as_sf() %>%
  bind_cols(all_site_walking %>% st_drop_geometry())

write_rds(highres_5min, "data/highres5min.rds")
