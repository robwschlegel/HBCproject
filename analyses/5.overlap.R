# overlap_workflow.R
# This code contains the workflow to quantify the overlap between predators and prey sensu Carroll et al. (2019)

# Setup -------------------------------------------------------------------

library(raster)
library(dplyr)
library(fuzzySim)
library(tidyverse)
library(geosphere)
library(sp)

# Data loading functions --------------------------------------------------

# Function for choosing the object name of data loaded from RData files
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Function to convert rasters to data.frames
rast_df <- function(rast){
  df_out <- as.data.frame(rast[[1]], xy = T) %>% 
    `colnames<-`(c("lon", "lat", "presence")) %>% 
    mutate(lon = round(lon, 5), lat = round(lat, 5),
           presence = as.numeric(presence)) %>% 
    na.omit()
  return(df_out)
}


# Species overlap functions -----------------------------------------------

#### overlap functions to accompany Carroll et al. (2019) 'A review of methods for quantifying predator-prey overlap,' Global Ecology and Biogeography 
## please see accompanying paper for detailed descriptions of metrics and their ecological interpretations

## area overlap - for binary data
# measures proportion of an area where two species co-occur (0-1)
area_overlapfn <- function(prey, pred, area){
  total_area <- sum(area, na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/total_area
}

## range overlap - for binary data
# measures the proportion of one species range where the other co-occurs (0-1)
range_overlapfn <- function(prey, pred, area){
  area_prey <- sum(area[prey > 0], na.rm = T)
  sum(area[pred > 0 & prey > 0], na.rm = T)/area_prey
}

## Schoener's D - for binary and biomass data
# density or probability of occurrence data
# measures how equally predator and prey share available resources (0-1)
schoeners_overlapfn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  1 - 0.5 * (sum(abs(p_prey-p_pred), na.rm = T))
}

## Bhattacharyya's coefficient - for binary and biomass data
# density or probability of occurrence data
# measures whether two species use space independently (0-1)
bhatta_coeffn <- function(prey, pred) {
  p_prey <- prey/sum(prey, na.rm = T)
  p_pred <- pred/sum(pred, na.rm = T)
  sum(sqrt(p_prey*p_pred), na.rm = T)
}


##To calculate the area (I ended up using this one only for the Arctic, but below it is found the way to
#calculate this in a general way. I do not know why at the moment of left_join in here now is not working
#the other way)

Arctic_AM <- loadRData(paste0("results/Rob_tests/Arctic_AM.RData"))

Arctic_AM <- Arctic_AM %>% mutate_at(1:2, round, 5)

# Load data (individual species) ---------------------------------------------------------------

d.predictions <- 'results/binary_final/NoData0Mask/'
# list the files
files <- list.files(d.predictions, pattern='tif$', full.names=TRUE )
files



# Data from Biomod2 output (extracted by depth mask)
Lsol_pres <- stack("results/binary_final/NoData0Mask/Lsol.tif")
Lsol_2050 <- stack("results/binary_final/NoData0Mask/Lsol_2050.tif")
Lsol_2100 <- stack("results/binary_final/NoData0Mask/Lsol_2100.tif")
#Lsol_pres_TSS <- loadRData("~/RES_project_HBC/HBCproject/Lsol/proj_presentEA/proj_presentEA_Lsol_ensemble_TSSbin.img")[[1]]

Goce_pres <- stack("results/binary_final/NoData0Mask/Goce.tif")
Goce_2050 <- stack("results/binary_final/NoData0Mask/Goce_2050.tif")
Goce_2100 <- stack("results/binary_final/NoData0Mask/Goce_2100.tif")

#Mdes_pres <- stack("results/binary_final/NoData0Mask/Mdes.tif")
#Mdes_2050 <- stack("C")
#Mdes_2100 <- stack("")

# Convert to data.frames
df_Lsol_pres <- rast_df(Lsol_pres) %>%
  dplyr::rename(prey = presence)
  
df_Goce_pres <- rast_df(Goce_pres) %>%
  dplyr::rename(pred = presence)

df_Lsol_2050 <- rast_df(Lsol_2050) %>%
  dplyr::rename(prey = presence)

df_Goce_2050 <- rast_df(Goce_2050) %>%
  dplyr::rename(pred = presence)

df_Lsol_2100 <- rast_df(Lsol_2100) %>%
  dplyr::rename(prey = presence)

df_Goce_2100 <- rast_df(Goce_2100) %>%
  dplyr::rename(pred = presence)


# Join data.frames
predprey <- left_join(df_Lsol_pres, df_Goce_pres, by = c("lon", "lat"))
# predprey <- predprey %>% mutate(pred = ifelse(is.na(pred), 0, pred)) %>%
#   mutate(prey = ifelse(is.na(prey), 0, prey)) ##to convert NA values to 0
ppLsolGoce <- left_join(df_Lsol_pres, df_Goce_pres, by = c("lon", "lat"))
ppLsolGoce2050 <- left_join(df_Lsol_2050, df_Goce_2050, by = c("lon", "lat"))
ppLsolGoce2100 <- left_join(df_Lsol_2100, df_Goce_2100, by = c("lon", "lat"))


#predpreyarea <- left_join(predprey, area, by = c("lon", "lat"))
#or
predpreyarea <- left_join(predprey, Arctic_AM[,c("lon", "lat", "sq_area")], by = c("lon", "lat"))

ppaLsolGoce <- left_join(ppLsolGoce, Arctic_AM[,c("lon", "lat", "sq_area")], by = c("lon", "lat"))
ppaLsolGoce2050 <- left_join(ppLsolGoce2050, Arctic_AM[,c("lon", "lat", "sq_area")], by = c("lon", "lat"))
ppaLsolGoce2100 <- left_join(ppLsolGoce2100, Arctic_AM[,c("lon", "lat", "sq_area")], by = c("lon", "lat"))


# Another method of extracting the raster info
# dat_Acla <- gridRecords(rst = project_present_Acla, pres.coords = df_Acla[c("lon", "lat")])
# dat_Aesc <- gridRecords(rst = project_present_Aesc, pres.coords = df_Aesc[c("lon", "lat")])


# Calculate overlap -------------------------------------------------------

### NB: Many of these functions require an 'area' object but it is not stated anywhere what this might be
# I assume that it comes from some of these other functions but it is not clear to me which ones that might be
# Perhaps one is meant to calculate the area of the study site beforehand and use that static value here

## area overlap - for binary data
# measures proportion of an area where two species co-occur
#area_overlapfn(prey = predprey$prey, pred = predprey$pred, area = ???)

area_overlapfn(prey = predpreyarea$prey, pred = predpreyarea$pred, area = predpreyarea$sq_area)

LsolGoce_area <- area_overlapfn(prey = ppaLsolGoce$prey, pred = ppaLsolGoce$pred, area = ppaLsolGoce$sq_area)
LsolGoce_area2050 <- area_overlapfn(prey = ppaLsolGoce2050$prey, pred = ppaLsolGoce2050$pred, area = ppaLsolGoce2050$sq_area)
LsolGoce_area2100 <- area_overlapfn(prey = ppaLsolGoce2100$prey, pred = ppaLsolGoce2100$pred, area = ppaLsolGoce2100$sq_area)

## range overlap - for binary data
# measures the proportion of one species range where the other co-occurs
range_overlapfn(prey = predpreyarea$prey, pred = predpreyarea$pred, area = predpreyarea$sq_area)

LsolGoce_range <- range_overlapfn(prey = ppaLsolGoce$prey, pred = ppaLsolGoce$pred, area = ppaLsolGoce$sq_area)
LsolGoce_range2050 <- range_overlapfn(prey = ppaLsolGoce2050$prey, pred = ppaLsolGoce2050$pred, area = ppaLsolGoce2050$sq_area)
LsolGoce_range2100 <- range_overlapfn(prey = ppaLsolGoce2100$prey, pred = ppaLsolGoce2100$pred, area = ppaLsolGoce2100$sq_area)

## Schoener's D
# density or probability of occurrence data
# measures how equally predator and prey share available resources
schoeners_overlapfn(prey = predpreyarea$prey, pred = predpreyarea$pred)

LsolGoce_Sch <- schoeners_overlapfn(prey = ppaLsolGoce$prey, pred = ppaLsolGoce$pred)
LsolGoce_Sch2050 <- schoeners_overlapfn(prey = ppaLsolGoce2050$prey, pred = ppaLsolGoce2050$pred)
LsolGoce_Sch2100 <- schoeners_overlapfn(prey = ppaLsolGoce2100$prey, pred = ppaLsolGoce2100$pred)

## Bhattacharyya's coefficient
# density or probability of occurrence data
# measures whether two species use space independently
bhatta_coeffn(prey = predpreyarea$prey, pred = predpreyarea$pred)

LsolGoce_Bh <- bhatta_coeffn(prey = ppaLsolGoce$prey, pred = ppaLsolGoce$pred)
LsolGoce_Bh2050 <- bhatta_coeffn(prey = ppaLsolGoce2050$prey, pred = ppaLsolGoce2050$pred)
LsolGoce_Bh2100 <- bhatta_coeffn(prey = ppaLsolGoce2100$prey, pred = ppaLsolGoce2100$pred)







## local index of collocation - for biomass
# estimates correlation of predator and prey densities
# loc_collocfn <- function(prey, pred) {
#   p_prey <- prey/sum(prey, na.rm = T)
#   p_pred <- pred/sum(pred, na.rm = T)
#   sum(p_prey*p_pred, na.rm = T)/(sqrt(sum(p_prey^2, na.rm = T)*sum(p_pred^2, na.rm = T)))
# }

## global index of collocation
# measures geographic distinctness by comparing centres of gravity and dispersion of sampled individuals
# glob_collocfn(prey_x = predpreyarea$lon, prey_y = predpreyarea$lat, prey = predpreyarea$prey,
#               pred_x = predpreyarea$lon, pred_y = predpreyarea$lat, pred = predpreyarea$pred)
# 
# ## AB ratio
# # measures predator production that can be attributed to spatial overlap with prey
# AB_overlapfn(prey = predpreyarea$prey, pred = predpreyarea$pred)


## asymmetrical alpha - for biomass
# measures pressure of predator on prey relative to underlying prey density
# asymmalpha_overlapfn <- function(prey, pred){
#   p_prey <- prey/sum(prey, na.rm = T)
#   p_pred <- pred/sum(pred, na.rm = T)
#   sum(p_pred*p_prey, na.rm = T)/sum(p_prey^2, na.rm = T)
# }

## biomass-weighted overlap (scaled to max) - for biomass
# measures amount of predator biomass interacting with prey relative to underlying prey biomass
# biomass_overlapfn <- function(prey, pred) {
#   sum((prey/max(prey, na.rm = T)) * (pred/max(pred, na.rm = T)), na.rm = T)/sum(prey/max(prey, na.rm = T), na.rm = T)
# }

## Hurlbert's overlap - for biomass
# measures interspecific encounter rate between predator and prey
# hurlbert_overlapfn <- function(prey, pred, area) {
#   area_occupied <- sum(area[pred > 0 | prey > 0], na.rm = T)
#   p_prey <- prey/sum(prey, na.rm = T)
#   p_pred <- pred/sum(pred, na.rm = T)
#   sum((p_pred*p_prey)/(area/area_occupied), na.rm = T)
# }

## global index of collocation - for biomass
# measures geographic distinctness by comparing centres of gravity and dispersion of sampled individuals
# glob_collocfn <- function(prey_x, prey_y, prey, pred_x, pred_y, pred){
#   prey_cgx <- sum(prey_x*prey, na.rm = T)/sum(prey, na.rm = T)
#   prey_cgy <- sum(prey_y*prey, na.rm = T)/sum(prey, na.rm = T)
#   prey_ix <- prey_x - prey_cgx
#   prey_iy <- prey_y - prey_cgy
#   prey_i <- sqrt(prey_ix^2 + prey_iy^2)
#   prey_inert <- sum(prey * (prey_i^2), na.rm = T)/sum(prey, na.rm = T)
#   pred_cgx <- sum(pred_x*pred, na.rm = T)/sum(pred, na.rm = T)
#   pred_cgy <- sum(pred_y*pred, na.rm = T)/sum(pred, na.rm = T)
#   pred_ix <- pred_x - pred_cgx
#   pred_iy <- pred_y - pred_cgy
#   pred_i <- sqrt(pred_ix^2 + pred_iy^2)
#   pred_inert <- sum(pred * (pred_i^2), na.rm = T)/sum(pred, na.rm = T)
#   GIC <- (((prey_cgx - pred_cgx)^2+(prey_cgy - pred_cgy)^2)/ (((prey_cgx-pred_cgx)^2+(prey_cgy-pred_cgy)^2)+prey_inert + pred_inert))
#   if(!is.na(GIC))
#     GIC <- 1-GIC
#   else GIC <- 1
#   GIC
# }

## AB ratio - for biomass
# measures predator production that can be attributed to spatial overlap with prey
# AB_overlapfn <- function(prey, pred) { 
#   mean((pred - mean(pred, na.rm = T)) * (prey - mean(prey, na.rm = T)), na.rm = T)/(mean(pred, na.rm = T) * mean(prey, na.rm = T)) 
# }


##calculate grid size

# grid_size.R
# Functions for finding the square area of pixels


# Setup -------------------------------------------------------------------

# Calculate the square kilometre surface area of each pixel
# sq_area <- function(df_single, pixel_width, pixel_height){
#   # Distance for longitude
#   lon_dist <- distm(c(df_single$lon-pixel_width, df_single$lat), c(df_single$lon+pixel_width, df_single$lat), fun = distHaversine)/1000
#   # Distance for latitude
#   lat_dist <- distm(c(df_single$lon, df_single$lat+pixel_height), c(df_single$lon, df_single$lat-pixel_height), fun = distHaversine)/1000
#   # Total area
#   sq_area <- data.frame(sq_area = lon_dist*lat_dist)
#   # Combine and exit
#   res <- cbind(df_single, sq_area)
#   return(res)
# }
# 
# # Run the function on a data.frame
# # round_diff: choose the decimal places to round the difference between pixel sizes
# # NB: This is just to make sure the grid cells are even, it's not used in the calculations
# grid_size <- function(df, round_diff = 6){
#   
#   # Find average size of pixels
#   unique_lon <- arrange(distinct(df["lon"]), lon) %>% 
#     mutate(diff = lon - lag(lon, default = first(lon))) %>% 
#     filter(diff != 0)
#   pixel_width <- mean(unique_lon$diff, na.rm = T)/2
#   pixel_width_range <- unique_lon %>% 
#     mutate(mean_diff = mean(diff, na.rm = T)) %>% 
#     summarise(pixel_width_range = round(mean(diff-mean_diff, na.rm = T), round_diff)) %>% 
#     as.numeric()
#   unique_lat <- arrange(distinct(df["lat"]), lat) %>% 
#     mutate(diff = lat - lag(lat, default = first(lat))) %>% 
#     filter(diff != 0)
#   pixel_height <- mean(unique_lat$diff, na.rm = T)/2
#   pixel_height_range <- unique_lat %>% 
#     mutate(mean_diff = mean(diff, na.rm = T)) %>% 
#     summarise(pixel_height_range = round(mean(diff-mean_diff, na.rm = T), round_diff)) %>% 
#     as.numeric()
#   
#   # Calculate square area for each grid cell
#   df_area <- df %>% 
#     mutate(plyr_idx = 1:n()) %>% 
#     plyr::ddply(c("plyr_idx"), sq_area, .parallel = F,
#                 pixel_height = pixel_height, pixel_width = pixel_width) %>% 
#     mutate(plyr_idx = NULL)
#   
#   # Report on how much variance there is between pixels
#   # These values should be very low
#   print(paste0("Height diff: ",pixel_height_range,"   Width diff: ",pixel_width_range))
#   return(df_area)
# }
# 
# 
# # Data --------------------------------------------------------------------
# 
# # Choose a species raster and convert to data.frame
# 
# e_layers <- 'data/Env_var/BO2.1pres/Fish_EA/'
# # list the files
# files <- list.files(e_layers, pattern='tif$', full.names=TRUE )
# files
# # 
# rast_sp_area <- stack("data/Env_var/BO2.1pres/Fish_EA/Benthic.Temperature.max.tif") 
# df_area <- as.data.frame(rast_sp_area, xy = T) %>% 
#   `colnames<-`(c("lon", "lat", "presence")) 
# 
# # Calculate area ----------------------------------------------------------
# 
# # Find square area of each pixel
# # system.time(
# #   area_Acla <- grid_size(df_Acla, round_diff = 20)
# # ) # 57 seconds on 7 cores
# 
# area <- grid_size(df_area, round_diff = 20)

## local index of collocation
# estimates correlation of predator and prey densities
# loc_collocfn(prey = predpreyarea$prey, pred = predpreyarea$pred)
# 
# ## asymmetrical alpha
# # measures pressure of predator on prey relative to underlying prey density
# asymmalpha_overlapfn(prey = predpreyarea$prey, pred = predpreyarea$pred)
# 
# ## biomass-weighted overlap (scaled to max)
# # measures amount of predator biomass interacting with prey relative to underlying prey biomass
# biomass_overlapfn(prey = predpreyarea$prey, pred = predpreyarea$pred)
# 
# ## Hurlbert's overlap
# # measures interspecific encounter rate between predator and prey
# hurlbert_overlapfn(prey = predpreyarea$prey, pred = predpreyarea$pred, area = predpreyarea$sq_area)