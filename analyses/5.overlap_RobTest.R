# overlap_workflow.R
# This code contains the workflow to quantify the overlap between predators and prey sensu Carroll et al. (2019)


# Setup -------------------------------------------------------------------

library(raster)
# library(dplyr) # This is loaded via tidyverse
library(fuzzySim)
library(tidyverse)
library(geosphere)
# library(sp) # This is loaded via raster


# Data loading functions --------------------------------------------------

# Function for choosing the object name of data loaded from RData files
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Function to convert rasters to data.frames
load_rast_sps <- function(rast_file){
  rast <- raster(rast_file)
  # NB: If the pathway changes it may be necessary to change the '4' to a different number
  sps_name <- gsub(".tif", "", lapply(str_split(rast_file, "/"), "[[", 4)[[1]])
  df_out <- as.data.frame(rast[[1]], xy = T) %>% 
    `colnames<-`(c("lon", "lat", "presence")) %>% 
    mutate(lon = round(lon, 5), lat = round(lat, 5),
           sps = sps_name,
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

# Load data (individual species) ---------------------------------------------------------------

## To calculate the area (I ended up using this one only for the Arctic, but below it is found the way to
# calculate this in a general way. I do not know why at the moment of left_join in here now is not working
# the other way)

Arctic_AM <- loadRData(paste0("results/Rob_tests/Arctic_AM.RData"))
# Arctic_AM <- loadRData("Arctic_AM.RData")

Arctic_AM <- Arctic_AM %>% mutate_at(1:2, round, 5)

d.predictions <- 'results/binary_final/NoData0Mask'
# list the files
files <- list.files(d.predictions, pattern = 'tif$', full.names = TRUE )
files

# Load, convert to data.frames, and join
df_Sps <- plyr::ldply(files, load_rast_sps, .parallel = F) %>%
  pivot_wider(names_from = "sps", values_from = "presence") %>%
  left_join(Arctic_AM, by = c("lon", "lat")) %>% 
  mutate(Amph_N = case_when(Goce %in% "1" ~ 1,
                            Gset %in% "1" ~ 1,
                            TRUE ~ NA_real_)) %>%
  mutate(Amph_N_2050 = case_when(Goce_2050 %in% "1" ~ 1,
                                 Gset_2050 %in% "1" ~ 1,
                                 TRUE ~ NA_real_)) %>%
  mutate(Amph_N_2100 = case_when(Goce_2100 %in% "1" ~ 1,
                                 Gset_2100 %in% "1" ~ 1,
                                 TRUE ~ NA_real_)) %>% #NA_real_ makes it numerical instead of character
  mutate(MdesMqua = case_when(Mdes %in% "1" ~ 1,
                              Mqua %in% "1" ~ 1,
                              TRUE ~ NA_real_)) %>%
  mutate(MdesMqua_2050 = case_when(Mdes_2050 %in% "1" ~ 1,
                                   Mqua_2050 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>%
  mutate(MdesMqua_2100 = case_when(Mdes_2100 %in% "1" ~ 1,
                                   Mqua_2100 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>%
  mutate(MscoGtri = case_when(Msco %in% "1" ~ 1,
                              Gtri %in% "1" ~ 1,
                              TRUE ~ NA_real_)) %>%
  mutate(MscoGtri_2050 = case_when(Msco_2050 %in% "1" ~ 1,
                                   Gtri_2050 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>%
  mutate(MscoGtri_2100 = case_when(Msco_2100 %in% "1" ~ 1,
                                   Gtri_2100 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>% 
  mutate(Kelp = case_when(Acla %in% "1" ~ 1,
                          Aesc %in% "1" ~ 1,
                          Slat %in% "1" ~ 1,
                          TRUE ~ NA_real_)) %>%
  mutate(Kelp_2050 = case_when(Acla_2050 %in% "1" ~ 1,
                               Aesc_2050 %in% "1" ~ 1,
                               Slat_2050 %in% "1" ~ 1,
                               TRUE ~ NA_real_)) %>%
  mutate(Kelp_2100 = case_when(Acla_2100 %in% "1" ~ 1,
                               Aesc_2100 %in% "1" ~ 1,
                               Slat_2100 %in% "1" ~ 1,
                               TRUE ~ NA_real_))


# Create predator prey comparison guide
# Kelp sps as prey: Acla, Aesc, Lsol, Slat, Zmar
# Amphipods that should be considered predators for Kelp but preys for Fish: Goce, Gset
# Sculpin fish that are predators of Amphipods: Gtri, Mdes, Msco, Mqua
# NB: There is probably a better way to do this...
pred_prey_df <- data.frame(idx = c(1:23),
                          pred = c("Goce", "Gset", "Goce", "Gset", "Goce", "Gset","Goce", "Gset","Goce", "Gset", "Gtri", "Mdes", "Msco", "Mqua", "Gtri", "Mdes", "Msco", "Mqua", "Zmar", "Lsol", "Kelp", "Amph_N", "Amph_N"),
                          prey = c("Acla", "Acla","Aesc", "Aesc", "Lsol", "Lsol", "Slat", "Slat", "Zmar", "Zmar", "Goce","Goce","Goce","Goce", "Gset", "Gset", "Gset", "Gset", "Amph_N", "Amph_N", "Amph_N", "MdesMqua", "MscoGtri"))

# Calculate overlap -------------------------------------------------------

### NB: Many of these functions require an 'area' object but it is not stated anywhere what this might be
# I assume that it comes from some of these other functions but it is not clear to me which ones that might be
# Perhaps one is meant to calculate the area of the study site beforehand and use that static value here

## area overlap - for binary data
# measures proportion of an area where two species co-occur
#area_overlapfn(prey = predprey$prey, pred = predprey$pred, area = ???)

# A function that expects a one row data.frame with a 'pred' and 'prey' column
pred_prey_funs <- function(df){
  
  # Get necessary columns
  prey_pres <- df_Sps[,df$prey]
  prey_2050 <- df_Sps[,paste0(df$prey,"_2050")]
  prey_2100 <- df_Sps[,paste0(df$prey,"_2100")]
  pred_pres <- df_Sps[,df$pred]
  pred_2050 <- df_Sps[,paste0(df$pred,"_2050")]
  pred_2100 <- df_Sps[,paste0(df$pred,"_2100")]
  sq_area <- df_Sps[, "sq_area"]
  
  # NB: There is probably a better way to do this...
  
  ## Area overlap
  area_overlap_pres <- area_overlapfn(prey = prey_pres, pred = pred_pres, area = sq_area)
  area_overlap_2050 <- area_overlapfn(prey = prey_2050, pred = pred_2050, area = sq_area)
  area_overlap_2100 <- area_overlapfn(prey = prey_2100, pred = pred_2100, area = sq_area)
  ## range overlap - for binary data
  # measures the proportion of one species range where the other co-occurs
  range_overlap_pres <- range_overlapfn(prey = prey_pres, pred = pred_pres, area = sq_area)
  range_overlap_2050 <- range_overlapfn(prey = prey_2050, pred = pred_2050, area = sq_area)
  range_overlap_2100 <- range_overlapfn(prey = prey_2100, pred = pred_2100, area = sq_area)
  ## Schoener's D
  # density or probability of occurrence data
  # measures how equally predator and prey share available resources
  schoeners_overlap_pres <- schoeners_overlapfn(prey = prey_pres, pred = pred_pres)
  schoeners_overlap_2050 <- schoeners_overlapfn(prey = prey_2050, pred = pred_2050)
  schoeners_overlap_2100 <- schoeners_overlapfn(prey = prey_2100, pred = pred_2100)
  ## Bhattacharyya's coefficient
  # density or probability of occurrence data
  # measures whether two species use space independently
  bhatta_coef_pres <- bhatta_coeffn(prey = prey_pres, pred = pred_pres)
  bhatta_coef_2050 <- bhatta_coeffn(prey = prey_2050, pred = pred_2050)
  bhatta_coef_2100 <- bhatta_coeffn(prey = prey_2100, pred = pred_2100)

  
  # NB: Add more functions as desired
  
  # Combine and exit
  res <- data.frame(pred = df$pred, prey = df$prey, 
                    area_overlap_pres, area_overlap_2050, area_overlap_2100,
                    range_overlap_pres, range_overlap_2050, range_overlap_2100,
                    schoeners_overlap_pres, schoeners_overlap_2050, schoeners_overlap_2100,
                    bhatta_coef_pres, bhatta_coef_2050, bhatta_coef_2100)
  return(res)
}

# Run the pred/prey functions on the index data.frame
pred_prey_res <- plyr::ddply(pred_prey_df, c("idx"), pred_prey_funs, .parallel = F)
write_csv(pred_prey_res, "results/EcologicalInteractions/pred_prey_res.csv")



##INVASIVE ANALOGUES 

# Load, convert to data.frames, and join
df_Sps <- plyr::ldply(files, load_rast_sps, .parallel = F) %>% 
  pivot_wider(names_from = "sps", values_from = "presence") %>% 
  left_join(Arctic_AM, by = c("lon", "lat")) %>%
  mutate(Amph_N = case_when(Goce %in% "1" ~ 1,
                            Gset %in% "1" ~ 1,
                            TRUE ~ NA_real_)) %>%
  mutate(Amph_N_2050 = case_when(Goce_2050 %in% "1" ~ 1,
                                 Gset_2050 %in% "1" ~ 1,
                                 TRUE ~ NA_real_)) %>%
  mutate(Amph_N_2100 = case_when(Goce_2100 %in% "1" ~ 1,
                                 Gset_2100 %in% "1" ~ 1,
                                 TRUE ~ NA_real_)) %>% #NA_real_ makes it numerical instead of character
  mutate(MdesMqua = case_when(Mdes %in% "1" ~ 1,
                            Mqua %in% "1" ~ 1,
                            TRUE ~ NA_real_)) %>%
  mutate(MdesMqua_2050 = case_when(Mdes_2050 %in% "1" ~ 1,
                                 Mqua_2050 %in% "1" ~ 1,
                                 TRUE ~ NA_real_)) %>%
  mutate(MdesMqua_2100 = case_when(Mdes_2100 %in% "1" ~ 1,
                                 Mqua_2100 %in% "1" ~ 1,
                                 TRUE ~ NA_real_)) %>%
  mutate(MscoGtri = case_when(Msco %in% "1" ~ 1,
                              Gtri %in% "1" ~ 1,
                              TRUE ~ NA_real_)) %>%
  mutate(MscoGtri_2050 = case_when(Msco_2050 %in% "1" ~ 1,
                                   Gtri_2050 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>%
  mutate(MscoGtri_2100 = case_when(Msco_2100 %in% "1" ~ 1,
                                   Gtri_2100 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>% 
  mutate(F_NN = case_when(Aatl %in% "1" ~ 1,
                              Aunc %in% "1" ~ 1,
                              TRUE ~ NA_real_)) %>%
  mutate(F_NN_2050 = case_when(Aatl_2050 %in% "1" ~ 1,
                                   Aunc_2050 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>%
  mutate(F_NN_2100 = case_when(Aatl_2100 %in% "1" ~ 1,
                                   Aunc_2100 %in% "1" ~ 1,
                                   TRUE ~ NA_real_)) %>% 
  mutate(Kelp = case_when(Acla %in% "1" ~ 1,
                          Aesc %in% "1" ~ 1,
                          Slat %in% "1" ~ 1,
                          TRUE ~ NA_real_)) %>%
  mutate(Kelp_2050 = case_when(Acla_2050 %in% "1" ~ 1,
                               Aesc_2050 %in% "1" ~ 1,
                               Slat_2050 %in% "1" ~ 1,
                               TRUE ~ NA_real_)) %>%
  mutate(Kelp_2100 = case_when(Acla_2100 %in% "1" ~ 1,
                               Aesc_2100 %in% "1" ~ 1,
                               Slat_2100 %in% "1" ~ 1,
                               TRUE ~ NA_real_))

#df_Sps <- ifelse(is.na(df_Sps), 0, df_Sps)


# Create predator prey comparison guide
# Kelp sps as prey: Acla, Aesc, Lsol, Slat, Zmar
# Amphipods that should be considered predators for Kelp but preys for Fish: Goce, Gset
# Sculpin fish that are predators of Amphipods: Gtri, Mdes, Msco, Mqua
# NB: There is probably a better way to do this...
pred_prey_df <- data.frame(idx = c(1:19),
                           pred = c("Lsol", "Zmar", "Acla", "Aesc", "Slat", "Goce","Gset", "Mdes","Mdes", "Mqua", "Mqua", "Gtri", "Gtri", "Msco", "Msco", "Amph_N", "MdesMqua", "MscoGtri", "Kelp"),
                           prey = c("Dcon", "Dcon","Dcon", "Dcon", "Dcon", "Gtig", "Gtig", "Aatl", "Aunc", "Aatl", "Aunc","Aatl", "Aunc","Aatl", "Aunc", "Gtig", "F_NN", "F_NN", "Dcon"))

# Calculate overlap -------------------------------------------------------

### NB: Many of these functions require an 'area' object but it is not stated anywhere what this might be
# I assume that it comes from some of these other functions but it is not clear to me which ones that might be
# Perhaps one is meant to calculate the area of the study site beforehand and use that static value here

## area overlap - for binary data
# measures proportion of an area where two species co-occur
#area_overlapfn(prey = predprey$prey, pred = predprey$pred, area = ???)

# A function that expects a one row data.frame with a 'pred' and 'prey' column
pred_prey_funs <- function(df){
  
  # Get necessary columns
  prey_pres <- df_Sps[,df$prey]
  prey_2050 <- df_Sps[,paste0(df$prey,"_2050")]
  prey_2100 <- df_Sps[,paste0(df$prey,"_2100")]
  pred_pres <- df_Sps[,df$pred]
  pred_2050 <- df_Sps[,paste0(df$pred,"_2050")]
  pred_2100 <- df_Sps[,paste0(df$pred,"_2100")]
  sq_area <- df_Sps[, "sq_area"]
  
  # NB: There is probably a better way to do this...
  
  ## Area overlap
  area_overlap_pres <- area_overlapfn(prey = prey_pres, pred = pred_pres, area = sq_area)
  area_overlap_2050 <- area_overlapfn(prey = prey_2050, pred = pred_2050, area = sq_area)
  area_overlap_2100 <- area_overlapfn(prey = prey_2100, pred = pred_2100, area = sq_area)
  ## range overlap - for binary data
  # measures the proportion of one species range where the other co-occurs
  range_overlap_pres <- range_overlapfn(prey = prey_pres, pred = pred_pres, area = sq_area)
  range_overlap_2050 <- range_overlapfn(prey = prey_2050, pred = pred_2050, area = sq_area)
  range_overlap_2100 <- range_overlapfn(prey = prey_2100, pred = pred_2100, area = sq_area)
  ## Schoener's D
  # density or probability of occurrence data
  # measures how equally predator and prey share available resources
  schoeners_overlap_pres <- schoeners_overlapfn(prey = prey_pres, pred = pred_pres)
  schoeners_overlap_2050 <- schoeners_overlapfn(prey = prey_2050, pred = pred_2050)
  schoeners_overlap_2100 <- schoeners_overlapfn(prey = prey_2100, pred = pred_2100)
  ## Bhattacharyya's coefficient
  # density or probability of occurrence data
  # measures whether two species use space independently
  bhatta_coef_pres <- bhatta_coeffn(prey = prey_pres, pred = pred_pres)
  bhatta_coef_2050 <- bhatta_coeffn(prey = prey_2050, pred = pred_2050)
  bhatta_coef_2100 <- bhatta_coeffn(prey = prey_2100, pred = pred_2100)
  
  
  # NB: Add more functions as desired
  
  # Combine and exit
  res <- data.frame(pred = df$pred, prey = df$prey, 
                    area_overlap_pres, area_overlap_2050, area_overlap_2100,
                    range_overlap_pres, range_overlap_2050, range_overlap_2100,
                    schoeners_overlap_pres, schoeners_overlap_2050, schoeners_overlap_2100,
                    bhatta_coef_pres, bhatta_coef_2050, bhatta_coef_2100)
  return(res)
}

# Run the pred/prey functions on the index data.frame
pred_prey_res <- plyr::ddply(pred_prey_df, c("idx"), pred_prey_funs, .parallel = F)
write_csv(pred_prey_res, "results/AnalogueSpecies/pred_prey_resAnalogue.csv")

# Calculate grid size -----------------------------------------------------

# grid_size.R
# Functions for finding the square area of pixels

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