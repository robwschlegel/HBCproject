library(biomod2)
library(ggplot2)
library(ggtext)
library(stringr)
library(tidyverse)
library(Rcpp)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
} ## To save the scripts later


# Load data ---------------------------------------------------------------


##Occurrence points

sps_files <- dir("data/rarefied_occ/clean_NoCurrVel/", full.names = T)
sps_names <- str_remove(dir("data/rarefied_occ/clean_NoCurrVel/", full.names = F), pattern = "_rarefied_points.csv")

# Choose a species
# Used for testing
sps_choice <- sps_files[7]
sps <- sps_names[7]



# Load biomod_models: ------------------------------------------------------


#sps_choice <- sps_files[9]
#sps <- sps_names[9]

# Load chosen biomod_model and print evaluation scores
# biomod_model <- loadRData(paste0(sps_names[8],"/",sps_names[8],".",sps_names[8],".models.out"))
# biomod_model # print summary



# Mean response curves ZOOBENTHOS, NATIVES ----------------------------------------------------

## Load model data per species
# NB: Due to BIOMOD2 internal structure, much of this needs to be run in the global environment...

# Goce

sps_choice <- sps_files[8]
sps <- sps_names[8]


biomod_model_Goce <- loadRData(paste0(sps_names[8],"/",sps_names[8],".",sps_names[8],".models.out"))
sp_name_Goce <- BIOMOD_LoadModels(biomod_model_Goce, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Goce <- response.plot2(models  = sp_name_Goce,
                              Data = get_formal_data(biomod_model_Goce, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Goce, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Goce, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Goce <- function(biomod_model_Goce, sp_dat_Goce){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Goce))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Goce %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Goce <- sps_response_data_Goce(biomod_model_Goce, sp_dat_Goce)


# Plot RC for Goce models

response_curve_species_mean_Goce <- curve_data_Goce %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Goce)

ggsave("figures/RC_Goce.png", response_curve_species_mean_Goce, height = 8, width = 12, dpi = 600)

rm(list = grep("Goce_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Goce",names(.GlobalEnv),value = TRUE))

# Gset

sps_choice <- sps_files[9]
sps <- sps_names[9]

# Load chosen biomod_model and print evaluation scores

biomod_model_Gset <- loadRData(paste0(sps_names[9],"/",sps_names[9],".",sps_names[9],".models.out"))
sp_name_Gset <- BIOMOD_LoadModels(biomod_model_Gset, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Gset <- response.plot2(models  = sp_name_Gset,
                              Data = get_formal_data(biomod_model_Gset, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Gset, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Gset, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Gset <- function(biomod_model_Gset, sp_dat_Gset){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Gset))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Gset %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Gset <- sps_response_data_Gset(biomod_model_Gset, sp_dat_Gset)


# Plot RC for Gset models

response_curve_species_mean_Gset <- curve_data_Gset %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Gset)

ggsave("figures/RC_Gset.png", response_curve_species_mean_Gset, height = 8, width = 12, dpi = 600)

rm(list = grep("Gset_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Gset",names(.GlobalEnv),value = TRUE))

# Lsax

sps_choice <- sps_files[13]
sps <- sps_names[13]

# Load chosen biomod_model and print evaluation scores

biomod_model_Lsax <- loadRData(paste0(sps_names[13],"/",sps_names[13],".",sps_names[13],".models.out"))
sp_name_Lsax <- BIOMOD_LoadModels(biomod_model_Lsax, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Lsax <- response.plot2(models  = sp_name_Lsax,
                              Data = get_formal_data(biomod_model_Lsax, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Lsax, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Lsax, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Lsax <- function(biomod_model_Lsax, sp_dat_Lsax){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Lsax))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Lsax %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Lsax <- sps_response_data_Lsax(biomod_model_Lsax, sp_dat_Lsax)


# Plot RC for Lsax models

response_curve_species_mean_Lsax <- curve_data_Lsax %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Lsax)

ggsave("figures/RC_Lsax.png", response_curve_species_mean_Lsax, height = 8, width = 12, dpi = 600)

rm(list = grep("Lsax_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Lsax",names(.GlobalEnv),value = TRUE))

# Mtru

sps_choice <- sps_files[19]
sps <- sps_names[19]

biomod_model_Mtru <- loadRData(paste0(sps_names[19],"/",sps_names[19],".",sps_names[19],".models.out"))
sp_name_Mtru <- BIOMOD_LoadModels(biomod_model_Mtru, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Mtru <- response.plot2(models  = sp_name_Mtru,
                              Data = get_formal_data(biomod_model_Mtru, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Mtru, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Mtru, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Mtru <- function(biomod_model_Mtru, sp_dat_Mtru){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Mtru))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Mtru %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Mtru <- sps_response_data_Mtru(biomod_model_Mtru, sp_dat_Mtru)


# Plot RC for Mtru models

response_curve_species_mean_Mtru <- curve_data_Mtru %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Mtru)

ggsave("figures/RC_Mtru.png", response_curve_species_mean_Mtru, height = 8, width = 12, dpi = 600)

rm(list = grep("Mtru_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Mtru",names(.GlobalEnv),value = TRUE))

# Sdro

sps_choice <- sps_files[20]
sps <- sps_names[20]


biomod_model_Sdro <- loadRData(paste0(sps_names[20],"/",sps_names[20],".",sps_names[20],".models.out"))
sp_name_Sdro <- BIOMOD_LoadModels(biomod_model_Sdro, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Sdro <- response.plot2(models  = sp_name_Sdro,
                              Data = get_formal_data(biomod_model_Sdro, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Sdro, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Sdro, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Sdro <- function(biomod_model_Sdro, sp_dat_Sdro){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Sdro))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Sdro %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Sdro <- sps_response_data_Sdro(biomod_model_Sdro, sp_dat_Sdro)


# Plot RC for Sdro models/species

response_curve_species_mean_Sdro <- curve_data_Sdro %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Sdro)

ggsave("figures/RC_Sdro.png", response_curve_species_mean_Sdro, height = 8, width = 12, dpi = 600)

rm(list = grep("Sdro_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Sdro",names(.GlobalEnv),value = TRUE))


# All ZB together now
curve_data_ALL_ZB <- rbind(curve_data_Goce, curve_data_Gset, curve_data_Lsax, curve_data_Mtru, curve_data_Sdro)

response_curve_species_mean_ALL_ZB <- curve_data_ALL_ZB %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'PuBu') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())

plot(response_curve_species_mean_ALL_ZB)

ggsave("figures/RC_ALL_ZB.png", response_curve_species_mean_ALL_ZB, height = 8, width = 12, dpi = 600)



# Mean response curves: FISH, NATIVES -------------------------------------

# Gacu

sps_choice <- sps_files[7]
sps <- sps_names[7]


biomod_model_Gacu <- loadRData(paste0(sps_names[7],"/",sps_names[7],".",sps_names[7],".models.out"))
sp_name_Gacu <- BIOMOD_LoadModels(biomod_model_Gacu, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Gacu <- response.plot2(models  = sp_name_Gacu,
                              Data = get_formal_data(biomod_model_Gacu, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Gacu, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Gacu, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Gacu <- function(biomod_model_Gacu, sp_dat_Gacu){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Gacu))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Gacu %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Mean" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Gacu <- sps_response_data_Gacu(biomod_model_Gacu, sp_dat_Gacu)


# Plot RC for Gacu models

response_curve_species_mean_Gacu <- curve_data_Gacu %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'RdBu') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Gacu)

ggsave("figures/RC_Gacu.png", response_curve_species_mean_Gacu, height = 8, width = 12, dpi = 600)

rm(list = grep("Gacu_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Gacu",names(.GlobalEnv),value = TRUE))

# Gtri

sps_choice <- sps_files[11]
sps <- sps_names[11]

# Load chosen biomod_model and print evaluation scores

biomod_model_Gtri <- loadRData(paste0(sps_names[11],"/",sps_names[11],".",sps_names[11],".models.out"))
sp_name_Gtri <- BIOMOD_LoadModels(biomod_model_Gtri, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Gtri <- response.plot2(models  = sp_name_Gtri,
                              Data = get_formal_data(biomod_model_Gtri, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Gtri, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Gtri, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Gtri <- function(biomod_model_Gtri, sp_dat_Gtri){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Gtri))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Gtri %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Mean" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Gtri <- sps_response_data_Gtri(biomod_model_Gtri, sp_dat_Gtri)


# Plot RC for Gtri models

response_curve_species_mean_Gtri <- curve_data_Gtri %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'RdBu') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Gtri)

ggsave("figures/RC_Gtri.png", response_curve_species_mean_Gtri, height = 8, width = 12, dpi = 600)

rm(list = grep("Gtri_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Gtri",names(.GlobalEnv),value = TRUE))

# Mdes

sps_choice <- sps_files[16]
sps <- sps_names[16]

# Load chosen biomod_model and print evaluation scores

biomod_model_Mdes <- loadRData(paste0(sps_names[16],"/",sps_names[16],".",sps_names[16],".models.out"))
sp_name_Mdes <- BIOMOD_LoadModels(biomod_model_Mdes, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Mdes <- response.plot2(models  = sp_name_Mdes,
                              Data = get_formal_data(biomod_model_Mdes, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Mdes, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Mdes, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Mdes <- function(biomod_model_Mdes, sp_dat_Mdes){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Mdes))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Mdes %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Lt.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Mdes <- sps_response_data_Mdes(biomod_model_Mdes, sp_dat_Mdes)


# Plot RC for Mdes models

response_curve_species_mean_Mdes <- curve_data_Mdes %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'RdBu') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Mdes)

ggsave("figures/RC_Mdes.png", response_curve_species_mean_Mdes, height = 8, width = 12, dpi = 600)

rm(list = grep("Mdes_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Mdes",names(.GlobalEnv),value = TRUE))

# Mqua

sps_choice <- sps_files[17]
sps <- sps_names[17]

biomod_model_Mqua <- loadRData(paste0(sps_names[17],"/",sps_names[17],".",sps_names[17],".models.out"))
sp_name_Mqua <- BIOMOD_LoadModels(biomod_model_Mqua, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Mqua <- response.plot2(models  = sp_name_Mqua,
                              Data = get_formal_data(biomod_model_Mqua, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Mqua, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Mqua, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Mqua <- function(biomod_model_Mqua, sp_dat_Mqua){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Mqua))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Mqua %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Lt.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)", 
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Mqua <- sps_response_data_Mqua(biomod_model_Mqua, sp_dat_Mqua)


# Plot RC for Mqua models

response_curve_species_mean_Mqua <- curve_data_Mqua %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'RdBu') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Mqua)

ggsave("figures/RC_Mqua.png", response_curve_species_mean_Mqua, height = 8, width = 12, dpi = 600)

rm(list = grep("Mqua_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Mqua",names(.GlobalEnv),value = TRUE))

# Msco

sps_choice <- sps_files[18]
sps <- sps_names[18]


biomod_model_Msco <- loadRData(paste0(sps_names[18],"/",sps_names[18],".",sps_names[18],".models.out"))
sp_name_Msco <- BIOMOD_LoadModels(biomod_model_Msco, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Msco <- response.plot2(models  = sp_name_Msco,
                              Data = get_formal_data(biomod_model_Msco, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Msco, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Msco, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Msco <- function(biomod_model_Msco, sp_dat_Msco){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Msco))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Msco %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Mean" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Msco <- sps_response_data_Msco(biomod_model_Msco, sp_dat_Msco)


# Plot RC for Msco models/species

response_curve_species_mean_Msco <- curve_data_Msco %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'RdBu') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Msco)

ggsave("figures/RC_Msco.png", response_curve_species_mean_Msco, height = 8, width = 12, dpi = 600)

rm(list = grep("Msco_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Msco",names(.GlobalEnv),value = TRUE))

# All Fish together now
curve_data_ALL_F <- rbind(curve_data_Gacu, curve_data_Gtri, curve_data_Mdes, curve_data_Mqua, curve_data_Msco)

response_curve_species_mean_ALL_F <- curve_data_ALL_F %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'Oranges') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())

plot(response_curve_species_mean_ALL_F)

ggsave("figures/RC_ALL_F.png", response_curve_species_mean_ALL_F, height = 8, width = 12, dpi = 600)



# Mean response curves, PHYTOBENTHOS, NATIVES -----------------------------

# Acla

sps_choice <- sps_files[2]
sps <- sps_names[2]


biomod_model_Acla <- loadRData(paste0(sps_names[2],"/",sps_names[2],".",sps_names[2],".models.out"))
sp_name_Acla <- BIOMOD_LoadModels(biomod_model_Acla, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Acla <- response.plot2(models  = sp_name_Acla,
                              Data = get_formal_data(biomod_model_Acla, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Acla, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Acla, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Acla <- function(biomod_model_Acla, sp_dat_Acla){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Acla))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Acla %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                 expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                 expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Acla <- sps_response_data_Acla(biomod_model_Acla, sp_dat_Acla)


# Plot RC for Acla models

response_curve_species_mean_Acla <- curve_data_Acla %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 7) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Acla)

ggsave("figures/RC_Acla.png", response_curve_species_mean_Acla, height = 8, width = 12, dpi = 600)

rm(list = grep("Acla_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Acla",names(.GlobalEnv),value = TRUE))

# Aesc

sps_choice <- sps_files[3]
sps <- sps_names[3]

# Load chosen biomod_model and print evaluation scores

biomod_model_Aesc <- loadRData(paste0(sps_names[3],"/",sps_names[3],".",sps_names[3],".models.out"))
sp_name_Aesc <- BIOMOD_LoadModels(biomod_model_Aesc, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Aesc <- response.plot2(models  = sp_name_Aesc,
                              Data = get_formal_data(biomod_model_Aesc, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Aesc, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Aesc, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Aesc <- function(biomod_model_Aesc, sp_dat_Aesc){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Aesc))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Aesc %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                  expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                  expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                  expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                  expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                  expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                  expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
                                 run_PA_model = paste0(run,"_",PA,"_",model))
           return(sp_res)
}

curve_data_Aesc <- sps_response_data_Aesc(biomod_model_Aesc, sp_dat_Aesc)


# Plot RC for Aesc models

response_curve_species_mean_Aesc <- curve_data_Aesc %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 7) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Aesc)

ggsave("figures/RC_Aesc.png", response_curve_species_mean_Aesc, height = 8, width = 12, dpi = 600)

rm(list = grep("Aesc_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Aesc",names(.GlobalEnv),value = TRUE))

# Lsol

sps_choice <- sps_files[14]
sps <- sps_names[14]

# Load chosen biomod_model and print evaluation scores

biomod_model_Lsol <- loadRData(paste0(sps_names[14],"/",sps_names[14],".",sps_names[14],".models.out"))
sp_name_Lsol <- BIOMOD_LoadModels(biomod_model_Lsol, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Lsol <- response.plot2(models  = sp_name_Lsol,
                              Data = get_formal_data(biomod_model_Lsol, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Lsol, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Lsol, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Lsol <- function(biomod_model_Lsol, sp_dat_Lsol){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Lsol))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Lsol %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Lt.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                 expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                 expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Lsol <- sps_response_data_Lsol(biomod_model_Lsol, sp_dat_Lsol)


# Plot RC for Lsol models

response_curve_species_mean_Lsol <- curve_data_Lsol %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 7) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Lsol)

ggsave("figures/RC_Lsol.png", response_curve_species_mean_Lsol, height = 8, width = 12, dpi = 600)

rm(list = grep("Lsol_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Lsol",names(.GlobalEnv),value = TRUE))

# Slat

sps_choice <- sps_files[21]
sps <- sps_names[21]

biomod_model_Slat <- loadRData(paste0(sps_names[21],"/",sps_names[21],".",sps_names[21],".models.out"))
sp_name_Slat <- BIOMOD_LoadModels(biomod_model_Slat, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Slat <- response.plot2(models  = sp_name_Slat,
                              Data = get_formal_data(biomod_model_Slat, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Slat, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Slat, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Slat <- function(biomod_model_Slat, sp_dat_Slat){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Slat))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Slat %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                 expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                 expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Slat <- sps_response_data_Slat(biomod_model_Slat, sp_dat_Slat)


# Plot RC for Slat models

response_curve_species_mean_Slat <- curve_data_Slat %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 7) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Slat)

ggsave("figures/RC_Slat.png", response_curve_species_mean_Slat, height = 8, width = 12, dpi = 600)

rm(list = grep("Slat_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Slat",names(.GlobalEnv),value = TRUE))

# Zmar

sps_choice <- sps_files[24]
sps <- sps_names[24]


biomod_model_Zmar <- loadRData(paste0(sps_names[24],"/",sps_names[24],".",sps_names[24],".models.out"))
sp_name_Zmar <- BIOMOD_LoadModels(biomod_model_Zmar, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Zmar <- response.plot2(models  = sp_name_Zmar,
                              Data = get_formal_data(biomod_model_Zmar, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Zmar, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Zmar, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Zmar <- function(biomod_model_Zmar, sp_dat_Zmar){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Zmar))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Zmar %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                 expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                 expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Zmar <- sps_response_data_Zmar(biomod_model_Zmar, sp_dat_Zmar)


# Plot RC for Zmar models/species

response_curve_species_mean_Zmar <- curve_data_Zmar %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 7) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Zmar)

ggsave("figures/RC_Zmar.png", response_curve_species_mean_Zmar, height = 8, width = 12, dpi = 600)

rm(list = grep("Zmar_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Zmar",names(.GlobalEnv),value = TRUE))

# All Phytobenthos together now
curve_data_ALL_PB <- rbind(curve_data_Acla, curve_data_Aesc, curve_data_Lsol, curve_data_Slat, curve_data_Zmar)

response_curve_species_mean_ALL_PB <- curve_data_ALL_PB %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'BuGn') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())

plot(response_curve_species_mean_ALL_PB)

ggsave("figures/RC_ALL_PB.png", response_curve_species_mean_ALL_PB, height = 8, width = 12, dpi = 600)



# Mean response curves, PHYTOBENTHOS NON-NATIVES ----------------------------

# Cfra

sps_choice <- sps_files[5]
sps <- sps_names[5]


biomod_model_Cfra <- loadRData(paste0(sps_names[5],"/",sps_names[5],".",sps_names[5],".models.out"))
sp_name_Cfra <- BIOMOD_LoadModels(biomod_model_Cfra, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Cfra <- response.plot2(models  = sp_name_Cfra,
                              Data = get_formal_data(biomod_model_Cfra, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Cfra, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Cfra, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Cfra <- function(biomod_model_Cfra, sp_dat_Cfra){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Cfra))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Cfra %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.mean" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                 expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                 expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Cfra <- sps_response_data_Cfra(biomod_model_Cfra, sp_dat_Cfra)


# Plot RC for Cfra models

response_curve_species_mean_Cfra <- curve_data_Cfra %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 2) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Cfra)

ggsave("figures/RC_Cfra.png", response_curve_species_mean_Cfra, height = 8, width = 12, dpi = 600)

rm(list = grep("Cfra_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Cfra",names(.GlobalEnv),value = TRUE))

# Dcon

sps_choice <- sps_files[6]
sps <- sps_names[6]

# Load chosen biomod_model and print evaluation scores

biomod_model_Dcon <- loadRData(paste0(sps_names[6],"/",sps_names[6],".",sps_names[6],".models.out"))
sp_name_Dcon <- BIOMOD_LoadModels(biomod_model_Dcon, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Dcon <- response.plot2(models  = sp_name_Dcon,
                              Data = get_formal_data(biomod_model_Dcon, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Dcon, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Dcon, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Dcon <- function(biomod_model_Dcon, sp_dat_Dcon){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Dcon))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Dcon %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                  expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                  expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                  expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                  expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                  expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                  expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
                                 run_PA_model = paste0(run,"_",PA,"_",model))
           return(sp_res)
}

curve_data_Dcon <- sps_response_data_Dcon(biomod_model_Dcon, sp_dat_Dcon)


# Plot RC for Dcon models

response_curve_species_mean_Dcon <- curve_data_Dcon %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 2) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Dcon)

ggsave("figures/RC_Dcon.png", response_curve_species_mean_Dcon, height = 8, width = 12, dpi = 600)

rm(list = grep("Dcon_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Dcon",names(.GlobalEnv),value = TRUE))

# Smut

sps_choice <- sps_files[22]
sps <- sps_names[22]

# Load chosen biomod_model and print evaluation scores

biomod_model_Smut <- loadRData(paste0(sps_names[22],"/",sps_names[22],".",sps_names[22],".models.out"))
sp_name_Smut <- BIOMOD_LoadModels(biomod_model_Smut, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Smut <- response.plot2(models  = sp_name_Smut,
                              Data = get_formal_data(biomod_model_Smut, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Smut, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Smut, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Smut <- function(biomod_model_Smut, sp_dat_Smut){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Smut))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Smut %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.mean" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                 expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                 expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Smut <- sps_response_data_Smut(biomod_model_Smut, sp_dat_Smut)


# Plot RC for Smut models

response_curve_species_mean_Smut <- curve_data_Smut %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 2) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Smut)

ggsave("figures/RC_Smut.png", response_curve_species_mean_Smut, height = 8, width = 12, dpi = 600)

rm(list = grep("Smut_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Smut",names(.GlobalEnv),value = TRUE))

# Upin

sps_choice <- sps_files[23]
sps <- sps_names[23]

biomod_model_Upin <- loadRData(paste0(sps_names[23],"/",sps_names[23],".",sps_names[23],".models.out"))
sp_name_Upin <- BIOMOD_LoadModels(biomod_model_Upin, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Upin <- response.plot2(models  = sp_name_Upin,
                              Data = get_formal_data(biomod_model_Upin, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Upin, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Upin, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Upin <- function(biomod_model_Upin, sp_dat_Upin){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Upin))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Upin %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.mean" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)",
                                 expl.name == "Surface.Iron.Mean" ~ "__F)__    Iron (umol.m-3)",
                                 expl.name == "parmean" ~ "__G)__    PAR (Einstein.m-2.day-1)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Upin <- sps_response_data_Upin(biomod_model_Upin, sp_dat_Upin)


# Plot RC for Upin models

response_curve_species_mean_Upin <- curve_data_Upin %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 2) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Upin)

ggsave("figures/RC_Upin.png", response_curve_species_mean_Upin, height = 8, width = 12, dpi = 600)

rm(list = grep("Upin_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Upin",names(.GlobalEnv),value = TRUE))

# All Phytobenthos Non-Native together now
curve_data_ALL_PB_NN <- rbind(curve_data_Cfra, curve_data_Dcon, curve_data_Smut, curve_data_Upin)

response_curve_species_mean_ALL_PB_NN <- curve_data_ALL_PB_NN %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'Greens') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())

plot(response_curve_species_mean_ALL_PB_NN)

ggsave("figures/RC_ALL_PB_NN.png", response_curve_species_mean_ALL_PB_NN, height = 8, width = 12, dpi = 600)



# Mean response curves, FISH NON-NATIVES ----------------------------------

# Aatl

sps_choice <- sps_files[1]
sps <- sps_names[1]


biomod_model_Aatl <- loadRData(paste0(sps_names[1],"/",sps_names[1],".",sps_names[1],".models.out"))
sp_name_Aatl <- BIOMOD_LoadModels(biomod_model_Aatl, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Aatl <- response.plot2(models  = sp_name_Aatl,
                              Data = get_formal_data(biomod_model_Aatl, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Aatl, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Aatl, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Aatl <- function(biomod_model_Aatl, sp_dat_Aatl){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Aatl))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Aatl %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Lt.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                 expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Aatl <- sps_response_data_Aatl(biomod_model_Aatl, sp_dat_Aatl)


# Plot RC for Aatl models

response_curve_species_mean_Aatl <- curve_data_Aatl %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'Spectral') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Aatl)

ggsave("figures/RC_Aatl.png", response_curve_species_mean_Aatl, height = 8, width = 12, dpi = 600)

rm(list = grep("Aatl_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Aatl",names(.GlobalEnv),value = TRUE))

# Aunc

sps_choice <- sps_files[4]
sps <- sps_names[4]

# Load chosen biomod_model and print evaluation scores

biomod_model_Aunc <- loadRData(paste0(sps_names[4],"/",sps_names[4],".",sps_names[4],".models.out"))
sp_name_Aunc <- BIOMOD_LoadModels(biomod_model_Aunc, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Aunc <- response.plot2(models  = sp_name_Aunc,
                              Data = get_formal_data(biomod_model_Aunc, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Aunc, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Aunc, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Aunc <- function(biomod_model_Aunc, sp_dat_Aunc){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Aunc))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Aunc %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Surface.Temperature.Mean" ~ "__A)__    Temperature (°C)",
                                  expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                  expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                  expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)",
                                  expl.name == "Surface.Nitrate.Mean" ~ "__E)__    Nitrate (umol.m-3)"),
                                 run_PA_model = paste0(run,"_",PA,"_",model))
           return(sp_res)
}

curve_data_Aunc <- sps_response_data_Aunc(biomod_model_Aunc, sp_dat_Aunc)


# Plot RC for Aunc models

response_curve_species_mean_Aunc <- curve_data_Aunc %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'Spectral') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Aunc)

ggsave("figures/RC_Aunc.png", response_curve_species_mean_Aunc, height = 8, width = 12, dpi = 600)

rm(list = grep("Aunc_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Aunc",names(.GlobalEnv),value = TRUE))

# All Fish Non-Native together now
curve_data_ALL_F_NN <- rbind(curve_data_Aatl, curve_data_Aunc)

response_curve_species_mean_ALL_F_NN <- curve_data_ALL_F_NN %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'Oranges') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())

plot(response_curve_species_mean_ALL_F_NN)

ggsave("figures/RC_ALL_F_NN.png", response_curve_species_mean_ALL_F_NN, height = 8, width = 12, dpi = 600)



# Mean response curves ZOOBENTHOS NON-NATIVE ------------------------------

# Gtig

sps_choice <- sps_files[10]
sps <- sps_names[10]


biomod_model_Gtig <- loadRData(paste0(sps_names[10],"/",sps_names[10],".",sps_names[10],".models.out"))
sp_name_Gtig <- BIOMOD_LoadModels(biomod_model_Gtig, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Gtig <- response.plot2(models  = sp_name_Gtig,
                              Data = get_formal_data(biomod_model_Gtig, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Gtig, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Gtig, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Gtig <- function(biomod_model_Gtig, sp_dat_Gtig){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Gtig))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Gtig %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Gtig <- sps_response_data_Gtig(biomod_model_Gtig, sp_dat_Gtig)


# Plot RC for Gtig models

response_curve_species_mean_Gtig <- curve_data_Gtig %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Gtig)

ggsave("figures/RC_Gtig.png", response_curve_species_mean_Gtig, height = 8, width = 12, dpi = 600)

rm(list = grep("Gtig_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Gtig",names(.GlobalEnv),value = TRUE))

# Llit

sps_choice <- sps_files[12]
sps <- sps_names[12]

# Load chosen biomod_model and print evaluation scores

biomod_model_Llit <- loadRData(paste0(sps_names[12],"/",sps_names[12],".",sps_names[12],".models.out"))
sp_name_Llit <- BIOMOD_LoadModels(biomod_model_Llit, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Llit <- response.plot2(models  = sp_name_Llit,
                              Data = get_formal_data(biomod_model_Llit, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Llit, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Llit, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Llit <- function(biomod_model_Llit, sp_dat_Llit){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Llit))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Llit %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                  expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                  expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                  expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
                                 run_PA_model = paste0(run,"_",PA,"_",model))
           return(sp_res)
}

curve_data_Llit <- sps_response_data_Llit(biomod_model_Llit, sp_dat_Llit)


# Plot RC for Llit models

response_curve_species_mean_Llit <- curve_data_Llit %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Llit)

ggsave("figures/RC_Llit.png", response_curve_species_mean_Llit, height = 8, width = 12, dpi = 600)

rm(list = grep("Llit_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Llit",names(.GlobalEnv),value = TRUE))

# Mare

sps_choice <- sps_files[15]
sps <- sps_names[15]

# Load chosen biomod_model and print evaluation scores

biomod_model_Mare <- loadRData(paste0(sps_names[15],"/",sps_names[15],".",sps_names[15],".models.out"))
sp_name_Mare <- BIOMOD_LoadModels(biomod_model_Mare, models = c('MARS', 'GLM', 'ANN', 'RF', 'GAM'))
sp_dat_Mare <- response.plot2(models  = sp_name_Mare,
                              Data = get_formal_data(biomod_model_Mare, 'expl.var'),
                              show.variables = get_formal_data(biomod_model_Mare, 'expl.var.names'),
                              data_species = get_formal_data(biomod_model_Mare, 'resp.var'),
                              do.bivariate = FALSE, fixed.var.metric = 'median', legend = FALSE, plot = FALSE)

# Function for extracting species response curve data

sps_response_data_Mare <- function(biomod_model_Mare, sp_dat_Mare){
  # Model scores
  model_scores <- t(data.frame(get_evaluations(biomod_model_Mare))) %>% 
    data.frame() %>% 
    mutate(idx = rownames(.),
           idx = str_replace(idx, "Testing.data", "Testing_data")) %>% 
    separate(idx, c("type", "model", "run", "PA"), sep = "[.]", remove = T) %>% 
    filter(type == "Testing_data") %>%
    `row.names<-`(NULL) %>% 
    dplyr::select(type, model, run, PA, everything())
  
  # Clean up and exit
  sp_res <- sp_dat_Mare %>% 
    separate(pred.name, c("species", "PA", "run", "model"), sep = "_", remove = T) %>% 
    left_join(model_scores, by = c("model", "run", "PA")) %>% 
    mutate(expl.name = case_when(expl.name == "Benthic.Temperature.max" ~ "__A)__    Temperature (°C)",
                                 expl.name == "Surface.Salinity.Mean" ~ "__B)__    Salinity (PSS)",
                                 expl.name == "Surface.Ice.thickness.Mean" ~ "__C)__    Ice thickness (m)",
                                 expl.name == "Surface.Chlorophyll.Mean" ~ "__D)__    Chlorophyll a (mg.m-3)"),
           run_PA_model = paste0(run,"_",PA,"_",model))
  return(sp_res)
}

curve_data_Mare <- sps_response_data_Mare(biomod_model_Mare, sp_dat_Mare)


# Plot RC for Mare models

response_curve_species_mean_Mare <- curve_data_Mare %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 3) +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())
plot(response_curve_species_mean_Mare)

ggsave("figures/RC_Mare.png", response_curve_species_mean_Mare, height = 8, width = 12, dpi = 600)

rm(list = grep("Mare_",names(.GlobalEnv),value = TRUE)); gc()
#rm(list = grep("_Mare",names(.GlobalEnv),value = TRUE))

# All Zoobenthos Non-Native together now
curve_data_ALL_ZB_NN <- rbind(curve_data_Gtig, curve_data_Llit, curve_data_Mare)

response_curve_species_mean_ALL_ZB_NN <- curve_data_ALL_ZB_NN %>%
  filter(TSS >= 0.7) %>% 
  ggplot(aes(x = expl.val, y = pred.val, colour = species)) +
  geom_point(alpha = 0.05) +
  geom_smooth(fill = "grey30") +
  facet_wrap(~expl.name, scales = 'free_x') +#
  labs(x = NULL, y = 'Probability of occurence', colour = 'Species') + 
  scale_color_brewer(type = 'qual', palette = 'Blues') +
  coord_cartesian(expand = F, ylim = c(0, 1)) +
  theme(legend.position = 'bottom',
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = ggtext::element_markdown(hjust = 0),
        strip.background = element_blank())

plot(response_curve_species_mean_ALL_ZB_NN)

ggsave("figures/RC_ALL_ZB_NN.png", response_curve_species_mean_ALL_ZB_NN, height = 8, width = 12, dpi = 600)
