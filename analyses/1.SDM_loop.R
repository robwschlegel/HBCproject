# analyses/1.SDM_loop.R


# Setup -------------------------------------------------------------------

# Libraries
#.libPaths(c("~/R-packages", .libPaths()))
library(biomod2)
library(raster)
library(sp)
library(FNN)
library(doParallel)
library(ggplot2)
library(usdm)
library(corrplot)
library(gridExtra)
library(dplyr)
library(FNN)
library(tidyverse)
#library(doParallel); registerDoParallel(cores = 4)

# Custom functions
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
} ## To save the scripts later



# Load Environmental Data -------------------------------------------------


# The environmental file pathways
var_filesNH <- dir("data/Env_var/loop/presentBO2.1/", full.names = T)
var_filesEA <- dir("data/Env_var/loop/presentBO2.1EA/", full.names = T)
var_2050_files <- dir("data/Env_var/loop/2050BO2.1/", full.names = T)
var_2100_files <- dir("data/Env_var/loop/2100BO2.1/", full.names = T)


# Occurrence points
sps_files <- dir("data/rarefied_occ/loop/", full.names = T)
sps_names <- str_remove(dir("data/rarefied_occ/loop/", full.names = F), pattern = "_rarefied_points.csv")

# Global coords from Jesi's data
global_coords <- as.data.frame(sp::read.asciigrid(var_filesNH[1]), xy = T)
global_coords$env_index <- 1:nrow(global_coords)


# The best variables per species
top_var <- read_csv("data/top_var.csv") %>% 
  dplyr::select(Code:var8) %>% 
  pivot_longer(cols = var1:var8) %>% 
  dplyr::select(-name) %>% 
  na.omit()

# Choose a species
# Used for testing
sps_choice <- sps_files[9]
sps <- sps_names[9]


# The full pipeline wrapped into a function JG: TO ACTIVATE ONCE IT WORKS
biomod_pipeline <- function(sps_choice){
  
  print(paste0("Began run on ",sps_choice))
  
  # 2: Load data ------------------------------------------------------------
  
  # Load the species
  sps <- read_csv(sps_choice) %>% 
    mutate(env_index = as.vector(knnx.index(as.matrix(global_coords[,c("s1", "s2")]),
                                            as.matrix(.[,6:7]), k = 1))) %>%
    left_join(global_coords, by = "env_index") %>% 
    dplyr::select(Sp, NEAR_X, NEAR_Y) %>%
    dplyr::rename(lon = NEAR_X, lat = NEAR_Y)
  
  # The species abbreviation
  sps_name <- sps$Sp[1]
  
  # Filter out the top variables
  top_var_sub <- top_var %>% 
    filter(Code == sps_name) %>% 
    mutate(value = paste0(value,".asc"))
  
  # Load the top variables for the species
  expl <- raster::stack(var_filesNH[which(sapply(str_split(var_filesNH, "/"), "[[", 5) %in% top_var_sub$value)])
  expl_EA <- raster::stack(var_filesEA[which(sapply(str_split(var_filesEA, "/"), "[[", 5) %in% top_var_sub$value)])
  expl_2050 <- raster::stack(var_2050_files[which(sapply(str_split(var_2050_files, "/"), "[[", 5) %in% top_var_sub$value)])
  expl_2100 <- raster::stack(var_2100_files[which(sapply(str_split(var_2100_files, "/"), "[[", 5) %in% top_var_sub$value)])
  
  
  # Crop layers to make test in smaller regions
  b <- as(extent(-90, -60, 50, 65), 'SpatialPolygons')
  crs(b) <- crs(expl)
  expl_baby <- crop(expl, b)
  expl_baby <- stack(expl_baby)
  
  # Set temp folder save locations
  dir.create(file.path(sps_name), showWarnings = FALSE)
  dir.create(file.path(sps_name,"/Temp"), showWarnings = FALSE)
  rasterOptions(tmpdir = paste0(sps_name,"/Temp"))
  
  
  # 3: Prep data ------------------------------------------------------------
  
  
  ##Formating Data
  
  biomod.data <- BIOMOD_FormatingData(
    resp.var= rep(1, nrow(sps)), 
    resp.xy = as.matrix(sps[,2:3]), 
    expl.var= expl_baby, #expl, 
    resp.name = sps_name, 
    PA.strategy = 'random',
    PA.nb.rep=1, #3,
    PA.nb.absences=5000
  )
  # biomod.data
  # plot(biomod.data)
  # biomod_data <- readRDS(paste0(sps_name,"/",sps_name,".base.Rds"))
  
  # Save the pre-model data for possible later use
  #saveRDS(biomod_data, file = paste0(sps_name,"/",sps_name,".base.Rds"))
  
  
  # 4: BIOMOD model runs -------------------------------------------------------
  
  
  # Model options
  #Print_Default_ModelingOptions() ##To see default options of each model and change if wanted
  biomod_option <- BIOMOD_ModelingOptions()
  
  
  # Run the model
  biomod_model <- BIOMOD_Modeling(
    biomod.data,
    models = c('RF', 'GLM'), # 'GAM', 'ANN', 'MARS'), ##not maxent because it uses pseudo-absence and final models are suitability and not probability. See Sillero and Barbosa 2021
    models.options = biomod_option,
    NbRunEval = 2, #5, 
    DataSplit = 70, 
    VarImport = 3, # Number of permutations to estimate variable importance
    models.eval.meth = c('TSS', 'ROC'), #'KAPPA', 'ACCURACY'),
    rescal.all.models = FALSE,
    do.full.models = FALSE,
    modeling.id = sps_name)
  #biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
  
  # Build the ensemble models
  biomod_ensemble <- BIOMOD_EnsembleModeling(
    modeling.output = biomod_model,
    chosen.models = 'all',  # defines models kept (useful for removing non-preferred models)
    em.by = 'all',
    eval.metric = c('TSS'),
    #eval.metric.quality.threshold = c(0.7), 
    models.eval.meth = c('TSS', 'ROC'), #'KAPPA', 'ACCURACY'),
    prob.mean = TRUE, # Mean probabilities across predictions
    prob.cv = F, #TRUE, # Coefficient of variation across predictions. Extent to which predictions agree (or diverge) between models
    prob.ci = F, # Confidence interval around prob.mean
    prob.ci.alpha = 0.05,
    VarImport = 5)
  
  
  # 5: BIOMOD Projections Eastern Arctic ---------------------------------------
  
  
  
  #Create projections in EA in the present
  biomod_projectionEA <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl_baby, #expl_EA 
    proj.name = 'presentEA',
    selected.models = 'all',
    binary.meth = 'TSS',
    output.format = '.img', ##if saved as img can then be read in ArcGIS (not if format is grd) #can be saved as .RData
    compress = "xz",
    build.clamping.mask = FALSE,
    do.stack = FALSE
  )
  
  #biomod_projectionEA <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".projection.out")
  
  # Clean out some space
  #rm(biomod_projectionEA); gc()
  
  # Create ensemble projections  
  
  biomod_ensemble_projectionEA <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projectionEA,
    binary.meth = 'TSS',
    output.format = '.img',
    do.stack = T)
  
  # Visualise
  plot(biomod_ensemble_projectionEA)
  
  
  # 6: BIOMOD Projections and forecasting in the future (2050 and 2100)-----------------------------------------------
  
  
  # Create projections in EA in 2050
  biomod_projectionEA2050 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl_2050, 
    proj.name = 'EA2050',
    selected.models = 'all',
    binary.meth = 'TSS',
    output.format = '.img', ##if saved as img can then be read in ArcGIS (not if format is grd) #can be saved as .RData
    compress = "xz",
    build.clamping.mask = FALSE,
    do.stack = FALSE
  )
  
  
  
  #plot(biomod_projection2050) ##It will plot each rep of each model run
  
  # Create ensemble projections   
  
  biomod_ensemble_projectionEA2050 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projectionEA2050,
    binary.meth = 'TSS',
    output.format = '.img',
    do.stack = T)
  
  # Visualise
  plot(biomod_ensemble_projectionEA2050)
  
  # Clean out 2050
  #rm(biomod_projectionEA2050, biomod_ensemble_projectionEA2050); gc()
  
  # Create projections in EA in 2100
  biomod_projectionEA2100 <- BIOMOD_Projection(
    modeling.output = biomod_model,
    new.env = expl_2100, 
    proj.name = 'EA2100',
    selected.models = 'all',
    binary.meth = 'TSS',
    output.format = '.img', ##if saved as img can then be read in ArcGIS (not if format is grd) #can be saved as .RData
    compress = "xz",
    build.clamping.mask = FALSE,
    do.stack = FALSE
  )
  
  #plot(biomod_projection) ##It will plot each rep of each model run
  
  # Create ensemble projections 2100  
  
  biomod_ensemble_projectionEA2100 <- BIOMOD_EnsembleForecasting(
    EM.output = biomod_ensemble,
    projection.output = biomod_projectionEA2100,
    binary.meth = 'TSS',
    output.format = '.img',
    do.stack = T)
  
  # Visualise
  plot(biomod_ensemble_projectionEA2100)
  
  # Clean out 2100
  #rm(biomod_projectionEA2100, biomod_ensemble_projectionEA2100); gc()
  
  # Flush local tmp drive. Better not to do this if running on mulitple cores
  unlink(paste0(normalizePath(tempdir()), "/", dir(tempdir())), recursive = TRUE)
  
  # Unlink Temp folder
  unlink(paste0(sps_name,"/Temp"), recursive = TRUE)
}

# 7: Run the pipeline -----------------------------------------------------

# Detect available cores automagically and set accordingly
# registerDoParallel(cores = detectCores()-1)

# Run one
registerDoParallel(cores = 1)
biomod_pipeline(sps_files[9])

#Run in sequential mode
lapply(sps_files, biomod_pipeline)

# Run them all
registerDoParallel(cores = 3)
plyr::l_ply(sps_files, biomod_pipeline, .parallel = TRUE)

if(require(snowfall)) {
  sfInit(parallel = TRUE, cpus = 3)
  sfExportAll()
  sfLibrary(biomod2)
  sf_out <- sfLapply(sps_files, biomod_pipeline)
  sfStop()
} else {
  for (Sp in sps_files){
    biomod_pipeline(Sp)
  }
}

# Load data/model results -------------------------------------------------

# biomod_model <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".models.out"))
# biomod_ensemble <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"ensemble.models.out"))
# biomod_projectionEA <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,".presentEA.projection.out")
# biomod_projectionEA2050 <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"EA2050.projection.out")
# biomod_projectionEA2100 <- loadRData(paste0(sps_name,"/",sps_name,".",sps_name,"EA2100.projection.out")


# 8: Evaluation scores -------------------------------------------------------

# Choose a species for the following code
sps_choice <- sps_names[9]

# Load chosen biomod_model and print evaluation scores
biomod_model <- loadRData(paste0(sps_choice,"/",sps_choice,".",sps_choice,".models.out"))
biomod_model # print summary

