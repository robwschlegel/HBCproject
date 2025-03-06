# analyses/1.SDM_loop.R


# Setup -------------------------------------------------------------------

# Libraries
library(biomod2)
library(ggplot2)
library(ggtext)
library(stringr)
library(ggtext)
library(raster)
library(FNN)
library(sdmpredictors)
library(sf)
library(sp)
library(tidyverse)

# Custom functions
## To save the scripts later
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
} 


# Load data (individual species) ---------------------------------------------------------------


##Occurrence points

sps_files <- dir("data/rarefied_occ/loop/clean_NoCurrVel/", full.names = T)
sps_names <- str_remove(dir("data/rarefied_occ/loop/clean_NoCurrVel/", full.names = F), pattern = "_rarefied_points.csv")

# Choose a species
# Used for testing
sps_choice <- sps_files[9]
sps <- sps_names[9]


# Load biomod_models ------------------------------------------------------

# Load chosen biomod_model and print evaluation scores
biomod_model <- loadRData(paste0(sps_names[9],"/",sps_names[9],".",sps_names[9],".models.out"))
biomod_model # print summary

#Get data
get_formal_data(biomod_model)

# Evaluation scores (individual species) -------------------------------------------------------


(Model_scores <- get_evaluations(biomod_model)) # get evaluation scores
apply(Model_scores, c(1,2,3), mean, na.rm = T)
 #dim(Model_scores)
 #dimnames(Model_scores)

#To extract only TSS values
scores_TSS <-
  as.numeric(Model_scores["TSS","Testing.data",,,])
scores_TSS

# Model evaluation by algorithm
models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS'), 
                    xlim = c(0.6,1), ylim = c(0.6,1)) + ggtitle("Algorithm") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

# Model evaluation by cross-validation
models_scores_graph(biomod_model, by = "cv_run", metrics = c('ROC','TSS'), 
                    xlim = c(0.6,1), ylim = c(0.6,1)) + ggtitle("Run") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

# Model evaluation by dataset
models_scores_graph(biomod_model, by = "data_set", metrics = c('ROC','TSS'), 
                    xlim = c(0.6,1), ylim = c(0.6,1)) + ggtitle("PA") +
  geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)

## Calculate mean of variable importance by algorithm
# JG: I have read that scores reported are raw in the table (to be easier to interpret, 
# it should be normalized on our own - sum to 1 across algorithms)
biomod_model <- loadRData(paste0(sps_names[9],"/",sps_names[9],".",sps_names[9],".models.out"))
biomod_model
(models_var_import <- get_variables_importance(biomod_model))
apply(models_var_import, c(1,2), mean, na.rm = T)
apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean) # Overall mean per variable



# Load ensemble.models ----------------------------------------------------


# Load if the model has already been run
sps <- sps_names[9]
biomod_ensemble <- loadRData(paste0(sps,"/",sps,".",sps,"ensemble.models.out"))
biomod_ensemble

# Evaluation scores ensembles(individual species) ----------------------------------------------



(models_scores_biomod_ensemble <- get_evaluations(biomod_ensemble))
get_variables_importance(biomod_ensemble)
(models_var_import <- get_variables_importance(biomod_ensemble))
apply(models_var_import, c(1,2), mean, na.rm = T)
apply(apply(models_var_import, c(1,2), mean, na.rm = T), 1, mean)



# GET MODELS' EVALUATION SCORES (all species in a loop) -------------------------------------------------

##From script 2_results in ArcticInvasion

.libPaths(c("~/R-packages", .libPaths()))
library(tidyverse)
library(ggpubr)
library(biomod2)
library(sp)
library(dtplyr)



biomod_res_table <- function(sps){
  
  # The full model results
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  model_scores <- models_scores_graph(biomod_model, plot = F,
                                      by = "models", metrics = c('ROC','TSS'), 
                                      xlim = c(0.6,1), ylim = c(0.6,1)) + ggtitle("Algorithm") +
    geom_hline(aes(yintercept = 0.7), colour = "red", size = 2)
   #ggsave(filename = paste0("results/model_scores/",sps,"_scores.png"), plot = model_scores)
  
  # Get the TSS scores, cutoffs for binary presence/absence, and specificity/sensitivity
  biomod_cutoff <- plyr::adply(get_evaluations(biomod_model), c(1,3,4,5)) %>% 
    dplyr::rename(test = X1, model = X2, run = X3, PA = X4, score = Testing.data) %>% 
    mutate(sps = sps) %>% 
    dplyr::select(sps, everything())
  
  return(biomod_cutoff)
}


# Run it all
all_res_table <- plyr::ldply(sps_names, biomod_res_table)
write_csv(all_res_table, "results/all_res_table_GsetSSTmax.csv")

# Create a table that shows the mean results
mean_res_table <- all_res_table %>%
  group_by(sps, test, model) %>%
  summarise_if(is.numeric, mean) %>%
  ungroup() %>%
  mutate_if(is.numeric, round, 2) %>%
  mutate(Cutoff = round(Cutoff))
write_csv(mean_res_table, "results/mean_res_table_GsetSSTmax.csv")


# VARIABLES USED ----------------------------------------------------------

# The wrapper functions
biomod_var_table <- function(sps){
  
  # The full model results
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Get list of variables used
  biomod_var <- as.data.frame(biomod_model@expl.var.names) %>% 
    `colnames<-`(c("var")) %>%
    mutate(sps = sps,
           var_count = 1:n()) %>% #sometimes loading plyr on top of dplyr can make n() or group_by not to work
    pivot_wider(names_from = var_count, values_from = var, names_prefix = "Var")
  return(biomod_var)
}

# Run it all
all_var_table <- plyr::ldply(sps_names, biomod_var_table)

# Save
write_csv(all_var_table, "results/all_var_table_GsetSSTmax.csv")


# COMPARE MODELS ----------------------------------------------------------
##I think this do not apply to ecological interactions paper and is not needed 
# 3: Create multi-model comparisons ---------------------------------------------
# library(raster)
# 
# # Function that loads an .Rds file and rounds it to the nearest 0.25 degree resolution
# readRDS_0.25 <- function(file_name, projection_name){
#   df <- readRDS(file_name) %>% 
#     na.omit() %>% 
#     mutate(x = plyr::round_any(x, 0.25), 
#            y = plyr::round_any(y, 0.25),
#            projection = projection_name) %>% 
#     group_by(model, projection, x, y) %>% 
#     summarise(z = round(mean(z, na.rm = T))) %>% 
#     ungroup()
# }
# 
# # Convenience function to process a raster into a long dataframe
# raster_to_long <- function(raster_file, model_name, projection_name){
#   res <- as.data.frame(raster(raster_file), xy = TRUE) %>%
#     na.omit() %>% 
#     `colnames<-`(c("x", "y", "z")) %>% 
#     mutate(x = plyr::round_any(x, 0.25), 
#            y = plyr::round_any(y, 0.25),
#            model = model_name,
#            projection = projection_name) %>% 
#     group_by(model, projection, x, y) %>% 
#     summarise(z = round(mean(z, na.rm = TRUE))) %>% 
#     ungroup()
#   return(res)
# }
# 
# # Load depth for result screening
# depth_long <- raster_to_long("data/present/depthclip.asc", "depth", "present") %>% 
#   dplyr::rename(depth = z) %>% 
#   dplyr::select(x, y, depth)
# 
# # Convenience function for multi-plotting
# single_plot <- function(df){
#   comp_multi_fig <- ggplot(data = df, aes(x = x, y = y)) +
#     geom_polygon(data = map_base, fill = "grey70",
#                  aes(x = lon, y = lat, group = group)) +
#     geom_tile(aes(fill = as.factor(z))) +
#     labs(x = NULL, y = NULL, 
#          fill = "Predicted habitat suitability") +
#     scale_fill_manual(values = c("red4")) +
#     facet_wrap(~model) +
#     coord_quickmap(expand = F) +
#     theme(legend.position = "bottom", 
#           title = element_text(size = 16),
#           strip.text = element_text(size = 12),
#           legend.text = element_blank(), 
#           axis.text = element_blank(),
#           axis.ticks = element_blank(), 
#           panel.grid = element_blank())
#   return(comp_multi_fig)
# }
# 
# comp_multi_plot <- function(df, sps){
#   
#   # Find species depth limit
#   sps_depth <- filter(sps_depths, Sps == sps)
#   
#   # Load species presence data
#   sps_data <- read_csv(sps_files[which(sps == sps_names)]) %>% 
#     `colnames<-`(c("Sps", "lon", "lat")) %>% 
#     mutate(lon =  plyr::round_any(lon, 0.25),
#            lat =  plyr::round_any(lat, 0.25))
#   
#   # Prep data for plotting
#   df_sub <- df %>% 
#     left_join(depth_long, by = c("x", "y")) %>% 
#     mutate(z = ifelse(depth > sps_depth$Depth, 0, z),
#            model = factor(model, levels = unique(model))) %>% 
#     filter(z != 0)
#   
#   # Create each panel
#   plot_ANN <- single_plot(filter(df_sub, model == "ANN")) +
#     geom_point(data = sps_data, aes(x = lon, y = lat), 
#                fill = "white", shape = 21, size = 1, stroke = 0.5) #+
#   # labs(title = sps_depth$Species)
#   plot_GAM <- single_plot(filter(df_sub, model == "GAM"))
#   plot_GLM <- single_plot(filter(df_sub, model == "GLM"))
#   plot_RF <- single_plot(filter(df_sub, model == "RF"))
#   plot_MaxEnt <- single_plot(filter(df_sub, model == "MaxEnt"))
#   plot_Ensemble <- single_plot(filter(df_sub, model == "Ensemble")) +
#     theme(panel.background = element_rect(colour = "black", size = 2))
#   
#   # Steek'em
#   plot_multi <- ggarrange(plot_ANN, plot_GAM, plot_GLM, plot_RF, plot_MaxEnt, plot_Ensemble, 
#                           ncol = 3, nrow = 2, common.legend = T, legend = "bottom", align = "hv")
#   plot_multi_title <- annotate_figure(p = plot_multi, fig.lab = df_sub$projection[1], fig.lab.face = "bold", fig.lab.size = 14,
#                                       top = text_grob(sps_depth$Species, color = "black", face = "italic", size = 16))
#   
#   # Save
#   ggsave(plot = plot_multi_title, width = 21, height = 8,
#          filename = paste0("graph/comparison_multi/",sps,"_",df$projection[1],".png"))
# }
# 
# # Wraper to run visuals for all species
# biomod_multi_visuals <- function(sps){
#   
#   # Create present figure
#   df_present <- rbind(readRDS_0.25(paste0("data/biomod/",sps,"_df_present.Rds"), "Present"),
#                       raster_to_long(paste0("data/maxent/",sps,"_avg_binary.tif"), "MaxEnt", "Present"),
#                       raster_to_long(paste0(sps,"/proj_present/proj_present_",sps,"_TSSbin.gri"), "Ensemble", "Present"))
#   comp_multi_present <- comp_multi_plot(df_present, sps)
#   rm(df_present, comp_multi_present); gc()
#   
#   # Create 2050 figure
#   df_2050 <- rbind(readRDS_0.25(paste0("data/biomod/",sps,"_df_2050.Rds"), "2050"),
#                    raster_to_long(paste0("data/maxent/",sps,"_2050_45_avg_binary.tif"), "MaxEnt", "2050"),
#                    raster_to_long(paste0(sps,"/proj_2050/proj_2050_",sps,"_TSSbin.gri"), "Ensemble", "2050"))
#   comp_multi_2050 <- comp_multi_plot(df_2050, sps)
#   rm(df_2050, comp_multi_2050); gc()
#   
#   # Create 2100 figure
#   df_2100 <- rbind(readRDS_0.25(paste0("data/biomod/",sps,"_df_2100.Rds"), "2100"),
#                    raster_to_long(paste0("data/maxent/",sps,"_2100_45_avg_binary.tif"), "MaxEnt", "2100"),
#                    raster_to_long(paste0(sps,"/proj_2100/proj_2100_",sps,"_TSSbin.gri"), "Ensemble", "2100"))
#   comp_multi_2100 <- comp_multi_plot(df_2100, sps)
#   rm(df_2100, comp_multi_2100); gc()
# }
# 
# # Run one
# # system.time(biomod_multi_visuals(sps_names[1])) # 223 seconds
# 
# # Run them all
# plyr::l_ply(sps_names, biomod_multi_visuals, .parallel = F)


# Boxplot_combined Evaluation Metrics---------------------------------------------------------------
#From ArcticKelp project with Rob

# The model values from the ensembles

#ZOOBENTHOS, NATIVES

model_stats_plot <- function(sps){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Create full species name
  if(sps == "Gset"){
    sps_title <- "Gammarus setosus"
  } else if(sps == "Goce"){
     sps_title <- "Gammarus oceanicus"
  # } else if(sps == "Lsax"){
  #   sps_title <- "Littorina saxatilis"
  # } else if(sps == "Mtru"){
  #   sps_title <- "Mya truncata"
  # } else if(sps == "Sdro"){
  #   sps_title <- "Strongylocentrotus droebachiensis"
  }
  
  # Extract raw results for geom_point()
  scores <- biomod2:::get_evaluations(biomod_model, as.data.frame = T) %>% 
    separate(Model.name, into = c("model", "run", "PA"), sep = "_") %>% 
    dplyr::select(model:Testing.data) %>% 
    pivot_wider(names_from = Eval.metric, values_from = Testing.data) 
  
  # Index of possible combos to catch 0 model psases
  scores_index <- data.frame(model = c("ANN", "GAM", "GLM", "MARS", "RF"),
                             total_count = 15)
  
  # Count of models that passed
  scores_passed <- scores %>% 
    group_by(model) %>% 
    mutate(total_count = n()) %>% 
    filter(TSS >= 0.7) %>% 
    group_by(model, total_count) %>% 
    summarise(pass_count = n(), .groups = "drop") %>% 
    right_join(scores_index, by = c("model", "total_count")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(pass_label = paste0(pass_count,"/",total_count)) %>% 
    dplyr::arrange(model)
  
  my_colors <- RColorBrewer::brewer.pal(9, 'Blues') [4:8]
  
  
  # Model evaluation by algorithm
  #model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
  #geom_point(data = scores, aes(x = ROC, y = TSS, colour = model)) +
  model_res <- ggplot(scores) +
    geom_label(data = scores_passed, label.size = 1.5, size = 4, show.legend = F,
               aes(x = model, y = 0.55, label = pass_label, colour = model)) +
    geom_boxplot(aes(x = model, y = TSS, fill = model)) +
    geom_hline(aes(yintercept = 0.7), colour = "black", size = 2) +
    labs(x = NULL, y = "TSS", title = sps_title, fill = "Model") +
    scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
    #scale_fill_brewer(palette = 'Blues', direction = 1, aesthetics = c("colour", "fill")) +
    scale_fill_manual(values = my_colors, aesthetics = c("colour", "fill")) +
    # guides(fill = guide_legend(override.aes = list(shape = 15))) +
     #coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # model_res
  
  return(model_res)
}

# Create plots
model_stats_Gset <- model_stats_plot("Gset")
model_stats_Goce <- model_stats_plot("Goce")
#model_stats_Lsax <- model_stats_plot("Lsax")
#model_stats_Mtru <- model_stats_plot("Mtru")
#model_stats_Sdro <- model_stats_plot("Sdro")

# Combine and save
# EvalMetr_ZB <- ggpubr::ggarrange(model_stats_Gset, model_stats_Goce, model_stats_Lsax, model_stats_Mtru, model_stats_Sdro, 
#                             ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", 
#                             labels = c("A)", "B)", "C)", "D)", "E)"))

EvalMetr_ZB <- ggpubr::ggarrange(model_stats_Gset, model_stats_Goce, 
                                 ncol = 2, nrow = 1, common.legend = TRUE, legend = "bottom", 
                                 labels = c("A)", "B)"))
plot(EvalMetr_ZB)

#ggsave("figures/EvaluationMetrics/EvalMetr_ZB_EcolInter_GsetSSTmax.png", EvalMetr_ZB, width = 12, height = 10.5,  dpi = 600)
#ggsave("figures/EvaluationMetrics/EvalMetr_ZB.jpg", fig_S4, width = 10, height = 10.5,  dpi = 600)

#ZOOBENTHOS, NON-NATIVES

model_stats_plot <- function(sps){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Create full species name
  if(sps == "Gtig"){
    sps_title <- "Gammarus tigrinus"
  # } else if(sps == "Llit"){
  #   sps_title <- "Littorina littorea"
  # } else if(sps == "Mare"){
  #   sps_title <- "Mya arenaria"
  } 
  
  # Extract raw results for geom_point()
  scores <- biomod2:::get_evaluations(biomod_model, as.data.frame = T) %>% 
    separate(Model.name, into = c("model", "run", "PA"), sep = "_") %>% 
    dplyr::select(model:Testing.data) %>% 
    pivot_wider(names_from = Eval.metric, values_from = Testing.data) 
  
  # Index of possible combos to catch 0 model psases
  scores_index <- data.frame(model = c("ANN", "GAM", "GLM", "MARS", "RF"),
                             total_count = 15)
  
  # Count of models that passed
  scores_passed <- scores %>% 
    group_by(model) %>% 
    mutate(total_count = n()) %>% 
    filter(TSS >= 0.7) %>% 
    group_by(model, total_count) %>% 
    summarise(pass_count = n(), .groups = "drop") %>% 
    right_join(scores_index, by = c("model", "total_count")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(pass_label = paste0(pass_count,"/",total_count)) %>% 
    dplyr::arrange(model)
  
  my_colors <- RColorBrewer::brewer.pal(9, 'Blues') [5:9]
  
  
  # Model evaluation by algorithm
  #model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
  #geom_point(data = scores, aes(x = ROC, y = TSS, colour = model)) +
  model_res <- ggplot(scores) +
    geom_label(data = scores_passed, label.size = 1.5, size = 4, show.legend = F,
               aes(x = model, y = 0.55, label = pass_label, colour = model)) +
    geom_boxplot(aes(x = model, y = TSS, fill = model)) +
    geom_hline(aes(yintercept = 0.7), colour = "black", size = 2) +
    labs(x = NULL, y = "TSS", title = sps_title, fill = "Model") +
    scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
    #scale_fill_brewer(palette = 'Blues', direction = 1, aesthetics = c("colour", "fill")) +
    scale_fill_manual(values = my_colors, aesthetics = c("colour", "fill")) +
    # guides(fill = guide_legend(override.aes = list(shape = 15))) +
    #coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # model_res
  
  return(model_res)
}

# Create plots
model_stats_Gtig <- model_stats_plot("Gtig")
# model_stats_Llit <- model_stats_plot("Llit")
# model_stats_Mare <- model_stats_plot("Mare")


# Combine and save
# EvalMetr_ZB_NN <- ggpubr::ggarrange(model_stats_Gtig, model_stats_Llit, model_stats_Mare,  
#                                  ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", 
#                                  labels = c("A)", "B)", "C)"))

EvalMetr_ZB_NN <- ggpubr::ggarrange(model_stats_Gtig,   
                                    ncol = 1, nrow = 1, common.legend = TRUE, legend = "bottom", 
                                    labels = c("A)"))
plot(EvalMetr_ZB_NN)

#ggsave("figures/EvaluationMetrics/EvalMetr_ZB_NN_EcolInt.png", EvalMetr_ZB_NN, width = 12, height = 10.5,  dpi = 600)
#ggsave("figures/EvaluationMetrics/EvalMetr_ZB.jpg", fig_S4, width = 10, height = 10.5,  dpi = 600)


#FISH, NATIVES

model_stats_plot <- function(sps){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Create full species name
  if(sps == "Gtri"){
    sps_title <- "Gymnacanthus tricuspis"
  } else if(sps == "Mdes"){
    sps_title <- "Myoxocephalus scorpioides"
  } else if(sps == "Mqua"){
    sps_title <- "Myoxocephalus quadricornis"
  } else if(sps == "Msco"){
    sps_title <- "Myoxocephalus scorpius"
  }
  
  # Extract raw results for geom_point()
  scores <- biomod2:::get_evaluations(biomod_model, as.data.frame = T) %>% 
    separate(Model.name, into = c("model", "run", "PA"), sep = "_") %>% 
    dplyr::select(model:Testing.data) %>% 
    pivot_wider(names_from = Eval.metric, values_from = Testing.data) 
  
  # Index of possible combos to catch 0 model psases
  scores_index <- data.frame(model = c("ANN", "GAM", "GLM", "MARS", "RF"),
                             total_count = 15)
  
  # Count of models that passed
  scores_passed <- scores %>% 
    group_by(model) %>% 
    mutate(total_count = n()) %>% 
    filter(TSS >= 0.7) %>% 
    group_by(model, total_count) %>% 
    summarise(pass_count = n(), .groups = "drop") %>% 
    right_join(scores_index, by = c("model", "total_count")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(pass_label = paste0(pass_count,"/",total_count)) %>% 
    dplyr::arrange(model)
  
  my_colors <- RColorBrewer::brewer.pal(9, 'Oranges') [4:8]
  
  
  # Model evaluation by algorithm
  #model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
  #geom_point(data = scores, aes(x = ROC, y = TSS, colour = model)) +
  model_res <- ggplot(scores) +
    geom_label(data = scores_passed, label.size = 1.5, size = 4, show.legend = F,
               aes(x = model, y = 0.55, label = pass_label, colour = model)) +
    geom_boxplot(aes(x = model, y = TSS, fill = model)) +
    geom_hline(aes(yintercept = 0.7), colour = "black", size = 2) +
    labs(x = NULL, y = "TSS", title = sps_title, fill = "Model") +
    scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
    #scale_fill_brewer(palette = 'Blues', direction = 1, aesthetics = c("colour", "fill")) +
    scale_fill_manual(values = my_colors, aesthetics = c("colour", "fill")) +
    # guides(fill = guide_legend(override.aes = list(shape = 15))) +
    #coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # model_res
  
  return(model_res)
}

# Create plots
#model_stats_Gacu <- model_stats_plot("Gacu")
model_stats_Gtri <- model_stats_plot("Gtri")
model_stats_Mdes <- model_stats_plot("Mdes")
model_stats_Mqua <- model_stats_plot("Mqua")
model_stats_Msco <- model_stats_plot("Msco")

# Combine and save
EvalMetr_F <- ggpubr::ggarrange(model_stats_Gtri, model_stats_Mdes, model_stats_Mqua, model_stats_Msco, 
                                 ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", 
                                 labels = c("A)", "B)", "C)", "D)", "E)"))
plot(EvalMetr_F)

#ggsave("figures/EvaluationMetrics/EvalMetr_F_noGacu.png", EvalMetr_F, width = 12, height = 10.5,  dpi = 600)
#ggsave("figures/EvaluationMetrics/EvalMetr_F.jpg", EvalMetrF, width = 10, height = 10.5,  dpi = 600)


##FISH, NON-NATIVES

model_stats_plot <- function(sps){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Create full species name
  if(sps == "Aatl"){
    sps_title <- "Artediellus atlanticus"
  } else if(sps == "Aunc"){
    sps_title <- "Artediellus uncinatus"
  } 
  
  # Extract raw results for geom_point()
  scores <- biomod2:::get_evaluations(biomod_model, as.data.frame = T) %>% 
    separate(Model.name, into = c("model", "run", "PA"), sep = "_") %>% 
    dplyr::select(model:Testing.data) %>% 
    pivot_wider(names_from = Eval.metric, values_from = Testing.data) 
  
  # Index of possible combos to catch 0 model psases
  scores_index <- data.frame(model = c("ANN", "GAM", "GLM", "MARS", "RF"),
                             total_count = 15)
  
  # Count of models that passed
  scores_passed <- scores %>% 
    group_by(model) %>% 
    mutate(total_count = n()) %>% 
    filter(TSS >= 0.7) %>% 
    group_by(model, total_count) %>% 
    summarise(pass_count = n(), .groups = "drop") %>% 
    right_join(scores_index, by = c("model", "total_count")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(pass_label = paste0(pass_count,"/",total_count)) %>% 
    dplyr::arrange(model)
  
  my_colors <- RColorBrewer::brewer.pal(9, 'Oranges') [4:8]
  
  
  # Model evaluation by algorithm
  #model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
  #geom_point(data = scores, aes(x = ROC, y = TSS, colour = model)) +
  model_res <- ggplot(scores) +
    geom_label(data = scores_passed, label.size = 1.5, size = 4, show.legend = F,
               aes(x = model, y = 0.55, label = pass_label, colour = model)) +
    geom_boxplot(aes(x = model, y = TSS, fill = model)) +
    geom_hline(aes(yintercept = 0.7), colour = "black", size = 2) +
    labs(x = NULL, y = "TSS", title = sps_title, fill = "Model") +
    scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
    #scale_fill_brewer(palette = 'Blues', direction = 1, aesthetics = c("colour", "fill")) +
    scale_fill_manual(values = my_colors, aesthetics = c("colour", "fill")) +
    # guides(fill = guide_legend(override.aes = list(shape = 15))) +
    #coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # model_res
  
  return(model_res)
}

# Create plots
model_stats_Aatl <- model_stats_plot("Aatl")
model_stats_Aunc <- model_stats_plot("Aunc")


# Combine and save
EvalMetr_F_NN <- ggpubr::ggarrange(model_stats_Aatl, model_stats_Aunc,  
                                ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", 
                                labels = c("A)", "B)"))
plot(EvalMetr_F_NN)

#ggsave("figures/EvaluationMetrics/EvalMetr_F_NN.png", EvalMetr_F_NN, width = 12, height = 10.5,  dpi = 600)
#ggsave("figures/EvaluationMetrics/EvalMetr_F.jpg", EvalMetrF, width = 10, height = 10.5,  dpi = 600)




#PHYTOBENTHOS, NATIVES

model_stats_plot <- function(sps){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Create full species name
  if(sps == "Acla"){
    sps_title <- "Agarum clathratum"
  } else if(sps == "Aesc"){
    sps_title <- "Alaria esculenta"
  } else if(sps == "Lsol"){
    sps_title <- "Laminaria solidungula"
  } else if(sps == "Slat"){
    sps_title <- "Saccharina latissima"
  } else if(sps == "Zmar"){
    sps_title <- "Zoostera marina"
  }
  
  # Extract raw results for geom_point()
  scores <- biomod2:::get_evaluations(biomod_model, as.data.frame = T) %>% 
    separate(Model.name, into = c("model", "run", "PA"), sep = "_") %>% 
    dplyr::select(model:Testing.data) %>% 
    pivot_wider(names_from = Eval.metric, values_from = Testing.data) 
  
  # Index of possible combos to catch 0 model psases
  scores_index <- data.frame(model = c("ANN", "GAM", "GLM", "MARS", "RF"),
                             total_count = 15)
  
  # Count of models that passed
  scores_passed <- scores %>% 
    group_by(model) %>% 
    mutate(total_count = n()) %>% 
    filter(TSS >= 0.7) %>% 
    group_by(model, total_count) %>% 
    summarise(pass_count = n(), .groups = "drop") %>% 
    right_join(scores_index, by = c("model", "total_count")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(pass_label = paste0(pass_count,"/",total_count)) %>% 
    dplyr::arrange(model)
  
  my_colors <- RColorBrewer::brewer.pal(9, 'Greens') [4:8]
  
  
  # Model evaluation by algorithm
  #model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
  #geom_point(data = scores, aes(x = ROC, y = TSS, colour = model)) +
  model_res <- ggplot(scores) +
    geom_label(data = scores_passed, label.size = 1.5, size = 4, show.legend = F, 
               aes(x = model, y = 0.55, label = pass_label, colour = model)) +
    geom_boxplot(aes(x = model, y = TSS, fill = model)) +
    geom_hline(aes(yintercept = 0.7), colour = "black", size = 2) +
    labs(x = NULL, y = "TSS", title = sps_title, fill = "Model") +
    scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
    #scale_fill_brewer(palette = 'Blues', direction = 1, aesthetics = c("colour", "fill")) +
    scale_fill_manual(values = my_colors, aesthetics = c("colour", "fill")) +
    # guides(fill = guide_legend(override.aes = list(shape = 15))) +
    #coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # model_res
  
  return(model_res)
}

# Create plots
model_stats_Acla <- model_stats_plot("Acla")
model_stats_Aesc <- model_stats_plot("Aesc")
model_stats_Lsol <- model_stats_plot("Lsol")
model_stats_Slat <- model_stats_plot("Slat")
model_stats_Zmar <- model_stats_plot("Zmar")

# Combine and save
EvalMetr_PB <- ggpubr::ggarrange(model_stats_Acla, model_stats_Aesc, model_stats_Lsol, model_stats_Slat, model_stats_Zmar, 
                                ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", 
                                labels = c("A)", "B)", "C)", "D)", "E)"))
plot(EvalMetr_PB)

#ggsave("figures/EvaluationMetrics/EvalMetr_PB.png", EvalMetr_PB, width = 12, height = 10.5,  dpi = 600)
#ggsave("figures/EvaluationMetrics/EvalMetr_ZB.jpg", fig_S4, width = 10, height = 10.5,  dpi = 600)


#PHYTOBENTHOS, NON-NATIVES

model_stats_plot <- function(sps){
  
  # Load chosen biomod_model and print evaluation scores
  biomod_model <- loadRData(paste0(sps,"/",sps,".",sps,".models.out"))
  
  # Create full species name
  if(sps == "Cfra"){
    sps_title <- "Codium fragile"
  } else if(sps == "Dcon"){
    sps_title <- "Dumontia contorta"
  } else if(sps == "Smut"){
    sps_title <- "Sargassum muticum"
  } else if(sps == "Upin"){
    sps_title <- "Undaria pinnatifida"
  } 
  
  # Extract raw results for geom_point()
  scores <- biomod2:::get_evaluations(biomod_model, as.data.frame = T) %>% 
    separate(Model.name, into = c("model", "run", "PA"), sep = "_") %>% 
    dplyr::select(model:Testing.data) %>% 
    pivot_wider(names_from = Eval.metric, values_from = Testing.data) 
  
  # Index of possible combos to catch 0 model psases
  scores_index <- data.frame(model = c("ANN", "GAM", "GLM", "MARS", "RF"),
                             total_count = 15)
  
  # Count of models that passed
  scores_passed <- scores %>% 
    group_by(model) %>% 
    mutate(total_count = n()) %>% 
    filter(TSS >= 0.7) %>% 
    group_by(model, total_count) %>% 
    summarise(pass_count = n(), .groups = "drop") %>% 
    right_join(scores_index, by = c("model", "total_count")) %>% 
    replace(is.na(.), 0) %>% 
    mutate(pass_label = paste0(pass_count,"/",total_count)) %>% 
    dplyr::arrange(model)
  
  my_colors <- RColorBrewer::brewer.pal(9, 'Greens') [5:9]
  
  
  # Model evaluation by algorithm
  #model_res <- biomod2::models_scores_graph(biomod_model, by = "models", metrics = c('ROC','TSS')) + 
  #geom_point(data = scores, aes(x = ROC, y = TSS, colour = model)) +
  model_res <- ggplot(scores) +
    geom_label(data = scores_passed, label.size = 1.5, size = 4, show.legend = F,
               aes(x = model, y = 0.55, label = pass_label, colour = model)) +
    geom_boxplot(aes(x = model, y = TSS, fill = model)) +
    geom_hline(aes(yintercept = 0.7), colour = "black", size = 2) +
    labs(x = NULL, y = "TSS", title = sps_title, fill = "Model") +
    scale_y_continuous(limits = c(0.5, 1.0), expand = c(0, 0)) +
    #scale_fill_brewer(palette = 'Blues', direction = 1, aesthetics = c("colour", "fill")) +
    scale_fill_manual(values = my_colors, aesthetics = c("colour", "fill")) +
    # guides(fill = guide_legend(override.aes = list(shape = 15))) +
    #coord_cartesian(xlim = c(0.6,1), ylim = c(0.3,1), expand = F) +
    theme(plot.title = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA))
  # model_res
  
  return(model_res)
}

# Create plots
model_stats_Cfra <- model_stats_plot("Cfra")
model_stats_Dcon <- model_stats_plot("Dcon")
model_stats_Smut <- model_stats_plot("Smut")
model_stats_Upin <- model_stats_plot("Upin")

# Combine and save
EvalMetr_PB_NN <- ggpubr::ggarrange(model_stats_Cfra, model_stats_Dcon, model_stats_Smut, model_stats_Upin,  
                                 ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom", 
                                 labels = c("A)", "B)", "C)", "D)"))
plot(EvalMetr_PB_NN)

ggsave("figures/EvaluationMetrics/EvalMetr_PB_NN.png", EvalMetr_PB_NN, width = 12, height = 10.5,  dpi = 600)
#ggsave("figures/EvaluationMetrics/EvalMetr_ZB.jpg", fig_S4, width = 10, height = 10.5,  dpi = 600)





# Agarum only for demo/talk
ggsave("figures/EvaluationMetrics/agarum.png", model_stats_Acla, width = 7, height = 5)

# Alaria only for talk
ggsave("talk/figure/fig_S4_alaria.png", model_stats_Aesc, width = 7, height = 5)

# Laminaria only for talk
ggsave("talk/figure/fig_S4_laminaria.png", model_stats_Lsol, width = 7, height = 5)

# Saccharina only for talk
ggsave("talk/figure/fig_S4_saccharina.png", model_stats_Slat, width = 7, height = 5)





# # To visualize species' modeled response to the given variables ##JG IT IS VISUALIZED IN THE NEXT STEP
# sp_name_GLM <- BIOMOD_LoadModels(biomod_model, models = 'GLM')
# sp_name_ANN <- BIOMOD_LoadModels(biomod_model, models = 'ANN')
# sp_name_RF <- BIOMOD_LoadModels(biomod_model, models = 'RF')
# sp_name_GAM <- BIOMOD_LoadModels(biomod_model, models = 'GAM')
# sp_name_MARS <- BIOMOD_LoadModels(biomod_model, models = 'MARS') 
# 
# # Evaluate individual models
# MARS_eval_strip <- biomod2::response.plot2(
#   models = sp_name_MARS,
#   Data = get_formal_data(biomod_model, 'expl.var'),
#   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
#   do.bivariate = F,
#   fixed.var.metric = 'mean',
#   legend = F,
#   display_title = F,
#   data_species = get_formal_data(biomod_model, 'resp.var')
# )
# 
# GLM_eval_strip <- biomod2::response.plot2(
#   models = sp_name_GLM,
#   Data = get_formal_data(biomod_model, 'expl.var'),
#   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
#   do.bivariate = F,
#   fixed.var.metric = 'mean',
#   legend = F,
#   display_title = F,
#   data_species = get_formal_data(biomod_model, 'resp.var')
# )
# 
# ANN_eval_strip <- biomod2::response.plot2(
#   models = sp_name_ANN,
#   Data = get_formal_data(biomod_model, 'expl.var'),
#   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
#   do.bivariate = F,
#   fixed.var.metric = 'mean',
#   legend = F,
#   display_title = F,
#   data_species = get_formal_data(biomod_model, 'resp.var')
# )
# 
# RF_eval_strip <- biomod2::response.plot2(
#   models = sp_name_RF,
#   Data = get_formal_data(biomod_model, 'expl.var'),
#   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
#   do.bivariate = F,
#   fixed.var.metric = 'mean',
#   legend = F,
#   display_title = F,
#   data_species = get_formal_data(biomod_model, 'resp.var')
# )
# 
# GAM_eval_strip <- biomod2::response.plot2(
#   models = sp_name_GAM,
#   Data = get_formal_data(biomod_model, 'expl.var'),
#   show.variables = get_formal_data(biomod_model, 'expl.var.names'),
#   do.bivariate = F,
#   fixed.var.metric = 'mean',
#   legend = F,
#   display_title = F,
#   data_species = get_formal_data(biomod_model, 'resp.var')
# )
# 

