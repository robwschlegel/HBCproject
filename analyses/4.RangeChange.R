# Species Range Change ----------------------------------------------------



## compute Species Range Change (SRC)
## load binary projections

sps_choice_bin_proj_current <- stack("Gset/proj_presentEA/proj_presentEA_Gset_ensemble_TSSbin.img") 

sps_choice_bin_proj_2050 <- stack("Gset/proj_EA2050/proj_EA2050_Gset_ensemble_TSSbin.img")

sps_choice_bin_proj_2100 <- stack("Gset/proj_EA2100/proj_EA2100_Gset_ensemble_TSSbin.img")


## SRC current -> 2050
SRC_current_2050 <- 
  BIOMOD_RangeSize(
    sps_choice_bin_proj_current,
    sps_choice_bin_proj_2050
  )

SRC_current_2050$Compt.By.Models

## SRC current -> 2100
SRC_current_2100 <- 
  BIOMOD_RangeSize(
    sps_choice_bin_proj_current,
    sps_choice_bin_proj_2100
  )

SRC_current_2100$Compt.By.Models

Sps_src_map <- 
  stack(
    SRC_current_2050$Diff.By.Pixel, 
    SRC_current_2100$Diff.By.Pixel
  )
names(Gtig_src_map) 

my.at <- seq(-2.5, 1.5, 1)
myColorkey <- 
  list(
    at = my.at, ## where the colors change
    labels = 
      list(
        labels = c("lost", "pres", "abs","gain"), ## labels
        at = my.at[-1] - 0.5 ## where to print labels
      )
  )

rasterVis::levelplot( 
  Gtig_src_map, 
  main = "Gtig range change",
  colorkey = myColorkey,
  col.regions=c('#f03b20', '#99d8c9', '#f0f0f0', '#2ca25f'),
  layout = c(2,2)
)