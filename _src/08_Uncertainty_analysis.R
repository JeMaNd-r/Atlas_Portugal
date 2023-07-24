#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Extract uncertainty               #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#


#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

# load number of occurrences per species and focal species names
speciesSub <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv")) %>% pull(SpeciesID)
speciesSub

covarsNames <- paste0("PC", 1:11)

#- - - - - - - - - - - - - - - - - - - - -
# load environmental variables
Env_clip <- terra::rast("_intermediates/EnvPredictor_PCA_1km_POR.tif")
Env_clip <- terra::subset(Env_clip, 1:11) #11 = >80%

#Env_clip_df <- terra::as.data.frame(Env_clip, xy=TRUE)

#uncertain_df <- Env_clip_df %>% dplyr::select(x, y)
uncertain_tif <- terra::rast(terra::ext(Env_clip), resolution=res(Env_clip)) 
crs(uncertain_tif) <- crs(Env_clip)
uncertain_tif

for(spID in speciesSub){try({
  
  print(paste0("Species: ", spID))

  # list files in species-specific BIOMOD folder
  temp_files <- list.files(paste0("_results/",Taxon_name, "/", stringr::str_replace(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
          
  #setwd(paste0(data_wd, "/results/", Taxon_name))

  myBiomodEnProj <- get(load(temp_files[stringr::str_detect(temp_files,"ensemble.RData")]))
       
  # save predictions as raster file
  temp_prediction <- myBiomodEnProj[,"pred"] #column with CoV
  temp_prediction <- as.numeric(temp_prediction)
  # add names of grid cell (only for those that have no NA in any layer)
  names(temp_prediction) <- rownames(Env_clip_df)
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_clip_df$x
  temp_prediction$y <- Env_clip_df$y
  temp_prediction <- temp_prediction %>% full_join(Env_clip_df %>% dplyr::select(x,y)) %>%
     rename("layer" = temp_prediction)
  temp_prediction$layer <- temp_prediction$layer / 1000
  temp_prediction[,spID] <- temp_prediction$layer
          
  temp_prediction <- terra::rast(temp_prediction[,c("x", "y", spID)])
  
  # add layer to stack
  uncertain_tif <- c(uncertain_tif, temp_prediction)
 
})}

uncertain_tif$Mean <- terra::app(uncertain_tif, mean, na.rm=TRUE)
uncertain_tif$SD <- terra::app(uncertain_tif, sd, na.rm=TRUE)

uncertain_tif

# save species' uncertainty map
terra::writeRaster(uncertain_tif, file=paste0("_results/", Taxon_name, "/SDM_Uncertainty_", Taxon_name, ".tif"))
uncertain_tif <- terra::rast(paste0("_results/", Taxon_name, "/SDM_Uncertainty_", Taxon_name, ".tif"))

# extract area with uncertainty lower than threshold
summary(uncertain_tif$Mean) #3rd Qu. E: 0.3, N: 0.445

uncertain_thresh <- stats::quantile(uncertain_tif$Mean, 0.9, na.rm=TRUE)
# 0.9-quantile E:0.326, N: 0.472

extent_df <- terra::as.data.frame(uncertain_tif, xy=TRUE) %>% filter(Mean<uncertain_thresh & !is.na(Mean)) %>% dplyr::select(x,y)
save(extent_df, file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, ".RData"))


