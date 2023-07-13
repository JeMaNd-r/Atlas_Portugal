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

Env_clip_df <- terra::as.data.frame(Env_clip, xy=TRUE)

uncertain_df <- Env_clip_df %>% dplyr::select(x, y)

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
          
  # add layer to stack
  uncertain_df <- full_join(uncertain_df, temp_prediction %>% dplyr::select(x,y, spID))
 
})}

uncertain_df$Mean <- rowMeans(uncertain_df %>% dplyr::select(-x, -y), na.rm=T)

# calculate sd of predictions
uncertain_df$SD <- apply(uncertain_df %>% dplyr::select(-x, -y, -Mean), 1, sd, na.rm = TRUE)

head(uncertain_df)

# save species' uncertainty map
save(uncertain_df, file=paste0("_results/", Taxon_name, "/SDM_Uncertainty_", Taxon_name, ".RData"))
load(file=paste0("_results/", Taxon_name, "/SDM_Uncertainty_", Taxon_name, ".RData")) #uncertain_df


# extract area with uncertainty lower than threshold
summary(uncertain_df$Mean)

extent_df <- uncertain_df %>% filter(Mean<0.1 & !is.na(Mean)) %>% dplyr::select(x,y)
save(extent_df, file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, ".RData"))


