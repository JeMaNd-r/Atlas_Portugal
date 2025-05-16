#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Calculate species richness under     #
#             current climate               #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

gc()
library(tidyverse)
library(raster)

Taxon_name <- "Crassiclitellata"

# load number of occurrences per species and focal species names
speciesSub <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))
if(nrow(speciesSub) != 0){
  speciesSub <- speciesSub %>% 
    full_join(read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))) %>%
    pull(species)
}else{
  speciesSub <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv")) %>% pull(species)
}
speciesSub

# load environmental variables
Env_clip <- terra::rast("_intermediates/EnvPredictor_PCA_1km_POR.tif")
Env_clip <- terra::subset(Env_clip, 1:11) #11 = >80%

# create directory to save species richness stacks
if(!dir.exists(paste0("_results/_Maps"))){ 
  dir.create(paste0("_results/_Maps"))
}

#- - - - - - - - - - - - - - - - - - - - - -
## Create maps and calculate richness ####
#- - - - - - - - - - - - - - - - - - - - - -
# create empty raster
species_rast <- terra::rast(extent=(Env_clip))
rm(Env_clip); gc()

# load model evaluation data (to extract threshold for species)
data_eval <- read.csv(paste0("_results/Model_evaluation_", Taxon_name, ".csv"))

# list all projections
species_rast <- list.files(paste0("_results/", Taxon_name, "/Projection"), full.names = TRUE) 
species_rast

# load into list
species_rast <- terra::rast(species_rast)
species_rast

# transform to binary
mean_thresh <- mean(data_eval$MaxTSS * 1000) # only used when no species-specific threshold available

for(spID in speciesSub){ try({
  print(spID)
  
  # Transform to binary
  temp_thresh <- data_eval[data_eval$Species==spID, "MaxTSS"] * 1000
  
  if(length(temp_thresh)==0) temp_thresh <- mean_thresh
  print(temp_thresh)
  
  temp_rast <- species_rast[[spID]]
  
  # ### ESM: Binary Projection based on max TSS of calibrated ESMs into new space                                                
  # my.ESM_EFproj_current_binary <- (my.ESM_EFproj_current > (my.ESM_thresholds$TSS.th*1000))*1
  
  ## from-to-becomes
  # classify the values into three groups 
  # all values >= temp_thresh become 1
  m <- c(0, temp_thresh, 0,
         temp_thresh, 1001, 1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  temp_rast <- terra::classify(temp_rast, rclmat, include.lowest=TRUE)
  
  species_rast[[spID]] <- temp_rast
  
  print("Binary transformation successful.")
  
  gc()
}, silent=TRUE)}
species_rast

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate richness ####
species_rast$Richness <- terra::app(species_rast, sum, na.rm=TRUE)

#- - - - - - - - - - - - - - - - - - - - - -
## Save species stack ####
terra::writeRaster(species_rast, file=paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))

#load(file=paste0(data_wd, "/_results/_Maps/SDM_stack_binary_", Taxon_name, ".RData")) #species_stack


# ## OLD - - - - - - - - - - - - - - - - - - - - - -
# # for loop through all species
# for(spID in speciesSub){ try({
#   
#   # list files in species-specific BIOMOD folder
#   temp_files <- list.files(paste0(data_wd, "/_results/", Taxon_name, "/", stringr::str_replace(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
#   
#   #setwd(paste0(data_wd, "/results/", Taxon_name))
#   
#   myBiomodEnProj <- terra::rast(temp_files[stringr::str_detect(temp_files,"ensemble.tif")])
#   
#   print(paste0(spID, " successfully loaded."))
#   
#   # extract prediction
#   temp_prediction <- myBiomodEnProj[[2]] #column with weighted mean
#   names(temp_prediction) <- spID
#   
#   #print(paste0("Saved binary prediction of ", spID))
#   
#   # add layer to stack
#   species_rast <- c(species_rast, temp_prediction)
#   
#   print(paste0("Added binary prediction of ", spID, " to the species stack"))
#   
#   rm(temp_thresh, best_pred)
# }, silent=T)}  
# 
# species_rast
# 
# #species_rast <- terra::app(species_rast, function(x){x/1000})
# 
