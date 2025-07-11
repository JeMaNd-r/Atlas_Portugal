#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Extract uncertainty               #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

gc()
library(tidyverse)
library(here)

library(raster)

for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes")){ #, "Bacteria", 
    
  print(Taxon_name)
  
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
  
  species10 <- read_csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))$species
  species100 <- read_csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))$species
  
  covarsNames <- paste0("PC", 1:11)
  
  #- - - - - - - - - - - - - - - - - - - - -
  # load environmental variables
  Env_clip <- terra::rast("_intermediates/EnvPredictor_PCA_1km_POR.tif")
  Env_clip <- terra::subset(Env_clip, 1:11) #11 = >80%
  
  #Env_clip_df <- terra::as.data.frame(Env_clip, xy=TRUE)
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Create maps and calculate richness ####
  #- - - - - - - - - - - - - - - - - - - - - -
  for(i in c(10, 100)){try({
    print(i)
    
    # list all projections
    uncertain_rast <- paste0("_results/", Taxon_name, "/Uncertainty/CV_", get(paste0("species", i)),".tif")
    
    # load into list
    uncertain_rast <- terra::rast(uncertain_rast)
    
    # calculate mean and SD
    uncertain_rast$Mean <- terra::app(uncertain_rast, mean, na.rm=TRUE)
    uncertain_rast$SD <- terra::app(uncertain_rast, sd, na.rm=TRUE)
    
    # save species' uncertainty map
    terra::writeRaster(uncertain_rast, file=paste0("_results/SDM_Uncertainty_", Taxon_name, "_", i, ".tif"), overwrite = TRUE)
    uncertain_rast <- terra::rast(paste0("_results/SDM_Uncertainty_", Taxon_name, "_", i, ".tif"))
    
    # extract area with uncertainty lower than threshold
    print(summary(uncertain_rast$Mean)) #3rd Qu. E: 421.2, N: 10-449.1 100-52.12, 
    
    uncertain_thresh <- stats::quantile(uncertain_rast$Mean, 0.9, na.rm=TRUE)
    print(uncertain_thresh)
    # 0.9-quantile E:452.75, N: 10-485.5714 100-57.45907, F: 10-488.0533  100-46.26587, Eu: 10-491.0213 100-48.25674, Pr: 10-495.4902 100-46.05416
    
    extent_df <- terra::as.data.frame(uncertain_rast, xy=TRUE) %>% filter(Mean<uncertain_thresh & !is.na(Mean)) %>% dplyr::select(x,y)
    save(extent_df, file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_", i, ".RData"))
  })}
  
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Combine 10 and 100 uncertainties ####
  #- - - - - - - - - - - - - - - - - - - - - -
  
  try(uncertain_10_tif <- terra::rast(paste0("_results/SDM_Uncertainty_", Taxon_name, "_10.tif")))
  try(uncertain_100_tif <- terra::rast(paste0("_results/SDM_Uncertainty_", Taxon_name, "_100.tif")))
  
  # rename mean and SD layers
  try(names(uncertain_10_tif)[c(length(names(uncertain_10_tif))-1,length(names(uncertain_10_tif)))] <- c("Mean_10", "SD_10"))
  try(names(uncertain_100_tif)[c(length(names(uncertain_100_tif))-1,length(names(uncertain_100_tif)))] <- c("Mean_100", "SD_100"))
  
  if(exists("uncertain_thresh_10") & exists("uncertain_thresh_100")){
    uncertain_tif <- c(uncertain_10_tif, uncertain_100_tif)
    uncertain_tif <- uncertain_tif[[sort(names(uncertain_tif))]]
  
    terra::writeRaster(uncertain_tif, file=paste0("_results/_Maps/SDM_Uncertainty_", Taxon_name, ".tif"), overwrite = TRUE)
  } else { 
    if(exists("uncertain_thresh_10")){
      uncertain_tif <- uncertain_10_tif
      uncertain_tif <- uncertain_tif[[sort(names(uncertain_tif))]]
      
      terra::writeRaster(uncertain_tif, file=paste0("_results/_Maps/SDM_Uncertainty_", Taxon_name, "_10.tif"), overwrite = TRUE)
    } else { 
      if(exists("uncertain_thresh_100")){
        uncertain_tif <- uncertain_100_tif
        uncertain_tif <- uncertain_tif[[sort(names(uncertain_tif))]]
        
        terra::writeRaster(uncertain_tif, file=paste0("_results/_Maps/SDM_Uncertainty_", Taxon_name, "_100.tif"), overwrite = TRUE)
  }}
  }
}
