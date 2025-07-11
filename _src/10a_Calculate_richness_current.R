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

for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes")){ #, "Bacteria", 
    
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
  
  # # load model evaluation data (to extract threshold for species)
  # if(file.exists(paste0("_results/Model_evaluation_", Taxon_name, ".csv"))){
  #   data_eval <- read.csv(paste0("_results/Model_evaluation_", Taxon_name, ".csv"))
  #   # transform to binary
  #   mean_thresh <- mean(data_eval$MaxTSS * 1000) # only used when no species-specific threshold available
  # } else {
  #   mean_thresh <- 500
  # }
  
  # load SDM thresholds
  if(file.exists(paste0("_results/Model_thresholds_SDM_", Taxon_name, ".csv"))){
    data_thres_sdm <- read_csv(paste0("_results/Model_thresholds_SDM_", Taxon_name, ".csv"))
  }
  
  # load ESM thresholds
  if(file.exists(paste0("_results/Model_thresholds_", Taxon_name, ".csv"))){
    data_thres <- read_csv(paste0("_results/Model_thresholds_", Taxon_name, ".csv"))
  }
  
  # list all projections
  species_rast <- list.files(paste0("_results/", Taxon_name, "/Ensembles"), full.names = TRUE) 
  species_rast
  
  # load into list
  species_rast <- terra::rast(species_rast)
  species_rast
  
  for(spID in speciesSub){ try({
    print(spID)
    
    if(exists("data_thres") & substr(spID, 1,9) %in% unique(data_thres$Species)){ 
      
      #take maximum threshold value of selected threshold metrices
      temp_thresh <- max(unlist(data_thres[data_thres$Species==substr(spID, 1,9),
                                           c("Kappa", "AUC", "SomersD", "TSS", "TSS.th", "MPA0.95", "Boyce.th.max")]),
                         na.rm=TRUE)* 1000 #max TSS value (from GLM and MAXENT full models)
      # temp_thresh <- max(unlist(data_thres[data_thres$Species==substr(spID, 1,9),
      #                                      c("TSS")]),
      #                    na.rm=TRUE)* 1000 #max TSS value (from GLM and MAXENT full models)
      print(temp_thresh)
      
    } else {
      if((exists("data_thres_sdm") & substr(spID, 1,9) %in% unique(data_thres_sdm$species))){
        # temp_thresh <- mean(c(data_eval[data_eval$Species==substr(spID, 1,9), "AUC"],
        #                       data_eval[data_eval$Species==substr(spID, 1,9), "MaxTSS"]))* 1000
        
        #take maximum threshold value of selected threshold metrices
        temp_thresh <- max(unlist(data_thres_sdm[data_thres_sdm$species==substr(spID, 1,9),
                                                 c("KAPPA", "ROC", "TSS", "MPA", "BOYCE")]),
                           na.rm=TRUE) #max TSS value (from GLM and MAXENT full models)
        # temp_thresh <- max(unlist(data_thres_sdm[data_thres_sdm$species==substr(spID, 1,9), 
        #                                          c("TSS")]), 
        #                    na.rm=TRUE) #max TSS value (from GLM and MAXENT full models)
        print(temp_thresh)
  
      } else {
      
      if(length(temp_thresh)==0 | temp_thresh == -Inf) temp_thresh <- 500
      print(temp_thresh)
    }}
    
    temp_rast <- species_rast[[spID]]
    
    # ### ESM: Binary Projection based on max TSS of calibrated ESMs into new space                                                
    # my.ESM_EFproj_current_binary <- (my.ESM_EFproj_current > (my.ESM_thresholds$TSS.th*1000))*1
    
    ## from-to-becomes
    # classify the values into three groups 
    # all values >= temp_thresh become 1
    m <- c(0, temp_thresh, 0,
           temp_thresh, 2000, 1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    temp_rast <- terra::classify(temp_rast, rclmat, include.lowest=TRUE)
    
    species_rast[[spID]] <- temp_rast
    
    print("Binary transformation successful.")
    
    gc()
  }, silent=TRUE)}
  species_rast
  
  plot(species_rast)
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Calculate richness ####
  species_rast$Richness <- terra::app(species_rast, sum, na.rm=TRUE)
  plot(species_rast$Richness)
  
  #- - - - - - - - - - - - - - - - - - - - - -
  ## Save species stack ####
  terra::writeRaster(species_rast, file=paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"), overwrite=TRUE)
  
  #load(file=paste0(data_wd, "/_results/_Maps/SDM_stack_binary_", Taxon_name, ".RData")) #species_stack
}

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
