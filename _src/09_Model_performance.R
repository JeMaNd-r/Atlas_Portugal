#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Extract model performance            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)
library(biomod2)

Taxon_name <- "Fungi"

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

#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -
## Summarize model output (evaluations into 1 table)
list_eval <- list.files(paste0("_results/", Taxon_name, "/SDMs"), full.names = TRUE)

data_eval <- lapply(list_eval, function(x) append(get(load(x)),
                                                c("Species" = substr(basename(x), 12, 11+nchar(speciesSub[2]))))) #get species ID
data_eval <- lapply(data_eval, function(x){
  x2 <- as_tibble(x$validation)
  # if ESM output, extract data differently than SDM [SDM output = 2 columns]
  if(ncol(x2)>3){
    x2 <- x2 %>% 
      pivot_wider(id_cols = "algo", names_from = metric.eval, 
                  values_from = "calibration") %>%
      rename(AUC = ROC, MaxTSS = TSS) %>%
      mutate(Model = "ensemble") %>% #because we will filter for this column later
      dplyr::select(-algo)
  }else{ 
    x2$KAPPA <- NA
    x2$Model <- rownames(x$validation)
  }
  
  x2$Species <- x$Species
  return(x2)
  })
data_eval <- do.call(rbind, data_eval)
data_eval <- data_eval %>% filter(Model == "ensemble")
data_eval

write_csv(data_eval %>% dplyr::select(-Model), 
          paste0("_results/Model_evaluation_", Taxon_name, ".csv"))


#- - - - - - - - - - - - - - - - - - - - -
## Extract thresholds from ESMs ####
#- - - - - - - - - - - - - - - - - - - - -

## Summarize model output (evaluations into 1 table)
list_eval <- list.files(paste0("_results/", Taxon_name, "/SDMs"), 
                        pattern = "ESM", full.names = TRUE)

data_eval <- lapply(list_eval, function(x) append(get(load(x)),
                                                  c("Species" = substr(basename(x), 12, 11+nchar(speciesSub[2]))))) #get species ID
data_eval <- lapply(data_eval, function(x){
  x2 <- as_tibble(x$threshold)
  x2$Species <- x$Species
  
  x2 <- x2 %>% dplyr::select(Species, model, everything())
  return(x2)
})
data_eval <- do.call(rbind, data_eval)
data_eval

write_csv(data_eval, 
          paste0("_results/Model_thresholds_", Taxon_name, ".csv"))


#- - - - - - - - - - - - - - - - - - - - -
## Extract thresholds from SDMs ####
#- - - - - - - - - - - - - - - - - - - - -

data_thres <- lapply(as.list(speciesSub), function(spID) {
  print(spID)
  
  # load observed values
  myBiomodData <- get(load(paste0("_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")))
  
  obs <- myBiomodData@data.species
  obs[is.na(obs)] <- 0 
  
  # load predicted values
  temp_prediction <- terra::rast(paste0("_results/", Taxon_name, "/Projection/", spID, ".tif"))
  fit <- unlist(terra::extract(temp_prediction, myBiomodData@coord, ID=FALSE))
  
  # threshold
  optimal_thresh <- do.call(rbind, lapply(list("POD", "FAR", "POFD", "SR", "ACCURACY", 
                                               "BIAS", "ROC", "TSS", "KAPPA", "OR", 
                                               "ORSS", "CSI", "ETS", "BOYCE", "MPA"),
                                          function(x) bm_FindOptimStat(
                                            metric.eval = x,
                                            obs = obs,
                                            fit = fit,
                                            nb.thresh = 1000,
                                            mpa.perc = 0.95
                                          )))
  optimal_thresh$species <- spID
  optimal_thresh <- optimal_thresh %>%
    dplyr::select(-sensitivity, -specificity, -best.stat) %>%
    pivot_wider(names_from = metric.eval, values_from = cutoff)
  
  return(optimal_thresh)
})

data_thres <- do.call(rbind, data_thres)
data_thres

write_csv(data_thres, paste0("_results/Model_thresholds_SDM_", Taxon_name, ".csv"))


