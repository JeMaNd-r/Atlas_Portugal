#- - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#                                                        #
#     Calculate ensemble models from single models       #
#              author: Romy Zeiss                        #
#                date: 2025-07-10                        #
#                                                        #
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# This script calculates the ensemble projection from the full GLM and
# full MaxEnt model projections predicted using Ensembles of Small Models
# (package ecospat). This step is only necessary because the wrong layers
# were saved from the ensemble projection output.

# load packages
library(dplyr)
library(readr)
library(terra)

for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes")){ #, "Bacteria", 

  print(Taxon_name)
  
  # load species names
  speciesSub <- read_csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv")) %>% pull(species)
  
  # create subfolder if not existing
  if(!dir.exists(paste0("_results/", Taxon_name, "/Ensembles"))){ 
    dir.create(paste0("_results/", Taxon_name, "/Ensembles"))
    dir.create(paste0("_results/", Taxon_name, "/OLD_SingleProjections"))
  }

  for(spID in speciesSub){
 
    print(spID)
    
    # next species if copy of projection files already exist (and therefore CV is overwritten)
    if(file.exists(paste0("_results/", Taxon_name, "/OLD_SingleProjections/GLM_", spID, ".tif"))) next
    
    # load model evaluation scores (as weights)
    mod_eval <- get(load(paste0("_results/", Taxon_name, "/SDMs/ESM_biomod_", spID, ".RData")))
    mod_eval <- mod_eval$validation
    
    # load GLM projection
    projections <- rast(paste0("_results/", Taxon_name, "/Uncertainty/CV_", spID, ".tif"))
    
    # load MaxEnt projection
    projections$MaxEnt <- rast(paste0("_results/", Taxon_name, "/Projection/", spID, ".tif"))
    
    # save projections into correct folder
    terra::writeRaster(projections[[1]], paste0("_results/", Taxon_name, "/OLD_SingleProjections/GLM_", spID, ".tif"))
    terra::writeRaster(projections[[2]], paste0("_results/", Taxon_name, "/OLD_SingleProjections/MaxEnt_", spID, ".tif"))
    
    # projections <- rast(paste0("_results/", Taxon_name, "/OLD_SingleProjections/GLM_", spID, ".tif"))
    # projections$MaxEnt <- rast(paste0("_results/", Taxon_name, "/OLD_SingleProjections/MaxEnt_", spID, ".tif"))
    
    # calculate weights based on evaluation score (raw TSS)
    eval_scores <- c(mod_eval["GLM", "MaxTSS"],
                     mod_eval["MAXENT", "MaxTSS"])
    
    #... if one of them NA or missing?
    
    weights <- eval_scores / sum(eval_scores)
    
    # weighted ensemble mean projection
    weighted_mean <- app(projections, fun = function(x) sum(x * weights, na.rm = TRUE))
    names(weighted_mean) <- spID
    
    # Unweighted mean and standard deviation across raw projections
    mean_unweighted <- app(projections, mean, na.rm = TRUE)
    sd_unweighted   <- app(projections, sd, na.rm = TRUE)
    
    # Coefficient of variation (unweighted)
    cv_unweighted <- sd_unweighted / mean_unweighted
    names(cv_unweighted) <- spID
    
    # save CV and projection
    terra::writeRaster(weighted_mean, paste0("_results/", Taxon_name, "/Ensembles/", spID, ".tif"))
    terra::writeRaster(cv_unweighted, paste0("_results/", Taxon_name, "/Uncertainty/CV_", spID, ".tif"), overwrite = TRUE)
}}

# 
# #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# ## Re-do maps for all (used wrong weights...) ####
# #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 
# for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes")){ #, "Bacteria", 
#   
#   print(Taxon_name)
#   
#   # load species names
#   speciesSub <- read_csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv")) %>% pull(species)
#   
#   for(spID in speciesSub){
#     
#     print(spID)
#     
#     # load model evaluation scores (as weights)
#     mod_eval <- get(load(paste0("_results/", Taxon_name, "/SDMs/ESM_biomod_", spID, ".RData")))
#     mod_eval <- mod_eval$validation
#    
#     projections <- rast(paste0("_results/", Taxon_name, "/OLD_SingleProjections/GLM_", spID, ".tif"))
#     projections$MaxEnt <- rast(paste0("_results/", Taxon_name, "/OLD_SingleProjections/MaxEnt_", spID, ".tif"))
#     
#     # calculate weights based on evaluation score (raw TSS)
#     eval_scores <- c(mod_eval["GLM", "MaxTSS"],
#                      mod_eval["MAXENT", "MaxTSS"])
#     
#     weights <- eval_scores / sum(eval_scores)
#     
#     # weighted ensemble mean projection
#     weighted_mean <- app(projections, fun = function(x) sum(x * weights, na.rm = TRUE))
#     names(weighted_mean) <- spID
#     
#     # Unweighted mean and standard deviation across raw projections
#     mean_unweighted <- app(projections, mean, na.rm = TRUE)
#     sd_unweighted   <- app(projections, sd, na.rm = TRUE)
#     
#     # Coefficient of variation (unweighted)
#     cv_unweighted <- sd_unweighted / mean_unweighted
#     names(cv_unweighted) <- spID
#     
#     # save CV and projection
#     terra::writeRaster(weighted_mean, paste0("_results/", Taxon_name, "/Ensembles/", spID, ".tif"), overwrite = TRUE)
#     terra::writeRaster(cv_unweighted, paste0("_results/", Taxon_name, "/Uncertainty/CV_", spID, ".tif"), overwrite = TRUE)
#   }}
