#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Predict ESMs in current climate      #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal") 

gc()
library(tidyverse)
library(here)

library(raster)

library(biomod2) # also to create pseudo-absences

library(parallel)
library(doParallel)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

# load number of occurrences per species and focal species names
speciesSub <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv")) %>% pull(SpeciesID)
speciesSub

covarsNames <- paste0("PC", 1:11)

# read fixed ecospat function
source("_src/FC_ecospat.ESM.Projection.R")

#- - - - - - - - - - - - - - - - -
## Predict in current climate at 5km ####

# load environmental variables
Env_clip <- terra::rast("_intermediates/EnvPredictor_PCA_1km_POR.tif")
Env_clip <- terra::subset(Env_clip, 1:11) #11 = >80%

#registerDoParallel(3)
# foreach(spID = speciesSub,
#         .export = c("Env_clip"),
#         .packages = c("tidyverse","biomod2")) %dopar% { try({   
for(spID in speciesSub){ try({     
          # list files in species-specific BIOMOD folder
          temp_files <- list.files(paste0("_results/", Taxon_name, "/ESM.BIOMOD.output_", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
          
          # load model output
          myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"ESM_Modeling.")]
          print(myBiomodModelOut)
          myBiomodModelOut <-get(load(myBiomodModelOut))
          
          # load ensemble model output
          myBiomodEM <- temp_files[stringr::str_detect(temp_files,"ESM_EnsembleModeling.")]
          myBiomodEM <- get(load(myBiomodEM))
          
          setwd(paste0("_results/", Taxon_name, "/ESM.BIOMOD.output_", stringr::str_replace(spID, "_", ".")))
          
          tmp <- proc.time()[3]
          ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
          # project single models (also needed for ensemble model)
          myBiomodProj <- ecospat::ecospat.ESM.Projection(ESM.modeling.output = myBiomodModelOut,
                                                     new.env = Env_clip,        #column/variable names have to perfectly match with training
                                                     name.env = "POR",
                                                     )
          
          # project ensemble of all models
          myBiomodEnProj <- ecospat::ecospat.ESM.EnsembleProjection(ESM.prediction.output = myBiomodProj,
                                                                    ESM.EnsembleModeling.output = myBiomodEM,
                                                                #... same arguments as above could be added but are not necessary when loading myBiomodProj
                                                                chosen.models = "all")
          
          temp_predict_time <- proc.time()[3] - tmp
          
          # # extracting the values for ensemble validation
          # myEnProjDF <- get_predictions(myBiomodEM) %>%
          #   filter(algo=="EMwmean") #for weighted probability mean
          
          # see the first few validations
          # note: the prediction scale of biomod is between 0 and 1000
          #head(myEnProjDF)
          
          # Get model evaluation values for later
          # myBiomodModelEval <- biomod2::get_evaluations(myBiomodEM) %>% 
          #   dplyr::select(algo:evaluation)
          ### Evaluation of the ensemble models based on the pooling procedure 
          myBiomodModelEval <- ecospat.ESM.EnsembleEvaluation(ESM.modeling.output= myBiomodModelOut,
                                                               ESM.EnsembleModeling.output = myBiomodEM,
                                                               metrics= c("AUC","MaxTSS"),
                                                               EachSmallModels = FALSE)
          myBiomodModelEval <- myBiomodModelEval$ESM.evaluations
          
          # # Calculate variable importance across all PA sets, eval runs and algorithms
          # # and extract only the one for weighed mean predictions (for later)
          # temp_varImp <- biomod2::get_variables_importance(myBiomodEM)  %>%
          #   filter(algo=="EMwmean") %>% 
          #   dplyr::select(expl.var:var.imp)
          # # average across 10 runs
          # temp_varImp <- temp_varImp %>% 
          #   rename(Predictors = expl.var,
          #          varImp = var.imp) %>%
          #   group_by(Predictors) %>%
          #   summarize(across(varImp, .fns=mean))
          
          # save predictions as raster file
          temp_prediction <- myBiomodEnProj@proj.out@val %>%
            filter(algo == "EMwmean") %>%
            dplyr::select(pred) %>%
            # add names of grid cell (only for those that have no NA in any layer) 
            mutate(grid_id = rownames(Env_clip_df),
                   x = Env_clip_df$x,
                   y = Env_clip_df$y) %>% 
            full_join(Env_clip_df %>% dplyr::select(x,y)) %>%
            rename(layer = pred) %>%
            mutate(layer = layer / 1000)
          
          biomod_list <- list(time_predict = temp_predict_time, 
                              validation = myBiomodModelEval, 
                              prediction = temp_prediction 
                              #varImp = temp_varImp
                              )
          save(biomod_list, file=paste0("SDM_biomod_", spID, ".RData"))
          
          rm(biomod_list, temp_predict_time, temp_prediction, 
             #temp_varImp, 
             myBiomodEnProj, myBiomodProj, myBiomodModelEval, myBiomodModelOut, myBiomodEM)
          
          setwd("../")
          setwd("../")
          setwd("../")
        })}

stopImplicitCluster()
