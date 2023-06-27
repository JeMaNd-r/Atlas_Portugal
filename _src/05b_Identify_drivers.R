#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#    Identify top predictors with MaxEnt    #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal")

options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) # expand Java memory

gc()
library(tidyverse)
library(here)

library(raster)

library(biomod2) # also to create pseudo-absences

library(dismo) # for MaxEnt and BRT
# download maxent.jar 3.3.3k, and place the file in the
# desired folder
# utils::download.file(url = "https://raw.githubusercontent.com/mrmaxent/Maxent/master/ArchivedReleases/3.3.3k/maxent.jar", 
#     destfile = paste0(system.file("java", package = "dismo"), 
#         "/maxent.jar"), mode = "wb")  ## wb for binary file, otherwise maxent.jar can not execute

library(rJava)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

Taxon_name <- "earthworms"

# load number of occurrences per species and focal species names
speciesSub <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))
speciesSub <- speciesSub %>% pull(SpeciesID)
speciesSub

#- - - - - - - - - - - - - - - - - - - - -
# note: we will load the datasets before each individual model

# load environmental variables (for projections)
Env_clip <- terra::rast("_intermediates/EnvPredictor_1km_POR.tif")

# remove correlated variables
env_vif <- read_csv(file="_results/VIF_predictors_1km_POR.csv")
env_vif

env_exclude_vif <- env_vif %>% filter(is.na(VIF)) %>% pull(Variables) 
Env_clip <- terra::subset(Env_clip, names(Env_clip)[!(names(Env_clip) %in% env_exclude_vif)])
names(Env_clip)

Env_clip_df <- terra::as.data.frame(Env_clip, xy=TRUE)

# define covariate names
covarsNames <- names(Env_clip)

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))

# create subfolder if not existing
if(!dir.exists(paste0("_results/_TopPredictor/", Taxon_name))){ 
  dir.create(paste0("_results/_TopPredictor/", Taxon_name))
}

#- - - - - - - - - - - - - - - - - - - - -
## Prepare model input ####

mySpeciesOcc <- read_csv(file=paste0("_intermediates/Occurrence_rasterized_1km_", Taxon_name, ".csv"))

for(spID in speciesSub) { try({

  myData <- mySpeciesOcc[!is.na(mySpeciesOcc[,spID]), c("x","y",spID)]
  myData_env <- terra::extract(x = Env_clip, y = myData[,c("x", "y")])
  
  myData <- myData %>% cbind(myData_env[,-1])
  
  save(myData, file=paste0("_results/_TopPredictor/", Taxon_name, "/MaxentData_noValid_", Taxon_name,"_", spID, ".RData"))

  rm(myData_env, myData)

})}


#- - - - - - - - - - - - - - - - - - - - -
## Modeling (no validation data) ####
# NOTE: This part might only work in R but not in RStudio.

# "We used five different regularization multipliers (0.5, 1, 2, 3 and 4)
# in combination with different features (L, LQ, H, LQH, LQHP) to find the 
# best parameters that maximizes the average AUC-ROC in CV."

modelName <- "MaxentData_noValid"

for(spID in speciesSub){ try({
  
  # load data
  load(paste0("_results/_TopPredictor/", Taxon_name, "/MaxentData_noValid_", Taxon_name,"_", spID, ".RData")) #myData
  
  tmp <- proc.time()[3]
  set.seed(32639)
   
  ## fit a maxent model with the tuned parameters
  maxent <- dismo::maxent(x = myData[, covarsNames],
                          p = myData[, spID],
                          removeDuplicates = FALSE, #remove occurrences that fall into same grid cell (not necessary)
                          path = paste0("_results/_TopPredictor/maxent_files/", spID), #wanna save files?
                          #args = c("responsecurves")
                          )
  
  temp_model_time <- proc.time()[3] - tmp
  
  
  tmp <- proc.time()[3]
  # create raster layer of predictions for whole environmental space
  #temp_prediction <- raster::predict(Env_clip, maxent)
  #temp_prediction <- data.frame(raster::rasterToPoints(temp_prediction))
  gc()
  #temp_prediction <- dismo::predict(maxent, Env_clip_df %>% dplyr::select(-x, -y)) # Java out of memory
  temp_prediction <- dismo::predict(maxent, Env_clip_df[,colnames(Env_clip_df) %in% covarsNames], type = c("cloglog"))
  temp_prediction <- as.numeric(temp_prediction)
  names(temp_prediction) <- rownames(Env_clip_df[complete.cases(Env_clip_df[,colnames(Env_clip_df) %in% covarsNames]),]) #add site names
  temp_prediction <- as.data.frame(temp_prediction)
  temp_prediction$x <- Env_clip_df[complete.cases(Env_clip_df[,colnames(Env_clip_df) %in% covarsNames]),]$x
  temp_prediction$y <- Env_clip_df[complete.cases(Env_clip_df[,colnames(Env_clip_df) %in% covarsNames]),]$y
  colnames(temp_prediction)[1] <- "layer"
  
  temp_runs <- 1
  
  temp_predict_time <- proc.time()[3] - tmp
  
  # varImp
  temp_varImp <- as.data.frame(maxent@results[str_detect(rownames(maxent@results),"permutation.importance"),])
  # Note: permutation importance = determine the importance of predictors 
  # calculated by permuting values of each predictor &  resulting reduction
  # in training AUC: large reduction = model is influenced by that predictor. 
  
  # extract predictor names
  temp_varImp$Predictor <- stringr::str_split_fixed(rownames(temp_varImp), "[.]perm", 2)[,1]
  colnames(temp_varImp) <- c("maxent", "Predictor")  
  
  maxent_list <- list(bg_data=modelName, time_model=temp_model_time, time_predict=temp_predict_time, runs=temp_runs, prediction=temp_prediction, varImp=temp_varImp)
  save(maxent_list, file=paste0("_results/_TopPredictor/", Taxon_name, "/SDM_maxent_noValid_", spID, ".RData"))
  
  rm(maxent, maxent_list, param_optim, temp_model_time, temp_predict_time, temp_runs, temp_prediction, temp_varImp)
  gc()
  
  rm(myData)

})}

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable importance (VI) ####

# create result data frame
var_imp <- data.frame("Predictor"= covarsNames,
                      "maxent"=NA, "Species"=NA)

for(spID in speciesSub){ try({
  
  print("=====================================")
  print(spID)
  
  temp_sdm <- get(load(file=paste0("_results/_TopPredictor/", Taxon_name, "/SDM_maxent_noValid_", spID, ".RData")))

  temp_vi <- temp_sdm[["varImp"]]
  
  # harmonizes predictor column structure
  temp_vi$Predictor <- as.character(temp_vi$Predictor)

  # round variable importance to 3 decimals
  temp_vi[,"maxent"] <- round(temp_vi[,"maxent"], 3)
      
  # add species name
  temp_vi$Species <- spID
  
  var_imp <- rbind(var_imp, temp_vi)

  rm(temp_vi, temp_sdm, temp_varImp)
    
  #var_imp
          
})}

rownames(var_imp) <- 1:nrow(var_imp)
var_imp <- var_imp[!is.na(var_imp$Species),] %>% unique()

var_imp
str(var_imp) # X species * Y variables

## Save
write_csv(var_imp, file=paste0("_results/Variable_importance_MaxEnt_", Taxon_name, ".csv"))
 
# save top 10 predictors
top_pred <- var_imp %>% group_by(Predictor) %>% 
  summarize(across(maxent, mean)) %>%
  rename(mean=maxent) %>%
  arrange(desc(mean)) %>% 
  #top_n(10) %>%
  mutate(Taxon = Taxon_name) %>%
  full_join(var_imp %>% group_by(Predictor) %>% filter(maxent>0) %>% count(name="No_species>0"),
            by="Predictor")
top_pred

write_csv(top_pred, file=paste0("_results/TopPredictor_MaxEnt_", Taxon_name, ".csv"))
