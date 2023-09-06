#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Species Distribution Models          #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# Note: dismo::maxent() might crash in RStudio


#setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal")

library(biomod2)

#options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8g")) # expand Java memory
gc()
library(tidyverse)

# read functions for ensemble of small models (ESM)
#devtools::source_url("https://raw.githubusercontent.com/cran/ecospat/master/R/ecospat.ESM.R")
library(ecospat)
#source("_src/F_ecospat.ESM.EnsembleModeling_fixed.R")

library(PresenceAbsence) # needed for ecospat
library(ENMeval)
library(gam)

library(raster)

library(parallel)
library(doParallel)

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal")
Taxon_name <- "nematodes"
species_csv <- paste0("ESM_", Taxon_name, ".csv") 
output_dir <- paste0(getwd(), "/_results")
input_dir <- getwd()

# Load focal species names
species_table <- read_csv(paste0(input_dir, "/_intermediates/", species_csv)) #for Computer

covarsNames <- paste0("PC", 1:11)

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))

#- - - - - - - - - - - - - - - - - - - - -
## Build models ####
# define parameters of the algorithms
myESMOption <- BIOMOD_ModelingOptions(
  GLM = list (type = "quadratic",
              interaction.level = 0,
              myFormula = NULL,
              test = "AIC",
              family = binomial(link = "logit") ),
  
  MAXENT = list(path_to_maxent.jar = paste0(input_dir, "/_results"), # change it to maxent directory
                memory_allocated = NULL, # use default from Java defined above
                visible = FALSE, 	# don't make maxEnt user interface visible
                linear = TRUE, 	# linear features allowed
                quadratic = TRUE, # quadratic allowed
                product = TRUE,	# product allowed
                threshold = TRUE,	# threshold allowed
                hinge = TRUE,	# hinge allowed
                lq2lqptthreshold = 80, # default
                l2lqthreshold = 10, # default
                hingethreshold = 15, # default
                beta_threshold = -1, # default
                beta_categorical = -1, # default
                beta_lqp = -1, # default
                beta_hinge = -1, # default
                betamultiplier = 1, # default
                defaultprevalence = 0.5 ), #default   
)


# registerDoParallel(3)
# foreach(spID = speciesSub, 
#         .export = c("Env_clip", "Env_clip_df", "form", "occ_points", "ecospat.ESM.EnsembleModeling_fixed"),
#         .packages = c("tidyverse","biomod2", "ecospat")) %dopar% { try({

for(spID in species_table$SpeciesID){
  try({ #some species may fail due to low number of presences
    
          load(paste0(input_dir, "/_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData
    
          # model fitting
          #tmp <- proc.time()[3]
          setwd(paste0(input_dir, "/_results/", Taxon_name))
          
          set.seed(32639)
          
          myEsmModelOut <- ecospat::ecospat.ESM.Modeling(data = myBiomodData,
                                        models = c("GLM", "CTA", "MAXENT"),
                                        models.options = myESMOption,
                                        Prevalence = NULL, 
                                        tune = FALSE, # TRUE: estimate optimal parameters for the models
                                        NbRunEval = 10,   # 10-fold crossvalidation evaluation run
                                        DataSplit = 80, # use subset of the data for training
                                        weighting.score = "TSS",
                                        cleanup = 2, #when to delete temporary unused files, in hours
                                        modeling.id = paste(spID,"_Modeling", sep = ""))
          
	        myEsmEM <- ecospat.ESM.EnsembleModeling(ESM.modeling.output=myEsmModelOut,
	                                                         weighting.score="TSS",
	                                                         threshold=NULL)

          #temp_model_time <- proc.time()[3] - tmp
          
          setwd(data_wd)
})}
stopImplicitCluster()


# #- - - - - - - - - - - - - - - - -
# ## Predict in current climate at 5km ####
# 
# # load environmental variables (for projections)
# Env_norm <- raster::stack(paste0(data_wd, "/results/EnvPredictor_5km_normalized.grd"))
# #Env_norm <- stack(Env_norm)
# 
# # as dataframe
# load(paste0(data_wd,"/results/EnvPredictor_5km_df_normalized.RData")) #Env_norm_df
# 
# registerDoParallel(no.cores)
# foreach(spID = speciesSub, 
#         .export = c("Env_norm", "Env_norm_df", "form"),
#         .packages = c("tidyverse","biomod2")) %dopar% { try({   
#           
#           # list files in species-specific BIOMOD folder
#           temp_files <- list.files(paste0(data_wd, "/", Taxon_name, "/BIOMOD_files/ESM.BIOMOD.output_", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
#           
#           setwd(paste0(data_wd, "/results/", Taxon_name))
#           
#           # load model output
#           myEsmModelOut  <- temp_files[stringr::str_detect(temp_files,"Modeling.out")]
#           print(myEsmModelOut )
#           myEsmModelOut  <-get(load(myEsmModelOut ))
#           
#           # load ensemble model output
#           myEsmEM <- temp_files[stringr::str_detect(temp_files,"Modelingensemble.models.out")]
#           myEsmEM <- get(load(myEsmEM))
#           
#           tmp <- proc.time()[3]
#           ## NOTE: because biomod output can hardly be stored in list file, we will do calculations based on model output now
#           # project single models (also needed for ensemble model)
#           myEsmProj <- biomod2::BIOMOD_Projection(modeling.output = myEsmModelOut,
#                                                      new.env = Env_norm_df[,colnames(Env_norm_df) %in% covarsNames],        #column/variable names have to perfectly match with training
#                                                      proj.name = "modeling",  #name of the new folder being created
#                                                      selected.models = "all", #use all models
#                                                      binary.meth = NULL,     #binary transformation according to criteria, or no transformation if NULL
#                                                      compress = TRUE,         #compression format of objects stored on hard drive
#                                                      build.clamping.mask = TRUE, #TRUE: clamping mask will be saved on hard drive different
#                                                      do.stack = TRUE,         #save output projections as rasterstack (if not too heavy)
#                                                      output.format = ".RData", #what format should projections have: RData, grd or img
#                                                      keep.in.memory = TRUE)  #FALSE: only story link to copy to projection file
#           
#           # project ensemble of all models
#           myEsmEnProj <- biomod2::BIOMOD_EnsembleForecasting(projection.output = myEsmProj,
#                                                                 EM.output = myEsmEM,
#                                                                 #... same arguments as above could be added but are not necessary when loading myBiomodProj
#                                                                 selected.models = "all")
#           
#           temp_predict_time <- proc.time()[3] - tmp
#           
#           # extracting the values for ensemble validation
#           myEnProjDF <- as.data.frame(get_predictions(myEsmEM)[,2]) #for weighted probability mean
#           
#           # see the first few validations
#           # note: the prediction scale of biomod is between 0 and 1000
#           #head(myEnProjDF)
#           
#           # Get model evaluation values for later
#           myEsmModelEval <- as.data.frame(biomod2::get_evaluations(myEsmEM))
#           
#           # Calculate variable importance across all PA sets, eval runs and algorithms
#           # and extract only the one for weighed mean predictions (for later)
#           temp_varImp <- biomod2::get_variables_importance(myEsmEM)[, , 2]
#           # average across 3 runs
#           temp_varImp <- temp_varImp %>% as.data.frame() %>%
#             mutate(mean_vi = as.numeric(rowMeans(temp_varImp, na.rm=T)),
#                    Predictor = rownames(temp_varImp))
#           colnames(temp_varImp)[colnames(temp_varImp) == "mean_vi"] <- "biomod"
#           temp_varImp <- temp_varImp[,c("biomod", "Predictor")]
#           
#           # save predictions as raster file
#           temp_prediction <- myEsmEnProj@proj@val[,2]
#           temp_prediction <- as.numeric(temp_prediction)
#           # add names of grid cell (only for those that have no NA in any layer)
#           names(temp_prediction) <- rownames(Env_norm_df)
#           temp_prediction <- as.data.frame(temp_prediction)
#           temp_prediction$x <- Env_norm_df$x
#           temp_prediction$y <- Env_norm_df$y
#           temp_prediction <- temp_prediction %>% full_join(Env_norm_df %>% dplyr::select(x,y)) %>%
#             rename("layer" = temp_prediction)
#           temp_prediction$layer <- temp_prediction$layer / 1000
#           
#           temp_runs <- 1
#           
#           biomod_list <- list(time_predict=temp_predict_time, validation=myEsmModelEval, prediction=temp_prediction, varImp=temp_varImp, evaluation=myEsmModelEval)
#           save(biomod_list, file=paste0("./SDMs/SDM_biomod_", spID, ".RData"))
#           
#           rm(biomod_list, temp_predict_time, temp_runs, temp_prediction, temp_varImp, myEsmEnProj, myEsmProj, myEsmModelEval, myEnProjDF, myEsmModelOut, myEsmEM)
#           
#           setwd(data_wd)
#           
#         })}
# stopImplicitCluster()
