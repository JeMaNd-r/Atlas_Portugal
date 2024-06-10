#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#    Species Distribution Models - full     #
#          author: Romy Zeiss               #
#            date: 2024-01-15               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

# Note: dismo::maxent() might crash in RStudio

setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal")
Taxon_name <- "Fungi"
species_csv <- paste0("SDM_", Taxon_name, ".csv") 
output_dir <- paste0(getwd(), "/_results")
input_dir <- getwd()

library(biomod2)

gc()
library(tidyverse)

library(parallel)
library(doParallel)

# # save history
# #devtools::install_github("leeper/rite")
# library("rite")

# Load number of occurrences per species and focal species names
species_table <- read_csv(paste0(input_dir, "/_intermediates/", species_csv)) #for Computer

# covariates
covarsNames <- paste0("PC", 1:11)

# create subfolder if not existing
if(!dir.exists(paste0("_results/", Taxon_name))){ 
  dir.create(paste0("_results/", Taxon_name))
  
  dir.create(paste0("_results/", Taxon_name, "/SDMs"))
  dir.create(paste0("_results/", Taxon_name, "/Uncertainty"))
  dir.create(paste0("_results/", Taxon_name, "/Projection"))
}

# define parameters of the algorithms
mySDMOption <- BIOMOD_ModelingOptions(
  GLM = list (type = "quadratic",
              interaction.level = 0,
              myFormula = NULL,
              test = "AIC",
              family = binomial(link = "logit") ),
  
  GAM = list (algo = "GAM_mgcv",
              myFormula = NULL,
              type = "s_smoother",
              interaction.level = 0,
              family =  binomial(link = "logit"),
              method = "GCV.Cp",
              optimizer = c("outer","newton"),
              select = FALSE,
              knots = NULL,
              paraPen = NULL,
              k = -1 ), 		#avoid error messages
  
  MARS = list(myFormula = NULL,
              nk = NULL, 		# maximum number of model terms, NULL: max(21, 2*nb_expl_var+1)
              penalty = 2, 	# default
              thresh = 0.001, 	# default
              nprune = 1+length(covarsNames), # max. number of terms including intercept
              pmethod = "backward" ), #pruning method
  
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
  
  GBM = list( distribution = "bernoulli",
              n.trees = 2500,	# default
              interaction.depth = 7, # default
              n.minobsinnode = 5, # default
              shrinkage = 0.001, # default, learning rate
              bag.fraction = 0.5, # default, proportion of observations used in selecting variables
              train.fraction = 0.8, # default 1, train.fraction * nrows(data) observations are used to fit the gbm 
              cv.folds = 10,	# default 3
              keep.data = FALSE, # default
              verbose = FALSE,	# default
              perf.method = "cv", # default
              n.cores = 1 ),	# default
  
  CTA = list(	method = "class", # default, response is factor
              parms = "default", # default
              cost = NULL ),	# default
  
  ANN = list(	NbCV = 10, 		# default, number CV
              size = NULL, 	# default, size parameter will be optimised by cross validation based on model AUC
              decay = NULL, 	# default, decay parameter will be optimised by cross validation
              rang = 0.1, 	# default, initial random weights on [-rang, rang] 
              maxit = 200 ), 	# default, maximum number of iterations
  
  SRE = list(quant = 0.025),	# default
  
  FDA = list(	method = "mars",	# default, regression method used in optimal scaling
              add_args = NULL ),# default
  
  RF = list(	do.classif = TRUE, # default classification random.forest computed, else regression random.forest 
             ntree = 500,	# default
             mtry = 10,		# number of variables randomly sampled as candidates at each split
             nodesize = 1,	# default 5, but 1 for classification, minimum size of terminal nodes
             maxnodes = NULL ) # default, maximum number of terminal nodes trees in the forest
)

# models to predict with
mymodels <- c("GLM","GBM","GAM","CTA","ANN", "SRE", "FDA","MARS","RF","MAXENT")

#- - - - - - - - - - - - - - - - - - - - -
## Perform models ####
#- - - - - - - - - - - - - - - - - - - - -

# load environmental variables (for projections)
Env_clip <- terra::rast("_intermediates/EnvPredictor_PCA_1km_POR.tif")
Env_clip <- terra::subset(Env_clip, 1:11) #11 = >80%
# wrap to make loop working
Env_clip <- terra::wrap(Env_clip)

setwd(paste0(input_dir, "/_results/", Taxon_name))
set.seed(32639)

for(spID in species_table$species){ 
  
  # save history
  #sink(paste0("./SDMs/ESM_biomod_", spID, ".txt"))
  print("#################################################")
  print(spID)

  tryCatch({
    load(paste0(input_dir, "/_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData
  }, error = function(e) {print(e); print("FAILED: Load data")}
  )
  
  #- - - - - - - - - - - - - - - - - - - - -
  ### Calibration of simple bivariate models
  tryCatch({
    myBiomodModelOut <- biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
                                                 models = mymodels,
                                                 bm.options = mySDMOption,
                                                 CV.nb.rep = 10,   # 10-fold crossvalidation evaluation run
                                                 CV.perc = 0.8, # use subset of the data for training
                                                 #weights = temp_weights$weight, # weight to observations, here based on year
                                                 metric.eval = "TSS",
                                                 var.import = 0,
                                                 scale.models = FALSE, #scale all predictions with binomial GLM?
                                                 CV.do.full.models = FALSE, # do evaluation & calibration with whole dataset
                                                 modeling.id = paste(spID,"_Modeling", sep = ""))
  }, error = function(e) {print(e); print("FAILED: Build Model")}
  )
  
  #- - - - - - - - - - - - - - - - - - - - -
  ### Ensemble models
  tryCatch({
    # ensemble modeling using mean probability
    myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                                   models.chosen = get_built_models(myBiomodModelOut),  # all working algorithms
                                                   em.by = "all",    #evaluated over evaluation data if given (it is not, see Prepare_input.R)
                                                   # note: evaluation not that important as we will calculate measures on independent data
                                                   metric.select = "TSS", # 'all' would takes same as above in BIOMOD_Modelling
                                                   metric.select.thresh = NULL, # since some species's auc are naturally low
                                                   em.algo = c("EMcv", "EMwmean"),
                                                   EMwmean.decay = "proportional", #the better a model (evaluation score), the higher weight
                                                   var.import = 0)    #number of permutations to estimate variable importance
  }, error = function(e) {print(e); print("FAILED: Ensemble Model")}
  )
  rm(list=ls(pattern="_Modeling"))
  
  #- - - - - - - - - - - - - - - - - - - - -
  ### Evaluation of the ensemble models based on the pooling procedure 
  myBiomodModelEval <- biomod2::get_evaluations(myBiomodEM) %>% 
    dplyr::select(algo:evaluation)
  #myBiomodModelEval$ESM.evaluations
  
  #- - - - - - - - - - - - - - - - - - - - -
  ### Projection of simple bivariate models into new space 
  current <- terra::unwrap(Env_clip)
  tryCatch({
    myBiomodProj <- biomod2::BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                               new.env = current,        #column/variable names have to perfectly match with training
                                               proj.name = "modeling",  #name of the new folder being created
                                               models.chosen = "all", #use all models
                                               metric.binary = NULL,     #binary transformation according to criteria, or no transformation if NULL
                                               compress = TRUE,         #compression format of objects stored on hard drive
                                               build.clamping.mask = TRUE, #TRUE: clamping mask will be saved on hard drive different
                                               do.stack = TRUE,         #save output projections as rasterstack (if not too heavy)
                                               #output.format = ".RData", #what format should projections have: RData, grd or img
                                               keep.in.memory = FALSE)  #FALSE: only story link to copy to projection file
  }, error = function(e) {print(e); print("FAILED: Project model")}
  )
  
  ### Projection of calibrated ESMs into new space 
  tryCatch({
    myBiomodEnProj <- biomod2::BIOMOD_EnsembleForecasting(bm.proj = myBiomodProj,
                                                          bm.em = myBiomodEM,
                                                          #... same arguments as above could be added but are not necessary when loading myBiomodProj
                                                          models.chosen = "all")
  }, error = function(e) {print(e); print("FAILED: Project ensemble")}
  )
  # extract prediction
  temp_prediction <- biomod2::get_predictions(myBiomodEnProj)[[2]] #column with weighted mean
  names(temp_prediction) <- spID
  
  # save prediction in different folder
  tryCatch({
    terra::writeRaster(temp_prediction, file=paste0("./Projection/", spID, ".tif"))
  }, error = function(e) {print(e); print("FAILED: Save projection")}
  )
  
  
  #- - - - - - - - - - - - - - - - - - - - -
  ### Uncertainty
  # extract prediction
  temp_prediction <- biomod2::get_predictions(myBiomodEnProj)[[1]] #column with cv
  names(temp_prediction) <- spID
  
  # save prediction in different folder
  tryCatch({
    terra::writeRaster(temp_prediction, file=paste0("./Uncertainty/CV_", spID, ".tif"))
  }, error = function(e) {print(e); print("FAILED: Save uncertainty")}
  )
  
  # save clamping mask outside folder
  tryCatch({
    temp_files <- list.files(paste0(output_dir, "/", Taxon_name, "/", stringr::str_replace(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
  }, error = function(e) {print(e); print("FAILED: List files [Clamping Mask]")}
  )
  tryCatch({
    myBiomod_Clamping <- terra::rast(temp_files[stringr::str_detect(temp_files,"ClampingMask.tif")])
    names(myBiomod_Clamping) <- spID
  
    terra::writeRaster(myBiomod_Clamping, file=paste0("./Uncertainty/ClampingMask_", spID, ".tif"))
  }, error = function(e) {print(e); print("FAILED: Save Clampling Mask")}
  )
  
  #- - - - - - - - - - - - - - - - - - - - -
  #### Note: varImp not calculated above (number of runs for varImp = 09)
  # ## get the variable contributions of ESMs
  # myBiomod_varImp <- biomod2::get_variables_importance(myBiomodEM)  %>%
  #   filter(algo=="EMwmean") %>% 
  #   dplyr::select(expl.var:var.imp)
  # # average across 10 runs
  # myBiomod_varImp <- myBiomod_varImp %>%
  #   rename(Predictors = expl.var,
  #          varImp = var.imp) %>%
  #   group_by(Predictors) %>%
  #   summarize(across(varImp, .fns=mean))
  
  # ## get the response plots of ESMs
  # my.ESM_responsePlot <- ecospat.ESM.responsePlot(my.ESM_EF,
  #                                                 my.ESM,
  #                                                 fixed.var.metric = 'mean')
   
  #- - - - - - - - - - - - - - - - - - - - -
  ## save model output
  tryCatch({
    biomod_list <- list(validation = myBiomodModelEval #, 
                      #varImp = myBiomod_varImp
                      )
    save(biomod_list, file=paste0("./SDMs/SDM_biomod_", spID, ".RData"))
  }, error = function(e) {print(e); print("FAILED: Save model output")}
  )
  
  #- - - - - - - - - - - - - - - - - - - - -
  ## delete folder & model files
  tryCatch({ 
    unlink(paste0(output_dir, "/", Taxon_name, "/", stringr::str_replace(spID, "_", ".")), recursive = TRUE)
  }, error = function(e) {print(e); print("FAILED: Remove taxa directory")}
  )
  
  gc()
  print("SUCCESSFUL")
  savehistory(file = paste0("./SDMs/.Rhistory_SDM_biomod_", spID))
  
}

# #- - - - - - - - - - - - - - - - - - - - -
# ## Extract names of working algorithms
# models <- tibble(data.frame("SpeciesID"="species",
#                             "GLM" = NA,
#                             "GBM" = NA,
#                             "GAM" = NA,
#                             "CTA" = NA,
#                             "ANN" = NA, 
#                             "SRE" = NA, 
#                             "FDA" = NA,
#                             "MARS" = NA,
#                             "RF" = NA,
#                             "MAXENT" = NA)[0,])
# 
# for(spID in species_table$species){try({
#   
#   # list files in species-specific BIOMOD folder
#   temp_files <- list.files(paste0(output_dir, "/", Taxon_name, "/", stringr::str_replace(spID, "_", ".")), full.names = TRUE)
#   
#   # load model output
#   myBiomodModelOut <- temp_files[stringr::str_detect(temp_files,"Modeling.models.out")]
#   print(myBiomodModelOut)
#   myBiomodModelOut <-get(load(myBiomodModelOut))
#   
#   # get working algorithms used in Ensemble model
#   temp_models <- get_built_models(myBiomodModelOut)
#   
#   # extract simple algorithm names
#   n_name <- nchar(spID)
#   temp_models <- substring(temp_models, n_name+11, n_name+17) %>%
#     stringr::str_replace("_", "") %>%
#     unique()
#   
#   # transform to wide format (each algorithm one column)
#   temp_models <- data.frame("model" = temp_models, "check" = 1) %>% 
#     mutate("SpeciesID" = spID) %>%
#     pivot_wider(id_cols = "SpeciesID", names_from = "model", values_from = "check")
#   
#   models <- models %>% full_join(temp_models)
# })}
# 
# models
# 
# # add number of presences and absences
# records <- read_csv(file=paste0(output_dir, "/Occurrence_rasterized_2km_BIOMOD_", Taxon_name, ".csv"))
# records_species <- records %>% group_by(SpeciesID) %>% summarize(across("occ", sum)) %>%
#   full_join(records %>% filter(occ==0) %>% group_by(SpeciesID) %>% count(name="Absences"))
# records_species
# 
# models <- models %>% full_join(records_species) %>%
#   rename("Presences" = occ) %>% 
#   arrange(SpeciesID)
# models
# 
# write_csv(models, file=paste0(output_dir, "/Models_BIOMOD_", Taxon_name, ".csv"))

