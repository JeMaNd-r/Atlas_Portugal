#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#     Ensembles of Small Models - full      #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

for(Taxon_name in c("Bacteria")){ #"Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes", 
  print(Taxon_name)
  
  # Note: dismo::maxent() might crash in RStudio
  
  setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal")
  species_csv <- paste0("ESM_", Taxon_name, ".csv") 
  output_dir <- paste0(getwd(), "/_results")
  input_dir <- getwd()
  
  library(biomod2)
  
  gc()
  #library(tidyverse)
  library(dplyr)
  library(readr)
  
  # read functions for ensemble of small models (ESM)
  #devtools::install_github("ecospat/ecospat/ecospat@dev")
  library(ecospat)
  
  library(PresenceAbsence) # needed for ecospat
  
  # library(raster)
  
  library(parallel)
  library(doParallel)
  
  # save history
  #devtools::install_github("leeper/rite")
  library("rite")
  
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
  
  # check if MAXENT working
  dismo::maxent()
  
  #- - - - - - - - - - - - - - - - - - - - -
  ## Perform models ####
  #- - - - - - - - - - - - - - - - - - - - -
  
  # load environmental variables (for projections)
  Env_clip <- terra::rast("_intermediates/EnvPredictor_PCA_1km_POR.tif")
  Env_clip <- terra::subset(Env_clip, 1:11) #11 = >80%
  # wrap to make loop working
  Env_clip <- terra::wrap(Env_clip)
  
  setwd(paste0(input_dir, "/_results/", Taxon_name))
  file.copy(from = paste0(input_dir, "/_results/maxent.jar"), 
            to = paste0(input_dir, "/_results/", Taxon_name),
            overwrite = FALSE,
            copy.mode = TRUE)
  set.seed(32639)
  
  if(Taxon_name == "Bacteria") species_table <- species_table[87:277,]
  
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
      my.ESM <- ecospat.ESM.Modeling(data = myBiomodData,
                                    models = c("GLM", "MAXENT"), #,"ANN"
                                    Prevalence = NULL,
                                    tune = TRUE, # TRUE: estimate optimal parameters for the models
                                    NbRunEval = 10,
                                    DataSplit = 80,
                                    weighting.score = c("TSS"),
                                    parallel = TRUE)  
    }, error = function(e) {print(e); print("FAILED: Build Model")}
    )
    
    #- - - - - - - - - - - - - - - - - - - - -
    ### Ensemble models
    tryCatch({
      my.ESM_EF <- ecospat.ESM.EnsembleModeling(my.ESM,
                                              weighting.score = c("TSS"),
                                              threshold=0)
    }, error = function(e) {print(e); print("FAILED: Ensemble Model")}
    )
    
    ### thresholds to produce binary maps
    my.ESM_thresholds <- ecospat.ESM.threshold(my.ESM_EF)
    
    #- - - - - - - - - - - - - - - - - - - - -
    # ### Evaluation of bivariate and ensemble models based on standard cross-validation
    # my.ESM_EF$ESM.evaluations
    # my.ESM_thresholds
    
    ### Evaluation of the ensemble models based on the pooling procedure 
    my.ESM_evaluations <- ecospat.ESM.EnsembleEvaluation(ESM.modeling.output = my.ESM,
                                                         ESM.EnsembleModeling.output = my.ESM_EF,
                                                         metrics = c("AUC","MaxTSS"),
                                                         EachSmallModels = FALSE)
    #my.ESM_evaluations$ESM.evaluations
    
    #- - - - - - - - - - - - - - - - - - - - -
    ### Projection of simple bivariate models into new space 
    current <- terra::unwrap(Env_clip)
    tryCatch({
      my.ESM_proj <- ecospat.ESM.Projection(ESM.modeling.output = my.ESM,
                                                new.env = current)
    }, error = function(e) {print(e); print("FAILED: Project model")}
    )
    
    ### Projection of calibrated ESMs into new space 
    tryCatch({
      my.ESM_EFproj <- ecospat.ESM.EnsembleProjection(ESM.prediction.output=my.ESM_proj,
                                                            ESM.EnsembleModeling.output=my.ESM_EF)
    }, error = function(e) {print(e); print("FAILED: Project ensemble")}
    )
    # extract prediction
    temp_prediction <- my.ESM_EFproj[[2]] #column with weighted mean
    names(temp_prediction) <- spID
    
    # save prediction in different folder
    tryCatch({
      terra::writeRaster(temp_prediction, file=paste0("./Projection/", spID, ".tif"), overwrite = TRUE)
    }, error = function(e) {print(e); print("FAILED: Save projection")}
    )
    
    
    #- - - - - - - - - - - - - - - - - - - - -
    ### Uncertainty
    # extract prediction
    temp_prediction <- my.ESM_EFproj[[1]] #column with cv
    names(temp_prediction) <- spID
    
    # save prediction in different folder
    tryCatch({
      terra::writeRaster(temp_prediction, file=paste0("./Uncertainty/CV_", spID, ".tif"), overwrite = TRUE)
    }, error = function(e) {print(e); print("FAILED: Save uncertainty")}
    )
    
    ## NO clamping mask for ESMs?
    # # save clamping mask outside folder
    # tryCatch({
    #   temp_files <- list.files(paste0(data_wd, "/_results/", Taxon_name, "/", stringr::str_replace(spID, "_", "."), "/proj_modeling"), full.names = TRUE)
    # }, error = function(e) {print(e); print("FAILED: List files [Clamping Mask]")}
    # )
    # my.ESM_Clamping <- terra::rast(temp_files[stringr::str_detect(temp_files,"ClampingMask.tif")])
    # names(my.ESM_Clamping) <- spID
    # tryCatch({
    #   terra::writeRaster(my.ESM_Clamping, file=paste0("./Uncertainty/ClampingMask_", spID, ".tif"))
    # }, error = function(e) {print(e); print("FAILED: Save Clampling Mask")}
    # )
    
    #- - - - - - - - - - - - - - - - - - - - -
    ## get the variable contributions of ESMs
    tryCatch({
      my.ESM_varImp <- ecospat.ESM.VarContrib(my.ESM, my.ESM_EF)                                                      
    }, error = function(e) {print(e); print("FAILED: Variable contribution")}
    )
    
    # ## get the response plots of ESMs
    # tryCatch({
    #   my.ESM_responsePlot <- ecospat.ESM.responsePlot(my.ESM_EF,
    #                                                 my.ESM,
    #                                                 fixed.var.metric = 'mean')
    # }, error = function(e) {print(e); print("FAILED: Response plot")}
    # )
  
    #- - - - - - - - - - - - - - - - - - - - -
    ## save model output
    tryCatch({
      biomod_list <- list(validation = my.ESM_evaluations$ESM.evaluations, 
                        varImp = my.ESM_varImp,
                        thresholds = ecospat.ESM.threshold(my.ESM_EF) #,
                        #varResponse = my.ESM_responsePlot
                        )
      save(biomod_list, file=paste0("./SDMs/ESM_biomod_", spID, ".RData"))
    }, error = function(e) {print(e); print("FAILED: Save model output")}
    )
    
    #- - - - - - - - - - - - - - - - - - - - -
    ## delete folder & model files
    tryCatch({ 
      unlink(paste0(output_dir, "/", Taxon_name, "/ESM.BIOMOD.output_", stringr::str_replace(spID, "_", ".")), recursive = TRUE)
    }, error = function(e) {print(e); print("FAILED: Remove taxa directory")}
    )
    
    gc()
    print("SUCCESSFUL")
    savehistory(file = paste0("./SDMs/.Rhistory_ESM_biomod_", spID))
    
  }
  
  # remove maxent.jar file from folder
  file.remove(paste0(input_dir, "/_results/", Taxon_name, "/maxent.jar"))
}  
