#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Build SDMs with BIOMOD           #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

arguments <- commandArgs(trailingOnly = TRUE)
Taxon_name <- arguments[1]
species_csv <- arguments[2]
output_dir <- arguments[3]
input_dir <- "/data/soilbio"

library(tidyverse)

library(raster)

library(biomod2) # also to create pseudo-absences

library(parallel)
library(doParallel)

#- - - - - - - - - - - - - - - - - - - - -
# load number of occurrences per species and focal species names
species_table <- read_csv(species_csv)
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) 
species <- species_table$SpeciesID[task]

# define formula for GLM (and biomod)
form <- paste0("occ ~ ", paste0(paste0("s(", covarsNames, ")"), collapse=" + "))

#- - - - - - - - - - - - - - - - - - - - -
## Build models ####
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
  
  MAXENT = list(path_to_maxent.jar = paste0(output_dir, "/_results"), # change it to maxent directory
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
## Run SDM ####  

load(paste0(input_dir, "/_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData

# model fitting
#tmp <- proc.time()[3]
setwd(paste0(input_dir, "/_results/", Taxon_name))

set.seed(32639)
myBiomodModelOut <- biomod2::BIOMOD_Modeling(bm.format = myBiomodData,
                                             models = mymodels,
                                             bm.options = mySDMOption,
                                             nb.rep = 10,   # 10-fold crossvalidation evaluation run
                                             data.split.perc = 80, # use subset of the data for training
                                             weights = temp_weights$weight, # weight to observations, here based on year
                                             metric.eval = "TSS",
                                             save.output = TRUE, #save output on hard drive?
                                             scale.models = FALSE, #scale all predictions with binomial GLM?
                                             do.full.models = FALSE, # do evaluation & calibration with whole dataset
                                             modeling.id = paste(spID,"_Modeling", sep = ""))

# ensemble modeling using mean probability
myBiomodEM <- biomod2::BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                               models.chosen = "all",  # all algorithms
                                               em.by = "all",    #evaluated over evaluation data if given (it is not, see Prepare_input.R)
                                               # note: evaluation not that important as we will calculate measures on independent data
                                               metric.select = "TSS", # 'all' would takes same as above in BIOMOD_Modelling
                                               metric.select.thresh = NULL, # since some species's auc are naturally low
                                               #prob.mean = FALSE, #estimate mean probabilities across predictions
                                               prob.cv = TRUE,   #estimate coefficient of variation across predictions
                                               #prob.ci = FALSE,  #estimate confidence interval around the prob.mean
                                               #prob.median = FALSE, #estimate the median of probabilities
                                               #committee.averaging = TRUE, #estimate committee averaging across predictions
                                               prob.mean.weight = TRUE, #estimate weighted sum of predictions
                                               EMwmean.decay = "proportional", #the better a model (evaluation score), the higher weight
                                               var.import = 5)    #number of permutations to estimate variable importance
#temp_model_time <- proc.time()[3] - tmp

setwd(input_dir)
