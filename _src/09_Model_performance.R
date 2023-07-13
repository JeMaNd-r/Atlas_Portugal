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

library(biomod2) # also to create pseudo-absences

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

# load number of occurrences per species and focal species names
speciesSub <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv")) %>% pull(SpeciesID)
speciesSub

#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

data_eval <- data.frame("SpeciesID"="species", "Kappa"=0, "TSS"=0, "ROC"=0)

# for loop through all species
for(spID in speciesSub){ try({
  
  ## Load probability maps 
  load(file=paste0("_results/", Taxon_name, "/SDM_biomod_", spID, ".RData")) #biomod_list
  temp_validation <- biomod_list$validation[,"validation"]
  
  data_eval <- rbind(data_eval, c(spID, temp_validation))
  
  print(paste0(spID, " successfully loaded."))
}, silent=T)}

data_eval <- data_eval %>% filter(SpeciesID!="species") 
data_eval$Kappa <- as.numeric(data_eval$Kappa)
data_eval$ROC <- as.numeric(data_eval$ROC)
data_eval$TSS <- as.numeric(data_eval$TSS)

data_eval; str(data_eval)
summary(data_eval)

write.csv(data_eval, paste0(data_wd, "/_results/Model_evaluation_", Taxon_name, ".csv"), row.names = F)


