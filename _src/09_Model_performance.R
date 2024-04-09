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

Taxon_name <- "Nematodes"

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
list_eval <- list.files(paste0(output_dir, "/", Taxon_name, "/SDMs"), full.names = TRUE)

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

if(Taxon_name == "Nematodes") data_eval$Species <- gsub(data_eval$Species, pattern="[.]", replacement="") 

write_csv(data_eval %>% dplyr::select(-Model), 
          paste0("_results/Model_evaluation_", Taxon_name, ".csv"))

#- - - - - - - - - - - - - - - - - - - - -
### similar step but in loop
# # for loop through all species
# for(spID in speciesSub){ try({
#   
#   file_name <- list.files(paste0("_results/", Taxon_name, "/SDMs"), full.names = TRUE,
#                           pattern = paste0("_biomod_", spID, ".RData"))
#     
#   ## Load probability maps 
#   load(file=file_name) #biomod_list
#   temp_validation <- biomod_list$validation %>% 
#     as_tibble() %>% 
#     mutate(Model = rownames(biomod_list$validation),
#            SpeciesID = spID) %>%
#     dplyr::select(SpeciesID, Model, AUC, MaxTSS)
#   
#   data_eval <- rbind(data_eval, temp_validation)
#   
#   print(paste0(spID, " successfully loaded."))
# }, silent=T)}
# 
# data_eval; str(data_eval)
# summary(data_eval)
# 
# write_csv(data_eval, paste0("_results/Model_evaluation_", Taxon_name, ".csv"))


