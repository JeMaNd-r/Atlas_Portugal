#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare input data for SDM/ESMs      #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

gc()
library(tidyverse)

library(biomod2) # also to create pseudo-absences

library(parallel)
library(doParallel)

# load environmental variables
Env_clip <- terra::rast("_intermediates/EnvPredictor_PCA_1km_POR.tif")
Env_clip <- terra::subset(Env_clip, 1:11) #11 = >80%

Taxon_name <- "nematodes"

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
mySpeciesOcc <- read_csv(file=paste0("_intermediates/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
speciesSub <- colnames(mySpeciesOcc %>% dplyr::select(-x, -y))
speciesSub

# create subfolder if not existing
if(!dir.exists(paste0("_intermediates/BIOMOD_data/", Taxon_name))){ 
  dir.create(paste0("_intermediates/BIOMOD_data/", Taxon_name))
}

# registerDoParallel(3)
# foreach(spID = speciesSub,
#         .export = c("Env_clip", "mySpeciesOcc"),
#         .packages = c("tidyverse","biomod2")) %dopar% { 
for(spID in speciesSub){
          try({
          
          # in case some spID run into memory allocation error, do only missing 
          # ones with 1 instead of 3 cores running in parallel
          #if(file.exists(paste0(data_wd, "/_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData"))==FALSE){
          
          myResp <- mySpeciesOcc[!is.na(mySpeciesOcc[,spID]), c("x","y",spID)]
          
          myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp %>% pull(spID) %>% as.numeric(),
                                                        expl.var = Env_clip,
                                                        resp.xy = myResp[,c("x", "y")],
                                                        resp.name = spID,
                                                        PA.nb.rep = 0,
                                                        #PA.nb.absences = 10000,
                                                        #PA.strategy = "random"
                                                        )
          
          # save data
          save(myBiomodData, file=paste0("_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData"))
          
          rm(myBiomodData, myResp, spID)
          gc()
        })}

#stopImplicitCluster()

#- - - - - - - - - - - - - - - - - - - - -
## Calculate number of records per species ####
records <- data.frame("x"=12, "y"=12,"occ"=1, "SpeciesID"="species")[0,]

for(spID in speciesSub){ try({
  
  print(paste0(spID, " will be added."))
  
  # load biomod data
  load(file=paste0("_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData

  # extract occurrences & pseudo-absences
  myData <- cbind(myBiomodData@data.species, myBiomodData@coord, myBiomodData@data.env.var)
  myData$SpeciesID <- spID
  myData <- myData %>% rename("occ" = "myBiomodData@data.species")
  myData[is.na(myData$occ),"occ"] <- 0
  
  records <- rbind(records, myData[,c("x", "y","occ", "SpeciesID")])
})
}

head(records)
nrow(records) # Crassiclitellata: 2,093, nematoda: 17,211
nrow(records %>% filter(occ==1)) # 171; N: 5630
nrow(records %>% filter(occ==0)) # 1,922; N: 11,581

records_species <- records %>% group_by(SpeciesID) %>% summarize(across("occ", sum)) %>%
  full_join(records %>% filter(occ==0) %>% group_by(SpeciesID) %>% count(name="Absences"))
records_species

records_species %>% filter(occ>=10) %>% count() #5 species
records_species %>% filter(occ>=100) %>% count() #0 species (max. 92)

write_csv(records, file=paste0("_results/Occurrence_rasterized_1km_BIOMOD_", Taxon_name, ".csv"))

## Prepare species list
species <- tibble(SpeciesID = speciesSub)
species #23; N: 44

write_csv(species, paste0("_intermediates/SDM_", Taxon_name, ".csv"))

