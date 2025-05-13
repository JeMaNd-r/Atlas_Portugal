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

Taxon_name <- "Nematodes"

#- - - - - - - - - - - - - - - - - - - - -
## Prepare data ####
if(Taxon_name == "Earthworms" | Taxon_name == "EarthGenus"){
  mySpeciesOcc <- read_csv(file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))
  mySpeciesOcc <- mySpeciesOcc %>% pivot_wider(id_cols = c(Sample_ID, x, y), names_from = SpeciesID, values_from = Presence, values_fn = max)
  
  data_env <- read_csv(file="_intermediates/SoilReCon_Data_4_23_LocationsEW_PCA.csv")
}else{
  mySpeciesOcc <- read_csv(file=paste0("_intermediates/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
}  

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
          
          if(Taxon_name == "Earthworm" | Taxon_name == "EarthGenus"){
            myExpl <- data_env %>% rename("Sample_ID" = ID) %>% dplyr::select(-x, -y) %>%
              full_join(mySpeciesOcc %>% dplyr::select("Sample_ID", "x", "y", all_of(spID)), by = "Sample_ID")
            myExpl <- myExpl[!is.na(myExpl[,spID]),]
            
            myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myExpl %>% pull(spID) %>% as.numeric(),
                                                          expl.var = myExpl %>% dplyr::select(PC1:PC11),
                                                          resp.xy = myExpl[,c("x", "y")],
                                                          resp.name = spID,
                                                          PA.nb.rep = 0
            )
          }else{
            myResp <- mySpeciesOcc[!is.na(mySpeciesOcc[,spID]), c("x","y",spID)]
            
            if(myResp[myResp[,spID]==0,] %>% nrow() < 
               myResp[myResp[,spID]==1,] %>% nrow()){
              myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp[myResp[,spID]!=0,] %>% pull(spID) %>% as.numeric(),
                                                            expl.var = Env_clip,
                                                            resp.xy = myResp[myResp[,spID]!=0,c("x", "y")],
                                                            resp.name = spID,
                                                            PA.nb.rep = 1,
                                                            PA.nb.absences = 1000, # not needed because true absence data available
                                                            PA.strategy = "random"
              )
            }else{
              myBiomodData <- biomod2::BIOMOD_FormatingData(resp.var = myResp %>% pull(spID) %>% as.numeric(),
                                                            expl.var = Env_clip,
                                                            resp.xy = myResp[,c("x", "y")],
                                                            resp.name = spID,
                                                            PA.nb.rep = 0
                                                            #PA.nb.absences = 10000, # not needed because true absence data available
                                                            #PA.strategy = "random"
                                                            )
            }
          }
          # save data
          save(myBiomodData, file=paste0("_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData"))
          
          rm(myBiomodData, myResp, spID)
          gc()
        })}

#stopImplicitCluster()


#- - - - - - - - - - - - - - - - - - - - -
## Calculate number of records per species ####
records <- lapply(as.list(speciesSub), function(x){
  spID <- x
  myBiomodData <- get(load(file=paste0("_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", x, ".RData")))
  
  # extract occurrences & pseudo-absences
  myData <- cbind(myBiomodData@data.species, myBiomodData@coord, myBiomodData@data.env.var)
  myData$SpeciesID <- spID
  myData <- myData %>% dplyr::rename("occ" = "myBiomodData@data.species")
  myData[is.na(myData$occ),"occ"] <- 0 #replace pseudo (NA) by 0
  
  myData
})

records <- do.call(rbind, records)
head(records)

## same but in for loop... takes >1 day for Bacteria
# for(spID in speciesSub){ try({
#   
#   print(paste0(spID, " will be added."))
#   
#   # load biomod data
#   load(file=paste0("_intermediates/BIOMOD_data/", Taxon_name, "/BiomodData_", Taxon_name,"_", spID, ".RData")) #myBiomodData
# 
#   # extract occurrences & pseudo-absences
#   myData <- cbind(myBiomodData@data.species, myBiomodData@coord, myBiomodData@data.env.var)
#   myData$SpeciesID <- spID
#   myData <- myData %>% rename("occ" = "myBiomodData@data.species")
#   myData[is.na(myData$occ),"occ"] <- 0 #replace pseudo (NA) by 0
#   
#   records <- rbind(records, myData[,c("x", "y","occ", "SpeciesID")])
# })
# }

nrow(records) # Crassiclitellata: 930, nematoda f: 23625; fungi s: 872009, g: 98564; bacteria: 16,641,669; protists: 109049, Eukaryotes: 56446 
nrow(records %>% filter(occ==1)) # 162; Nf: 5218; Fs: 61483, g: 14713; B: 2,383,587; P: 20455; Eu: 8979
nrow(records %>% filter(occ==0)) # 768; Nf: 18417; Fs: 810526, g: 83851; B: 14,258,082; P: 88594, Eu: 47467

records_species <- records %>% group_by(SpeciesID) %>% summarize(across("occ", sum)) %>%
  full_join(records %>% filter(occ==0) %>% group_by(SpeciesID) %>% count(name="Absences")) 
records_species

records_species %>% filter(occ>=10) %>% count() # C: 5 species/ 4 genera, N: 29f / 17g, Fs: 1023, g:157, B: 25653, P: 133, EU: 82
records_species %>% filter(occ>=100) %>% count() # C: 0 species (max. possible occ=92, max. occ=50) / 0 genera, N: 22f / 8g, Fs: 166/g:50, B: 7997, P: 82, Eu: 35

write_csv(records, file=paste0("_results/Occurrence_rasterized_1km_BIOMOD_", Taxon_name, ".csv"))
records <- read_csv(file=paste0("_results/Occurrence_rasterized_1km_BIOMOD_", Taxon_name, ".csv"))

# merge with species list
species_list <- read_csv(paste0("_results/Species_list_", Taxon_name, ".csv"))

species_list <- species_list %>%
  rename(SpeciesID = taxaID) %>%
  full_join(records_species, by = "SpeciesID")
write_csv(species_list, file=paste0("_results/Species_list_", Taxon_name, ".csv"))
species_list <- read_csv(file=paste0("_results/Species_list_", Taxon_name, ".csv"))

# save species IDs for species with >100 BIOMOD presences
write_csv(species_list %>% filter(occ>=100) %>% 
            dplyr::select(SpeciesID) %>% dplyr::rename(species = SpeciesID) %>% unique(),
          file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))

# save IDS for species with >10 and <100 presences
write_csv(species_list %>% filter(occ>=10 & 
                                       occ<100) %>% 
            dplyr::select(SpeciesID) %>% dplyr::rename(species = SpeciesID) %>% unique(),
          file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))


# # check sampling locations
# terra::plot(Env_clip[[1]])
# terra::plot(terra::vect(myResp, geom = c("x", "y")), add = TRUE, pch = 0)
# terra::plot(terra::vect(myData, geom = c("x", "y")), add = TRUE, pch = 15)

