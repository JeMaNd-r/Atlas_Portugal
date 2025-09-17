#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare input for Zonation           #
#          author: Romy Zeiss               #
#            date: 2023-11-02               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal")

gc()
library(tidyverse)
library(terra)
library(sf)

# define names of taxonomic groups
taxaNames <- c("Crassiclitellata", "Fungi", "Nematodes", "Protists", "Eukaryotes", "Bacteria") #
taxaNames

#- - - - - - - - - - - - - - - - - - - - - -
### Save each layer as individual tif but with threshold probability of 500/0.5 ####
#- - - - - - - - - - - - - - - - - - - - - -

dir.create(paste0("_results/Zonation"))
#dir.create(paste0("_results/Zonation/NoData"))
for(i in taxaNames) dir.create(paste0("_results/Zonation/", i))

names_list <- list()
for(Taxon_name in taxaNames){
  print(Taxon_name)
  
  # load number of occurrences per species and focal species names
  try(species100 <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv")))
  if(nrow(species100) != 0){
    species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))
    speciesSub <- species100 %>% 
      full_join(species10) %>%
      pull(species)
  }else{
    speciesSub <- species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv")) %>% pull(species)
  }
  speciesSub
  try(species100 <- species100 %>% pull(species))
  try(species10 <- species10 %>% pull(species))
  
  for(i in c(10, 100)){try({
    print(i)
    
    # list all projections
    projection_rast <- paste0("_results/", Taxon_name, "/Ensembles/", get(paste0("species", i)),".tif")
    projection_rast <- projection_rast[projection_rast %in% paste0("_results/", Taxon_name, "/Ensembles/", list.files(paste0("_results/", Taxon_name, "/Ensembles")))]
    temp_species <- get(paste0("species", i))[str_detect( projection_rast, str_c(get(paste0("species", i)), collapse = "|"))]
    
    # load into list
    projection_rast <- terra::rast(projection_rast)
    names(projection_rast) <- paste0(i, "_", temp_species)
    
    projection_rast[projection_rast<500 & !is.na(projection_rast)] <- 0
    
    terra::writeRaster(projection_rast, paste0("_results/Zonation/", Taxon_name, "/", i, "_", temp_species, ".tif"), overwrite = TRUE)
    
    # number of cells with probability >500
    n_probability <- terra::ncell(projection_rast) - unlist(lapply(terra::cells(projection_rast, 0), length))
    
    temp_list <- data.frame("Taxon" = Taxon_name, 
                            "ID" = paste0(i, "_", temp_species), 
                            "SpeciesID" = temp_species, 
                            "No_cells_probability" = n_probability)
    
    names_list <- c(names_list, list(temp_list))
    
    ## Uncertainty
    uncertainty_rast <- paste0("_results/", Taxon_name, "/Ensembles/", temp_species,".tif")
    uncertainty_rast <- terra::rast(uncertainty_rast)
    names(uncertainty_rast) <- paste0(i, "_", temp_species)
    
    # standardize uncertainty
    unc_min <- global(uncertainty_rast, "min", na.rm = TRUE)[,1]
    unc_max <- global(uncertainty_rast, "max", na.rm = TRUE)[,1]
    
    uncertainty_rast <- (uncertainty_rast - unc_min) / (unc_max - unc_min)
    
    # invert (high uncertainty -> low certainty)
    certainty_rast <- 1 - uncertainty_rast
    terra::writeRaster(certainty_rast, paste0("_results/Zonation/", Taxon_name, "/", i, "_", temp_species, "_certainty.tif"), overwrite = TRUE)
  })
  }
}

# save name coding
names_list <- do.call(rbind, names_list)
colnames(names_list) <- c("Taxon", "ID", "SpeciesID", "No_cells_probability")
names_list <- as_tibble(names_list)
names_list

write_csv(names_list,  paste0("_results/Zonation/TaxaNames_legend.csv"))

# define weight and group for each layer
names_list <- read_csv(paste0("_results/Zonation/TaxaNames_legend.csv"))
write_delim(names_list %>%
              filter(SpeciesID != "Richness") %>%
              mutate("weight" = "1",
                     "group" = as.numeric(as.factor(Taxon)),
                     "uncertainty" = paste0("../data/", Taxon, "/", ID, "_certainty.tif"),
                     "filename" = paste0("../data/", Taxon, "/", ID, ".tif")) %>%
              dplyr::select("weight", "uncertainty", "group","filename"),
            paste0("_results/Zonation/features.txt"))

write_delim(names_list %>%
              filter(SpeciesID != "Richness") %>%
              filter(str_detect(ID, "10_")) %>%
              mutate("weight" = "1",
                     "group" = as.numeric(as.factor(Taxon)),
                     "uncertainty" = paste0("../data/", Taxon, "/", ID, "_certainty.tif"),
                     "filename" = paste0("../data/", Taxon, "/", ID, ".tif")) %>%
              dplyr::select("weight", "uncertainty", "group","filename"),
            paste0("_results/Zonation/features_10.txt"))

write_delim(names_list %>%
              filter(SpeciesID != "Richness") %>%
              filter(str_detect(ID, "100_")) %>%
              mutate("weight" = "1",
                     "group" = as.numeric(as.factor(Taxon)),
                     "uncertainty" = paste0("../data/", Taxon, "/", ID, "_certainty.tif"),
                     "filename" = paste0("../data/", Taxon, "/", ID, ".tif")) %>%
              dplyr::select("weight", "uncertainty", "group","filename"),
            paste0("_results/Zonation/features_100.txt"))

#- - - - - - - - - - - - - - - - - - - - - -
## Prepare protection layers ####

# load one extent
env_por <- terra::rast( "_intermediates/EnvPredictor_1km_POR_normalized.tif")
r <- env_por[[1]]
r$all <- 1
r <- terra::mask(r, env_por[[1]])
r <- terra::subset(r, "all")
rm(env_por)

## load all files provided by protectedplanet.net 
shp0 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_0/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp1 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_1/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp2 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_2/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")

shp0_por <- terra::crop(shp0, ext(r))
shp1_por <- terra::crop(shp1, ext(r))
shp2_por <- terra::crop(shp2, ext(r))

shp_por <- rbind(shp0_por, shp1_por, shp2_por)
terra::plot(shp_por)

# ## merge all files
# shp_all <- rbind(shp0, shp1)
# rm(shp0, shp1)
# 
# shp_all <- rbind(shp_all, shp2)
# rm(shp2)
# 
# shp_all
# 
# ## crop to Portuguese extent
# shp_por <- terra::crop(shp_all, r, snap = "out")
# rm(shp_all)
# gc()

# save
terra::writeVector(shp_por, filename = paste0("_intermediates/Protection_POR.shp"), overwrite = TRUE)

# load protection for portugal
#shp_por <- sf::st_read(paste0("_intermediates/Protection_POR.shp"))
shp_por <- terra::vect(paste0("_intermediates/Protection_POR.shp"))

# split multipolygon by type of protection area
# type of protected area: protect_raw$DESIG_ENG
# IUCN category: IUCN_CAT

shp_por_merged <- terra::aggregate(shp_por, by = "IUCN_CAT", fun = "mean") 

# create empty stack to store PA types as layers
protect_stack <- terra::rast()

# for(i in unique(shp_por$IUCN_CAT)){
#   temp_shp <- shp_por[shp_por$IUCN_CAT == i,]
#   
#   # calculate percent of grid cells covered by PA type i
#   temp_cover <- sum(do.call(c, exactextractr::coverage_fraction(r, temp_shp)))
#   names(temp_cover) <- i
#   
#   # replace values higher than 1 by 1
#   temp_cover[temp_cover > 1] <- 1
#   temp_cover[temp_cover < 0] <- 0
#   
#   protect_stack <- c(protect_stack, temp_cover)
#   
#   print(paste0(i, " ready."))
# }

for(i in unique(shp_por_merged$IUCN_CAT)){
  temp_shp <- shp_por_merged[shp_por_merged$IUCN_CAT == i,]
 
  temp_rast <- rasterize(temp_shp, r,
                                field = 1,
                                touches = TRUE,
                                background = 0,
                                cover=TRUE)
  names(temp_rast) <- i
  
  protect_stack <- c(protect_stack, temp_rast)

  print(paste0(i, " ready."))
}

protect_stack

# calculate sum for IUCN categories Ia-VI
protect_stack[protect_stack==0,] <- NA
protect_stack$PA_coverage <- mean(protect_stack, na.rm=TRUE) #not assigned and not reported are overlapping
protect_stack$IUCN_coverage <- sum(subset(protect_stack, names(protect_stack) %in% c("Ia", "Ib", "II", "III", "IV", "V", "VI")), na.rm=TRUE)
protect_stack$PA <- ifelse(values(protect_stack$PA_coverage) > 0.5 & !is.na(values(protect_stack$PA_coverage)), 1, 0)
protect_stack$IUCN_PA <- ifelse(values(protect_stack$IUCN_coverage) > 0.5 & !is.na(values(protect_stack$IUCN_coverage)), 1, 0)

# # add areas with no IUCN category reported
# protect_stack$nonIUCN <- sum(do.call(c, exactextractr::coverage_fraction(r, shp_por[shp_por$PA_DEF==1,])))>1
# protect_stack$nonIUCN[protect_stack$nonIUCN > 1] <- 1
# protect_stack$nonIUCN[protect_stack$nonIUCN < 0] <- 0

protect_stack <- mask(protect_stack, r)

# save
terra::writeRaster(protect_stack, filename = paste0("_intermediates/Protection_POR_coverage.tif"), overwrite = TRUE)

# split layers
protect_stack <- terra::rast(paste0("_intermediates/Protection_POR_coverage.tif"))
#try(dir.create(paste0("_results/Zonation/protection")))

# save layers for Zonation
for(ly in names(protect_stack)[names(protect_stack) %in% c("PA", "IUCN_PA")]){
  protect_stack[[ly]]
  terra::writeRaster(protect_stack[[ly]], 
                     filename = paste0("_results/Zonation/POR_", 
                                       ly, 
                                       ".tif"), 
                     datatype = "INT1U",  # 8-bit unsigned integer
                     overwrite = TRUE)
  print(paste0(ly, " ready."))
}

# features_pa_category <- data.frame("ID" = list.files(paste0("_results/Zonation/protection"), pattern = "^POR_"))
# features_pa_category <- features_pa_category %>%
#   filter(!(str_detect(ID, "PA"))) %>%
#   mutate("weight" = "1", 
#          "group" = "Protection") %>%
#   dplyr::select(weight, ID, group)
# 
# write_delim(features_pa_category,
#             paste0("_results/Zonation/features_pa_category.txt"))

#- - - - - - - - - - - - - - - - - - - - - -
## Threat layers ####

dir.create(paste0("_results/Zonation/threats"))

degr_drivers <- terra::rast(paste0("../Soil_degradation_drivers/Degradation_POR_binary.tif"))

# save layers for Zonation
for(ly in names(degr_drivers)){
  terra::writeRaster(degr_drivers[[ly]], 
                     filename = paste0("_results/Zonation/threats/POR_", 
                                       ly, 
                                       ".tif"))
  print(paste0(ly, " ready."))
}

costs <- data.frame("cost_file" = list.files(paste0("_results/Zonation/threats/"), pattern = "^POR_"))
costs <- costs %>%
  mutate("weight" = c("3", "2", "1", "2", "1", "2"),
         "cost_file" = paste0("data/", cost_file)) %>%
  dplyr::select(weight, cost_file)

write_delim(costs,
            paste0("_results/Zonation/costs.txt"))





#- - - - - - - - - - - - - - - - - - - - - -
## Test if all same extent ####

rast1 <- terra::rast(paste0("_results/Zonation/Fungi/Fungi_3.tif"))
rast2 <- terra::rast(paste0("_results/Zonation/Fungi/Fungi_Rich.tif"))
rast3 <- terra::rast(paste0("_results/Zonation/Nematodes/Nemat_5.tif"))
rast4 <- terra::rast(paste0("_results/Zonation/Nematodes/Nemat_20.tif"))
rast5 <- terra::rast(paste0("_intermediates/Protection_POR_coverage.tif"))

terra::compareGeom(rast1, rast2)
terra::compareGeom(rast1, rast3)
terra::compareGeom(rast3, rast4)
terra::compareGeom(rast1, rast5)

par(mfrow = c(2,2))
terra::plot(rast1)
terra::plot(rast2)
terra::plot(rast3)
terra::plot(rast5$PA_coverage)
