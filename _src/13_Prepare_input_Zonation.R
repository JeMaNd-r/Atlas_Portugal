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
dir.create(paste0("_results/Zonation/data"))
for(i in taxaNames) dir.create(paste0("_results/Zonation/data/", i))

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

    terra::writeRaster(projection_rast, paste0("_results/Zonation/data/", Taxon_name, "/", i, "_", temp_species, ".tif"), overwrite = TRUE)

    # number of cells with probability >500
    n_probability <- terra::ncell(projection_rast) - unlist(lapply(terra::cells(projection_rast, 0), length))

    temp_list <- data.frame("Taxon" = Taxon_name,
                            "ID" = paste0(i, "_", temp_species),
                            "SpeciesID" = temp_species,
                            "No_cells_probability" = n_probability)

    names_list <- c(names_list, list(temp_list))
    
    ## Uncertainty
    uncertainty_rast <- paste0("_results/", Taxon_name, "/Uncertainty/CV_", temp_species,".tif")
    uncertainty_rast <- terra::rast(uncertainty_rast)
    names(uncertainty_rast) <- paste0(i, "_", temp_species)
    
    # standardize uncertainty
    unc_min <- global(uncertainty_rast, "min", na.rm = TRUE)[,1]
    unc_max <- global(uncertainty_rast, "max", na.rm = TRUE)[,1]
    
    uncertainty_rast <- (uncertainty_rast - unc_min) / (unc_max - unc_min)
    
    # invert (high uncertainty -> low certainty)
    certainty_rast <- 1 - uncertainty_rast
    terra::writeRaster(certainty_rast, paste0("_results/Zonation/data/", Taxon_name, "/", i, "_", temp_species, "_certainty.tif"), overwrite = TRUE)
  })
  }
}

# save name coding
names_list <- do.call(rbind, names_list)
colnames(names_list) <- c("Taxon", "ID", "SpeciesID", "No_cells_probability")
names_list <- as_tibble(names_list)
names_list

# numerical groups
names_list <- names_list %>% mutate("wgrp" = as.numeric(as.factor(Taxon)))

write_csv(names_list,  paste0("_results/Zonation/TaxaNames_legend.csv"))

# define weight and group for each layer
names_list <- read_csv(paste0("_results/Zonation/TaxaNames_legend.csv"))
write_delim(names_list %>%
              filter(SpeciesID != "Richness") %>%
              mutate("weight" = "1",
                     "group" = wgrp,
                     "uncertainty" = paste0("../data/", Taxon, "/", ID, "_certainty.tif"),
                     "filename" = paste0("../data/", Taxon, "/", ID, ".tif")) %>%
              dplyr::select("weight", "uncertainty", "wgrp","filename", "group"),
            paste0("_results/Zonation/features.txt"))

write_delim(names_list %>%
              filter(SpeciesID != "Richness") %>%
              filter(str_detect(ID, "10_")) %>%
              mutate("weight" = "1",
                     "group" = wgrp,
                     "uncertainty" = paste0("../data/", Taxon, "/", ID, "_certainty.tif"),
                     "filename" = paste0("../data/", Taxon, "/", ID, ".tif")) %>%
              dplyr::select("weight", "uncertainty", "wgrp","filename", "group"),
            paste0("_results/Zonation/features_10.txt"))

write_delim(names_list %>%
              filter(SpeciesID != "Richness") %>%
              filter(str_detect(ID, "100_")) %>%
              mutate("weight" = "1",
                     "group" = wgrp,
                     "uncertainty" = paste0("../data/", Taxon, "/", ID, "_certainty.tif"),
                     "filename" = paste0("../data/", Taxon, "/", ID, ".tif")) %>%
              dplyr::select("weight", "uncertainty", "wgrp","filename", "group"),
            paste0("_results/Zonation/features_100.txt"))

# specify group weights to weight all taxa equally across groups
names_list <- read_csv(paste0("_results/Zonation/TaxaNames_legend.csv"))

f_group_weights <- function(x){
  x$max <- max(x$n)
  x$weight <- x$max / x$n
  x <- x %>%
    dplyr::select(wgrp, weight)
}

group_weights <- f_group_weights(names_list %>% 
                                   group_by(wgrp) %>% 
                                   count() )

group_weights_10 <- f_group_weights(names_list %>% 
          filter(str_detect(ID, "10_")) %>% 
          group_by(wgrp) %>% count())

group_weights_100 <- f_group_weights(names_list %>% 
          filter(str_detect(ID, "100_")) %>% 
          group_by(wgrp) %>% count())

group_weights
group_weights_10
group_weights_100

write_delim(group_weights,
            col_names = FALSE,
            paste0("_results/Zonation/group_weights.txt"))
write_delim(group_weights_10,
            col_names = FALSE,
            paste0("_results/Zonation/group_weights_10.txt"))
write_delim(group_weights_100,
            col_names = FALSE,
            paste0("_results/Zonation/group_weights_100.txt"))

#- - - - - - - - - - - - - - - - - - - - - -
## Prepare protection layers ####

# load one extent
env_por <- terra::rast( "_intermediates/EnvPredictor_1km_POR_normalized.tif")
r <- env_por[[1]]
r$all <- 1
r <- terra::mask(r, env_por[[1]])
r <- terra::subset(r, "all")
rm(env_por)

terra::writeRaster(r, file = "_results/Zonation/Analysis_area.tif",
                   datatype = "INT1U",
                   overwrite = TRUE)

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

# calculate layer with distance to PA to rank areas closer to PA higher
protect_layers <- list.files("_results/Zonation/", pattern = "POR_", full.names = TRUE)
protect_layers <- protect_layers[!(str_detect(protect_layers, "distance"))]
protect_layers <- terra::rast(protect_layers)

protect_layers[protect_layers==0,] <- NA #disance is calculated for all NAs

protect_dist_layers <- lapply(1:nlyr(protect_layers), function(lyr){
  temp_dist <- terra::distance(protect_layers[[lyr]]) #distance given in meters
  # # rescale 0–1
  # temp_min <- as.numeric(global(temp_dist, "min", na.rm=TRUE)[1])
  # temp_max <- as.numeric(global(temp_dist, "max", na.rm=TRUE)[1])
  # temp_dist[!is.na(temp_dist)] <- (temp_dist[!is.na(temp_dist),] - temp_min) / (temp_max - temp_min)
  # temp_dist[!is.na(temp_dist)]  <- (1 - temp_dist[!is.na(temp_dist)]) 
  temp_dist
})

protect_dist_layers <- do.call(c, protect_dist_layers)

terra::writeRaster(protect_dist_layers, file = "_intermediates/POR_protection_distance.tif", overwrite = TRUE)
protect_dist_layers <- terra::rast("_intermediates/POR_protection_distance.tif")

# # rescale to categories 
# max_distance <- max(global(do.call(c,protect_dist_layers), max, na.rm=TRUE))
# m_dist <- matrix(c(
#   5000,    max_distance,  1,
#   4000,  5000,  1,
#   3000,  4000,  1.2,
#   2000,  3000,  1.4,
#   1000,  2000,  1.6,
#   0.1,  1000,  1.8,
#   # 0.1, 5000, 1.7,
#   -0.1, 0.1, 2
# ), ncol=3, byrow=TRUE)
# 
# protect_dist_layers <- lapply(1:nlyr(protect_layers), function(lyr)
#   classify(protect_dist_layers[[lyr]], m_dist)
#   )
# 
# protect_dist_layers <- do.call(c, protect_dist_layers)

# save distance layers for Zonation
for(ly in names(protect_dist_layers)){
  protect_dist_layers[[ly]]
  terra::writeRaster(protect_dist_layers[[ly]], 
                     filename = paste0("_results/Zonation/POR_distance_", 
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

# dir.create(paste0("_results/Zonation/threats"))

degr_drivers <- terra::rast(paste0("../Soil_degradation_drivers/Degradation_POR_binary.tif"))

degr_drivers$erosion <- 1 * degr_drivers$erosion 
degr_drivers$climvelo <- 3 * degr_drivers$climvelo
degr_drivers$compaction <- 2 * degr_drivers$compaction
degr_drivers$invasion <- 2 * degr_drivers$invasion
degr_drivers$lu_instability <- 1 * degr_drivers$lu_instability
degr_drivers$salinization <- 2 * degr_drivers$salinization

degr_drivers_sum <- sum(degr_drivers, na.rm = TRUE)
terra::plot(degr_drivers_sum)

#degr_drivers_sum[is.na(degr_drivers_sum)] <- 0

terra::writeRaster(degr_drivers_sum, file = "_results/Zonation/Degradation_sum.tif", overwrite = TRUE)
terra::writeRaster(degr_drivers_sum * 10, file = "_results/Zonation/Degradation_sumx10.tif", overwrite = TRUE)

degr_drivers_max <- max(degr_drivers, na.rm = TRUE)
terra::plot(degr_drivers_max)
#degr_drivers_max[is.na(degr_drivers_max)] <- 0
terra::writeRaster(degr_drivers_max, file = "_results/Zonation/Degradation_max.tif", overwrite = TRUE)
terra::writeRaster(degr_drivers_max/3, file = "_results/Zonation/Degradation_maxdiv3.tif", overwrite = TRUE)

degr_drivers_subset <- terra::subset(degr_drivers, c("erosion", "compaction", "salinization"))
degr_drivers_subset_max <- max(degr_drivers_subset, na.rm = TRUE)
terra::plot(degr_drivers_subset_max)
#degr_drivers_subset_max[is.na(degr_drivers_subset_max)] <- 0

terra::writeRaster(degr_drivers_subset_max, file = "_results/Zonation/Degradation_maxSubset.tif", overwrite = TRUE)

degr_drivers_subset_sum <- sum(degr_drivers_subset, na.rm = TRUE)
terra::plot(degr_drivers_subset_sum)
#degr_drivers_subset_sum[is.na(degr_drivers_subset_sum)] <- 0
terra::writeRaster(degr_drivers_subset_sum, file = "_results/Zonation/Degradation_sumSubset.tif", overwrite = TRUE)

terra::plot(c(degr_drivers_sum, degr_drivers_max, degr_drivers_subset_max), nc=3)


# # save layers for Zonation
# for(ly in names(degr_drivers)){
#   terra::writeRaster(degr_drivers[[ly]], 
#                      filename = paste0("_results/Zonation/threats/POR_", 
#                                        ly, 
#                                        ".tif"))
#   print(paste0(ly, " ready."))
# }
# 
# costs <- data.frame("cost_file" = list.files(paste0("_results/Zonation/threats/"), pattern = "^POR_"))
# costs <- costs %>%
#   mutate("weight" = c("3", "2", "1", "2", "1", "2"),
#          "cost_file" = paste0("data/", cost_file)) %>%
#   dplyr::select(weight, cost_file)
# 
# write_delim(costs,
#             paste0("_results/Zonation/costs.txt"))

#- - - - - - - - - - - - - - - - - - - - - -
## Write scripts for Zonation analysis ####

zonation_dir <- "D:/EIE_Macroecology/_students/Romy/SoilBioPrio_Zonation/POR"

approaches <- c("target", "complement", "prevent")
species_numbers <- c("", "_10", "_100")
group_weights <- c("", "_grW")
protected_areas <- c("", "_allPA")
degradation_weights <- c("_max", "_maxSubset")

zonation_scenarios <- c(
  apply(expand.grid(
  approaches[1], species_numbers, group_weights), 1, paste, collapse = ""),
  apply(expand.grid(
    approaches[2], species_numbers, group_weights, protected_areas), 1, paste, collapse = ""), #c("", "_dist") = distance of PAs considering for expanding
  apply(expand.grid(
    approaches[3], species_numbers, group_weights, protected_areas, degradation_weights), 1, paste, collapse = "")
)
zonation_scenarios

# create directories with correct files
lapply(paste0(zonation_dir, "/", zonation_scenarios), function(x) dir.create(x))

# prepare settings files
for(zonation_scenario in zonation_scenarios){
  approach <- na.omit(str_extract(zonation_scenario, approaches))[1]
  
  # write settings file
  sink(paste0(zonation_dir, "/", zonation_scenario, "/settings_", approach,".z5"))
    if(str_detect(zonation_scenario, "_10")){
      if(str_detect(zonation_scenario, "_100")){
        cat(paste0("feature list file = features_100.txt"), "\n")
        if(str_detect(zonation_scenario, "grW")) cat(paste0("weight groups file = group_weights_100.txt"), "\n")
      } else {
        cat(paste0("feature list file = features_10.txt"), "\n")
        if(str_detect(zonation_scenario, "grW")) cat(paste0("weight groups file = group_weights_10.txt"), "\n")
      }
    } else {
      cat(paste0("feature list file = features.txt"), "\n")
      if(str_detect(zonation_scenario, "grW")) cat(paste0("weight groups file = group_weights.txt"), "\n")
    }
    cat("zero mode = \"all\"", "\n")
    if(approach != approaches[1]){
      if(str_detect(zonation_scenario, "allPA")) { 
         cat("hierarchic mask layer = ../data/POR_PA.tif", "\n")
      } else {
        cat("hierarchic mask layer = ../data/POR_IUCN_PA.tif", "\n")
      }
    }
    if(approach == approaches[3]){
      if(str_detect(zonation_scenario, "maxSubset")){
        cat("cost layer = ../data/Degradation_maxSubset.tif", "\n")
      } else {
        cat("cost layer = ../data/Degradation_max.tif", "\n")
      }
    }
    cat("analysis area mask layer = ../data/Analysis_area.tif", "\n")
  sink()
  
  print(paste0("Successfully wrote ", zonation_dir, "/", zonation_scenario, "/settings_", approach,".z5"))
}
sink()

# copy species layers
file.copy(paste0("_results/Zonation/data"),
          paste0(zonation_dir), recursive = TRUE)
# copy analysis file
file.copy(paste0("_results/Zonation/Analysis_area.tif"),
          paste0(zonation_dir, "/data/Analysis_area.tif"))

temp_files <- list.files(paste0("_results/Zonation/"), full.names = TRUE)
temp_files <- temp_files[str_detect(temp_files, "POR_") | str_detect(temp_files, "Degradation_")]
temp_files
file.copy(temp_files,
          paste0(zonation_dir, "/data"))

# copy feature files in correct directories
for(zonation_scenario in zonation_scenarios){

  if(str_detect(zonation_scenario, "_10")){
    if(str_detect(zonation_scenario, "_100")){
      file.copy(paste0("_results/Zonation/features_100.txt"),
                paste0(zonation_dir, "/", zonation_scenario, "/features_100.txt"))
      if(str_detect(zonation_scenario, "grW")) 
        file.copy(paste0("_results/Zonation/group_weights_100.txt"),
                  paste0(zonation_dir, "/", zonation_scenario, "/group_weights_100.txt"))
      
    } else {
      file.copy(paste0("_results/Zonation/features_10.txt"),
                paste0(zonation_dir, "/", zonation_scenario, "/features_10.txt"))
      if(str_detect(zonation_scenario, "grW")) 
        file.copy(paste0("_results/Zonation/group_weights_10.txt"),
                  paste0(zonation_dir, "/", zonation_scenario, "/group_weights_10.txt"))
    }
    
  } else {
    file.copy(paste0("_results/Zonation/features.txt"),
              paste0(zonation_dir, "/", zonation_scenario, "/features.txt"))
    if(str_detect(zonation_scenario, "grW")) 
      file.copy(paste0("_results/Zonation/group_weights.txt"),
                paste0(zonation_dir, "/", zonation_scenario, "/group_weights.txt"))
  }
  
  print(list.files(paste0(zonation_dir, "/", zonation_scenario)))
  print(paste0("Successfully copied files to ", zonation_dir, "/", zonation_scenario))
}

## create script to run all analysis
# zonation settings tried
# " -wga " #without group weights and with area analysis
# "-wWga" #with group weights and with area analysis
# "-wWgxXa" #with group weights and cost layer and cost prio and area analysis
# "-wWxXha" #same plus hierarchical mask, and without groups

# example call:
#"C:/Program Files/Zonation5/z5w.exe" --mode=CAZ2 -wWga --area=1 D:/EIE_Macroecology/_students/Romy/SoilBioPrio_Zonation/POR/targeted_groupWeights/settings_targeted.z5 D:/EIE_Macroecology/_students/Romy/SoilBioPrio_Zonation/POR/targeted_groupWeights 

zonation_exe   <- '"C:/Program Files/Zonation5/z5w.exe"'
zonation_options <- "--mode=CAZ2 -SCENARIO --area=1"

# Output file
script_output <- file.path(zonation_dir, "/Zonation_script.txt")

sink(script_output)

for (zonation_scenario in zonation_scenarios) {
  approach <- na.omit(str_extract(zonation_scenario, approaches))[1]

  # construct Zonation call
  temp_options <- "-w"
  if(str_detect(zonation_scenario, "_grW")) temp_options <- paste0(temp_options, "W")
  temp_options <- paste0(temp_options, "g")
  if(str_detect(zonation_scenario, "prevent")) temp_options <- paste0(temp_options, "xX")
  if(str_detect(zonation_scenario, "complement")) temp_options <- paste0(temp_options, "h")
  if(str_detect(zonation_scenario, "prevent")) temp_options <- paste0(temp_options, "h")
  temp_options <- paste0(temp_options, "a")
  
  temp_options <- str_replace(zonation_options, "-SCENARIO", temp_options)
  
  cat(":: Scenario: ", zonation_scenario, "\n")
  cat(
    paste(
      zonation_exe,
      temp_options,
      file.path(zonation_dir, zonation_scenario, paste0("settings_", approach, ".z5")),
      file.path(zonation_dir, zonation_scenario)
    ),
    "\n\n",     # blank line between scenarios
    sep = ""
  )
}

sink()  # very important to close the sink!
cat("Commands written to:", script_output, "\n")




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
