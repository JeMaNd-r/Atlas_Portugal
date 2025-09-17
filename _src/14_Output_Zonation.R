#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Check output from Zonation           #
#          author: Romy Zeiss               #
#            date: 2025-09-17               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

zonation_dir <- "../SoilBioPrio_Zonation/POR"
list.files(zonation_dir, include.dirs = TRUE)

library(terra)
library(tidyverse)

approaches <- c("targeted", "complementing", "prevention")
species_numbers <- c("", "_10", "_100")
protected_areas <- c("", "_allPA")

# load priority maps
priority_rast <- terra::rast()
for(approach in approaches){
  for(species_number in species_numbers){ 
    for(protected_area in protected_areas) try({
    temp_rast <- terra::rast(paste0(zonation_dir, "/", approach, protected_area, species_number, "/rankmap.tif"))
    names(temp_rast) <- paste0(approach, protected_area, species_number)
    
    priority_rast <- c(priority_rast, temp_rast)
    
    print(paste0("Priority map of ", paste0(approach, protected_area, species_number), " added."))
  }, silent = TRUE)}
}
priority_rast

# mask to extent to get rid of strange orange background
env_por <- terra::rast( "_intermediates/EnvPredictor_1km_POR_normalized.tif")
#priority_rast <- terra::crop(priority_rast, env_por[[1]])
priority_rast <- terra::mask(priority_rast, env_por[[1]])

terra::plot(priority_rast)
terra::plot(terra::diff(priority_rast))

#- - - - - - - - - - - - - - - - - - - - - -#
## Classify top areas ####

priority_rast_c <- priority_rast

# top 30, 10, 5, and 1% of area (if no protected areas)
m <- matrix(c(
  0,    0.7,  0,
  0.7,  0.9,  1,
  0.9,  0.95, 2,
  0.95, 0.99, 3,
  0.99, 1.0,  4
), ncol=3, byrow=TRUE)

# maps without protected areas
for (i in names(priority_rast_c)[names(priority_rast_c) %in% paste0("targeted", species_numbers)]) {
  
  priority_rast_c[[i]] <- classify(priority_rast_c[[i]], m)
  
  temp_name <- names(priority_rast_c[[i]])
  levels(priority_rast_c[[i]]) <- data.frame(
    ID    = 0:4,
    class = c(">30%", "30%","10%","5%","1%")
  )
  names(priority_rast_c[[i]]) <- temp_name
}

# top 30, 10, 5, and 1% of area WITH protected areas
protect_stack <- terra::rast(paste0("_intermediates/Protection_POR_coverage.tif"))

coverage_all <- terra::ncell(protect_stack)
coverage_pa <- unlist(global(protect_stack$PA==1, fun = "sum", na.rm = TRUE))
coverage_iucnpa <- unlist(global(protect_stack$IUCN_PA==1, fun = "sum", na.rm = TRUE))
print(paste0("Number of cells in total: ", coverage_all)); print(paste0("Number of cells under any protection: ", coverage_pa)); print(paste0("Number of cells under IUCN protection: ", coverage_iucnpa))

percentage_pa <- coverage_pa / coverage_all
percentage_iucnpa <- coverage_iucnpa / coverage_all
print(paste0("Percentage of area coverage by any PA: ", round(percentage_pa*100,2), "% and by IUCN PAs: ", round(percentage_iucnpa*100,2), "%."))

m_pa <- matrix(c(
  0, 1-percentage_pa-0.3,  0,
  1-percentage_pa-0.3,  1-percentage_pa-0.1,  1,
  1-percentage_pa-0.1,  1-percentage_pa-0.05, 2,
  1-percentage_pa-0.05, 1-percentage_pa-0.01, 3,
  1-percentage_pa-0.01, 1-percentage_pa+0.00001, 4,
  1-percentage_pa+0.00001, 1.0, 5
), ncol=3, byrow=TRUE)
m_pa; percentage_pa

m_iucnpa <- matrix(c(
  0, 1-percentage_iucnpa-0.3,  0,
  1-percentage_iucnpa-0.3,  1-percentage_iucnpa-0.1,  1,
  1-percentage_iucnpa-0.1,  1-percentage_iucnpa-0.05, 2,
  1-percentage_iucnpa-0.05, 1-percentage_iucnpa-0.01, 3,
  1-percentage_iucnpa-0.01, 1-percentage_iucnpa+0.00001, 4,
  1-percentage_iucnpa+0.00001, 1.0, 5
), ncol=3, byrow=TRUE)
m_iucnpa; percentage_iucnpa

# maps without protected areas
layer_names <- apply(expand.grid("complementing", protected_areas, species_numbers), 1, paste, collapse = "")
for (i in names(priority_rast_c)[names(priority_rast_c) %in% layer_names]) {
  
  if(str_detect(i, "allPA")){
    priority_rast_c[[i]] <- classify(priority_rast_c[[i]], m_pa)
  } else {
    priority_rast_c[[i]] <- classify(priority_rast_c[[i]], m_iucnpa)
  }
  
  temp_name <- names(priority_rast_c[[i]])
  levels(priority_rast_c[[i]]) <- data.frame(
    ID    = 0:5,
    class = c(">30%", "30%","10%","5%","1%", "PA")
  )
  names(priority_rast_c[[i]]) <- temp_name
}

terra::plot(priority_rast_c)
terra::plot(terra::diff(priority_rast>=0.95))

# save priority maps
terra::writeRaster(priority_rast, "_results/_Maps/Zonation_priorities_raw.tif", overwrite=TRUE)
terra::writeRaster(priority_rast_c, "_results/_Maps/Zonation_priorities.tif", overwrite=TRUE)

# plot classified priorities (top 30, 5 etc. %)
pdf("_figures/Zonation_priorities.pdf")
terra::plot(priority_rast_c)
dev.off()

# plot difference between top 5% areas
pdf("_figures/Zonation_priorities_diff_top5.pdf")
terra::plot(terra::diff(priority_rast>=0.95))
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Compare with richness maps ####
features <- read_delim(paste0(zonation_dir, "/features.txt"))
features_10 <- read_delim(paste0(zonation_dir, "/features_10.txt"))
features_100 <- read_delim(paste0(zonation_dir, "/features_100.txt"))

richness <- lapply(unique(features$group), function(x){
  temp_features <- features %>% filter(group == x)
  temp_rast <- terra::rast(paste0(zonation_dir, "/data/", temp_features$filename))
  temp_rast <- terra::app(temp_rast, sum)/1000
  names(temp_rast) <- str_extract(temp_features$filename[1], "(?<=\\.\\./data/)[^/]+")
  temp_rast
})
richness <- do.call(c, richness)

terra::plot(richness)
terra::writeRaster(richness, "_results/_Maps/Richness_allTaxa_features.tif", overwrite=TRUE)

pdf("_figures/Richness_allTaxa_features.pdf")
terra::plot(richness)
dev.off()
