#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Point occurrences to grid            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

gc()
library(tidyverse)
library(terra)

#- - - - - - - - - - - - - - - - - - - - -
## load grid
r <- terra::rast("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")
r <- r[[1]]
r

Taxon_name <- "Nematodes"
# Note: don't do this for "Earthworms". Because of their low occurrence, we will use point-data instead of the gridded data.

## Load occurrence data
occ <- read.csv(paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# spatial thinning to 5km² scale
occ_points <- data.frame(x=0, y=0)[0,]
occ_abund <- data.frame(x=0, y=0)[0,]

for(spID in unique(occ$SpeciesID)){
  
  print(paste0(spID, " is now tried to add to raster stack."))
  
  # ignore errors
  try({
    
    # load occurrences for specific species
    temp_occ <- occ[occ$SpeciesID==spID & !is.na(occ$SpeciesID), c("x", "y", "Presence", "Abundance", "SpeciesID")] 
    
    # make occurrences as SpatVector object
    temp_occ <- terra::vect(temp_occ, crs="+proj=longlat +datum=WGS84", geom=c("x", "y"))
    
    # make points to raster
    grid_points <- terra::rasterize(temp_occ, r, "Presence", fun=max)
    names(grid_points) <- spID
    
    # same but with abundance data
    grid_abund <- terra::rasterize(temp_occ, r, "Abundance", fun=mean)
    names(grid_abund) <- spID
    
    # crop to raster extent
    grid_points <- terra::mask(grid_points, r)
    grid_abund <- terra::mask(grid_abund, r)
    
    # add to point data frame
    temp_points <- as.data.frame(grid_points, xy=TRUE, row.names = FALSE)
    temp_abund <- as.data.frame(grid_abund, xy=TRUE, row.names = FALSE)
    
    occ_points <- dplyr::full_join(occ_points, temp_points, by=c("x", "y"))
    occ_abund <- dplyr::full_join(occ_abund, temp_abund, by=c("x", "y"))
    
    print(paste0("It has ", nrow(occ[occ$SpeciesID==spID & occ$Presence==1,] %>% dplyr::select(x, y) %>% unique()), " records."))
    print("######################################################")
    
  }) # end of try loop
}

# save point data frame
write_csv(occ_points, file=paste0("_intermediates/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
write_csv(occ_abund, file=paste0("_intermediates/Abundance_rasterized_1km_", Taxon_name, ".csv"))

# nrow(occ_points) #92?
# ggplot()+
#   geom_tile(data = terra::as.data.frame(r, xy=TRUE), aes(x=x, y=y), fill="white", color="grey")+
#   geom_point(data = occ_points, aes(x=x, y=y), size=4)+
#   geom_point(data = occ, aes(x=x, y=y, color="raw"))


