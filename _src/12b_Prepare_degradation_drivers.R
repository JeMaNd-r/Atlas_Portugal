#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare degradation drivers          #
#          author: Romy Zeiss               #
#            date: 2025-05-20               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

## Workflow:
# Load datasets
# Harmonize data
# Define thresholds
# Save binary raster
# fill NA with 0 (as rasters have different extents or only report non NA values)

library(tidyverse)
library(terra)

# load mask (North of Portugal)
env_por <- terra::rast( "_intermediates/EnvPredictor_1km_POR_normalized.tif")
env_por

# function to crop and mask to extend
f_crop <- function(x, input_mask) {
  temp_x <- terra::project(x, crs(input_mask))
  temp_x <- resample(x=temp_x, y=input_mask, method = "near")
  temp_x <- terra::crop(x = temp_x, y = input_mask)
  temp_x <- terra::mask(x = temp_x, mask = input_mask)
  return(temp_x)
}


#- - - - - - - - - - - - - - - - - - - - - -
## ESDAC degradation indicators ####
#- - - - - - - - - - - - - - - - - - - - - -
degr_ind <- terra::rast( "../Soil_degradation_drivers/ESDAC_Degradation/multiband.tif")

degr_compaction <- degr_ind$multiband_8
degr_salinization <- degr_ind$multiband_12
degr_erosion <- degr_ind$multiband_17

degr_sum <- ifel(degr_ind$multiband_20 == 126, NA, degr_ind$multiband_20) #126 = sealing

par(mfrow = c(2,2))
terra::plot(degr_compaction)
terra::plot(degr_salinization)
terra::plot(degr_erosion)
terra::plot(degr_sum) 

degr_compaction <- f_crop(degr_compaction, env_por[[1]])
degr_salinization <- f_crop(degr_salinization, env_por[[1]])
degr_erosion <- f_crop(degr_erosion, env_por[[1]])
degr_sum <- f_crop(degr_sum, env_por[[1]])

par(mfrow = c(2,2))
terra::plot(degr_compaction)
terra::plot(degr_salinization)
terra::plot(degr_erosion)
terra::plot(degr_sum) 

# save
terra::writeRaster(degr_compaction, "../Soil_degradation_drivers/ESDAC_Degradation/Degradation_compaction_POR.tif")
terra::writeRaster(degr_salinization, "../Soil_degradation_drivers/ESDAC_Degradation/Degradation_salinization_POR.tif")
terra::writeRaster(degr_erosion, "../Soil_degradation_drivers/ESDAC_Degradation/Degradation_erosion_POR.tif")
terra::writeRaster(degr_sum, "../Soil_degradation_drivers/ESDAC_Degradation/Degradation_sum_POR.tif")

#- - - - - - - - - - - - - - - - - - - - - -
## Climate velocity ####
#- - - - - - - - - - - - - - - - - - - - - -
climvelo <- terra::rast( "../Soil_degradation_drivers/NCC_RE-RUNv2/data/RevisedTimelines/VEL_RCP262050.tif")

climvelo <- terra::project(climvelo, crs(env_por[[1]]))

climvelo <- f_crop(climvelo, env_por[[1]])
terra::plot(climvelo)
terra::writeRaster(climvelo, "../Soil_degradation_drivers/NCC_RE-RUNv2/ClimVelo_POR.tif")

climvelo <- terra::rast("../Soil_degradation_drivers/NCC_RE-RUNv2/ClimVelo_POR.tif")
summary(climvelo) # Mean :1.116,  Median :1.032
climvelo <- climvelo > 1.032
terra::plot(climvelo)
terra::writeRaster(climvelo, "../Soil_degradation_drivers/NCC_RE-RUNv2/ClimVelo_POR_binary.tif")

#- - - - - - - - - - - - - - - - - - - - - -
## Erosion ####
#- - - - - - - - - - - - - - - - - - - - - -
erosion <- terra::rast( "../Soil_degradation_drivers/ESDAC_Erosion/2050/RCP26.tif")

erosion <- f_crop(erosion, env_por[[1]])
terra::plot(erosion)

terra::writeRaster(erosion, "../Soil_degradation_drivers/ESDAC_Erosion/Erosion_POR.tif")

erosion <- terra::rast("../Soil_degradation_drivers/ESDAC_Erosion/Erosion_POR.tif")
erosion <- erosion > 2 #2 tonnes per ha & year; EUSO dashboard thresholds
terra::plot(erosion)
terra::writeRaster(erosion, "../Soil_degradation_drivers/ESDAC_Erosion/Erosion_POR_binary.tif")


#- - - - - - - - - - - - - - - - - - - - - -
## Land-use instability ####
#- - - - - - - - - - - - - - - - - - - - - -
lu_instability <- terra::rast( "../Soil_degradation_drivers/NCC_RE-RUNv2/data/RevisedTimelines/INS_SSP52050.tif")

lu_instability <- terra::project(lu_instability, crs(env_por[[1]]))

lu_instability <- f_crop(lu_instability, env_por[[1]])
terra::plot(lu_instability)
terra::writeRaster(lu_instability, "../Soil_degradation_drivers/NCC_RE-RUNv2/LUInstability_POR.tif")

lu_instability <- terra::rast("../Soil_degradation_drivers/NCC_RE-RUNv2/LUInstability_POR.tif")
summary(lu_instability) # Mean :0.283, Median :0.120
lu_instability <- lu_instability > 0.12 
terra::plot(lu_instability)
terra::writeRaster(lu_instability, "../Soil_degradation_drivers/NCC_RE-RUNv2/LUInstability_POR_binary.tif")


#- - - - - - - - - - - - - - - - - - - - - -
## Compaction ####
#- - - - - - - - - - - - - - - - - - - - - -
compaction <- terra::vect( "../Soil_degradation_drivers/ESDAC_Compaction/soil_compaction_pack/data/natural_soil_susc_EU27_laea1052.shp")

compaction <- terra::project(compaction, crs(env_por[[1]]))
compaction <- terra::rasterize(compaction, env_por[[1]], field = "Evaluation")

# # replace cells that could not be evaluated with NA
# compaction[compaction == 9] <- NA #9 = urban areas

compaction <- terra::crop(x = compaction, y = env_por[[1]])
compaction <- terra::mask(x = compaction, mask = env_por[[1]])

terra::plot(compaction)
terra::writeRaster(compaction, "../Soil_degradation_drivers/ESDAC_Compaction/Compaction_POR.tif")

compaction <- terra::rast("../Soil_degradation_drivers/ESDAC_Compaction/Compaction_POR.tif")
compaction <- compaction > 1.75 #EUSO threshold
terra::plot(compaction)
terra::writeRaster(compaction, "../Soil_degradation_drivers/ESDAC_Compaction/Compaction_POR_binary.tif")


#- - - - - - - - - - - - - - - - - - - - - -
## Invasion ####
#- - - - - - - - - - - - - - - - - - - - - -
 
## GeoJSON from EASIN portal
# invasion <- sf::read_sf("../Soil_degradation_drivers/EASIN_Invasion/easin-layer-export-2025-05-21.json")
# invasion <- terra::vect(invasion)
# 
# invasion <- terra::project(invasion, crs(env_por[[1]]))
# invasion <- terra::rasterize(invasion, env_por[[1]])
# 
# invasion <- terra::crop(x = invasion, y = env_por[[1]])
# invasion <- terra::mask(x = invasion, mask = env_por[[1]])
# 
# terra::plot(invasion)
# terra::writeRaster(invasion, "../Soil_degradation_drivers/EASIN_Invasion/Invasion_POR.tif")

## Data packages: SHP
# files <- list.files("../Soil_degradation_drivers/EASIN_Invasion/Art.24.EASIN.Reporting.2018.PT/27.05.2019/EasinSpeciesDistribution/SpeciesDistribution/SHP",
#                     full.names = TRUE)
# i = 1
# unzip(files[i], exdir = substr(files[i], 1, nchar(files[i])-4))
# invasion <- terra::vect(list.files(substr(files[i], 1, nchar(files[i])-4), pattern = "shp", full.names = TRUE))
# 
# invasion <- terra::project(invasion, crs(env_por[[1]]))
# invasion <- terra::rasterize(invasion, env_por[[1]], field = "FID")
# 
# # replace cells that could not be evaluated with NA
# invasion[invasion == 9] <- NA
# 
# invasion <- terra::crop(x = invasion, y = env_por[[1]])
# invasion <- terra::mask(x = invasion, mask = env_por[[1]])
# 
# terra::plot(invasion)
# terra::writeRaster(invasion, "../Soil_degradation_drivers/EASIN_Invasion/Invasion_POR.tif")

## Data packages: GeoJSON

files <- list.files("../Soil_degradation_drivers/EASIN_Invasion/Art.24.EASIN.Reporting.2018.PT/27.05.2019/EasinSpeciesDistribution/SpeciesDistribution/GeoJSON",
                    full.names = TRUE,
                    pattern = "EASIN")

invasion <- lapply(files, function(x){
  temp_raster <- terra::vect(x)
  temp_raster$occ <- 1
  temp_raster <- terra::rasterize(temp_raster, env_por[[1]], field = "occ")
  
  temp_raster <- terra::crop(x = temp_raster, y = env_por[[1]])
  temp_raster <- terra::mask(x = temp_raster, mask = env_por[[1]])
})

invasion <- do.call(c, invasion)
invasion <- sum(invasion, na.rm=TRUE)

terra::plot(invasion)
terra::writeRaster(invasion, "../Soil_degradation_drivers/EASIN_Invasion/Invasion_POR.tif")

invasion <- terra::rast("../Soil_degradation_drivers/EASIN_Invasion/Invasion_POR.tif")
invasion <- invasion > 0 #EUSO threshold
terra::plot(invasion)
terra::writeRaster(invasion, "../Soil_degradation_drivers/EASIN_Invasion/Invasion_POR_binary.tif")


#- - - - - - - - - - - - - - - - - - - - - -
## Salinization ####
#- - - - - - - - - - - - - - - - - - - - - -
salinization <- terra::vect( "../Soil_degradation_drivers/ESDAC_Salinization/data/sal_alk_eu27_laea1052.shp")

salinization <- terra::project(salinization, crs(env_por[[1]]))
salinization <- terra::rasterize(salinization, env_por[[1]], field = "Sal_AlkCL")

salinization <- terra::crop(x = salinization, y = env_por[[1]])
salinization <- terra::mask(x = salinization, mask = env_por[[1]])

terra::plot(salinization)
terra::writeRaster(salinization, "../Soil_degradation_drivers/ESDAC_Salinization/Salinization_POR.tif")

salinization <- terra::rast("../Soil_degradation_drivers/ESDAC_Salinization/Salinization_POR.tif")
salinization <- salinization > 0 #EUSO threshold
terra::plot(salinization)
terra::writeRaster(salinization, "../Soil_degradation_drivers/ESDAC_Salinization/Salinization_POR_binary.tif")



#- - - - - - - - - - - - - - - - - - - - - -
## Align all layers ####
#- - - - - - - - - - - - - - - - - - - - - -

degr_drivers <- terra::rast(ext(env_por[[1]]), resolution = res(env_por[[1]]))
crs(degr_drivers) <- crs(env_por[[1]])

degr_drivers$erosion <- terra::rast("../Soil_degradation_drivers/ESDAC_Erosion/Erosion_POR_binary.tif")
degr_drivers$climvelo <- terra::rast("../Soil_degradation_drivers/NCC_RE-RUNv2/ClimVelo_POR_binary.tif")
degr_drivers$compaction <- terra::rast("../Soil_degradation_drivers/ESDAC_Compaction/Compaction_POR_binary.tif")
degr_drivers$lu_instability <- terra::rast("../Soil_degradation_drivers/NCC_RE-RUNv2/LUInstability_POR_binary.tif")
degr_drivers$invasion <- terra::rast("../Soil_degradation_drivers/EASIN_Invasion/Invasion_POR_binary.tif")
degr_drivers$salinization <- terra::rast("../Soil_degradation_drivers/ESDAC_Salinization/Salinization_POR_binary.tif")

terra::plot(degr_drivers)

# fill missing areas
degr_drivers[is.na(degr_drivers)] <- 0

degr_drivers <- terra::mask(degr_drivers, env_por[[1]])
terra::writeRaster(degr_drivers, "../Soil_degradation_drivers/Degradation_POR_binary.tif")

pdf("../Soil_degradation_drivers/Degradation_POR_binary.pdf")
terra::plot(degr_drivers)
dev.off()
