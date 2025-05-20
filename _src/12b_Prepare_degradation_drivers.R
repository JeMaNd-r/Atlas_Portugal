#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Prepare degradation drivers          #
#          author: Romy Zeiss               #
#            date: 2025-05-20               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

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
#... error with latitude...

climvelo <- f_crop(climvelo, env_por[[1]])
terra::plot(erosion)
terra::writeRaster(erosion, "../Soil_degradation_drivers/ESDAC_Erosion/Erosion_POR.tif")


#- - - - - - - - - - - - - - - - - - - - - -
## Erosion ####
#- - - - - - - - - - - - - - - - - - - - - -
erosion <- terra::rast( "../Soil_degradation_drivers/ESDAC_Erosion/2050/RCP26.tif")

erosion <- f_crop(erosion, env_por[[1]])
terra::plot(erosion)

terra::writeRaster(erosion, "../Soil_degradation_drivers/ESDAC_Erosion/Erosion_POR.tif")


#- - - - - - - - - - - - - - - - - - - - - -
## Land-use instability ####
#- - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - -
## Compaction ####
#- - - - - - - - - - - - - - - - - - - - - -
compaction <- terra::vect( "../Soil_degradation_drivers/ESDAC_Compaction/soil_compaction_pack/data/natural_soil_susc_EU27_laea1052.shp")

compaction <- terra::project(compaction, crs(env_por[[1]]))
compaction <- terra::rasterize(compaction, env_por[[1]], field = "Evaluation")

# replace cells that could not be evaluated with NA
compaction[compaction == 9] <- NA

compaction <- terra::crop(x = compaction, y = env_por[[1]])
compaction <- terra::mask(x = compaction, mask = env_por[[1]])

terra::plot(compaction)
terra::writeRaster(compaction, "../Soil_degradation_drivers/ESDAC_Compaction/Compaction_POR.tif")


#- - - - - - - - - - - - - - - - - - - - - -
## Invasion ####
#- - - - - - - - - - - - - - - - - - - - - -


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





