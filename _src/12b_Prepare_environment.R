#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#          Prepare environment              #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(tidyverse)
library(terra)

#- - - - - - - - - - - - - - - - - - - - - -
## Load data ####
# load mask (North of Portugal)
temp_mask <- terra::vect( "T:/01_On_Going_Projects/SoilRecon/Modeling_SDMs/mask.shp")
temp_mask

# load environmental stack
env <- terra::rast("D:/00_datasets/EnvPredictor_1km_all_v2.grd")
names(env)

env <- terra::subset(env, c("MAP_Chelsa.EU", "MAPseas_Chelsa.EU", "MAT_Chelsa.EU", "MATseas_Chelsa.EU",
                            "Aridity", "Snow", "Agriculture", "Dist_Urban", "Forest_Coni", "Forest_Deci",
                            "NDVI", "Pastures", "Pop_Dens", "Shrubland", "Aspect", "Dist_Coast", "Dist_River",
                            "Elev", "Slope", "CEC", "Clay.Silt", "Cu", "Hg", "Moisture","N", "P", "pH", "SOC",
                            "SoilT"))
env

#- - - - - - - - - - - - - - - - - - - - - -
## Mask ####
env_masked <- terra::mask(x = env, mask = temp_mask)
env_masked

env_masked <- terra::crop(x = env_masked, y = temp_mask)

#terra::plot(env_masked[[1]])

terra::writeRaster(env_masked, "D:/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km.tif")

pdf("D:/_students/Romy/Atlas_Portugal/_figures/Predictors_1km_POR.pdf")
terra::plot(env_masked, maxnl = 30)
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Scale predictors ####

env_masked <- terra::rast("D:/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR.tif")

# save scaling parameters for future climate preparation (later)
scale_values <- data.frame("predictor"=names(env_masked))
scale_values$scale_mean <- terra::global(env_masked, fun=mean, na.rm=TRUE)
scale_values$scale_sd <-  terra::global(env_masked, fun=sd, na.rm=TRUE)
scale_values

# save scale parameters
write.csv(scale_values, file="D:/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_scaling.csv",
          row.names = FALSE)

# scale rasterStack
#env_norm <- terra::scale(env_masked, center=TRUE, scale=TRUE)
env_norm <- env_masked

for(temp_name in names(env_masked)){
  temp_raster <- env_masked[[temp_name]]
  
  temp_mean <- scale_values[scale_values$predictor==temp_name, ] %>% pull(scale_mean) %>% as.numeric()
  temp_sd <- scale_values[scale_values$predictor==temp_name, ] %>% pull(scale_sd) %>% as.numeric()
  
  temp_raster <- terra::app(temp_raster, 
                              fun = function(x, na.rm=T){(x - temp_mean) / temp_sd},
                              na.rm=TRUE)
  names(temp_raster) <- temp_name
  env_norm[[temp_name]] <- temp_raster
  
  print(paste0("Stacked file ", names(temp_raster)))
}

env_norm

# save Env_norm
raster::writeRaster(env_norm, file="D:/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif", overwrite=T)

pdf("D:/_students/Romy/Atlas_Portugal/_figures/Predictors_1km_POR_scaled.pdf")
terra::plot(env_norm, maxnl = 30)
dev.off()






