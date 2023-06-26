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

#- - - - - - - - - - - - - - - - - - - - - -
## Calculate variable inflation factor (VIF) ####
# VIF is the extent of correlation between one predictor and all others.
# The lower VIF, the better we can tell what predictor contributed (most) to the model

# load env. stack
env_norm <- raster::stack("D:/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")

## VIF basen on raw data (explanatory raster stack)
env_vif <- usdm::vif(env_norm)

# which predictors should be excluded?
vif_cor <- usdm::vifcor(env_norm, th=0.8)  #th = threshold correlation for exclusion
# how: first find a pair of variables which has the maximum linear correlation 
# (greater than th), and exclude one of them which has greater VIF. The 
# procedure is repeated until no variable with a high correlation coefficient 
# (grater than threshold) with other variables remains.

vif_step <- usdm::vifstep(env_norm, th=10) #VIF >10

# merge both data.frames
env_vif <- env_vif %>% rename("VIF_raw" = VIF) %>% full_join(vif_step@results) %>%
  full_join(as.data.frame(vif_cor@corMatrix) %>% mutate("Variables"=rownames(vif_cor@corMatrix)))

env_vif

write_csv(env_vif, file="_results/VIF_predictors_1km_POR.csv")

rm(env_vif, vif_cor)

#- - - - - - - - - - - - - - - - - - - - - -
## Perform PCA ####

# library(devtools)
# install_github("bleutner/RStoolbox")
library(RStoolbox) # for PCA of rasters
library(vegan) #for variance partitioning
library(raster)

env_vif <- read_csv("_results/VIF_predictors_1km_POR.csv")

# load env. stack
env_norm <- raster::stack("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")

env_exclude_vif <- env_vif %>% filter(is.na(VIF)) %>% pull(Variables) 
env_norm <- raster::subset(env_norm, names(env_norm)[!(names(env_norm) %in% env_exclude_vif)])
names(env_norm)

# run PCA
Env_pca <- RStoolbox::rasterPCA(img = env_norm, 
                                nSample = NULL, #use all pixels, not subset
                                #nComp = nlayers(env_norm), #= default
                                spca = FALSE, # don't standardize as variables are already standardized
                                maskCheck = TRUE)
save(Env_pca, file="_intermediates/EnvPredictor_PCA_1km_POR.RData")

#Env_pca
str(print(summary(Env_pca$model)))

# get cumulative importance
vars <- Env_pca$model$sdev^2
vars <- vars/sum(vars)
cumsum(vars)

sink("_intermediates/PCA_EnvPredictor_1km_POR.txt")
print(summary(Env_pca$model))
print(loadings(Env_pca$model))
sink()

#load(paste0(data_wd, "/_intermediates/EnvPredictor_PCA.RData")) #Env_pca
Env_pca

Env_tif <- raster::stack(Env_pca$map)
Env_tif

# save PC axes
Env_tif <- terra::rast(Env_tif)
terra::writeRaster(Env_tif, "_intermediates/EnvPredictor_PCA_1km_POR.tif") #, overwrite=TRUE

# Env_df <- terra::as.data.frame(Env_clip) 
# ggplot()+
#   geom_point(data=Env_df %>% sample_n(100),
#              aes(x=PCA_1, y=PCA_2))

## Plot PCA ####
# # Extract the necessary components
# loadings <- Env_pca$model$loadings
# sdev <- Env_pca$model$sdev
# 
# # Get the variable names
# variable_names <- names(Env_clip)
# 
# # Set the dimension names of the loadings matrix
# dimnames(loadings) <- list(variable_names, colnames(loadings))
# 
# # Flatten the raster data into a matrix
# raster_data_matrix <- as.matrix(Env_clip)
# 
# # Calculate the scores by matrix multiplication
# scores <- raster_data_matrix %*% t(loadings)
# 
# # Plot the biplot
# biplot(loadings, scores, xlim = c(-max(sdev), max(sdev)), ylim = c(-max(sdev), max(sdev)), main = "Biplot")




