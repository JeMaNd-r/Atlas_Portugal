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

## Define input data directory ####
data_wd <- "D:/EIE_Macroecology/_students/Romy/Atlas_Portugal"

# define names of taxonomic groups
taxaNames <- c("Nematodes") 
taxaNames

#- - - - - - - - - - - - - - - - - - - - - -
### Get minimum extent for all taxa ####
#- - - - - - - - - - - - - - - - - - - - - -

uncertain_list <- as.list(taxaNames)

# load uncertainty extent for all taxa
uncertain_list <- lapply(uncertain_list, function(x){
  #load(file=paste0(data_wd, "/_results/SDM_Uncertainty_extent_", x, "_10.RData")) #extent_df
  
  # load uncertainty
  uncertain_tif <- terra::rast(paste0(data_wd, "/_results/_Maps/SDM_Uncertainty_", x, ".tif"))
  
  # load species names
  species10 <- read_csv(file=paste0(data_wd, "/_intermediates/ESM_", x, ".csv"))$species
  species100 <- read_csv(file=paste0(data_wd, "/_intermediates/SDM_", x, ".csv"))$species
  
  # save threshold for uncertainty
  uncertain_thresh_10 <- as.numeric(stats::quantile(uncertain_tif$Mean_10, 0.9, na.rm=TRUE))
  uncertain_thresh_100 <- as.numeric(stats::quantile(uncertain_tif$Mean_100, 0.9, na.rm=TRUE))
  
  # transform to binary based on uncertainty
  uncertain_10$extent <- ifelse(values(uncertain_tif$Mean_10) >=uncertain_thresh_10, 0, 1)
  uncertain_100$extent <- ifelse(values(uncertain_tif$Mean_100) >=uncertain_thresh_100, 0, 1)
  
  uncertain_tif$extent <- ifelse(values(uncertain_10$extent) != 0 & values(uncertain_100$extent) != 0, 1, 0)
  
  uncertain_tif
})

names(uncertain_list) <- taxaNames

#- - - - - - - - - - - - - - - - - - - - - -
### Load species richness & crop to same extent ####
#- - - - - - - - - - - - - - - - - - - - - -

richness_list <- as.list(taxaNames)

richness_list <- lapply(names(uncertain_list), function(x){
  species_stack <- terra::rast(paste0(data_wd, "/_results/_Maps/SDM_stack_binary_", x, ".tif"))
  species_stack <- terra::mask(species_stack, mask=uncertain_list[[x]]$extent, maskvalues = 0)
  
  species_stack
})
names(richness_list) <- taxaNames

richness_stacks <- do.call(c, richness_list)
 
# #- - - - - - - - - - - - - - - - - - - - - -
# ### Crop all layers to same extent ####
# #- - - - - - - - - - - - - - - - - - - - - -
# 
# # define minimum extent
# common_extent <- terra::ext(do.call(c, richness_stacks)[[1]])
# 
# # crop all rasters
# richness_list <- lapply(richness_list, function(x) crop(x, common_extent))
# 
# # crop to rasters that needs to be validated (European ones)
# eu_mask <- terra::rast("D:/EIE_Macroecology/_students/Romy/DATA_SoilBioPrio/_results/_Maps/Zonation/Collembola/Colle_1.tif")
# #eu_ext <- terra::ext(eu_mask)
# 
# # crop all rasters
# richness_list <- lapply(richness_list, function(x) terra::project(x, eu_mask, 
#                                                                   method = "med",
#                                                                   mask=TRUE))
# richness_list

#- - - - - - - - - - - - - - - - - - - - - -
### Save each layer as individual tif ####
#- - - - - - - - - - - - - - - - - - - - - -

dir.create(paste0(data_wd, "/_results/Zonation"))
#dir.create(paste0(data_wd, "/_results/Zonation/NoData"))
for(i in taxaNames) dir.create(paste0(data_wd, "/_results/Zonation/", i))

names_list <- lapply(names(richness_list), function(x){
  names_list <- c()
  
  for(i in 1:nlyr(richness_list[[x]])){ try({
    taxa_name <- substr(x, 1, 5)
    names_list <- rbind(names_list, 
                        cbind(taxa_name,
                              paste0(taxa_name, "_", i),
                              names(richness_list[[x]][[i]]),
                              sum(values(richness_list[[x]][[i]]), na.rm=TRUE))
                        ) 
    
    # save but last layer (Richness) not species -> rename in names_list
    if(names(richness_list[[x]][[i]]) == "Richness"){
      names_list[i,] <- c(taxa_name,
                          paste0(taxa_name, "_Rich"),
                          "Richness",
                          sum(values(richness_list[[x]][[i]])>0 & !is.na(values(richness_list[[x]][[i]])), na.rm=TRUE))
      
      # save
      writeRaster(richness_list[[x]][[i]], 
                  paste0(data_wd, "/_results/Zonation/", x, "/", taxa_name, "_Rich.tif"),
                  names = paste0(taxa_name, "_Rich")
                  )
    }else{
      writeRaster(richness_list[[x]][[i]], 
                  paste0(data_wd, "/_results/Zonation/", x, "/", taxa_name, "_", i, ".tif"),
                  names = paste0(taxa_name, "_", i)
                  )
    }
  })}
  names_list
})

# save name coding
names_list <- do.call(rbind, names_list)
colnames(names_list) <- c("Taxon", "ID", "SpeciesID", "No_cells_presence")
names_list <- as.tibble(names_list)
names_list

write_csv(names_list,  paste0(data_wd, "/_results/Zonation/TaxaNames_legend.csv"))

# define weigth for each layer
names_list <- read_csv(paste0(data_wd, "/_results/Zonation/TaxaNames_legend.csv"))
write_delim(names_list %>%
              filter(SpeciesID != "Richness") %>%
              mutate("weight" = "1.0", 
                     "filename" = paste0("data/", ID, ".tif"),
                     "group" = as.numeric(as.factor(Taxon))) %>%
              dplyr::select(weight, filename, group),
            paste0(data_wd, "/_results/Zonation/features.txt"))


#- - - - - - - - - - - - - - - - - - - - - -
## Prepare threat & protection layers ####

# load one cropped species layer
r <- terra::rast(paste0(data_wd, "/_results/Zonation/", taxaNames[1], "/", substr(taxaNames[1], 1, 5), "_", 1, ".tif"))
r$all <- 1
r <- terra::subset(r, "all")

## load all files provided by protectedplanet.net 
shp0 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_0/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp1 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_1/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")
shp2 <- terra::vect("D:/EIE_Macroecology/_students/Romy/_Shapefiles/WDPA_WDOECM_Nov2022_Public_all_shp/WDPA_WDOECM_Nov2022_Public_all_shp_2/WDPA_WDOECM_Nov2022_Public_all_shp-polygons.shp")

## merge all files
shp_all <- rbind(shp0, shp1)
rm(shp0, shp1)

shp_all <- rbind(shp_all, shp2)
rm(shp2)

shp_all

## crop to Portuguese extent ####
shp_por <- terra::crop(shp_all, r)
rm(shp_all)
gc()

# save
terra::writeVector(shp_por, filename = paste0(data_wd, "/_intermediates/Protection_POR.shp"))

# load protection for portugal
shp_por <- sf::st_read(paste0(data_wd, "/_intermediates/Protection_POR.shp"))

# create empty stack to store PA types as layers
protect_stack <- terra::rast()

# split multipolygon by type of protection area
# type of protected area: protect_raw$DESIG_ENG
# IUCN category: IUCN_CAT
for(i in unique(shp_por$IUCN_CAT)){
  temp_shp <- shp_por[shp_por$IUCN_CAT == i,]
  
  # calculate percent of grid cells covered by PA type i
  temp_cover <- exactextractr::coverage_fraction(r, temp_shp)[[1]]
  names(temp_cover) <- i
  
  # replace values higher than 1 by 1
  temp_cover[temp_cover > 1] <- 1
  temp_cover[temp_cover < 0] <- 0
  
  protect_stack <- c(protect_stack, temp_cover)
  
  print(paste0(i, " ready."))
}

protect_stack

# add areas with no IUCN category reported
protect_stack$nonIUCN <- exactextractr::coverage_fraction(r, shp_por)[[1]]
protect_stack$nonIUCN[protect_stack$nonIUCN > 1] <- 1
protect_stack$nonIUCN[protect_stack$nonIUCN < 0] <- 0

# calculate sum for IUCN categories Ia-VI
protect_stack$PA_coverage <- sum(protect_stack)
protect_stack$IUCN_coverage <- sum(subset(protect_stack, names(protect_stack) %in% c("Ia", "Ib", "II", "III", "IV", "V", "VI")))
protect_stack$PA <- ifelse(values(protect_stack$PA_coverage) > 0.5, 1, 0)
protect_stack$IUCN_PA <- ifelse(values(protect_stack$IUCN_coverage) > 0.5, 1, 0)

# save
terra::writeRaster(protect_stack, filename = paste0(data_wd, "/_intermediates/Protection_POR_coverage.tif"))

# save layers for Zonation
for(ly in names(protect_stack)[names(protect_stack) %in% c("Ia", "Ib", "II", "III", "IV", "V", "VI")]){
  terra::writeRaster(protect_stack[[ly]], 
                     filename = paste0(data_wd, "/_intermediates/POR_", 
                                       ly, 
                                       ".tif"))
  print(paste0(ly, " ready."))
}

features_pa <- data.frame("ID" = list.files(paste0(data_wd, "/_intermediates/"), pattern = "^POR_"))
features_pa <- features_pa %>%
  mutate("weight" = "1.0", 
         "group" = "Protection") %>%
  dplyr::select(weight, ID, group)
         
write_delim(features_pa,
            paste0(data_wd, "/_results/Zonation/features_pa.txt"))

#- - - - - - - - - - - - - - - - - - - - - -
## Test if all same extent ####

rast1 <- terra::rast(paste0(data_wd, "/_results/Zonation/Earthworms/Earth_3.tif"))
rast2 <- terra::rast(paste0(data_wd, "/_results/Zonation/Earthworms/Earth_Rich.tif"))
rast3 <- terra::rast(paste0(data_wd, "/_results/Zonation/Nematodes/Nemat_5.tif"))
rast4 <- terra::rast(paste0(data_wd, "/_results/Zonation/Nematodes/Nemat_20.tif"))

terra::compareGeom(rast1, rast2)
terra::compareGeom(rast1, rast3)
terra::compareGeom(rast3, rast4)

par(mfrow = c(2,2))
terra::plot(rast1)
terra::plot(rast2)
terra::plot(rast3)
terra::plot(rast4)
