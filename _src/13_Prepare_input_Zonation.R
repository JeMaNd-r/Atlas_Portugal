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

## Define input data directory ####
data_wd <- "D:/EIE_Macroecology/_students/Romy/Atlas_Portugal"

# define names of taxonomic groups
taxaNames <- c("earthworms", "nematodes") 
taxaNames

#- - - - - - - - - - - - - - - - - - - - - -
### Get minimum extent for all taxa ####
#- - - - - - - - - - - - - - - - - - - - - -

uncertain_list <- as.list(taxaNames)

# load uncertainty extent for all taxa
uncertain_list <- lapply(uncertain_list, function(x){
  load(file=paste0(data_wd, "/_results/SDM_Uncertainty_extent_", x, ".RData")) #extent_df
  
  # load uncertainty
  uncertain_tif <- terra::rast(paste0(data_wd, "/_results/", x, "/SDM_Uncertainty_", x, ".tif"))
  
  # save threshold for uncertainty
  uncertain_thresh <- as.numeric(stats::quantile(uncertain_tif$Mean, 0.9, na.rm=TRUE))
  uncertain_tif$extent <- ifelse(values(uncertain_tif$Mean) >=uncertain_thresh, 0, 1)

  uncertain_tif$extent
})

names(uncertain_list) <- taxaNames

#- - - - - - - - - - - - - - - - - - - - - -
### Load species richness ####
#- - - - - - - - - - - - - - - - - - - - - -

richness_list <- as.list(taxaNames)

richness_list <- lapply(names(uncertain_list), function(x){
  load(paste0(data_wd, "/_results/", x, "/SDM_stack_bestPrediction_binary_", x, ".RData")) #species_stack
  species_stack <- terra::rast(species_stack, crs=terra::crs(uncertain_list[[x]]))
  species_stack <- terra::mask(species_stack, mask=uncertain_list[[x]], maskvalues = 0)
  
  species_stack
})
names(richness_list) <- taxaNames

richness_stack <- do.call(c, richness_list)

#- - - - - - - - - - - - - - - - - - - - - -
### Crop all layers to same extent ####
#- - - - - - - - - - - - - - - - - - - - - -

# define minimum extent
common_extent <- terra::ext(richness_stack)

# crop all rasters
richness_list <- lapply(richness_list, function(x) crop(x, common_extent))

# crop to rasters that needs to be validated (European ones)
eu_mask <- terra::rast("D:/EIE_Macroecology/_students/Romy/DATA_SoilBioPrio/_results/_Maps/Zonation/Collembola/Colle_1.tif")
#eu_ext <- terra::ext(eu_mask)

# crop all rasters
richness_list <- lapply(richness_list, function(x) terra::project(x, eu_mask, 
                                                                  method = "med",
                                                                  mask=TRUE))
richness_list

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
                              names(richness_list[[x]][[i]]))
                        ) 
    
    # last layer not species -> rename
    #...todo
    
    if(names(richness_list[[x]][[i]]) == "Richness"){
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

names_list <- do.call(rbind, names_list)
colnames(names_list) <- c("Taxon", "ID", "SpeciesID")
names_list <- as.tibble(names_list)
names_list

write_csv(names_list,  paste0(data_wd, "/_results/Zonation/TaxaNames_legend_validation.csv"))

names_list <- read_csv(paste0(data_wd, "/_results/Zonation/TaxaNames_legend_validation.csv"))
write_delim(names_list %>% 
              mutate("weight" = "1.0", 
                     "filename" = paste0("validation/", ID, ".tif"),
                     "group" = as.numeric(as.factor(Taxon))) %>%
              dplyr::select(weight, filename, group),
            paste0(data_wd, "/_results/Zonation/features_validation.txt"))
# note: manually rename last taxa in each group to Rich OR delete from feature table

#- - - - - - - - - - - - - - - - - - - - - -
## Test if all same extent

rast1 <- terra::rast(paste0(data_wd, "/_results/Zonation/earthworms/earth_3.tif"))
rast2 <- terra::rast(paste0(data_wd, "/_results/Zonation/earthworms/earth_Rich.tif"))
rast3 <- terra::rast(paste0(data_wd, "/_results/Zonation/nematodes/nemat_5.tif"))
rast4 <- terra::rast(paste0(data_wd, "/_results/Zonation/nematodes/nemat_20.tif"))

terra::compareGeom(rast1, rast2)
terra::compareGeom(rast1, rast3)
terra::compareGeom(rast3, rast4)

par(mfrow = c(2,2))
terra::plot(rast1)
terra::plot(rast2)
terra::plot(rast3)
terra::plot(rast4)
