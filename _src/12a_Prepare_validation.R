#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#        Prepare data & analysis            #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(tidyverse)

#- - - - - - - - - - - - - - - - - - - - - -
## Load data ####
# load raw data (if using presence-absence)
# we don't need rarefied abundances for presence-absence only models
data_18s_raw <- read_delim("_data/euk.even5000.txt")
data_18s_raw

data_16s_raw <- read_delim("_data/even13500.txt")
data_16s_raw

# Note: rarefied data is better when using abundance-based models 

# separate taxonomy column
data_18s <- data_18s_raw %>% 
  separate(taxonomy, sep="; ", 
           into=c("Division", "Phylum", "Class", "Order", "Family", "Genus", "Species")) 
data_18s

data_16s <- data_16s_raw %>% 
  separate(taxonomy, sep="; ", 
           into=c("Division", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
data_16s

# remove "d__" from Division column
data_18s <- data_18s %>% 
  mutate_at(vars(c("Division", "Phylum", "Class", "Order", "Family", "Genus", "Species")), 
            ~ str_replace(., "[:alpha:]{1}__", ""))
data_16s <- data_16s %>% 
  mutate_at(vars(c("Division", "Phylum", "Class", "Order", "Family", "Genus", "Species")), 
            ~ str_replace(., "[:alpha:]{1}__", ""))

#### This part can be done in a later stage....
# #- - - - - - - - - - - - - - - - - - - - - -
# ## Split data into taxonomic groups ####
# 
# ## Separate Eukaryota taxa by division-specific taxonomic levels
# # get taxonomic backbone
# data_taxonomy <- data_18s %>% dplyr::select(Division, Phylum, Class, Order) %>% 
#   unique() %>% 
#   arrange(Division, Phylum, Class, Order) %>%
#   full_join(data_16s %>% dplyr::select(Division, Phylum, Class, Order) %>% 
#               unique())
# data_taxonomy
# #write_csv(data_taxonomy, file=paste0(data_wd, "/_intermediates/Validation/VALID_Taxonomy_occurrences.csv"))
# data_taxonomy <- read_csv(file=paste0(data_wd, "/_intermediates/Validation/VALID_Taxonomy_occurrences.csv"))
# 
# # define division-specific levels
# levels_taxa <- tibble("Division" = c("Archaea", "Bacteria", "Eukaryota", "Unassigned"),
#                       "Level" = c("Order", # only 1 OTU for Archaea
#                                   "Division", # resolution large enough and already >30 phylum
#                                   "Phylum",  # but then Phylum- and Class-specific separation
#                                   "Division")) #only 1 OTU with no taxonomic information 
# 
# 
# levels_euka <- tibble("Phylum" = c("Amoebozoa", "Annelida", "Aphelidea", "Apicomplexa",
#                                    "Apusomonadidae", "Arthropoda" ,"Mollusca",
#                                    "Nematozoa"),
#                       "Level" = c("Phylum", # classes not of interest here
#                                   "Phylum", # only Clitellata as class
#                                   "Phylum", # only 1 class & order
#                                   "Phylum", # classes not of interest here
#                                   "Phylum",
#                                   "Class", #note: Ellipura have to be split into Collembola & Protura, and Insecta maybe too to order level
#                                   "Phylum", # only 1 class (Gastropoda)
#                                   "Phylum")) 
# 
# # add to table
# redefine_taxon <- function(data){ 
#   data <- data %>% mutate("Taxon" = Division)
#   
#   for(i in 1:nrow(data)){
#     temp_level <- levels_taxa[levels_taxa$Division==data[i,"Division"] %>% pull(), "Level"] %>% pull()
#     if(!is.na(data[i,temp_level])) data[i,"Taxon"] <- data[i, temp_level]
#   }
#   
#   for(i in 1:nrow(data)){
#     if(data[i,"Division"]=="Eukaryota" &
#        !is.na(data[i,"Phylum"])){
#       data[i,"Taxon"] <- data[i,"Phylum"]
#     }
#     
#     if(data[i, "Taxon"] %in% levels_euka$Phylum){
#       temp_level <- levels_euka[levels_euka$Phylum==data[i,"Taxon"] %>% pull(), "Level"] %>% pull()
#       if(!is.na(data[i,temp_level])) data[i,"Taxon"] <- data[i, temp_level]
#     }
#   }
#   
#   return(data)
# }
# 
# data_18s <- redefine_taxon(data_18s)
# data_16s <- redefine_taxon(data_16s)

#- - - - - - - - - - - - - - - - - - - - - -
## Transform to long format ####

data_18s <- data_18s %>% 
  rename("OTU" = `#OTU ID`) %>%
  pivot_longer(cols=`1`:`999`, names_to = "Sample", values_to = "Abundance")
data_18s

data_16s <- data_16s %>% 
  rename("OTU" = `#OTU ID`) %>%
  pivot_longer(cols=`1`:`424`, names_to = "Sample", values_to = "Abundance") #note: for raw datasets, you need to change ...cols=`001`:...
data_16s

#- - - - - - - - - - - - - - - - - - - - - -
## Merge data ####

data_full <- data_18s %>% mutate("OTU_ID" = paste0("18_", OTU)) %>%
  full_join(data_16s %>% mutate("OTU_ID" = paste0("16_", OTU)))
data_full

#- - - - - - - - - - - - - - - - - - - - - -
## Clean data ####
# check taxa with less than 10 occurrences
data_full %>% filter(Abundance>0) %>% group_by(OTU_ID) %>% count() %>% nrow() # 43,987 OTUs
otus_remove <- data_full %>% filter(Abundance>0) %>% group_by(Division, OTU_ID) %>% count() %>% filter(n<10) 
otus_remove # 19,453 OTUs

otus_remove %>% ungroup() %>% group_by(Division) %>% count()
# Division       n
# Bacteria   12278
# Eukaryota   6787
# Unassigned   388

# remove taxa with less than 10 occurrences
data_full <- data_full %>% filter(!(OTU_ID %in% otus_remove$OTU_ID)) 
#data_full %>% filter(Abundance_raw>0) %>% group_by(OTU_ID) %>% count() %>% nrow() 
#-> nrow of unique OTUs should be 34,613

#- - - - - - - - - - - - - - - - - - - - - -
## Clean data ####

# # separate taxa by Division
# list_val <- data_full %>% split(.$Division)
# 
# # check for NAs
# lapply(list_val, function(x) x %>% summarize(across(everything(), function(x) sum(is.na(x)))))

#- - - - - - - - - - - - - - - - - - - - - -
## Save occurrence data ####
write_csv(data_full, file="_intermediates/Occurrences_clean.csv")

rm(list=ls())

#- - - - - - - - - - - - - - - - - - - - - -
## Add locations ####

data_full <- read_csv(file="_intermediates/Occurrences_clean.csv")
data_xy <- read_csv(file="_data/SoilReCon_Data_4_23_Locations.csv")

occ <- data_full %>% mutate(Sample=as.double(Sample)) %>% 
  full_join(data_xy, by=c("Sample"= "SampleID")) %>%
  rename(x=Longitude, y=Latitude) %>%
  mutate(Presence = ifelse(Abundance>0, 1, 0)) %>%
  dplyr::select(x,y,OTU_ID, Abundance, Presence)
occ

write_csv(occ, file="_intermediates/Occurrences_filtered.csv")

rm(occ, data_full, data_xy)

#- - - - - - - - - - - - - - - - - - - - - -
## Prepare data for SDM ####

# load grid
r <- terra::rast("D:/00_datasets/Grids/grid_1k_0p008.tif")

# load data
occ <- read_csv(file="_intermediates/Occurrences_filtered.csv")

# spatial thinning to 5km² scale
occ_points <- data.frame(x=0, y=0)[0,]
occ_abund <- data.frame(x=0, y=0)[0,]

for(spID in unique(occ$OTU_ID)){
  
  print(paste0(spID, " is now tried to add to raster stack."))
  
  # ignore errors
  try({
    
    # load occurrences for specific species
    temp_occ <- occ[occ$OTU_ID==spID & !is.na(occ$OTU_ID), c("x", "y", "Presence", "Abundance", "OTU_ID")] 
    
    # make occurrences as SpatVector object
    temp_occ <- terra::vect(temp_occ, crs="+proj=longlat +datum=WGS84", geom=c("x", "y"))
    
    # make points to raster
    grid_points <- terra::rasterize(temp_occ, r, "Presence", fun=sum)
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
    
    print(paste0("It has ", nrow(occ[occ$OTU_ID==spID & occ$Presence==1,] %>% dplyr::select(x, y) %>% unique()), " records."))
    print("######################################################")
    
  }) # end of try loop
}

# save point data frame
write_csv(occ_points, file="/_intermediates/Occurrence_rasterized_1km.csv")


# calculate richness per site



# SDMs

