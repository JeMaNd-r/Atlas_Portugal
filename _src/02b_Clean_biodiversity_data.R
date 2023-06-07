#- - - - - - - - - - - - - - - - -#
#        Data cleaning            #
#                                 #
#     author: Romy Zeiss          #
#       date: 2023-06-06          #
#- - - - - - - - - - - - - - - - -#

gc()
library(tidyverse)
library(terra) # to get spatial extent of environmental variables
library(CoordinateCleaner) # to check for spatial issues

# load dataset
recon_raw <- read_csv(file="_data/SoilReCon_earthworms_clean.csv")

# remove non-species level and NA species
recon_raw <- recon_raw[recon_raw$Species!="Octolasion sp." & !is.na(recon_raw$Species),]

# get absences (= sites where species are not present)
recon_long <- recon_raw %>%
  dplyr::select(Sample_ID, Species, Abundance) %>%
  pivot_wider(id_cols = Species, names_from = Sample_ID, 
              values_from = Abundance, values_fill = 0)
recon_long

recon <- recon_raw %>% 
  dplyr::select(Sample_ID, POINT_X, POINT_Y) %>%
  full_join(recon_long %>% 
              pivot_longer(cols="1":"462", 
                           names_to = "Sample_ID",
                           values_to = "Abundance") %>%
              mutate(Sample_ID = as.double(Sample_ID)), multiple="all") %>%
  full_join(recon_raw %>% dplyr::select(Sample_ID , LU_type), multiple="all") %>%
  full_join(recon_raw %>% dplyr::select(Species, 'Ecological grouping'), multiple="all") %>%
  arrange(Sample_ID) %>%
  unique
recon

# filter relevant columns
recon <- recon %>% dplyr::select(Species, POINT_X, POINT_Y, Abundance) %>%
  mutate("SpeciesID" = paste0(substr(word(recon$Species), 1, 5), "_", substr(word(recon$Species,2), 1, 4)),
         #"Year" = 2021, "Datasource" = "SoilReCon"
         ) %>%
  mutate("SpeciesID" = ifelse(Species == "Dendrodrilus rubidus tenuis", "Dendr_ruTe",
                              ifelse(Species == "Dendrodrilus rubidus subrubicundus", "Dendr_ruSu", SpeciesID))) %>%
  rename("x" = POINT_X,
         "y" = POINT_Y) %>%
  dplyr::select(x, y, SpeciesID, Abundance)
recon

# Save
write_csv(recon, "_intermediates/Occurrence_raw_earthworms.csv")

# - - - - - - - - - - - - - - - - - - -

data_raw <- read_csv("_intermediates/Occurrence_raw_earthworms.csv")

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$x, data_raw$y, data_raw$SpeciesID),] #nrow=83,260
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
r <- terra::rast("D:/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")[[1]]
extent_Portugal <- terra::ext(r)
data <- data %>% filter(extent_Portugal[1] <= x &  x <= extent_Portugal[2]) %>% 
  filter(extent_Portugal[3] <= y &  y <= extent_Portugal[4])
data # nrow=78,055

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Portugal", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner ####
# flag problems with coordinates
dat_cl <- data.frame(data)
flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "x", lat = "y",
                                              species = "SpeciesID", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
                                              country_ref = rnaturalearth::ne_countries("small"), 
                                              country_refcol = "iso_a3", urban_ref = NULL)
sum(flags$.summary) #those NOT flagged!
# 0 records flagged

# save flagged coordinates
#write_csv(flags %>% filter(!.summary), file="/_results/FlaggedRecords_earthworms.csv")

# remove flagged records from the clean data (i.e., only keep non-flagged ones)
dat_cl <- dat_cl[flags$.summary, ]

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))
df_cleaning

# - - - - - - - - - - - - - - - - - - -
## Save clean data ####
write_csv(dat_cl, file="_intermediates/Occurrences_clean_earthworms.csv")

# save updated number of records during cleaning process
write_csv(df_cleaning, file="_results/NoRecords_cleaning_earthworms.csv")

