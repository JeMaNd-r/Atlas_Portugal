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

#- - - - - - - - - - - - - - - - -
## Earthworms ####
#- - - - - - - - - - - - - - - - -
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
  mutate(Genus = str_split_fixed(Species, " ", 2)[,1]) %>%
  unique
recon


# check taxonomic levels
recon %>% filter(Abundance!=0) %>% count(Species)
recon %>% filter(Abundance!=0) %>% count(Genus)

# filter relevant columns
recon <- recon %>% dplyr::select(Genus, POINT_X, POINT_Y, Abundance) %>%
  mutate("SpeciesID" = recon$Genus,
         #"Year" = 2021, "Datasource" = "SoilReCon"
         ) %>%
  rename("x" = POINT_X,
         "y" = POINT_Y) %>%
  dplyr::select(x, y, SpeciesID, Abundance)%>%
  mutate(Presence = ifelse(Abundance>0, 1, 0))
recon

# Save
write_csv(recon, "_intermediates/Occurrence_raw_Earthworms.csv")

# - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Earthworms"
data_raw <- read_csv(paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$x, data_raw$y, data_raw$SpeciesID),] #nrow=83,260
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
r <- terra::rast("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")[[1]]
extent_Portugal <- terra::ext(r)
data <- data %>% filter(extent_Portugal[1] <= x &  x <= extent_Portugal[2]) %>% 
  filter(extent_Portugal[3] <= y &  y <= extent_Portugal[4])
data # nrow=78,055

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Portugal", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner 
# flag problems with coordinates
dat_cl <- data.frame(data)
flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "x", lat = "y",
                                              species = "SpeciesID", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
                                              country_ref = rnaturalearth::ne_countries("small"), 
                                              country_refcol = "iso_a3", urban_ref = NULL)
sum(flags$.summary) #those NOT flagged!
# 0 records flagged

# save flagged coordinates
#write_csv(flags %>% filter(!.summary), file=paste0("/_results/FlaggedRecords_", Taxon_name, ".csv"))

# remove flagged records from the clean data (i.e., only keep non-flagged ones)
dat_cl <- dat_cl[flags$.summary, ]

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))
df_cleaning

# - - - - - - - - - - - - - - - - - - -
### Save clean data 
write_csv(dat_cl, file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# save updated number of records during cleaning process
write_csv(df_cleaning, file=paste0("_results/NoRecords_cleaning_", Taxon_name, ".csv"))

#- - - - - - - - - - - - - - - - -
## Nematodes ####
#- - - - - - - - - - - - - - - - -
# load dataset
recon_raw <- read_delim(file="_data/SoilReCon_nematodes.csv")
recon_raw

recon_raw %>% dplyr::select(taxonomy_f, taxonomy_g) %>% unique() #45

recon <- recon_raw %>% 
  mutate_all(as.character) %>%
  pivot_longer(cols=`1`:`462`, names_to = "SampleID", values_to = "Abundance") %>%
  mutate(Abundance = as.double(Abundance),
         SampleID = as.double(SampleID)) %>%
  filter(!is.na(ID)) %>% #remove empty samples
  arrange(SampleID)
recon #17,864

data_xy <- read_csv(file="_data/SoilReCon_Data_4_23_Locations.csv")

recon <- recon %>%
  full_join(data_xy, by="SampleID") %>%
  rename(x=Longitude, y=Latitude) %>%
  mutate(Presence = ifelse(Abundance>0, 1, 0)) %>%
  dplyr::select(-SampleID)
recon

# fix spelling mistake
recon <- recon %>%
  mutate(taxonomy_f = ifelse(str_detect(taxonomy_f,"Tylenchidae"), "Tylenchidae", taxonomy_f))

# check taxonomic levels
recon %>% filter(Abundance!=0) %>% count(taxonomy_g) #19 (-2)
recon %>% filter(Abundance!=0) %>% count(taxonomy_f) %>% print(n=100) #35 (-4)

# filter relevant columns & rows
recon <- recon  %>%
  # filter(!is.na(taxonomy_g)) %>% #NemaGenus
  # mutate(SpeciesID = substr(taxonomy_g, 1, 9)) %>% #NemaGenus
  mutate(SpeciesID = substr(taxonomy_f, 1, 8)) %>% #Nematodes (family level)
  dplyr::select(x,y, SpeciesID, Abundance, Presence)

Taxon_name <- "Nematodes" #NemaGenus or Nematodes

# Save
write_csv(recon, "_intermediates/Occurrence_raw_", Taxon_name,".csv")

# - - - - - - - - - - - - - - - - - - -
data_raw <- read_csv(paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$x, data_raw$y, data_raw$SpeciesID),] #nrow=17,864
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
r <- terra::rast("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")[[1]]
extent_Portugal <- terra::ext(r)
data <- data %>% filter(extent_Portugal[1] <= x &  x <= extent_Portugal[2]) %>% 
  filter(extent_Portugal[3] <= y &  y <= extent_Portugal[4])
data # nrow=17,864

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Portugal", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner 
# flag problems with coordinates
dat_cl <- data.frame(data)
# flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "x", lat = "y",
#                                               species = "SpeciesID", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
#                                               country_ref = rnaturalearth::ne_countries("small"), 
#                                               country_refcol = "iso_a3", urban_ref = NULL)
# sum(flags$.summary) #those NOT flagged!
# # 44 records flagged
# 
# # save flagged coordinates
# #write_csv(flags %>% filter(!.summary), file=paste0("_results/FlaggedRecords_", Taxon_name, ".csv"))
# 
# # remove flagged records from the clean data (i.e., only keep non-flagged ones)
# #dat_cl <- dat_cl[flags$.summary, ]
# # NOTE: will not be removed as they are 44 taxa from the same sample.

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))
df_cleaning

# - - - - - - - - - - - - - - - - - - -
## Save clean data
write_csv(dat_cl, file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# save updated number of records during cleaning process
write_csv(df_cleaning, file=paste0("_results/NoRecords_cleaning_", Taxon_name, ".csv"))


#- - - - - - - - - - - - - - - - -
## 16s data ####
#- - - - - - - - - - - - - - - - -
# load dataset
recon_raw <- read_delim(file="_data/even13500.txt")
recon_raw

# split taxonomy into multiple columns
recon_raw <- recon_raw %>% separate_wider_regex( col = taxonomy,
                                    patterns = c(D = "d_+\\w+",
                                                 "; ",
                                                 P = "p_+\\w+",
                                                 "; ",
                                                 C = "c_+\\w+",
                                                 "; ",
                                                 O = "o_+\\w+",
                                                 "; ",
                                                 F = "f_+\\w+",
                                                 "; ",
                                                 G = "g_+\\w+",
                                                 "; ",
                                                 S = "s_+\\w+"), 
                                    too_few = "align_start")  %>%  
  mutate(D = gsub(x = D, pattern = "d__", replacement = ""),
         P = gsub(x = P, pattern = "p__", replacement = ""),
         C = gsub(x = C, pattern = "c__", replacement = ""),
         O = gsub(x = O, pattern = "o__", replacement = ""),
         F = gsub(x = F, pattern = "f__", replacement = ""),
         G = gsub(x = G, pattern = "g__", replacement = ""),
         S = gsub(x = S, pattern = "s__", replacement = ""))
unique(recon_raw %>% dplyr::select(D:S))

#- - - - - - - - - - - - - - - - -
### Bacteria ####
#- - - - - - - - - - - - - - - - -
# filter per taxon group
Taxon_name <- "Bacteria"

recon <- recon_raw %>% filter(D == Taxon_name)

# check taxonomic level
#View(recon %>% dplyr::select(D:S) %>% unique() %>% arrange(D, P, C, O, F, G, S))

# add ID column for taxa
species_list <- recon %>% 
  dplyr::select(D:S, `#OTU ID`) %>% 
  unique() %>% 
  arrange(D, P, C, O, F, G, S, `#OTU ID`) %>% #removing OTU & S column here and above will give ID per genus
  mutate(taxaID = paste0(substr(Taxon_name, 1, 1), sprintf("%05d", 1:n())))
species_list

write_csv(species_list, paste0("_results/Species_list_", Taxon_name, ".csv"))

recon <- recon %>% 
  full_join(species_list %>% dplyr::select("#OTU ID", taxaID), by = "#OTU ID") %>%
  rename("SpeciesID" = taxaID) %>%
  mutate_all(as.character) %>%
  pivot_longer(cols=`1`:`424`, names_to = "SampleID", values_to = "Abundance") %>%
  mutate(Abundance = as.double(Abundance),
         SampleID = as.double(SampleID)) %>%
  filter(!is.na(SpeciesID)) %>% #remove empty samples
  arrange(SampleID)
recon # 12,246,228

data_xy <- read_csv(file="_data/SoilReCon_Data_4_23_Locations.csv")

recon <- recon %>%
  full_join(data_xy, by="SampleID") %>%
  rename(x=Longitude, y=Latitude) %>%
  mutate(Presence = ifelse(Abundance>0, 1, 0)) %>%
  dplyr::select(x,y,SpeciesID, Abundance, Presence)
recon

# Save
write_csv(recon, paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# - - - - - - - - - - - - - - - - - - -
data_raw <- read_csv(paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$x, data_raw$y, data_raw$SpeciesID),] #nrow=12,214,584
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
r <- terra::rast("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")[[1]]
extent_Portugal <- terra::ext(r)
data <- data %>% filter(extent_Portugal[1] <= x &  x <= extent_Portugal[2]) %>% 
  filter(extent_Portugal[3] <= y &  y <= extent_Portugal[4])
data # nrow=12,214,584

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Portugal", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner
# flag problems with coordinates
dat_cl <- data.frame(data)
flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "x", lat = "y",
                                              species = "SpeciesID", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
                                              country_ref = rnaturalearth::ne_countries("small"), 
                                              country_refcol = "iso_a3", urban_ref = NULL)
sum(flags$.summary) #those NOT flagged!
#  records flagged

# save flagged coordinates
#write_csv(flags %>% filter(!.summary), file=paste0("_results/FlaggedRecords_", Taxon_name, ".csv"))

# remove flagged records from the clean data (i.e., only keep non-flagged ones)
#dat_cl <- dat_cl[flags$.summary, ]
# NOTE: will not be removed as they are 44 taxa from the same sample.

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))
df_cleaning

# - - - - - - - - - - - - - - - - - - -
## Save clean data
write_csv(dat_cl, file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# save updated number of records during cleaning process
write_csv(df_cleaning, file=paste0("_results/NoRecords_cleaning_", Taxon_name, ".csv"))

#- - - - - - - - - - - - - - - - -
## 18s data ####
#- - - - - - - - - - - - - - - - -
# load dataset
recon_raw <- read_delim(file="_data/euk.even5000.txt")
recon_raw

# split taxonomy into multiple columns
recon_raw <- recon_raw %>% separate_wider_regex( col = taxonomy,
                                                 patterns = c(D = "d_+\\w+",
                                                              "; ",
                                                              P = "p_+\\w+",
                                                              "; ",
                                                              C = "c_+\\w+",
                                                              "; ",
                                                              O = "o_+\\w+",
                                                              "; ",
                                                              F = "f_+\\w+",
                                                              "; ",
                                                              G = "g_+\\w+",
                                                              "; ",
                                                              S = "s_+\\w+"), 
                                                 too_few = "align_start")  %>%  
  mutate(D = gsub(x = D, pattern = "d__", replacement = ""),
         P = gsub(x = P, pattern = "p__", replacement = ""),
         C = gsub(x = C, pattern = "c__", replacement = ""),
         O = gsub(x = O, pattern = "o__", replacement = ""),
         F = gsub(x = F, pattern = "f__", replacement = ""),
         G = gsub(x = G, pattern = "g__", replacement = ""),
         S = gsub(x = S, pattern = "s__", replacement = ""))
unique(recon_raw %>% dplyr::select(D:S))

# - - - - - - - - - - - - - - - - - - -
### Fungi ####
# - - - - - - - - - - - - - - - - - - -

# filter per taxon group
Taxon_name <- "Fungi"

recon <- recon_raw %>% filter(str_detect(P, "mycota"))

# check taxonomic level
View(recon %>% dplyr::select(D:S) %>% unique() %>% arrange(D, P, C, O, F, G, S))

# add ID column for taxa
species_list <- recon %>% 
  dplyr::select(D:S, `#OTU ID`) %>% 
  unique() %>% 
  arrange(D, P, C, O, F, G, S, `#OTU ID`) %>% #removing OTU & S column here and above would give ID per genus
  mutate(taxaID = paste0(substr(Taxon_name, 1, 1), sprintf("%05d", 1:n())))
species_list

write_csv(species_list, paste0("_results/Species_list_", Taxon_name, ".csv"))

recon <- recon %>% 
  full_join(species_list %>% dplyr::select("#OTU ID", taxaID), by = "#OTU ID") %>%
  rename("SpeciesID" = taxaID) %>%
  mutate_all(as.character) %>%
  pivot_longer(cols=`1`:`999`, names_to = "SampleID", values_to = "Abundance") %>%
  mutate(Abundance = as.double(Abundance),
         SampleID = as.double(SampleID)) %>%
  filter(!is.na(SpeciesID)) %>% #remove empty samples
  arrange(SampleID)
recon #844,602

data_xy <- read_csv(file="_data/SoilReCon_Data_4_23_Locations.csv")

recon <- recon %>%
  full_join(data_xy, by="SampleID") %>%
  rename(x=Longitude, y=Latitude) %>%
  mutate(Presence = ifelse(Abundance>0, 1, 0)) %>%
  dplyr::select(x,y,SpeciesID, Abundance, Presence)
recon

# Save
write_csv(recon, paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# - - - - - - - - - - - - - - - - - - -
data_raw <- read_csv(paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$x, data_raw$y, data_raw$SpeciesID),] #nrow=17,864
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
r <- terra::rast("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")[[1]]
extent_Portugal <- terra::ext(r)
data <- data %>% filter(extent_Portugal[1] <= x &  x <= extent_Portugal[2]) %>% 
  filter(extent_Portugal[3] <= y &  y <= extent_Portugal[4])
data # nrow=17,864

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Portugal", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner: not needed as same sites as above
dat_cl <- data.frame(data)

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))
df_cleaning

# - - - - - - - - - - - - - - - - - - -
## Save clean data
write_csv(dat_cl, file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# save updated number of records during cleaning process
write_csv(df_cleaning, file=paste0("_results/NoRecords_cleaning_", Taxon_name, ".csv"))


# - - - - - - - - - - - - - - - - - - -
### Other Eukaryotes ####
# - - - - - - - - - - - - - - - - - - -
View(unique(recon_raw %>% dplyr::select(D:S)))

write_csv(recon_raw %>% dplyr::select(D, P) %>% unique() %>% arrange(P),
  file = "_intermediates/Species_list_raw_Eukaryota.csv")

# filter per taxon group
Taxon_name <- "Protists"
recon <- recon_raw %>% 
  filter(str_detect(P, "mycota"))

# check taxonomic level
View(recon %>% dplyr::select(D:S) %>% unique() %>% arrange(D, P, C, O, F, G, S))

# add ID column for taxa
species_list <- recon %>% 
  dplyr::select(D:S, `#OTU ID`) %>% 
  unique() %>% 
  arrange(D, P, C, O, F, G, S, `#OTU ID`) %>% #removing OTU & S column here and above would give ID per genus
  mutate(taxaID = paste0(substr(Taxon_name, 1, 1), sprintf("%05d", 1:n())))
species_list

write_csv(species_list, paste0("_results/Species_list_", Taxon_name, ".csv"))

recon <- recon %>% 
  full_join(species_list %>% dplyr::select("#OTU ID", taxaID), by = "#OTU ID") %>%
  rename("SpeciesID" = taxaID) %>%
  mutate_all(as.character) %>%
  pivot_longer(cols=`1`:`999`, names_to = "SampleID", values_to = "Abundance") %>%
  mutate(Abundance = as.double(Abundance),
         SampleID = as.double(SampleID)) %>%
  filter(!is.na(SpeciesID)) %>% #remove empty samples
  arrange(SampleID)
recon #1,940,178

data_xy <- read_csv(file="_data/SoilReCon_Data_4_23_Locations.csv")

recon <- recon %>%
  full_join(data_xy, by="SampleID") %>%
  rename(x=Longitude, y=Latitude) %>%
  mutate(Presence = ifelse(Abundance>0, 1, 0)) %>%
  dplyr::select(x,y,SpeciesID, Abundance, Presence)
recon

# Save
write_csv(recon, paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# - - - - - - - - - - - - - - - - - - -
data_raw <- read_csv(paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$x, data_raw$y, data_raw$SpeciesID),] #nrow=17,864
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
r <- terra::rast("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")[[1]]
extent_Portugal <- terra::ext(r)
data <- data %>% filter(extent_Portugal[1] <= x &  x <= extent_Portugal[2]) %>% 
  filter(extent_Portugal[3] <= y &  y <= extent_Portugal[4])
data # nrow= 1,935,099

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Portugal", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner: not needed as same sites as above
dat_cl <- data.frame(data)

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))
df_cleaning

# - - - - - - - - - - - - - - - - - - -
## Save clean data
write_csv(dat_cl, file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# save updated number of records during cleaning process
write_csv(df_cleaning, file=paste0("_results/NoRecords_cleaning_", Taxon_name, ".csv"))

# - - - - - - - - - - - - - - - - - - -
### Other Eukaryotes ####
# - - - - - - - - - - - - - - - - - - -
View(unique(recon_raw %>% dplyr::select(D:S)))

recon_raw %>% dplyr::select(D, P) %>% unique() %>% print(n=100)

# filter per taxon group
Taxon_name <- "Eukaroytes"
recon <- recon_raw %>% 
  filter(!str_detect(P, "mycota")) %>% # no fungi
  filter() %>% # no protists
  filter(P != "Nematozoa") # no nematodes (earthworms anyway not present)

# check taxonomic level
View(recon %>% dplyr::select(D:S) %>% unique() %>% arrange(D, P, C, O, F, G, S))

# add ID column for taxa
species_list <- recon %>% 
  dplyr::select(D:S, `#OTU ID`) %>% 
  unique() %>% 
  arrange(D, P, C, O, F, G, S, `#OTU ID`) %>% #removing OTU & S column here and above would give ID per genus
  mutate(taxaID = paste0(substr(Taxon_name, 1, 1), sprintf("%05d", 1:n())))
species_list

write_csv(species_list, paste0("_results/Species_list_", Taxon_name, ".csv"))

recon <- recon %>% 
  full_join(species_list %>% dplyr::select("#OTU ID", taxaID), by = "#OTU ID") %>%
  rename("SpeciesID" = taxaID) %>%
  mutate_all(as.character) %>%
  pivot_longer(cols=`1`:`999`, names_to = "SampleID", values_to = "Abundance") %>%
  mutate(Abundance = as.double(Abundance),
         SampleID = as.double(SampleID)) %>%
  filter(!is.na(SpeciesID)) %>% #remove empty samples
  arrange(SampleID)
recon #1,940,178

data_xy <- read_csv(file="_data/SoilReCon_Data_4_23_Locations.csv")

recon <- recon %>%
  full_join(data_xy, by="SampleID") %>%
  rename(x=Longitude, y=Latitude) %>%
  mutate(Presence = ifelse(Abundance>0, 1, 0)) %>%
  dplyr::select(x,y,SpeciesID, Abundance, Presence)
recon

# Save
write_csv(recon, paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# - - - - - - - - - - - - - - - - - - -
data_raw <- read_csv(paste0("_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$x, data_raw$y, data_raw$SpeciesID),] #nrow=17,864
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
r <- terra::rast("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal/_intermediates/EnvPredictor_1km_POR_normalized.tif")[[1]]
extent_Portugal <- terra::ext(r)
data <- data %>% filter(extent_Portugal[1] <= x &  x <= extent_Portugal[2]) %>% 
  filter(extent_Portugal[3] <= y &  y <= extent_Portugal[4])
data # nrow= 1,935,099

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Portugal", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner: not needed as same sites as above
dat_cl <- data.frame(data)

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))
df_cleaning

# - - - - - - - - - - - - - - - - - - -
## Save clean data
write_csv(dat_cl, file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# save updated number of records during cleaning process
write_csv(df_cleaning, file=paste0("_results/NoRecords_cleaning_", Taxon_name, ".csv"))