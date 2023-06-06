#- - - - - - - - - - - - - - - - -#
#        Data cleaning            #
#                                 #
#     author: Romy Zeiss          #
#       date: 2023-06-06          #
#- - - - - - - - - - - - - - - - -#

gc()
library(tidyverse)

# load dataset
data_raw <- read_csv(file=paste0(data_wd, "/_intermediates/Occurrence_raw_", Taxon_name, ".csv"))

# - - - - - - - - - - - - - - - - - - -
# create a table to see how many records get removed.
df_cleaning <- tibble::tibble(CleaningStep="merged_RawData", NumberRecords=nrow(data_raw))

# remove data without coordinates or taxa names
data <- data_raw[complete.cases(data_raw$longitude, data_raw$latitude, data_raw$species),] #nrow=83,260
data

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_coordinates", NumberRecords=nrow(data))

# remove records outside of Europe
extent_Europe <- c(-23, 60, 31, 75)
data <- data %>% filter(extent_Europe[1] <= longitude &  longitude <= extent_Europe[2]) %>% 
  filter(extent_Europe[3] <= latitude &  latitude <= extent_Europe[4])
data # nrow=78,055

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_Europe", NumberRecords=nrow(data))

# remove species with only sp in name
data <- data %>% filter(!stringr::str_detect(data$species, "[[:blank:]]sp"))
data # nrow=78,036

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_species", NumberRecords=nrow(data))
df_cleaning

data$OBJECTID <- 1:nrow(data) 

# - - - - - - - - - - - - - - - - - - -
## CoordinateCleaner ####
# flag problems with coordinates
dat_cl <- data.frame(data)
flags <- CoordinateCleaner::clean_coordinates(x = dat_cl, lon = "longitude", lat = "latitude",
                                              species = "species", tests = c("capitals", "centroids", "equal", "gbif", "zeros", "seas"), #normally: test "countries"
                                              country_ref = rnaturalearth::ne_countries("small"), 
                                              country_refcol = "iso_a3", urban_ref = NULL)
sum(flags$.summary) #those NOT flagged!
# Crassiclitellata: 98732 out of 98732
# Collembola: 219367 of 236221 records

# save flagged coordinates
write_csv(flags %>% filter(!.summary), file=paste0(data_wd, "/_results/FlaggedRecords_", Taxon_name, ".csv"))

# remove flagged records from the clean data (i.e., only keep non-flagged ones)
dat_cl <- dat_cl[flags$.summary, ]

df_cleaning <- df_cleaning %>% add_row(CleaningStep="merged_CoordinateCleaner", NumberRecords=nrow(dat_cl))

## Plot flagged records
world.inp <- map_data("world")

pdf(paste0(data_wd, "/_figures/CoordinateCleaner_flagged_records_", Taxon_name, ".pdf"), width=20)
ggplot() + 
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") + 
  xlim(min(data$longitude, na.rm = T), max(data$longitude, na.rm = T)) + 
  ylim(min(data$latitude, na.rm = T), max(data$latitude, na.rm = T)) + 
  geom_point(data = data_raw, aes(x = longitude, y = latitude), colour = "darkred", size = 1, show.legend = T) +
  geom_point(data = dat_cl, aes(x = longitude, y = latitude), colour = "darkgreen", size = 1, show.legend = T) + 
  coord_fixed() + 
  scale_color_manual(name='CoordinateCleaner',
                     values=c('RawRecords = red'='darkred', 'CleanRecords = green'='darkgreen'))+ 
  theme_bw() + theme(axis.title = element_blank())
dev.off()

# - - - - - - - - - - - - - - - - - - -
## Some structuring ####

# create x and y column
dat_cl <- dat_cl %>% mutate("x"=longitude, "y"=latitude)

# remove NA in coordinates (should not change nrow() as we filtered for complete.cases before)
dat_cl <- dat_cl[complete.cases(dat_cl$x),]
dat_cl <- dat_cl[complete.cases(dat_cl$y),]

# - - - - - - - - - - - - - - - - - - -
## Save clean data ####
write_csv(dat_cl, file=paste0(data_wd, "/_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

# load number of records during cleaning process
df_cleaning_gbif <- read_csv(file=paste0(data_wd, "/_intermediates/GBIF_NoRecords_cleaning_", Taxon_name, ".csv"))

df_cleaning <- df_cleaning_gbif %>% full_join(df_cleaning)
df_cleaning

rm(df_cleaning_gbif)

# save updated number of records during cleaning process
write_csv(df_cleaning, file=paste0(data_wd, "/_results/NoRecords_cleaning_", Taxon_name, ".csv"))

