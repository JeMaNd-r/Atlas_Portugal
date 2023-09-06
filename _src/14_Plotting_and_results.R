#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Plotting and other results           #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/_students/Romy/SoilBiodiversity")

gc()
library(tidyverse)
library(here)

library(raster)

library(biomod2) # also to create pseudo-absences

library(ggplot2) # for plotting the curves
library(ggpubr)

library(sf)
library(ggsn) #for scale and north arrow

library(parallel)
library(doParallel)

# plotting
library(gridExtra)
library(colorRamps) #color gradient in Supplementary figure

#write("TMPDIR = 'D:/00_datasets/Trash'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

# change temporary directory for files
#raster::rasterOptions(tmpdir = "D:/00_datasets/Trash")

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "earthworms"
speciesSub <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv")) %>% pull(SpeciesID)
speciesSub

# geographic extent of Europe/Portugal
#extent_Europe <- c(-23, 60, 31, 75)
extent_Portugal <- c()

# load background map
world.inp <- map_data("world")

# function to extract legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 

#- - - - - - - - - - - - - - - - - - - - -
## Occurrences ####
#- - - - - - - - - - - - - - - - - - - - -

# # load raw data
# occ_raw <- read.csv(file=paste0(data_wd, "/_intermediates/Earthworm_occurrence_GBIF-sWorm-Edapho-SoilReCon-JM.csv"))
# occ_raw <- occ_raw %>% rename("species"=?..species)
# 
# # load summary number of records during processing
# occ_process <- read.csv(file=paste0(data_wd, "/_results/NoRecords_summary_Crassiclitellata.csv"))
# occ_process <- occ_process %>%
#   filter(!str_detect(Subset,"species")) %>%
#   filter(!str_detect(Subset,"[_]")) %>%
#   filter(!str_detect(Subset,"merged")) %>%
#   filter(!str_detect(Subset,"total")) 
# 
# ## barplot number of occurrences per datasource
# pdf(paste0(data_wd, "/_figures/OccurrenceRaw_perDatasource_barplot.pdf"), width=10)
# ggplot(data=occ_process, 
#        aes(x=reorder(Subset, NumberRecords), y=NumberRecords, fill=ProcessingStep))+
#   geom_bar(stat="identity", position="dodge")+
#   geom_text(label=occ_process$NumberRecords, hjust=-0.5, position=position_dodge(width=1))+
#   xlab("")+ ylab("Number of occurrence records")+
#   scale_y_continuous(expand=c(0,0), limits=c(0,65000))+
#   coord_flip()+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position =c(0.1,0.8),legend.direction = "horizontal",
#         axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), strip.text = element_text(size=5),
#         panel.border = element_blank(),
#         panel.background = element_rect(fill="grey95"))
# dev.off()
# 
# # load cleaned occurrence records 
# occ_clean <- read.csv(file=paste0(data_wd, "/_intermediates/Occurrences_", Taxon_name, ".csv"))
# 
# # load matrix containing information on number of occurrence records in grid
# occ_points <- read.csv(file=paste0(data_wd, "/_results/Occurrence_rasterized_2km_", Taxon_name, ".csv"))
# occ_points <- occ_points %>% rename("x"=?..x)
# 
# ## plot raw occurrences colored by year
# png(paste0(data_wd, "/_figures/OccurrencesGridded_", Taxon_name, "_perYear.png"), height=4000, width=4000, res=400)
# ggplot()+
#   geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white")+
#   
#   geom_point(data=occ_points %>% arrange(year), 
#              aes(x=x, y=y, col=year),cex=0.3)+ theme_bw()+
#   xlim(-10, 30) +
#   ylim(35, 70) +
#   scale_color_steps2(breaks=c(1970, 1980, 1990, 2000, 2010), midpoint=1995, 
#                      high="#10a53dFF", mid="#ffcf20FF", low="#541352FF")+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position =c(0.15,0.85),legend.direction = "horizontal",
#         axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), strip.text = element_text(size=5),
#         panel.border = element_blank(),
#         panel.background = element_rect(fill="grey95"))
# dev.off()
# 
# # calculate number of records per datasource
# full_join(occ_raw %>% group_by(datasource) %>% count(name="raw"), 
#           occ_clean %>% group_by(datasource) %>% count(name="clean"))
# 
# # load number removed records during cleaning
# read.csv(file=paste0(data_wd, "/_results/NoRecords_cleaning_", Taxon_name, ".csv"))
# 
# # count occurrences per species & data source
# count_data <- occ_raw %>% 
#   full_join(speciesNames[,c("Species", "SpeciesID")], by=c("species"="Species")) %>%
#   group_by(datasource, SpeciesID) %>% count() %>%
#   pivot_wider(names_from = datasource, values_from = n) %>%
#   full_join(occ_clean %>% group_by(datasource, SpeciesID) %>% count() %>%
#               pivot_wider(names_from = datasource, values_from = n), suffix = c("_raw", "_clean"), by="SpeciesID")
# count_data$RawOcc <- rowSums(count_data[,2:7], na.rm=T)
# count_data$CleanOcc <- rowSums(count_data[,8:13], na.rm=T)
# 
# # add number of records after
# count_data <- count_data %>% 
#   full_join(speciesNames[,c("SpeciesID", #"Acc_name", "Species_final", "Species", #"NumCells_1km", 
#                              "NumCells_2km")], by=c("SpeciesID")) %>%
#   filter(!is.na(SpeciesID))
# 
# # replace NA with 0
# #count_data[is.na(count_data)] <- 0
# 
# # add info if species was analysed or not
# count_data$Included <- FALSE
# count_data[count_data$NumCells_2km >=10, "Included"] <- TRUE 
# 
# # sort by included or not, and have a look
# count_data <- count_data %>% arrange(desc(Included), SpeciesID) %>%
#   dplyr::select(SpeciesID, RawOcc, CleanOcc, NumCells_2km, everything()) %>%
#   unique()
# count_data
# 
# # save
# write.csv(count_data, file=paste0(data_wd, "/_results/NoRecords_perSpecies_full_", Taxon_name, ".csv"), row.names = F)
# 
# count_data <- read.csv(file=paste0(data_wd, "/_results/NoRecords_perSpecies_full_", Taxon_name, ".csv"))
# 
# # look at speciesID with most records
# count_data %>% arrange(desc(NumCells_2km), SpeciesID)
# # Eisen_tetr, Aporr_cali, Aporr_rose, Lumbr_terr, Lumbr_rube...
# 
# count_data2 <- count_data %>% full_join(speciesNames %>%
#                                           dplyr::select(-NumCells_2km, -Group_name), by=c("SpeciesID"))
# count_data2
# 
# ## plot count occurrences
# ggplot(count_data2 %>% filter(!is.na(SpeciesID), Included=TRUE), aes(x=RawOcc-CleanOcc, y=SpeciesID, fill=Ecogroup))+
#   geom_bar(stat = "identity")
# 
# ## plot total species' occurrences per data source
# plotOccClean <- ggplot()+
#   geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white")+
#   xlim(-10, 30) +
#   ylim(35, 70) +
#   
#   geom_point(data=occ_clean %>% mutate(datasource=ifelse(datasource=="jean", "jerome", datasource)), 
#              aes(x=longitude, y=latitude, color=datasource), 
#              cex=0.5, shape=19)+
#   scale_color_manual(values=c("#824351", "#C18746", "#FFCF20", "#BAC22F", "#65B039"))+ 
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position =c(0.28,0.95), legend.direction = "horizontal",
#         axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         legend.text = element_text(size=10), #legend.key.size = unit(3, 'cm'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(), strip.text = element_text(size=5),
#         panel.border = element_blank(),
#         panel.background = element_rect(fill="grey95"))+
#   guides(color = guide_legend(override.aes = list(size = 5))) #makes legend icons bigger
# 
# plotOccClean
# 
# #png(paste0(data_wd, "/_figures/OccurrencesClean_", Taxon_name, "_perDatasource.png"),height=4000, width=4000, res=400); plotOccClean; dev.off()
# rm(plotOccClean)
# 
# occ_clean %>% group_by(datasource) %>% count() 
# # edapho 13054, gbif 54883, jean 24860, jerome 732, soilrecon 175, sworm 5028
# 
# ## plot in Germany
# german.inp <- map_data("world", "Germany")
# 
# # plot total species' occurrences in Germany
# plotOccRawGER <- ggplot()+ #, alpha=`Number of Records`
#   #geom_polygon(data=bg.map)+
#   geom_map(data=world.inp, map = world.inp, aes(map_id = region), fill="grey90") +
#   geom_map(data=german.inp, map = german.inp, aes(map_id = region), fill="white")+
#   xlim(-10, 30) +
#   ylim(35, 70) +
#   
#   geom_point(data=occ_clean, aes(x=longitude, y=latitude, color=datasource, shape="."), cex=0.4)+
#   #scale_x_continuous(limits=c(5, 17))+ 
#   #scale_y_continuous(limits=c(46,57))+
#   theme_bw()+
#   theme(legend.position = "bottom")+
#   theme(panel.background = element_rect(fill = "grey70",
#                                         colour = "grey70",
#                                         size = 0.5, 
#                                         linetype = "solid"))+
#   guides(color = guide_legend(override.aes = list(size = 3))) #makes legend icons bigger
# 
# plotOccRawGER
# 
# # pdf(paste0(data_wd, "/_figures/OccurrencesRaw_", Taxon_name, "_perDatasource_GER.pdf")); plotOccRawGER; dev.off()
# 
# rm(plotOccRawGER, occ_clean)
# 
# # Calculate individuals species' occurrence based on BIOMOD input
# occ_points_species <- read.csv(file=paste0(data_wd, "/_results/Occurrence_rasterized_2km_BIOMOD_", Taxon_name, ".csv"))
# head(occ_points_species)
# 
# # number of (true) records per species
# View(occ_points_species %>% group_by(SpeciesID) %>% filter(occ>=1) %>% count())
# View(occ_points_species %>% group_by(SpeciesID) %>% filter(occ>=1 & SpeciesID %in% species100) %>% count())
# occ_points_species %>% group_by(SpeciesID) %>% filter(occ>=1 & SpeciesID %in% species100) %>% count() %>%
#   ungroup() %>% summarize(across(everything(), list(mean, median, sd)))
# 
# # summarize individual species occurrence
# occ_points_species <- occ_points_species %>%
#   mutate("Latitude"=round(x,0), "Longitude"=round(y,0)) %>%
#   group_by(Latitude, Longitude, SpeciesID) %>%
#   summarize( "Records"= sum(occ), .groups="keep") %>%
#   filter(Records > 0)
# 
# # occ_points_species <- occ_points %>%
# #   pivot_longer(cols=speciesNames$SpeciesID[speciesNames$SpeciesID %in% colnames(occ_points)],
# #                names_to = "SpeciesID") %>%
# #   mutate("Latitude"=round(x,0), "Longitude"=round(y,0)) %>%
# #   group_by(Latitude, Longitude, SpeciesID) %>%
# #   filter(!is.na(value)) %>%
# #   summarize("Number of Records"= n(), .groups="keep") %>%
# #   filter("Number of Records" > 0)
# # 
# # # only keep species that will be analyzed (i.e., present in at least 5 grid cells)
# # occ_points_species <- occ_points_species[occ_points_species$SpeciesID %in%
# #                                            speciesNames[speciesNames$NumCells_2km >=5, "SpeciesID"],]
# 
# ## plot some of the individual species' occurrences
# plotOccSpecies <- ggplot()+
#   geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white") +
#   xlim(-10, 30) +
#   ylim(35, 70) +
# 
#   geom_point(data=occ_points_species,
#              aes(x=Latitude, y=Longitude, color=Records, group=SpeciesID),
#              cex=0.07, pch=15)+
#   facet_wrap(vars(SpeciesID))+
#   scale_color_gradientn(    # provide any number of colors
#     colors = c("black", "blue", "orange"),
#     values = scales::rescale(c(1, 5, 20, 30, 50, 100, 300)),
#     breaks = c(5, 20, 50, 100, 200))+
# 
#   # add number of grid cells in which the species is present
#   geom_text(data=occ_points_species %>% group_by(SpeciesID) %>% summarize("n"=sum(Records)),
#             aes(x=-5, y=68, label=paste0("n=", n)), color="black",
#             inherit.aes=FALSE, parse=FALSE, cex=2, hjust=0)+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position =c(0.6, 0.05),legend.direction = "horizontal",
#         axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_rect(fill="grey95"))
# plotOccSpecies
# 
# # png(paste0(data_wd, "/_figures/OccurrencesGridded_", Taxon_name, "_perSpecies.png"), height=4000, width=4000, res=400); plotOccSpecies; dev.off()
# 
# rm(plotOccSpecies, occ_points_species)
# 
# ## Barplot: occurrence per year
# summary(occ_points$year)
# 
# ggplot(data=occ_points, aes(x=year))+
#   geom_bar()+theme_bw()
# 
# 
# ## Calculate number of occurrences per country
# library(rworldmap)
# m <- rworldmap::getMap()
# 
# occ_points_sp <- occ_points
# occ_points_sp <- occ_clean
# 
# coordinates(occ_points_sp) <- ~x+y
# proj4string(occ_points_sp) = proj4string(m)
# 
# # extract country names
# occ_country <- droplevels(over(occ_points_sp,m)$NAME)
# unique(occ_country)
# 
# table(occ_country)

- - - - - - - - - - - - - - - - - - - - -
# Variable importance (MaxEnt) ####
- - - - - - - - - - - - - - - - - - - - -

## Visualize and get top 10
var_imp <- read.csv(file=paste0("_results/Variable_importance_MaxEnt_", Taxon_name, ".csv"))
var_imp

var_imp[var_imp$Predictor=="MAP_Chelsa.EU","Predictor"] <- "MAP"
var_imp[var_imp$Predictor=="MAPseas_Chelsa.EU","Predictor"] <- "MAP_Seas"
var_imp[var_imp$Predictor=="MAT_Chelsa.EU","Predictor"] <- "MAT"
var_imp[var_imp$Predictor=="MATseas_Chelsa.EU","Predictor"] <- "MAT_Seas"

# load predictor table to get classification of variables
# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0("_data/METADATA_Predictors.csv"))

# transform to long format and add variable categories
var_imp <- var_imp %>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")

# add category for clay.silt
var_imp[var_imp$Predictor=="Clay.Silt","Category"] <- "Soil"
var_imp[var_imp$Predictor=="Snow","Category"] <- "Climate"

# summarize mean & SD
View(var_imp %>%
       dplyr::select(-Species, -Category) %>% group_by(Predictor) %>% summarize_all(c(mean, median, sd)))

# plot Permutation importance
plotVarImp <- ggplot()+
  geom_boxplot(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent,mean)),cex=0.2, outlier.size=0.2)+
  geom_point(data=var_imp %>%
               group_by(Predictor, Category) %>%
               summarize(across(maxent, mean)),
             aes(x=maxent, y=reorder(Predictor, maxent, mean)), color="black", fill="grey30",
             size=3, alpha=1) +
  # stat_summary(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent,median)),
  #              geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75, linetype="dashed")+
  geom_jitter(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent, mean),
                                color=Category
                                ), alpha=0.3, height=0.2) +
  xlab("Variable importance (Permutation importance)")+
  ylab("Predictor")+
  geom_hline(yintercept=length(unique(var_imp$Predictor))-9.5, lty=2)+
  theme_bw()+
  theme(axis.text = element_text(size=15), axis.title = element_text(size=15), axis.title.y=element_blank(),
        legend.text = element_text(size=15), legend.title = element_blank(), legend.position = "none")
plotVarImp

png(paste0("_figures/VariableImportance_MaxEnt_", Taxon_name, ".png")); plotVarImp; dev.off()

## [Supplementary Figure 4] stacked barplot all species
var_imp$Predictor <- factor(var_imp$Predictor, levels=c("MAP", "MAP_Seas", "MAT", "Snow",
                                                        "Aspect", "Dist_Coast", "Elev", "Dist_River", "Slope",
                                                        "Agriculture","Dist_Urban", "Forest_Coni", "Forest_Deci", "NDVI", "Pastures", "Pop_Dens", "Shrubland",
                                                        "CEC", "Clay.Silt", "Cu", "Hg","Moisture", "N", "P", "pH", "SOC"))

plotAllVI <- ggplot(var_imp %>%
                      full_join(var_imp %>% filter(Category=="Climate") %>%
                                  dplyr::select(Species, Category, maxent)  %>%
                                  group_by(Species, Category) %>% summarize_all(sum) %>%
                                  group_by(Species) %>% top_n(1, maxent) %>% unique() %>%
                                  dplyr::select(-Category) %>% rename("SumClimate"=maxent), by="Species"),
                    aes(fill=Predictor, alpha=Predictor, y=maxent, x=reorder(Species, SumClimate))) +
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  xlab("Species")+
  scale_y_continuous(expand = c(0, 0))+
  scale_alpha_manual(values=c("MAP"=0.75, "MAP_Seas"=0.55, "MAT"=0.35, "Snow" = 0.25,
                              "Aspect"=0.85, "Dist_Coast"=0.65, "Elev"=0.55,"Dist_Coast"=0.45, "Elev"=0.35, "Dist_River"=0.25, "Slope"=0.15,
                              "Agriculture"=0.85,"Dist_Urban"=0.65, "Forest_Coni"=0.45, "Forest_Deci"=0.25, "NDVI"=0.75, "Pastures"=0.55,  "Pop_Dens"=0.35, "Shrubland"=0.15,
                              "CEC"=0.75,"Clay.Silt"=0.65, "Cu"=0.55, "Hg"=0.45,"Moisture"=0.75, "N"=0.65, "P"=0.55, "pH"=0.45, "SOC"=0.35))+
  scale_fill_manual(values=c("MAP"="#F8766D", "MAP_Seas"="#F8766D", "MAT"="#F8766D", "Snow" = "#F8766D", 
                             "Aspect"="#00BFC4", "Dist_Coast"="#00BFC4", "Elev"="#00BFC4","Dist_Coast"="#00BFC4", "Elev"="#00BFC4", "Dist_River"="#00BFC4", "Slope"="#00BFC4",
                             "Agriculture"="#7CAE00","Dist_Urban"="#7CAE00", "Forest_Coni"="#7CAE00", "Forest_Deci"="#7CAE00", "NDVI"="#698B22", "Pastures"="#698B22",  "Pop_Dens"="#698B22", "Shrubland"="#698B22",
                             "CEC"="#C77CFF","Clay.Silt"="#C77CFF", "Cu"="#C77CFF", "Hg"="#C77CFF","Moisture"="#BF3EFF", "N"="#BF3EFF", "P"="#BF3EFF", "pH"="#BF3EFF", "SOC"="#BF3EFF"))+
  theme_bw()+
  theme(legend.position = "bottom", axis.text = element_text(size=15), axis.title = element_text(size=15),
        legend.text = element_text(size=15), legend.title = element_blank())
plotAllVI

png(paste0("_figures/VariableImportance_maxent_species_", Taxon_name, ".png"), height=800, width=600); plotAllVI; dev.off()


# #- - - - - - - - - - - - - - - - - - - - -
# ## Variable importance (BioMod) ####
# #- - - - - - - - - - - - - - - - - - - - -
# 
# var_imp <- read.csv(file=paste0(data_wd, "/_results/Variable_importance_biomod_", Taxon_name, ".csv"))
# var_imp
# 
# # load predictor table to get classification of variables
# pred_tab <- readr::read_csv(file=paste0(data_wd, "/data_environment/METADATA_Predictors.csv"))
# 
# # transform to long format and add variable categories
# var_imp <- var_imp %>%
#   left_join(pred_tab %>% dplyr::select(Predictor, Category, Long_name) %>%
#               mutate("Long_name"=str_replace_all(Long_name, "_", " ")), by="Predictor")
# 
# # add category for clay.silt
# var_imp[var_imp$Predictor=="Clay.Silt","Category"] <- "Soil"
# var_imp[var_imp$Predictor=="Clay.Silt","Long_name"] <- "Clay and silt content"
# 
# # summarize mean and sd 
# View(var_imp %>% filter(Species %in% unique(speciesNames[speciesNames$NumCells_2km_biomod>=100,"SpeciesID"])) %>% 
#        dplyr::select(-Species, -Category, -Long_name) %>% group_by(Predictor) %>% summarize_all(c(mean, sd)))
# 
# # plot VIF
# plotVarImp <- ggplot(data=var_imp, aes(x=biomod, y=reorder(Long_name, biomod), fill=Category))+
#   geom_boxplot(cex=0.2, outlier.size=1.5, show.legend = F)+
#   guides(fill=(guide_legend(override.aes=list(alpha = 1))))+
#   geom_point(alpha = 0, shape=21, show.legend = T)+
#   geom_jitter(height=0.2, alpha=0.3, show.legend = F)+
#   ylab("")+
#   theme_bw()+
#   theme(axis.text.y = element_text(size = 25), axis.title = element_blank(), axis.text.x = element_text(size=15),
#         legend.position=c(0.75, 0.15), legend.text = element_text(size=20), legend.title = element_blank(),
#         panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
# plotVarImp
# 
# png(paste0(data_wd, "/_figures/VariableImportance_biomod_", Taxon_name, ".png"), height=600, width=700); plotVarImp; dev.off()
# 
# # # plot barplot with top 10
# # plotTopVI <- var_imp %>% dplyr::select(biomod, Predictor, Category) %>%
# #   group_by(Predictor, Category) %>% summarize_all(mean, na.rm=T) %>% arrange(desc(biomod)) %>%
# #   ggplot(aes(x=reorder(Predictor, biomod), y=biomod, fill=Category)) + 
# #   geom_segment(aes(x=reorder(Predictor, biomod), xend=reorder(Predictor, biomod), y=0, yend=biomod), color="black") +
# #   geom_point(aes(color=Category), size=4, alpha=1) +
# #   coord_flip() +
# #   xlab("Predictors")+ylab("Mean variable importance")+
# #   theme_bw()+theme(aspect.ratio=1/1)
# # plotTopVI
# # 
# # png(paste0(data_wd, "/_figures/VariableImportance_biomod_top10_", Taxon_name, ".png")); plotTopVI; dev.off()
# 
# # mean varImp
# var_imp %>% group_by(Predictor) %>% dplyr::select(-Species, -Category, -Long_name) %>% summarize_all(mean)
# 
# # plot varImp of each species
# var_imp$Predictor <- factor(var_imp$Predictor, levels=c("MAP_Seas", "MAT",
#                                                         "Dist_Coast", "Elev",
#                                                         "Agriculture", "Pop_Dens",
#                                                         "CEC", "Clay.Silt", "P", "pH"))
# plotAllVI <- ggplot(var_imp, aes(fill=Predictor, alpha=Predictor, y=biomod, x=reorder(Species, biomod))) + 
#   geom_bar(position="stack", stat="identity")+
#   coord_flip()+
#   xlab("Species")+
#   scale_y_continuous(expand = c(0, 0))+
#   scale_alpha_manual(values=c("MAP_Seas"=0.75, "MAT"=0.5, "Dist_Coast"=0.75, "Elev"=0.5,
#                               "Agriculture"=0.75, "Pop_Dens"=0.5, "CEC"=0.75,"Clay.Silt"=0.55, "P"=0.35, "pH"=0.15))+
#   scale_fill_manual(values=c("MAP_Seas"="#F8766D", "MAT"="#F8766D", "Dist_Coast"="#00BFC4", "Elev"="#00BFC4",
#                              "Agriculture"="#7CAE00", "Pop_Dens"="#7CAE00", "CEC"="#C77CFF","Clay.Silt"="#C77CFF", "P"="#C77CFF", "pH"="#C77CFF"))+
#   theme_bw()+
#   theme(legend.position = "bottom", axis.text = element_text(size=15), axis.title = element_text(size=15),
#         legend.text = element_text(size=15), legend.title = element_blank())
# plotAllVI
# 
# png(paste0(data_wd, "/_figures/VariableImportance_biomod_species_", Taxon_name, ".png"), height=800, width=600); plotAllVI; dev.off()


#- - - - - - - - - - - - - - - - - - - - - -
## Load UNCERTAINTY (relevant for all plots) ####
#- - - - - - - - - - - - - - - - - - - - - -

# load uncertainty extent for all maps
load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, ".RData")) #extent_df

# load uncertainty
uncertain_tif <- terra::rast(paste0("_results/", Taxon_name, "/SDM_Uncertainty_", Taxon_name, ".tif"))

# save threshold for uncertainty
uncertain_thresh <- stats::quantile(uncertain_tif$Mean, 0.9, na.rm=TRUE)

#- - - - - - - - - - - - - - - - - - - - - -
## Uncertainty ####
#- - - - - - - - - - - - - - - - - - - - - -

terra::plot(uncertain_tif$Mean)

png(file=paste0("_figures/Uncertainty_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  #xlim(-10, -5) +
  #ylim(40, 45) +

  geom_tile(data=uncertain_tif %>% terra::as.data.frame(xy=TRUE) %>% dplyr::select(x,y), aes(x=x, y=y))+
  #geom_sf(uncertain_tif %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
  coord_map()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80", color="white") +
  geom_tile(data=uncertain_tif %>% terra::as.data.frame(xy=TRUE) %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
  
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  theme_bw()+

  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="bottom",legend.direction = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()


# # While scalebar in WGS84 makes no sence:
# #https://stackoverflow.com/questions/41371170/ggsn-global-map-scale-bar-looks-wrong/41373569#41373569
#
# uncertain_st <- sf::st_multipoint(uncertain_df %>% dplyr::select(x,y,Mean) %>% as.matrix(), dim="XYZ")
# uncertain_st <- raster::rasterFromXYZ(uncertain_df %>% dplyr::select(x,y,Mean))
# uncertain_st <- raster::rasterToPolygons(uncertain_st)
# uncertain_st <- sf::st_as_sf(uncertain_st)

# png(file=paste0(data_wd, "/_figures/Uncertainty_", Taxon_name, ".png"), width=1000, height=1000)
# ggplot()+
#   geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
#   xlim(-10, 30) +
#   ylim(35, 70) +
#
#   geom_sf(data=uncertain_st, aes(color=Mean))+
#   scale_color_viridis_c(option="E")+
#   north(x.min=-10, x.max=30, y.min=35, y.max=70, symbol=10,
#         location = "topleft")+
#   scalebar(x.min=-10, x.max=30, y.min=35, y.max=70,
#            dist=500, dist_unit="km", st.bottom = TRUE,
#           height=0.05,location = "bottomleft",
#            transform = TRUE, model="WGS84")+
#   ggtitle("Coefficient of variation averaged across SDMs")+
#
#
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position =c(0.2, 0.95), legend.direction = "horizontal",
#         axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())
# dev.off()

summary(uncertain_tif$Mean) #3rd Qu. E: 0.3, N: 0.445

png(file=paste0("_figures/Uncertainty_", round(uncertain_thresh,3), "_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_tile(data=uncertain_tif %>% terra::as.data.frame(xy=TRUE), aes(x=x, y=y, fill=Mean))+
  coord_map()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  geom_tile(data=uncertain_tif %>% terra::as.data.frame(xy=TRUE) %>% filter(Mean<uncertain_thresh), aes(x=x, y=y, fill=Mean))+
  
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  geom_tile(data=uncertain_tif %>% terra::as.data.frame(xy=TRUE) %>% filter(Mean>=uncertain_thresh), aes(x=x, y=y), fill="linen")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position =c(0.2, 0.85),legend.direction = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()


- - - - - - - - - - - - - - - - - - - - - -
## Map species uncertainty maps ####

names(uncertain_tif)
max(uncertain_tif)
min(uncertain_tif)

#pdf(file=paste0("_figures/Uncertainty_allSpecies_", Taxon_name, ".pdf"))
png(file=paste0("_figures/Uncertainty_allSpecies_", Taxon_name, ".png"),width=3000, height=3000)
terra::plot(terra::subset(uncertain_tif, 1:39), maxnl=100,
            plg=list( # parameters for drawing legend
              title = "",
              title.cex = 2, # Legend title size
              cex = 4 # Legend text size
            ), 
            cex.main = 6, # Title text size
            range = c(0, 1),
            axes = FALSE,
            legend = "bottom")
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Species richness (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

load(file=paste0("_results/", Taxon_name, "/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack

# Calculate area with 19 species
summary(species_stack$Richness)
sd(species_stack$Richness)

species_stack %>% filter(Richness==5 & !is.na(Richness)) %>% count()
species_stack %>% filter(Richness==4 & !is.na(Richness)) %>% count()
#(species_stack %>% filter(Richness>=10 & !is.na(Richness)) %>% count())/nrow(species_stack)
#(species_stack %>% filter(Richness>=15 & !is.na(Richness)) %>% count())/nrow(species_stack)

# extract most prominent species
View(as.data.frame(colSums(species_stack, na.rm=T)) %>% arrange(colSums(species_stack, na.rm=T)))

# species richness

png(file=paste0("_figures/SpeciesRichness_cert", round(uncertain_thresh, 3), "_", Taxon_name, ".png"), width=1000, height=1000)
ggplot()+
  geom_tile(data=extent_df %>% mutate(x=round(x, 5), y=round(y,5)), aes(x,y))+
  coord_map()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  geom_tile(data=extent_df %>% mutate(x=round(x, 5), y=round(y,5)) %>% 
              inner_join(species_stack %>% mutate(x=round(x, 5), y=round(y,5)) %>% filter(Richness>0), by=c("x","y")),
            aes(x=x, y=y, fill=Richness))+
  
  ggtitle("Species richness (number of species)")+
  scale_fill_viridis_c()+
  geom_tile(data=extent_df  %>% mutate(x=round(x, 5), y=round(y,5)) %>% 
              inner_join(species_stack  %>% mutate(x=round(x, 5), y=round(y,5)) %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="bottom",legend.direction = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
dev.off()

while (!is.null(dev.list()))  dev.off()

# # Species richness >=18
# ggplot()+
#   geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
#   xlim(-10, 30) +
#   ylim(35, 70) +
# 
#   geom_tile(data=extent_df %>% inner_join(species_stack %>% filter(Richness>0 & Richness>=18), by=c("x","y")),
#             aes(x=x, y=y, fill=Richness))+
#   ggtitle("Species richness (number of species)")+
#   scale_fill_viridis_c()+
#   geom_tile(data=extent_df %>% inner_join(species_stack %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
#   theme_bw()+
#   theme(axis.title = element_blank(), legend.title = element_blank(),
#         legend.position =c(0.2, 0.85),legend.direction = "horizontal",
#         axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
#         legend.text = element_text(size=30), legend.key.size = unit(2, 'cm'),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())

#- - - - - - - - - - - - - - - - - - - - - -
## Species distributions (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

# map binary species distributions
plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
  print(s-2)
  temp_data <- extent_df %>% mutate(x=round(x, 5), y=round(y,5)) %>% 
    inner_join(species_stack[!is.na(species_stack[,s]),] %>% mutate(x=round(x, 5), y=round(y,5)) )
  ggplot()+
    geom_tile(data=temp_data,
              aes(x=x, y=y))+
    coord_map()+
    ggtitle(colnames(species_stack)[s])+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    geom_tile(data=temp_data,
              aes(x=x, y=y, fill=as.factor(temp_data[,s])))+
    
    scale_fill_manual(values=c("1"="#440154","0"="grey60","NA"="lightgrey"))+
    theme_bw()+
    guides(fill = guide_legend(# title.hjust = 1, # adjust title if needed
      label.position = "bottom",
      label.hjust = 0.5))+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.9,0.15), legend.direction = "horizontal",
          legend.key.size = unit(2, 'cm'),
          legend.text = element_text(size=40),
          title = element_text(size=50),
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
})
})

require(gridExtra)
#pdf(file=paste0("_figures/SpeciesDist_bestBinary_", Taxon_name, ".pdf"))
png(file=paste0("_figures/SpeciesDist_bestBinary_", Taxon_name, ".png"), width=3000, height=3300)
do.call(grid.arrange, plots)
dev.off()

while (!is.null(dev.list()))  dev.off()


# #- - - - - - - - - - - - - - - - - - - - -
# ## Model performance ####
# #- - - - - - - - - - - - - - - - - - - - -
# 
# data_eval <- read.csv(paste0("_results/Model_evaluation_", Taxon_name, ".csv"))
# 
# data_eval %>% dplyr::select(-SpeciesID) %>% summarize_all(mean)
# data_eval %>% dplyr::select(-SpeciesID) %>% summarize_all(sd)
# 
# mod_eval <- read.csv(file=paste0(data_wd, "/_results/ModelEvaluation_", Taxon_name, ".csv"))
# 
# # point plot with lables, tss over roc
# pdf(paste0(data_wd, "/_figures/Model_performance_", Taxon_name, "_tss-roc_perSpecies.pdf"))
# ggplot(mod_eval %>% filter(!is.na(species)), aes(x=tss, y=roc, color=model))+
#   geom_text(label=mod_eval[!is.na(mod_eval$species),"model"], nudge_x = 0, nudge_y = 0, check_overlap = F, cex=1)+
#   facet_wrap(vars(species))+
#   xlim(0,1)+
#   theme_bw()
# dev.off()
# 
# # boxplot, tss per algorithm
# pdf(paste0(data_wd, "/_figures/Model_performance", Taxon_name, "_boxplot.pdf"))
# ggplot(mod_eval %>% filter(!is.na(species)), aes(x=tss, y=model))+
#   geom_boxplot()+
#   xlim(0,1)+
#   theme_bw()
# dev.off()

# #- - - - - - - - - - - - - - - - - - - -
# ## Extract GBIF taxon keys ####
# #- - - - - - - - - - - - - - - - - - - -
# 
# dat <- read.csv(file=paste0(data_wd, "/_intermediates/Occurrences_GBIF_Crassiclitellata.csv")) #dat
# dat %>% dplyr::select(species, speciesKey) %>% unique() %>% arrange(species)
# 
# 
# #- - - - - - - - - - - - - - - - - - - -
# ## Some more numbers ####
# #- - - - - - - - - - - - - - - - - - - -
# 
# # relationship between life form (ecological group) and number of records
# bar_groups <- 
#   ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
#          mutate(Ecogroup = ifelse(Ecogroup=="Epi-Endogeic", "Epi-endogeic", Ecogroup)), aes(x=Ecogroup, y=NumCells_2km_biomod, color=SpeciesID))+
#   geom_bar(stat="identity")+
#   #facet_wrap(vars(Ecogroup))+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
# ggsave(filename=paste0(data_wd, "/_figures/NumberOccurrences_ecogroup_", Taxon_name, ".png"), bar_groups)
# 
# bar_groups2 <- 
#   ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
#            mutate(Ecogroup = ifelse(Ecogroup=="Epi-Endogeic", "Epi-endogeic", Ecogroup)), aes(x=Ecogroup, y=Records, color=SpeciesID))+
#   geom_bar(stat="identity")+
#   #facet_wrap(vars(Ecogroup))+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
# ggsave(filename=paste0(data_wd, "/_figures/NumberRecords_ecogroup_", Taxon_name, ".png"), bar_groups2)
# 
# #- - - - - - - - - - - - - - - - - - - -
# # records vs. threat level
# speciesStatus <- read.csv(paste0(here::here(), "/doc/Species_status_", Taxon_name, ".csv"))
# 
# bar_status <- 
#   ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
#            right_join(speciesStatus, by="SpeciesID") %>%
#            mutate(RedList_Germany=ifelse(RedList_Germany=="", NA, RedList_Germany)), aes(x=RedList_Germany, y=NumCells_2km_biomod, color=SpeciesID))+
#   geom_bar(stat="identity")+
#   #facet_wrap(vars(Ecogroup))+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
# ggsave(filename=paste0(data_wd, "/_figures/NumberOccurrences_status_", Taxon_name, ".png"), bar_status)
# 
# bar_status2 <- 
#   ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
#            right_join(speciesStatus, by="SpeciesID") %>%
#            mutate(RedList_Germany=ifelse(RedList_Germany=="", NA, RedList_Germany)), aes(x=RedList_Germany, y=Records, color=SpeciesID))+
#   geom_bar(stat="identity")+
#   #facet_wrap(vars(Ecogroup))+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
# ggsave(filename=paste0(data_wd, "/_figures/NumberRecords_status_", Taxon_name, ".png"), bar_status2)

