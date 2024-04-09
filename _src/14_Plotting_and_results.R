#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Plotting and other results           #
#          author: Romy Zeiss               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

#setwd("D:/EIE_Macroecology/_students/Romy/Atlas_Portugal")

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

# covariates
covarsNames <- paste0("PCA_", 1:11)

# define geographical extent for plotting
extent_portugal <- c(-9, -6, 40.5, 42.3)

# load background map
world.inp <- subset(map_data("world"))#, region == "Portugal")

# function to extract legend
g_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
#- - - - - - - - - - - - - - - - - - - - -
## Running progress ####
#- - - - - - - - - - - - - - - - - - - - -
# extract information of how many taxa are processed, will be processed and are plotted

df_overview <- data.frame(Taxon = "Taxon", 
           No_taxa = 1,
           No_taxa10 = 1, No_taxa100 = 1,
           Data_clean = TRUE, Data_gridded = TRUE, 
           Data_biomod = 1,
           No_varImp = 1,
           No_SDMs = 1, No_ESMs = 1,
           No_predict = 1, No_uncertain = 1,
           No_uncertain_clamp = 1, No_uncertain_cv = 1,
           Plot_occ_clean = TRUE, Plot_occ_grid = TRUE,
           Plot_varImp = TRUE, 
           Plot_uncertain = TRUE,
           Plot_predict = TRUE)[0,]

for(Taxon_name in c("Earthworms", "EarthGenus", "Nematodes", "NemaGenus", "Fungi", "Bacteria", "Eukaryotes", "Protists")){
  print(Taxon_name)
  temp_overview <- data.frame(test=1)
  
  temp_overview[,"Taxon"] <- Taxon_name
  
  # number of taxa
  try({
    species100 <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))
    if(nrow(species100) != 0){
      species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))
      speciesSub <- species100 %>% 
        pull(species) %>%
        c(species10 %>% pull(species))
    }else{
      species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))
      speciesSub <- species10 %>% pull(species)
    }
    temp_overview[,"No_taxa"] <- length(speciesSub)
    temp_overview[,"No_taxa10"] <- length(species10 %>% pull(species))
    temp_overview[,"No_taxa100"] <- length(species100 %>% pull(species))
  })
    
  # data prepared?
  temp_overview[,"Data_clean"] <- file.exists(paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))
  temp_overview[,"Data_gridded"] <- file.exists(paste0("_intermediates/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
  temp_overview[,"Data_biomod"] <- length(list.files(paste0("_intermediates/BIOMOD_data/", Taxon_name)))
  
  # how many varImp ready?
  temp_overview[,"No_varImp"] <- length(list.files(paste0("_results/_TopPredictor/", Taxon_name), 
                                                  pattern = "SDM_"))
  
  # how many SDMs and ESMs ready?
  temp_overview[,"No_SDMs"] <- length(list.files(paste0("_results/", Taxon_name, "/SDMs"), 
                                                  pattern = "SDM_"))
  temp_overview[,"No_ESMs"] <- length(list.files(paste0("_results/", Taxon_name, "/SDMs"), 
                                                pattern = "ESM_"))
  
  # how many predictions & uncertainty ready?
  temp_overview[,"No_predict"] <- length(list.files(paste0("_results/", Taxon_name, "/Projection")))
  temp_overview[,"No_uncertain_clamp"] <- length(list.files(paste0("_results/", Taxon_name, "/Uncertainty"),
                                                           pattern = "ClampingMask_"))
  temp_overview[,"No_uncertain_cv"] <- length(list.files(paste0("_results/", Taxon_name, "/Uncertainty"),
                                                        pattern = "CV_"))
  temp_overview[,"No_uncertain"] <- as.numeric(temp_overview["No_uncertain_clamp"]) + as.numeric(temp_overview["No_uncertain_cv"])

  # visualization ready?
  temp_overview[,"Plot_occ_clean"] <- file.exists(paste0("_figures/OccurrencesCleaned_", Taxon_name, ".pdf"))
  temp_overview[,"Plot_occ_grid"] <- file.exists(paste0("_figures/OccurrencesGridded_", Taxon_name, ".pdf"))
  temp_overview[,"Plot_varImp"] <- file.exists(paste0("_figures/VariableImportance_MaxEnt_", Taxon_name, ".pdf"))
  temp_overview[,"Plot_uncertain"] <- file.exists(paste0("_figures/Uncertainty_", Taxon_name, ".pdf"))
  temp_overview[,"Plot_predict"] <- file.exists(paste0("_figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"))
  
  # add to summary df
  df_overview <- full_join(df_overview, temp_overview %>% dplyr::select(-test))
  
}

df_overview


#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Nematodes"

# load number of occurrences per species and focal species names
species100 <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))
if(nrow(species100) != 0){
  species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))
  speciesSub <-  species100 %>% 
    pull(species) %>%
    c(species10 %>% pull(species))
}else{
  species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv")) %>% pull(species)
  speciesSub <- species10
}

species10
species100
speciesSub


#- - - - - - - - - - - - - - - - - - - - -
## Occurrences ####
#- - - - - - - - - - - - - - - - - - - - -
### Cleaned occurrences ####
#- - - - - - - - - - - - - - - - - - - - -
# load cleaned occurrence records 
occ_clean <- read_csv(file=paste0("_intermediates/Occurrences_clean_", Taxon_name, ".csv"))

## plot raw occurrences colored by year
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white", color = "grey80")+
  
  geom_point(data=occ_clean, 
             aes(x=x, y=y, color=as.factor(Presence), size=Abundance))+ 
  theme_bw()+
  coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]))+
  #scale_color_binned()+
  scale_color_manual(values = c("1"="forestgreen", "0"="orange"))+
  facet_wrap(vars(SpeciesID))+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        #legend.position = c(0.85,0.2),legend.direction = "vertical",
        #legend.box = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), #strip.text = element_text(size=5),
        panel.border = element_blank(),
        panel.background = element_rect(fill="grey80"))
ggsave(paste0("_figures/OccurrencesCleaned_", Taxon_name, ".pdf"),
       last_plot(),
       height = 8, width = 12)


# ## plot maps
# temp_files <- list.files(paste0("_results/_TopPredictor/", Taxon_name))
# temp_files <- temp_files[stringr::str_detect(temp_files, "SDM_maxent_[:graph:]*.RData")]
# 
# plots <- lapply(c(1:length(temp_files)), function(m) {try({
#   print(temp_files[m]); print(m)
#   temp_pred <- get(load(file=paste0("_results/_TopPredictor/",Taxon_name, "/", temp_files[m])))[["prediction"]]
#   #print(m)
#   spID <- str_split_fixed(temp_files[m], "_", 4)
#   spID <- str_split_fixed(spID[,4], "[.]", 2)[,1]
#   ggplot(data=temp_pred, aes(x=x, y=y, fill=layer))+
#     geom_tile()+
#     ggtitle(spID)+
#     scale_fill_viridis_c(limits = c(0,1))+
#     theme_bw()+
#     theme(axis.title = element_blank(), legend.title = element_blank(),
#           legend.position = "none")
# })})
# 
# lapply(plots, class) # if error, remove that species
# 
# pdf(file=paste0("_figures/DistributionMap_MaxEnt_", Taxon_name, ".pdf"), height = 10, width = 10)
# #png(file=paste0("_figures/DistributionMap_MaxEnt_", Taxon_name, ".png"),width=3000, height=3000)
# do.call(gridExtra::grid.arrange, plots)
# dev.off()


#- - - - - - - - - - - - - - - - - - - - -
### Cleaned, gridded Occurrences ####
#- - - - - - - - - - - - - - - - - - - - -
# load matrix containing information on number of occurrence records in grid
occ_points <- read_csv(file=paste0("_intermediates/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
#occ_points <- occ_points %>% rename("x"=?..x)

## plot raw occurrences colored by year
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white", color = "grey80")+
  
  geom_point(data=occ_points %>% mutate(rowsum = rowSums(.[3:ncol(.)])), 
             aes(x=x, y=y, color=rowsum, size=rowsum))+ 
  theme_bw()+
  geom_text(aes(x=mean(c(extent_portugal[1], extent_portugal[2])), 
                y=extent_portugal[4]), 
            label = Taxon_name, cex = 10)+
  coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]))+
  scale_color_binned()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position =c(0.85,0.2),legend.direction = "vertical",
        legend.box = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), strip.text = element_text(size=5),
        panel.border = element_blank(),
        panel.background = element_rect(fill="grey80"))
ggsave(paste0("_figures/OccurrencesGridded_", Taxon_name, ".pdf"),
       last_plot(),
       height = 8, width = 12)

#- - - - - - - - - - - - - - - - - - - - -
### BIOMOD Occurrences per species ####
#- - - - - - - - - - - - - - - - - - - - -
# Calculate individuals species' occurrence based on BIOMOD input
occ_points_species <- read_csv(file=paste0("_results/Occurrence_rasterized_1km_BIOMOD_", Taxon_name, ".csv"))
head(occ_points_species)

# number of (true) records per species
#View(occ_points_species %>% group_by(SpeciesID) %>% filter(occ>=1) %>% count())
#View(occ_points_species %>% group_by(SpeciesID) %>% filter(occ>=1 & SpeciesID %in% species100) %>% count())

## plot some of the individual species' occurrences
# Note: too large vector for Bacteria
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white", color = "grey80") +
  coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]))+

  geom_point(data=occ_points_species,
             aes(x=x, y=y, group=SpeciesID, color = as.factor(occ)),
             cex=0.5, pch=19, alpha = 0.7)+
  scale_color_manual(values = c("1"="forestgreen", "0"="orange"))+
  facet_wrap(vars(SpeciesID))+

  # add number of grid cells in which the species is present
  geom_text(data=occ_points_species %>% group_by(SpeciesID) %>% summarize("n"=sum(occ)),
            aes(x=-5, y=68, label=paste0("n=", n)), color="black",
            inherit.aes=FALSE, parse=FALSE, cex=2, hjust=0)+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = "bottom",legend.direction = "horizontal",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill="grey80"))

ggsave(paste0("_figures/OccurrencesGridded_", Taxon_name, "_perSpecies.pdf"),
       last_plot(),
       height = 12, width = 12)

rm(occ_points_species)

#- - - - - - - - - - - - - - - - - - - - -
## Variable importance (MaxEnt) ####
#- - - - - - - - - - - - - - - - - - - - -

## Visualize and get top 10
var_imp <- read_csv(file=paste0("_results/Variable_importance_MaxEnt_", Taxon_name, ".csv"))
var_imp

# load predictor table to get classification of variables
# load the predictor table containing the individual file names
pred_tab <- readr::read_csv(file=paste0("_data/METADATA_Predictors.csv"))

# transform to long format and add variable categories
var_imp <- var_imp %>%
  left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor") %>%
  filter(!is.na(Predictor))

# add category for clay.silt
var_imp[var_imp$Predictor=="Clay.Silt","Category"] <- "Soil"

# summarize mean & SD
# View(var_imp %>% 
#        dplyr::select(-Species, -Category) %>% group_by(Predictor) %>% summarize_all(list(mean=mean, sd=sd)))

# plot Permutation importance
ggplot()+
  geom_boxplot(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent, median)),cex=0.2, outlier.size=0.2)+
  # stat_summary(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent,median)),
  #              geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75, linetype="dashed")+
  geom_jitter(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent, median), 
                                color=Category), height=0.2, alpha = 0.5) +
  geom_point(data=var_imp %>% 
               group_by(Predictor, Category) %>% summarize(across(maxent, mean)), 
             aes(x=maxent, y=reorder(Predictor, maxent, median)), color="black", fill="grey30",
             size=3, alpha=1) +
  xlab("Variable importance (Permutation importance)")+
  ylab("Predictor")+
  geom_hline(yintercept=length(unique(var_imp$Predictor))-9.5, lty=2)+
  theme_bw()+
  theme(axis.text = element_text(size=15), 
        axis.title = element_text(size=15), axis.title.y=element_blank(),
        legend.text = element_text(size=15), legend.title = element_blank(), 
        legend.position = "none")
ggsave(paste0("_figures/VariableImportance_MaxEnt_", Taxon_name, ".pdf"), 
       last_plot(),
       height = 8, width = 8)

## stacked barplot all species
var_imp$Predictor <- factor(var_imp$Predictor, levels=c("MAP_Chelsa.EU", "MAPseas_Chelsa.EU", "MAT_Chelsa.EU", "MATseas_Chelsa.EU", "Snow",
                                                        "Aspect", "Dist_Coast", "Elev", "Dist_River", "Slope",
                                                        "Agriculture","Dist_Urban", "Forest_Coni", "Forest_Deci", "NDVI", "Pastures", "Pop_Dens", "Shrubland",
                                                        "CEC", "Clay.Silt", "Cu", "Hg","Moisture", "N", "P", "pH", "SOC"))

ggplot(var_imp %>%
          full_join(var_imp %>% filter(Category=="Climate") %>%
                      dplyr::select(Species, Category, maxent)  %>%
                      group_by(Species, Category) %>% summarize_all(sum) %>%
                      group_by(Species) %>% top_n(1, maxent) %>% unique() %>%
                      dplyr::select(-Category) %>% rename("SumClimate"=maxent), by="Species"),
        aes(fill=Predictor, alpha=Predictor, y=maxent, x=reorder(Species, SumClimate))) +
  geom_bar(position="stack", stat="identity")+
  coord_flip()+
  xlab("Taxa")+
  scale_y_continuous(expand = c(0, 0))+
  scale_alpha_manual(values=c("MAP_Chelsa.EU"=0.75, "MAPseas_Chelsa.EU"=0.55, "MAT_Chelsa.EU"=0.35, "Snow" = 0.20,
                              "Aspect"=0.85, "Dist_Coast"=0.65, "Elev"=0.55,"Dist_Coast"=0.45, "Elev"=0.35, "Dist_River"=0.25, "Slope"=0.15,
                              "Agriculture"=0.85,"Dist_Urban"=0.65, "Forest_Coni"=0.45, "Forest_Deci"=0.25, "NDVI"=0.75, "Pastures"=0.55,  "Pop_Dens"=0.35, "Shrubland"=0.15,
                              "CEC"=0.75,"Clay.Silt"=0.65, "Cu"=0.55, "Hg"=0.45,"Moisture"=0.75, "N"=0.65, "P"=0.55, "pH"=0.45, "SOC"=0.35))+
  scale_fill_manual(values=c("MAP_Chelsa.EU"="#F8766D", "MAPseas_Chelsa.EU"="#F8766D", "MAT_Chelsa.EU"="#F8766D", "Snow" = "#F8766D",
                             "Aspect"="#00BFC4", "Dist_Coast"="#00BFC4", "Elev"="#00BFC4","Dist_Coast"="#00BFC4", "Elev"="#00BFC4", "Dist_River"="#00BFC4", "Slope"="#00BFC4",
                             "Agriculture"="#7CAE00","Dist_Urban"="#7CAE00", "Forest_Coni"="#7CAE00", "Forest_Deci"="#7CAE00", "NDVI"="#698B22", "Pastures"="#698B22",  "Pop_Dens"="#698B22", "Shrubland"="#698B22",
                             "CEC"="#C77CFF","Clay.Silt"="#C77CFF", "Cu"="#C77CFF", "Hg"="#C77CFF","Moisture"="#BF3EFF", "N"="#BF3EFF", "P"="#BF3EFF", "pH"="#BF3EFF", "SOC"="#BF3EFF"))+
  theme_bw()+
  theme(legend.position = "bottom", axis.text = element_text(size=15), axis.title = element_text(size=15),
        legend.text = element_text(size=15), legend.title = element_blank())

ggsave(paste0("_figures/VariableImportance_MaxEnt_species_", Taxon_name, ".pdf"), height=11, width=9)
 
# ## plot maps
# temp_files <- list.files(paste0("_results/_TopPredictor/", Taxon_name))
# temp_files <- temp_files[stringr::str_detect(temp_files, "SDM_maxent_[:graph:]*.RData")]
# 
# plots <- lapply(c(1:length(temp_files)), function(m) {try({
#   print(temp_files[m]); print(m)
#   temp_pred <- get(load(file=paste0("_results/_TopPredictor/",Taxon_name, "/", temp_files[m])))[["prediction"]]
#   #print(m)
#   spID <- str_split_fixed(temp_files[m], "_", 4)
#   spID <- str_split_fixed(spID[,4], "[.]", 2)[,1]
#   ggplot(data=temp_pred, aes(x=x, y=y, fill=layer))+
#     geom_tile()+
#     ggtitle(spID)+
#     scale_fill_viridis_c(limits = c(0,1))+
#     theme_bw()+
#     theme(axis.title = element_blank(), legend.title = element_blank(),
#           legend.position = "none")
# })})
# 
# lapply(plots, class) # if error, remove that species
# 
# pdf(file=paste0("_figures/DistributionMaps_MaxEnt_", Taxon_name, ".pdf"), height = 10, width = 10)
# #png(file=paste0("_figures/DistributionMaps_MaxEnt_", Taxon_name, ".png"),width=3000, height=3000)
# do.call(gridExtra::grid.arrange, plots)
# dev.off()


# #- - - - - - - - - - - - - - - - - - - - -
# ## Calculate variable importance for richness (lm) ####
# #- - - - - - - - - - - - - - - - - - - - -
# 
# # load environmental variables (for projections) and species stack as dataframe
# load(paste0(data_wd, "/_results/EnvPredictor_5km_df_clipped.RData")) #Env_clip_df
# load(file=paste0(data_wd, "/_results/_Maps/SDM_stack_bestPrediction_binary_", Taxon_name, ".RData")) #species_stack
# 
# data_stack <- species_stack %>% full_join(Env_clip_df)
# 
# lm1 <- lm(data=data_stack, Richness~MAT+Dist_Coast+MAP_Seas+CEC+Elev+P+Pop_Dens+Agriculture+pH+Clay.Silt)
# summary(lm1)
# confint(lm1)
# 
# lm_varImp <- data.frame("t_value"=summary(lm1)[["coefficients"]][,"t value"])
# lm_varImp$Predictor <- rownames(lm_varImp)
# lm_varImp <- lm_varImp %>% filter(Predictor != "(Intercept)")
# lm_varImp$t_abs <- abs(lm_varImp$t_value)
# lm_varImp$Direction <- factor(sign(lm_varImp$t_value), 1:(-1), c("positive", "neutral", "negative"))
# 
# # transform to long format and add variable categories
# lm_varImp <- lm_varImp%>%
#   left_join(pred_tab %>% dplyr::select(Predictor, Category), by="Predictor")
# 
# # add category for clay.silt
# lm_varImp[lm_varImp$Predictor=="Clay.Silt","Category"] <- "Soil"
# 
# plotTopVI <- lm_varImp %>% dplyr::select(t_abs, Predictor, Category, Direction) %>% arrange(desc(t_abs)) %>%
#   ggplot(aes(x=reorder(Predictor, t_abs), y=t_abs, fill=Category)) + 
#   geom_segment(aes(x=reorder(Predictor, t_abs), xend=reorder(Predictor, t_abs), y=0, yend=t_abs, lty=Direction), color="black") +
#   geom_point(aes(color=Category), size=4, alpha=1) +
#   coord_flip() +
#   xlab("Predictors")+ylab("Variable importance (SR)")+
#   theme_bw()+theme(aspect.ratio=1/1)
# plotTopVI
# 
# png(paste0(data_wd, "/_figures/VariableImportance_biomod_lm_", Taxon_name, ".png")); plotTopVI; dev.off()
# 
# # save model summary
# sink(paste0(data_wd, "/_results/Summary_lm1_Crassiclitellata_varImp.txt"))
# print(summary(lm1))
# sink()

#- - - - - - - - - - - - - - - - - - - - - -
## Load UNCERTAINTY (relevant for all plots) ####
#- - - - - - - - - - - - - - - - - - - - - -

# load uncertainty extent for all maps
load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, ".RData")) #extent_df

# load uncertainty
uncertain_tif <- terra::rast(paste0("_results/SDM_Uncertainty_", Taxon_name, ".tif")) 
uncertain_df <- terra::as.data.frame(uncertain_tif, xy=TRUE)

#- - - - - - - - - - - - - - - - - - - - - -
## Uncertainty ####
#- - - - - - - - - - - - - - - - - - - - - -
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  geom_tile(data=uncertain_df %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
  ggtitle("Coefficient of variation averaged across SDMs")+
  coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]))+
  
  scale_fill_viridis_c(option="E")+
  theme_bw()+  
  
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position ="right",legend.direction = "vertical",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggsave(file=paste0("_figures/Uncertainty_", Taxon_name, ".pdf"), height = 5, width = 8)

# uncertainty threshold
if(Taxon_name == "Earthworms") temp_thresh <- 452.75
if(Taxon_name == "Nematodes") temp_thresh <- 227.9611
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]))+
  
  geom_tile(data=uncertain_df %>% filter(Mean<temp_thresh), aes(x=x, y=y, fill=Mean))+
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  geom_tile(data=uncertain_df %>% filter(Mean>=temp_thresh), aes(x=x, y=y), fill="linen")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = "right",legend.direction = "vertical",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggsave(file=paste0("_figures/Uncertainty_", Taxon_name, "_", temp_thresh, ".pdf"), height = 5, width = 8)

#- - - - - - - - - - - - - - - - - - - - - -
## Map species uncertainty maps ####
# to plot uncertainty of species with temp_thresh, add commented text and
# change limits of scale_fill to c(0, temp_thresh) in both plots and legend

plots <- lapply(3:(ncol(uncertain_df)-2), function(s) {try({
  print(s-2)
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                    ylim = c(extent_portugal[3], extent_portugal[4]))+
    
    geom_tile(data=uncertain_df[!is.na(uncertain_df[,s]) & uncertain_df[,s]>0,], 
              aes(x=x, y=y, 
                  fill=uncertain_df[!is.na(uncertain_df[,s]) & uncertain_df[,s]>0,s]))+
    ggtitle(colnames(uncertain_df)[s])+
    annotate(geom="text", x=-3, y=68, label=colnames(uncertain_df)[s], color="black", size=15)+
    scale_fill_gradient2(midpoint = temp_thresh, limits = c(0,1000),
                         low = "blue",
                         high = "red")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "none",
          axis.line = element_blank(), axis.text = element_blank(), 
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  })
})

legend <- g_legend(ggplot(data=uncertain_df[!is.na(uncertain_df[,3]) & uncertain_df[,3]>0,], 
                          aes(x=x, y=y, fill=uncertain_df[!is.na(uncertain_df[,3]),3]))+
                     geom_tile()+
                     scale_fill_gradient2(midpoint = temp_thresh, limits = c(0,1000),
                                          low = "blue",
                                          high = "red")+
                     theme_bw()+
                     theme(axis.title = element_blank(), legend.title = element_blank(),
                           legend.position = c(0.5, 0.5), legend.direction = "horizontal",
                           legend.text = element_text(size=20), legend.key.size = unit(1.5, 'cm'),
                           axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
                           panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank(),
                           panel.border = element_blank(),
                           panel.background = element_blank())) 

plots2 <- c(plots, list(legend))

pdf(file=paste0("_figures/Uncertainty_allSpecies_", Taxon_name, ".pdf"), height = 10, width = 10)
#png(file=paste0("_figures/Uncertainty_allSpecies_", Taxon_name, ".png"),width=3000, height=3000)
do.call(grid.arrange, plots2)
dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Species richness (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif")) 
species_stack <- terra::as.data.frame(species_stack, xy = TRUE)

summary(species_stack$Richness)
sd(species_stack$Richness)

# species richness
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]))+
  
  geom_tile(data=extent_df %>% inner_join(species_stack %>% filter(Richness>0), by=c("x","y")), 
            aes(x=x, y=y, fill=Richness))+
  ggtitle("Species richness (number of species)")+
  scale_fill_viridis_c()+
  geom_tile(data=extent_df %>% inner_join(species_stack %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = "right",legend.direction = "vertical",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=30), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggsave(file=paste0("_figures/SpeciesRichness_", Taxon_name, ".pdf"), 
       last_plot(),
       height = 5, width = 8)

#- - - - - - - - - - - - - - - - - - - - - -
## Species distributions (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

# map binary species distributions
plots <- lapply(3:(ncol(species_stack)-1), function(s) {try({
  print(s-2)
  temp_data <- extent_df %>% inner_join(species_stack[!is.na(species_stack[,s]),])
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                    ylim = c(extent_portugal[3], extent_portugal[4]))+
    
    geom_tile(data=temp_data, 
              aes(x=x, y=y, fill=as.factor(temp_data[,s])))+
    ggtitle(colnames(species_stack)[s])+
    scale_fill_manual(values=c("1"="#440154","0"="grey60","NA"="lightgrey"))+
    theme_bw()+
    guides(fill = guide_legend(# title.hjust = 1, # adjust title if needed
      label.position = "bottom",
      label.hjust = 0.5))+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = c(0.9, 0.1), legend.direction = "horizontal",
          legend.key.size = unit(0.2, 'cm'),
          legend.text = element_text(size=30),
          title = element_text(size=30),
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
})
})

pdf(file=paste0("_figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"))
#png(file=paste0("_figures/DistributionMap_bestBinary_", Taxon_name, ".png"),width=3000, height=3300)
do.call(grid.arrange, plots)
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

data_eval <- read_csv(paste0("_results/Model_evaluation_", Taxon_name, ".csv"))

data_eval %>% dplyr::select(-Species) %>% summarize_all(mean)
data_eval %>% dplyr::select(-Species) %>% summarize_all(sd)

# point plot with lables, tss over roc
ggplot(data_eval, 
       aes(x=MaxTSS, y=AUC, group=Species))+
  #geom_point()+
  geom_text(aes(label = Species), nudge_x = 0, nudge_y = 0, check_overlap = F, cex=2)+
  #facet_wrap(vars(Species))+
  xlim(0,1)+
  theme_bw()+
  theme(legend.position = "none")
ggsave(paste0("_figures/Model_performance_", Taxon_name, ".pdf"))

# boxplot, tss per algorithm
ggplot(data_eval, aes(x=MaxTSS))+
  geom_boxplot()+
  geom_jitter(aes(y=0), height = 0.1, color = "grey80")+
  xlim(0,1)+
  coord_flip()+
  theme_bw()
ggsave(paste0("_figures/Model_performance_", Taxon_name, "_boxplot.pdf"),
       height = 4, width = 2)


#- - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - -
# UPDATING to be continued from here ####
#- - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - -
## Protection of ranges ####
#- - - - - - - - - - - - - - - - - - - - -

cover_df_current <- read.csv(file=paste0(data_wd, "/_results/ProtectionStatus_current_", Taxon_name, ".csv"))
cover_df_current$IUCNcat <- factor(cover_df_current$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))
cover_sr_current <- read.csv(file=paste0(data_wd, "/_results/ProtectionStatus_SR_future_", Taxon_name, ".csv")) 
cover_sr_current$IUCNcat <- factor(cover_sr_current$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# boxplot of percent area covered by PA overall
a <- ggplot(data=cover_sr_current %>%  filter(IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=reorder(IUCNcat, IUCNcat), y=current_mean, fill=IUCNcat))+
  geom_violin(width=1.4, alpha=0.7)+
  ggtitle("Current")+
  # geom_boxplot(width=0.1, color="black", fill="white", alpha=1)+
  stat_summary(fun = "mean",geom = "point",color = "black", size=3.5, show.legend = F)+
  geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                fun.args = list(mult = 1),
                position =  position_dodge(width = 0.9),
                width=0.1) +
  #geom_jitter(alpha=0.6, width=0.2)+
  theme_bw()+ 
  xlab("Type of protected area (IUCN categories)")+ ylab("Number of species")+
  scale_fill_manual(values=c("olivedrab1","olivedrab3", "olivedrab4", "darkolivegreen", "goldenrod3","goldenrod1", "gold1",
                             "lemonchiffon2","gainsboro", "grey","white" ))+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position="left", axis.text.x=element_text(angle=30, hjust=1),
        legend.text = element_text(size=5), legend.title = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())

# bar chart of percent area covered by PA per species
b <- ggplot(data=cover_df_current %>% group_by(SpeciesID, IUCNcat) %>% 
              summarize(across(c(coverage, coverage_km2), mean)) %>%
              filter(IUCNcat != "Presence" & IUCNcat != "Unprotected" &
                       IUCNcat != " Protected" & IUCNcat != "Outside.PA" &
                       IUCNcat != "ProtectedAreas"& SpeciesID!="_SD") %>%
              mutate(IUCNcat = factor(IUCNcat, levels = c("Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported")))%>% 
                filter(!is.na(IUCNcat)), 
            
            aes(y=coverage, x=reorder(SpeciesID, desc(SpeciesID)), fill=reorder(IUCNcat, desc(IUCNcat))))+
  geom_bar(position="stack", stat="identity")+
  theme_bw()+
  xlab("Species")+ ylab("Proportion of range covered by protected area network")+
  coord_flip()+
  #scale_fill_viridis_d()+
  scale_fill_manual(values=rev(c("olivedrab1","olivedrab3", "olivedrab4", "darkolivegreen", "goldenrod3","goldenrod1", "gold1",
                                 "lemonchiffon2","gainsboro", "grey" )))+
  scale_y_continuous(expand=c(0,0), limits=c(0,0.8))+
  geom_vline(xintercept=19.5, lty=2)+
  theme(legend.position="none",legend.text = element_text(size=5), legend.title = element_text(size=5),
        axis.text.y = element_text(size=15),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

# png(paste0(data_wd, "/_figures/ProtectionStatus_current_", Taxon_name, ".png"))
# grid.arrange(a, b, heights=c(1,1))
# dev.off()

## Load future protection
cover_df_future <- read.csv(file=paste0(data_wd, "/_results/ProtectionStatus_future_", Taxon_name, ".csv"))

cover_sr_future <- read.csv(file=paste0(data_wd, "/_results/ProtectionStatus_SR_future_", Taxon_name, ".csv"))

cover_df_future$SpeciesID <- substr(cover_df_future$SpeciesID, 1, 10)

cover_df <- cover_df_future %>% full_join(cover_df_current %>% mutate("SSP"="current"))

cover_df$IUCNcat <- factor(cover_df$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))
cover_sr_future$IUCNcat <- factor(cover_sr_future$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# save for Supplementary [Appendix S13]
write.csv(cover_df, paste0(data_wd, "/_results/ProtectionStatus_current+future_", Taxon_name, ".csv"),
          row.names = F)

# boxplot of percent area covered by PA overall
boxplot_template <- function(data, col_name){
  data$layer <- c(as.numeric(data[,paste0(col_name, "_mean")] %>% unlist()))
  ggplot(data=data %>%  filter(IUCNcat!="Presence" & IUCNcat!="Protected"), aes(x=reorder(IUCNcat, IUCNcat), y=layer, fill=IUCNcat))+
    geom_violin(width=1.4, alpha=0.7)+
    ggtitle(col_name)+
    # geom_boxplot(width=0.1, color="black", fill="white", alpha=1)+
    stat_summary(fun = "mean",geom = "point",color = "black", size=3.5, show.legend = FALSE)+
    geom_errorbar(stat = "summary", fun.data = "mean_sdl", 
                  fun.args = list(mult = 1),
                  position =  position_dodge(width = 0.9),
                  width=0.1) +#geom_jitter(alpha=0.6, width=0.2)+
    theme_bw()+ 
    xlab("Type of protected area (IUCN categories)")+ ylab("Number of species")+
    scale_fill_manual(values=c("olivedrab1","olivedrab3", "olivedrab4", "darkolivegreen", "goldenrod3","goldenrod1", "gold1",
                               "lemonchiffon2","gainsboro", "grey","white" ))+
    scale_y_continuous(expand=c(0,0))+
    theme(legend.position="right", axis.text.x=element_text(angle=30, hjust=1),
          panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
}

#a <- boxplot_template(cover_sr_future, "current")
#c <- boxplot_template(cover_sr_future, "ssp126")
#d <- boxplot_template(cover_sr_future, "ssp370")
#e <- boxplot_template(cover_sr_future, "ssp585")

cover_sr_matrix <- cover_sr_future %>% pivot_longer(cols=ssp126_mean:ssp585_mean, names_to = "SSP", values_to = "SR") %>%
  mutate("name"="Species richness")
cover_sr_matrix$SR_change <- cover_sr_matrix$SR - cover_sr_matrix$current_mean

cover_sr_matrix %>% dplyr::select(-x, -y, -name) %>% filter(SSP!="current" & IUCNcat!="Presence") %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_sr_matrix %>% dplyr::select(-x, -y, -name)%>% filter(SSP!="current" & IUCNcat!="Presence") %>% group_by(IUCNcat) %>% summarize_all(sd)
cover_sr_matrix %>% dplyr::select(-x, -y, -name) %>% filter(SSP!="current" & IUCNcat!="Presence") %>% group_by(IUCNcat) %>% summarize_all(max)
cover_sr_matrix %>% dplyr::select(-x, -y, -name) %>% filter(SSP!="current" & IUCNcat!="Presence") %>% group_by(IUCNcat) %>% summarize_all(min)


# change in protection in species richness
c <- ggplot(cover_sr_matrix %>% dplyr::select(-x, -y) %>% 
              filter(SSP!="current" & IUCNcat!="Presence" & IUCNcat!="Protected") %>%
              group_by(SSP, IUCNcat, name) %>% summarize(across(everything(), mean)), aes(x=SSP, y=name))+
  geom_tile(aes(fill=SR_change))+
  scale_fill_gradient2(high="dodgerblue3", low="brown3", mid="white", name="No. species gained (+) or lost(-)        ")+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  facet_wrap(vars(IUCNcat), nrow=1, strip.position="bottom")+
  theme(axis.text.x = element_text(angle=30, hjust=1),
        legend.position="top", legend.text = element_text(size=10), legend.title = element_text(size=5),
        axis.text.y = element_text(size=10), axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),strip.background = element_blank(),  panel.spacing = unit(0, "lines"), 
        panel.border = element_rect(color="grey"), strip.text = element_text(angle=30, hjust=1))

# calculate change in protection
cover_matrix <- cover_df %>% 
  mutate(SpeciesID=substr(SpeciesID, 1, 10)) %>%
  full_join(cover_df %>% filter(SSP=="current") %>% dplyr::select(-SSP), by=c("SpeciesID", "IUCNcat"),
            suffix = c("", "_current"))
cover_matrix$coverage_change <- (cover_matrix$coverage_km2 - cover_matrix$coverage_km2_current) / cover_matrix$coverage_km2_current

#cover_matrix <- cover_matrix %>% pivot_wider(id_cols=c(SpeciesID, SSP), names_from=IUCNcat, values_from=coverage_change)

#temp_matrix <- cover_matrix %>% filter(SSP=="ssp126") %>% dplyr::select(-SSP) %>% as.data.frame()
#rownames(temp_matrix) <- temp_matrix$SpeciesID

# plot heatmap-like matrix to show change
cover_matrix$IUCNcat <- factor(cover_matrix$IUCNcat, level=c("Presence", "ProtectedAreas","Ia", "Ib", "II", "III", "IV", "V", "VI", "Not.Applicable", "Not.Assigned", "Not.Reported", "Unprotected", "Protected", "Outside.PA"))

# NOTE: "Unprotected" includes Not.assigned, Not.applicable and Not.reported. 
# We therefore use Outside.PA here, which excludes any protected areas.

f <- ggplot(cover_matrix %>%
              mutate(coverage_change=ifelse(coverage_change>1, 1.1, coverage_change)) %>%
              filter(IUCNcat!="Unprotected" & SSP!="current" & IUCNcat!="Presence" & IUCNcat!="Protected" & SpeciesID!="_SD"), 
              aes(x=SSP, y=reorder(SpeciesID, desc(SpeciesID))))+
  geom_tile(aes(fill=coverage_change*100))+
  scale_fill_gradient2(high="dodgerblue3", low="brown3", mid="white", name="Change in coverage [%]    ",
                       na.value="grey85")+
  theme_bw()+
  ylab("")+
  xlab("")+
  scale_x_discrete(expand = c(0, 0))+
  scale_y_discrete(expand = c(0, 0))+
  facet_wrap(vars(IUCNcat), nrow=1, strip.position="top")+
  theme(axis.text.x = element_text(angle=30, hjust=1),
        legend.position="top", legend.text = element_text(size=10), legend.title = element_text(size=15),
        axis.text.y = element_text(size=10), axis.text.x.bottom = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank(),strip.background = element_blank(),  panel.spacing = unit(0, "lines"), 
        panel.border = element_rect(color="grey"), strip.text = element_blank())

require(gridExtra)
png(paste0(data_wd, "/_figures/ProtectionStatus_heatmap_", Taxon_name, ".png"), height=1000, width=1000)
grid.arrange(a, b, c, f, 
             layout_matrix = rbind(c(1,1,1,3,3),
                                   c(1,1,1,3,3),
                                   #c(1,1,1,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4),
                                   c(2,2,2,4,4)))
dev.off()

#- - - - - - - - - - - - - - - - - - - - 
## Protection: some numbers ####
#- - - - - - - - - - - - - - - - - - - -

protect_df <- extent_df %>% inner_join(protect_df)
sum(protect_df$protected); sum(protect_df$protected) / nrow(extent_df) 

# summarize future protection by IUCN categories across species
cover_matrix %>%
  filter(IUCNcat!="Unprotected" & SSP!="current" & IUCNcat!="Presence" & 
           IUCNcat!="Protected" & SpeciesID!="_SD" & SpeciesID!="_Mean") %>% group_by(IUCNcat) %>% 
  summarize(across(everything(), mean, na.rm=T))
cover_matrix %>%
  filter(IUCNcat!="Unprotected" & SSP!="current" & IUCNcat!="Presence" & 
           IUCNcat!="Protected" & SpeciesID!="_SD" & SpeciesID!="_Mean") %>% group_by(IUCNcat) %>% 
  summarize(across(everything(), sd, na.rm=T))

cover_matrix %>%
  filter(IUCNcat!="Unprotected" & SSP!="current" & IUCNcat!="Presence" & 
           IUCNcat!="Protected" & SpeciesID!="_SD" & SpeciesID!="_Mean") %>%
  dplyr::select(coverage_change) %>%
  summarize(across(everything(), mean, na.rm=T))
cover_matrix %>%
  filter(IUCNcat!="Unprotected" & SSP!="current" & IUCNcat!="Presence" & 
           IUCNcat!="Protected" & SpeciesID!="_SD" & SpeciesID!="_Mean") %>%
  dplyr::select(coverage_change) %>%
  summarize(across(everything(), median, na.rm=T))
cover_matrix %>%
  filter(IUCNcat!="Unprotected" & SSP!="current" & IUCNcat!="Presence" & 
           IUCNcat!="Protected" & SpeciesID!="_SD" & SpeciesID!="_Mean") %>%
  dplyr::select(coverage_change) %>%
  summarize(across(everything(), sd, na.rm=T))

# current mean protected area across species
cover_df %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% filter(IUCNcat=="Protected" & SSP=="current") %>%
  dplyr::select(coverage, coverage_km2, sumCell) %>%
  summarize_all(list("mean"=mean, "sd"=sd))

# current protected area per species
cover_df %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% filter(IUCNcat=="Protected" & SSP=="current") %>% arrange(coverage)

# current mean and SD of protected area across species
cover_df %>% filter(SpeciesID=="_Mean" | SpeciesID=="_SD") %>% arrange(SpeciesID, coverage_km2)

# current mean and SD of protected area across species and scenarios
cover_df %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD" & SSP=="current") %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_df %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD" & SSP=="current") %>% dplyr::select(-SpeciesID) %>% group_by(IUCNcat) %>% summarize_all(sd)

# categories Ia and Ib coverage
min(cover_df[cover_df$IUCNcat=="Ia" & cover_df$SSP=="current",]$coverage)
max(cover_df[cover_df$IUCNcat=="Ia" & cover_df$SSP=="current",]$coverage)
min(cover_df[cover_df$IUCNcat=="Ib" & cover_df$SSP=="current",]$coverage)
max(cover_df[cover_df$IUCNcat=="Ib" & cover_df$SSP=="current",]$coverage)

# SR
cover_sr_current %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_sr_current %>% group_by(IUCNcat) %>% summarize_all(sd)

# SR (future)
cover_sr_matrix %>% filter(SSP!="current") %>% group_by(IUCNcat) %>% summarize_all(mean)
cover_sr_matrix %>% filter(SSP!="current") %>% group_by(IUCNcat) %>% summarize_all(sd)

# Range area
# protected area
cover_matrix %>% filter(IUCNcat=="Protected" & SSP=="current" & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% 
  dplyr::select(-SpeciesID, -IUCNcat, -SSP) %>% summarize_all(mean)
cover_matrix %>% filter(IUCNcat=="Protected" & SSP=="current"  & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% 
  dplyr::select(-SpeciesID, -IUCNcat, -SSP) %>% summarize_all(sd)
cover_matrix %>% filter(IUCNcat=="Protected" & SSP=="current"  & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% 
  dplyr::select(-SpeciesID, -IUCNcat, -SSP) %>% summarize_all(median)

# protected area in future
cover_matrix %>% filter(IUCNcat=="Protected" & SSP!="current" & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% 
  dplyr::select(-SpeciesID, -IUCNcat, -SSP) %>% summarize_all(mean)
cover_matrix %>% filter(IUCNcat=="Protected" & SSP!="current"  & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% 
  dplyr::select(-SpeciesID, -IUCNcat, -SSP) %>% summarize_all(sd)
cover_matrix %>% filter(IUCNcat=="Protected" & SSP!="current"  & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% 
  dplyr::select(-SpeciesID, -IUCNcat, -SSP) %>% summarize_all(min)
cover_matrix %>% filter(IUCNcat=="Protected" & SSP!="current"  & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% 
  dplyr::select(-SpeciesID, -IUCNcat, -SSP) %>% summarize_all(max)

# absolute change in area covered by protected areas
# cover_matrix %>% filter(IUCNcat=="Protected" & SSP!="current" & SpeciesID!="_Mean" & SpeciesID!="_SD") %>%
#   dplyr::select(coverage_change) %>% summarize_all(.funs=list("mean"=function(x)mean(abs(x)),
#                                                               "sd"=function(x)sd(abs(x)),
#                                                               "median"=function(x)median(abs(x)),
#                                                               "min"=function(x)min(x),
#                                                               "max"=function(x)max(x)))

cover_matrix %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID) %>% 
  group_by(IUCNcat) %>% dplyr::select(coverage_change) %>% summarize_all(function(x)mean(x, na.rm=T)) %>% arrange(coverage_change)
cover_matrix %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID) %>% 
  group_by(IUCNcat) %>% dplyr::select(coverage_change) %>% summarize_all(function(x)sd(x, na.rm=T))

# # absolute (no minus)
# cover_matrix %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID) %>% 
#   group_by(IUCNcat) %>% dplyr::select(coverage_change) %>% summarize_all(function(x)mean(abs(x), na.rm=T)) %>% arrange(coverage_change)
# cover_matrix %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID) %>% 
#   group_by(IUCNcat) %>% dplyr::select(coverage_change) %>% summarize_all(function(x)sd(abs(x), na.rm=T)) %>% arrange(coverage_change)
# cover_matrix %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID) %>% 
#   group_by(IUCNcat) %>% summarize_all(max, na.rm=T) %>% arrange(coverage_change)
# cover_matrix %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID) %>% 
#   group_by(IUCNcat) %>% summarize_all(min, na.rm=T) %>% arrange(coverage_change)
# cover_matrix %>% filter(SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID) %>% 
#   group_by(IUCNcat) %>% summarize_all(sd, na.rm=T) %>% arrange(coverage_change)

# more pronounced changes under SSPs
cover_matrix %>% filter(IUCNcat!="Presence" & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID, -IUCNcat) %>%
  ungroup() %>%
  group_by(SSP) %>% summarize_all(function(x) {mean(abs(x), na.rm=T)})
cover_matrix %>% filter(IUCNcat!="Presence" & SpeciesID!="_Mean" & SpeciesID!="_SD") %>% dplyr::select(-SpeciesID, -IUCNcat) %>%
  ungroup() %>%
  group_by(SSP) %>% summarize_all(function(x) {sd(abs(x), na.rm=T)})

# Dendr_illy change across SSPs
View(cover_matrix %>% filter(SpeciesID=="Dendr_illy" & IUCNcat!="Presence") %>% 
  group_by(IUCNcat, SSP))
cover_matrix %>% filter(SpeciesID=="Dendr_illy" & IUCNcat!="Presence") %>% ungroup() %>%
  dplyr::select(-IUCNcat, -SpeciesID) %>% group_by(SSP) %>% summarize_all(mean, na.rm=T)
cover_matrix %>% filter(SpeciesID=="Dendr_illy" & IUCNcat!="Presence")  %>% ungroup() %>%
  dplyr::select(-IUCNcat, -SpeciesID)%>% group_by(SSP) %>% summarize_all(sd, na.rm=T)
cover_matrix %>% filter(SpeciesID=="Dendr_illy" & IUCNcat!="Presence")  %>% ungroup() %>%
  dplyr::select(-SSP, -SpeciesID)%>% group_by(IUCNcat) %>% summarize_all(mean, na.rm=T)
cover_matrix %>% filter(SpeciesID=="Dendr_illy" & IUCNcat!="Presence")  %>% ungroup() %>%
  dplyr::select(-SSP, -SpeciesID)%>% group_by(IUCNcat) %>% summarize_all(mean, na.rm=T)

# # save mean coverage in km2 per IUCN protected area type
# write.csv(cover_df %>% ungroup() %>%
#             group_by(SpeciesID, IUCNcat, SSP) %>% 
#             dplyr::select(SpeciesID, IUCNcat, coverage_km2, SSP) %>% 
#             pivot_wider(names_from=IUCNcat, values_from=coverage_km2), 
#           file=paste0(data_wd, "/_results/ProtectionStatus_current_coveragePerCategory_", Taxon_name, ".csv"), row.names=F)

#- - - - - - - - - - - - - - - - - - - - 
## Map: BD high, PA high and BD low, PA high etc. ####
#- - - - - - - - - - - - - - - - - - - - 

protect_df$Protected <- round(rowSums(protect_df %>% dplyr::select(-x,-y,-Not.Reported, -Not.Assigned, -Not.Applicable), na.rm=T),2)
protect_df[protect_df$Protected>=1 & !is.na(protect_df$Protected), "Protected"] <- 1

biplot_df <- right_join(protect_df %>% dplyr::select(x,y,Protected), 
                        species_stack %>% dplyr::select(x,y,Richness))

# area (of grid) that is currently protected
sum(biplot_df$Protected)
sum(biplot_df$Protected)*5
sum(biplot_df$Protected) / nrow(biplot_df) * 100

biplot_df <- biplot_df %>% 
  full_join(data.frame("Richness"=c(NA,0:19),
                       "Richness_cont"=c(0,0,rep(1,4), rep(5,5), rep(10,5), rep(15,5)))) %>%
  mutate("Earthworm_richness"=as.factor(Richness_cont)) %>%
  filter(Richness>0)

biplot_df$Protection <- "0"
biplot_df[biplot_df$Protected>0 & !is.na(biplot_df$Protected), "Protection"] <- "1" 
biplot_df[biplot_df$Protected>=0.5 & !is.na(biplot_df$Protected), "Protection"] <- "2"
biplot_df[biplot_df$Protected>=0.75 & !is.na(biplot_df$Protected), "Protection"] <- "3"

biplot_df$Protection <- factor(biplot_df$Protection, levels=c(0,1,2,3))
biplot_df$Earthworm_richness <- factor(biplot_df$Earthworm_richness, levels=c(1,5,10,15))

head(biplot_df)

# define fill scale
biplot_input <- biscale::bi_class(biplot_df, x=Protection, y=Earthworm_richness, dim=4)

# load map
world.inp <- map_data("world")

biplot <- ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), col="grey90", fill="white") +
  xlim(-10, 30) +
  ylim(35, 70) +
  
  geom_tile(data=biplot_input, aes(x=x, y=y, fill=bi_class))+
  bi_scale_fill(pal = "BlueGold", dim=4, flip_axes = T)+
  theme_gray()+
  theme(legend.position="none", axis.text.x=element_text(angle=45, hjust=1),
        panel.background = element_rect(fill="#e8e8e8"))


legend <- bi_legend(pal = "BlueGold",
                    dim = 4,
                    xlab = "Protection",
                    ylab = "Earthworm diversity",
                    size = 12,
                    flip_axes = T)
#legend

finalPlot <- cowplot::ggdraw(biplot) +
  cowplot::draw_plot(legend, 0.03, 0.8, 0.28, 0.15)

png(paste0(data_wd, "/_figures/Protection_vs_richness_", Taxon_name, ".png"), width=1500, height=1300)
finalPlot +
  annotate("text", -Inf, Inf, hjust = -0.5, vjust = 1.5, label = "A", size = 16/.pt,
           fontface = "bold")
dev.off()

#- - - - - - - - - - - - - - - - - - - -
## Extract GBIF taxon keys ####
#- - - - - - - - - - - - - - - - - - - -

dat <- read.csv(file=paste0(data_wd, "/_intermediates/Occurrences_GBIF_Crassiclitellata.csv")) #dat
dat %>% dplyr::select(species, speciesKey) %>% unique() %>% arrange(species)


#- - - - - - - - - - - - - - - - - - - -
## Some more numbers ####
#- - - - - - - - - - - - - - - - - - - -

# relationship between life form (ecological group) and number of records
bar_groups <- 
  ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
         mutate(Ecogroup = ifelse(Ecogroup=="Epi-Endogeic", "Epi-endogeic", Ecogroup)), aes(x=Ecogroup, y=NumCells_2km_biomod, color=SpeciesID))+
  geom_bar(stat="identity")+
  #facet_wrap(vars(Ecogroup))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
ggsave(filename=paste0(data_wd, "/_figures/NumberOccurrences_ecogroup_", Taxon_name, ".png"), bar_groups)

bar_groups2 <- 
  ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
           mutate(Ecogroup = ifelse(Ecogroup=="Epi-Endogeic", "Epi-endogeic", Ecogroup)), aes(x=Ecogroup, y=Records, color=SpeciesID))+
  geom_bar(stat="identity")+
  #facet_wrap(vars(Ecogroup))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
ggsave(filename=paste0(data_wd, "/_figures/NumberRecords_ecogroup_", Taxon_name, ".png"), bar_groups2)

#- - - - - - - - - - - - - - - - - - - -
# records vs. threat level
speciesStatus <- read.csv(paste0(here::here(), "/doc/Species_status_", Taxon_name, ".csv"))

bar_status <- 
  ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
           right_join(speciesStatus, by="SpeciesID") %>%
           mutate(RedList_Germany=ifelse(RedList_Germany=="", NA, RedList_Germany)), aes(x=RedList_Germany, y=NumCells_2km_biomod, color=SpeciesID))+
  geom_bar(stat="identity")+
  #facet_wrap(vars(Ecogroup))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
ggsave(filename=paste0(data_wd, "/_figures/NumberOccurrences_status_", Taxon_name, ".png"), bar_status)

bar_status2 <- 
  ggplot(data=speciesNames %>% filter(SpeciesID %in% species10) %>%
           right_join(speciesStatus, by="SpeciesID") %>%
           mutate(RedList_Germany=ifelse(RedList_Germany=="", NA, RedList_Germany)), aes(x=RedList_Germany, y=Records, color=SpeciesID))+
  geom_bar(stat="identity")+
  #facet_wrap(vars(Ecogroup))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")
ggsave(filename=paste0(data_wd, "/_figures/NumberRecords_status_", Taxon_name, ".png"), bar_status2)


#- - - - - - - - - - - - - - - - - - - -
## DRAFT (old) plots ####
#- - - - - - - - - - - - - - - - - - - -
# 
# # bar chart of percent area covered by PA per species
# f <- ggplot(data=cover_df %>% filter(IUCNcat!="Presence" & IUCNcat!="Unprotected" & IUCNcat!="Outside.PA" & IUCNcat!="Protected"), 
#             aes(y=coverage, x=SSP, fill=IUCNcat))+
#   geom_bar(position="stack", stat="identity")+
#   theme_bw()+
#   xlab("Species")+ ylab("Proportion of range covered by protected area network")+
#   coord_flip()+
#   facet_wrap(vars(SpeciesID))+
#   #scale_fill_viridis_d()+
#   scale_fill_manual(values=c("darkgoldenrod4","darkgoldenrod3", "darkgoldenrod2", "darkgoldenrod1", "goldenrod2","goldenrod1", "gold1",
#                              "lightgoldenrod1","palegoldenrod", "lemonchiffon2","gainsboro" ))+
#   scale_y_continuous(expand=c(0,0), limits=c(0,0.65))+
#   #scale_x_discrete(labels=SpeciesID)+
#   geom_vline(xintercept=1.5, lty=2)+
#   theme(legend.position="none")

# # calculate raw species richness
# #### needs to be fixed ######
# occ_rich <- occ_points %>% 
#   group_by(Latitude = round(x,0), Longitude=round(y,0)) %>%
#   summarise_at(vars(colnames(occ_points %>% dplyr::select(-x, -y))), mean, na.rm=T)
# occ_rich$Richness <- apply(occ_rich > 0, 1, sum, na.rm=T)
# 
# # plot total species' occurrences
# plotOcc <- ggplot()+ 
#   geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
#   xlim(min(extent_Europe[1], na.rm = T), max(extent_Europe[2], na.rm = T)) +
#   ylim(min(extent_Europe[3], na.rm = T), max(extent_Europe[4], na.rm = T)) +
#   
#   geom_point(data=occ_rich %>%
#                dplyr::select(c(Latitude, Longitude, Richness)),
#              aes(x=Latitude, y=Longitude, color=Richness))+ #, alpha=`Number of Species`
#   
#   scale_color_gradient2(5,    # provide any number of colors
#                         low = "black", high="orange", mid= "blue",
#                         midpoint = 10,
#                         #values = scales::rescale(c(1, 2, 3, 5, 10, 30)), 
#                         breaks = c(1, 2, 5, 10, 20, 30), limits=c(0,30))+
#   theme_bw()+
#   theme(legend.position = "bottom", legend.text = element_text(size=8), legend.key.width = unit(2, "cm"))
# 
# plotOcc
# 


