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

library(terra)

library(biomod2) # also to create pseudo-absences

library(ggplot2) # for plotting the curves
library(ggpubr)

library(sf)
library(ggsn) #for scale and north arrow

library(eurostat) # to get NUTS regions
library(giscoR)  # to get NUTS regions

library(parallel)
library(doParallel)

# plotting
library(gridExtra)
library(colorRamps) #color gradient in Supplementary figure

# NMDS
library(vegan)

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


# load NUTS II regions of Germany in R as sf
nuts_sf <- eurostat::get_eurostat_geospatial(nuts_level = "2", year = "2021", crs = "4326", resolution = "01") #crs: "4326" - WGS84, "3035" - ETRS89 / ETRS-LAEA, "3857" - Pseudo-Mercator

# extract Portuguese sf
por_sf <- nuts_sf[nuts_sf$CNTR_CODE=="PT",]
por_sf <- por_sf[por_sf$id == "PT11",] #only continental Portugal
por_sf <- sf::st_crop(por_sf, sf::st_bbox(c(xmin=extent_portugal[1], xmax=extent_portugal[2], ymin=extent_portugal[3], ymax=extent_portugal[4])))

# define names of taxonomic groups
taxaNames <- c("Crassiclitellata", "Fungi", "Nematodes", "Protists", "Eukaryotes", "Bacteria") #
taxaNames

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
           Plot_eval = TRUE,
           Plot_uncertain = TRUE,
           Plot_predict = TRUE)[0,]

for(temp_taxon in c("Crassiclitellata", "Nematodes", "Fungi", "Bacteria", "Eukaryotes", "Protists")){ #"Earthworms", "EarthGenus", "NemaGenus", 
  print(temp_taxon)
  temp_overview <- data.frame(test=1)
  
  temp_overview[,"Taxon"] <- temp_taxon
  
  # number of taxa
  try({
    species100 <- read.csv(file=paste0("_intermediates/SDM_", temp_taxon, ".csv"))
    if(nrow(species100) != 0){
      species10 <- read.csv(file=paste0("_intermediates/ESM_", temp_taxon, ".csv"))
      speciesSub <- species100 %>% 
        pull(species) %>%
        c(species10 %>% pull(species))
    }else{
      species10 <- read.csv(file=paste0("_intermediates/ESM_", temp_taxon, ".csv"))
      speciesSub <- species10 %>% pull(species)
    }
    temp_overview[,"No_taxa"] <- length(speciesSub)
    temp_overview[,"No_taxa10"] <- length(species10 %>% pull(species))
    temp_overview[,"No_taxa100"] <- length(species100 %>% pull(species))
  })
    
  # data prepared?
  temp_overview[,"Data_clean"] <- file.exists(paste0("_intermediates/Occurrences_clean_", temp_taxon, ".csv"))
  temp_overview[,"Data_gridded"] <- file.exists(paste0("_intermediates/Occurrence_rasterized_1km_", temp_taxon, ".csv"))
  temp_overview[,"Data_biomod"] <- length(list.files(paste0("_intermediates/BIOMOD_data/", temp_taxon)))
  
  # how many varImp ready?
  temp_overview[,"No_varImp"] <- length(list.files(paste0("_results/_TopPredictor/", temp_taxon), 
                                                  pattern = "SDM_"))
  
  # how many SDMs and ESMs ready?
  temp_overview[,"No_SDMs"] <- length(list.files(paste0("_results/", temp_taxon, "/SDMs"), 
                                                  pattern = "SDM_"))
  temp_overview[,"No_ESMs"] <- length(list.files(paste0("_results/", temp_taxon, "/SDMs"), 
                                                pattern = "ESM_"))
  
  # how many predictions & uncertainty ready?
  temp_overview[,"No_predict"] <- length(list.files(paste0("_results/", temp_taxon, "/Projection")))
  temp_overview[,"No_uncertain_clamp"] <- length(list.files(paste0("_results/", temp_taxon, "/Uncertainty"),
                                                           pattern = "ClampingMask_"))
  temp_overview[,"No_uncertain_cv"] <- length(list.files(paste0("_results/", temp_taxon, "/Uncertainty"),
                                                        pattern = "CV_"))
  temp_overview[,"No_uncertain"] <- as.numeric(temp_overview["No_uncertain_clamp"]) + as.numeric(temp_overview["No_uncertain_cv"])

  # visualization ready?
  temp_overview[,"Plot_occ_clean"] <- file.exists(paste0("_figures/OccurrencesCleaned_", temp_taxon, ".pdf"))
  temp_overview[,"Plot_occ_grid"] <- file.exists(paste0("_figures/OccurrencesGridded_", temp_taxon, ".pdf"))
  temp_overview[,"Plot_varImp"] <- file.exists(paste0("_figures/VariableImportance_MaxEnt_", temp_taxon, ".pdf"))
  temp_overview[,"Plot_eval"] <- file.exists(paste0("_figures/Model_performance_", temp_taxon, ".pdf"))
  temp_overview[,"Plot_uncertain"] <- file.exists(paste0("_figures/Uncertainty_", temp_taxon, ".pdf"))
  temp_overview[,"Plot_predict"] <- file.exists(paste0("_figures/DistributionMap_bestBinary_", temp_taxon, ".pdf"))
  
  # add to summary df
  df_overview <- full_join(df_overview, temp_overview %>% dplyr::select(-test))

}

df_overview


# number of taxa passing model evaluation threshold
thresholds_auc <- c(0.7, 0.6, 0.5) #values <0.5 indicate performance worse than random (Elith et al. 2006), ESMs only include AUC>0.5 in ensemble (Breiner et al. 2015)

df_overview[paste0(rep(c("No_taxa10_AUC", "No_taxa100_AUC"), each=length(thresholds_auc)), rep(thresholds_auc, 2))] <- NA


for(temp_taxon in c("Crassiclitellata", "Nematodes", "Fungi", "Bacteria", "Eukaryotes", "Protists")){ #"Earthworms", "EarthGenus", "NemaGenus", 
  print(temp_taxon)
  
  # number of taxa
  try({
    species100 <- read.csv(file=paste0("_intermediates/SDM_", temp_taxon, ".csv"))
    if(nrow(species100) != 0){
      species10 <- read.csv(file=paste0("_intermediates/ESM_", temp_taxon, ".csv"))
      speciesSub <- species100 %>% 
        pull(species) %>%
        c(species10 %>% pull(species))
    }else{
      species10 <- read.csv(file=paste0("_intermediates/ESM_", temp_taxon, ".csv"))
      speciesSub <- species10 %>% pull(species)
    }
    })
  
  if(file.exists(paste0("_results/Model_evaluation_", Taxon_name, ".csv"))){
    mod_eval <- read_csv(paste0("_results/Model_evaluation_", temp_taxon, ".csv"))
    for(thresh_auc in thresholds_auc){
      # taxa that pass evaluation threshold of AUC>thresh_auc
      df_overview[df_overview$Taxon==temp_taxon,paste0("No_taxa10_AUC", thresh_auc)] <- mod_eval %>% filter(Species %in% species10$species) %>% filter(AUC>=thresh_auc) %>% count() %>% pull(n)
      df_overview[df_overview$Taxon==temp_taxon,paste0("No_taxa100_AUC", thresh_auc)] <- mod_eval %>% filter(Species %in% species100$species) %>% filter(AUC>=thresh_auc) %>% count() %>% pull(n)
        
    }
    rm(mod_eval)
  }
}
df_overview

#- - - - - - - - - - - - - - - - - - - - -
Taxon_name <- "Protists"

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

if(Taxon_name == "Fungi") species10 <- species10[species10 != "F02025"]

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
### BIOMOD Occurrences per species (all taxa) ####
#- - - - - - - - - - - - - - - - - - - - -

# define names of taxonomic groups
taxaNames <- c("Crassiclitellata", "Fungi", "Nematodes", "Protists", "Eukaryotes", "Bacteria") #
taxaNames

# Calculate individuals species' occurrence based on BIOMOD input
occ_list <- lapply(taxaNames, function(tax){
  read_csv(file=paste0("_results/Occurrence_rasterized_1km_BIOMOD_", tax, ".csv")) %>%
    mutate(taxon = tax)})
occ_list <- do.call(rbind, occ_list)

taxa_list <- lapply(taxaNames, function(tax){
  read_csv(file=paste0("_intermediates/SDM_", tax, ".csv")) %>%
    mutate(taxon = tax, records = 100) %>%
    full_join(read_csv(file=paste0("_intermediates/ESM_", tax, ".csv")) %>%
                mutate(taxon = tax, records = 10))})
taxa_list <- do.call(rbind, taxa_list)
taxa_list %>% count(records, taxon)

# plot all
p_occ_bm <- ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white", color = "grey80") +
  geom_sf(data = por_sf)+
  coord_sf(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]), expand = FALSE)+
  
  geom_sf(data=st_as_sf(occ_list,
                           coords = c("x", "y"),   # columns with coordinates
                           crs = 4326) %>% #WGS84
            filter(occ >= 1) %>% dplyr::select(taxon, geometry, occ) %>% unique(),
             aes(group=taxon),
             cex=0.5, pch=19, alpha = 0.7, color = "forestgreen")+
  #scale_color_manual(values = c("1"="forestgreen", "0"="orange"))+
  facet_wrap(vars(taxon))+
  
  # add number of grid cells in which the species is present
  geom_text(data=occ_list %>% group_by(taxon) %>% summarize("n"=sum(occ)),
            aes(x=-6.5, y=41, label=paste0("n=", n)), color="black",
            inherit.aes=FALSE, parse=FALSE, cex=2, hjust=0)+
  geom_text(data=occ_list %>% group_by(taxon) %>% 
              filter(occ >= 1) %>%
              dplyr::select(SpeciesID, taxon) %>% unique() %>% 
              mutate(occ = 1) %>% summarize("n"=sum(occ)),
            aes(x=-6.5, y=40.9, label=paste0("taxa=", n)), color="black",
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
p_occ_bm
ggsave(paste0("_figures/OccurrencesBIOMOD_allGroups_perGroup.pdf"),
       p_occ_bm,
       height = 12, width = 12)

p_occ_bm_sum <- ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white", color = "grey80") +
  geom_sf(data = por_sf)+
  coord_sf(xlim = c(extent_portugal[1], extent_portugal[2]),
           ylim = c(extent_portugal[3], extent_portugal[4]), expand = FALSE)+
  
  geom_sf(data=sf::st_as_sf(occ_list,
                            coords = c("x", "y"),   # columns with coordinates
                            crs = 4326) %>% filter(occ >=1) %>% 
            dplyr::select(geometry, taxon, occ) %>% unique() %>%
               group_by(geometry) %>% 
               summarize(across(occ, sum)),
             aes(color = occ), 
             pch=19,  cex = 3)+
  #scale_color_manual(values = c("1"="forestgreen", "0"="orange"))+
  scale_color_steps(name = "Number of taxonomic groups",
                       low = "#d9b3c6",
                       high = "#800040",
                    limits= c(0, 5))+
  # add number of grid cells in which the species is present
  theme_bw()+
  theme(axis.title = element_blank(), #legend.title = element_blank(),
        legend.title.position = "top",
        legend.position = "inside",legend.direction = "horizontal",
        legend.position.inside = c(0.84, 0.17),
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=10), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()#,
        #panel.background = element_rect(fill="grey80")
        )
p_occ_bm_sum
ggsave(paste0("_figures/OccurrencesBIOMOD_allGroups_sum.pdf"),
       p_occ_bm_sum,
       height = 12, width = 12)


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
  theme(legend.position = "bottom", 
        axis.text.x = element_text(size=15), axis.text.y = element_text(size=5),
        axis.title = element_text(size=15),
        legend.text = element_text(size=15), legend.title = element_blank())

ggsave(paste0("_figures/VariableImportance_MaxEnt_species_", Taxon_name, ".pdf"), height=11, width=9)
 

#- - - - - - - - - - - - - - - - - - - - -
## Variable importance (MaxEnt, all taxa) ####
#- - - - - - - - - - - - - - - - - - - - -

## Visualize and get top 10
var_imp <- lapply(taxaNames, function(x){
  y <- read_csv(paste0("_results/Variable_importance_MaxEnt_", x, ".csv"))
  
  # save top 10 predictors
  z <- y %>% group_by(Predictor) %>% 
    summarize(across(maxent, median)) %>%
    rename(median=maxent) %>%
    arrange(desc(median)) %>% 
    top_n(10) %>% # top 10 per taxon group
    mutate(Taxon = x) %>%
    left_join(y %>% group_by(Predictor) %>% filter(maxent>0) %>% count(name="No_species>0"),
              by="Predictor")
  z
})
var_imp <- do.call(rbind, var_imp)
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

write_csv(var_imp, paste0("_figures/VariableImportance_MaxEnt_top10_allTaxa.csv"))
var_imp <- read_csv(paste0("_figures/VariableImportance_MaxEnt_top10_allTaxa.csv"))

# summarize mean & SD
# View(var_imp %>% 
#        dplyr::select(-Species, -Category) %>% group_by(Predictor) %>% summarize_all(list(mean=mean, sd=sd)))

# plot Permutation importance
ggplot()+
  geom_bar(data=var_imp, 
           aes(x=median, y=reorder(Predictor, median), fill = Category),
           stat = "identity",
           position = "dodge")+
  # stat_summary(data=var_imp, aes(x=maxent, y=reorder(Predictor, maxent,median)),
  #              geom = "errorbar", fun.min = mean, fun = mean, fun.max = mean, width = .75, linetype="dashed")+
  xlab("Variable importance (Permutation importance)")+
  ylab("Predictor")+
  scale_x_continuous(expand = c(0,0))+
  facet_wrap(vars(Taxon))+
  theme_bw()+
  theme(axis.text = element_text(size=12), 
        axis.title = element_text(size=20), axis.title.y=element_blank(),
        strip.text = element_text(size = 20),
        legend.text = element_text(size=12), legend.title = element_blank(), 
        panel.grid.major.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.555,0.635))
ggsave(paste0("_figures/VariableImportance_MaxEnt_top10_allTaxa.png"), 
       last_plot(),
       height = 7, width = 9)


#- - - - - - - - - - - - - - - - - - - - -
## Variable importance (NMDS) ####
#- - - - - - - - - - - - - - - - - - - - -

var_imp <- lapply(taxaNames, function(x){
  y <- read_csv(paste0("_results/Variable_importance_MaxEnt_", x, ".csv"))
  
  # save top 10 predictors
  z <- y %>%
    arrange(desc(maxent)) %>% 
    mutate(Taxon = x)
  z
})
var_imp <- do.call(rbind, var_imp)
var_imp

var_imp <- var_imp %>% 
  pivot_wider(names_from=Predictor, values_from = maxent)

var_imp <- data.frame(row.names = var_imp$Species,
                      var_imp %>% dplyr::select(-Species))

# run NMDS
set.seed(293)
var_nmds <- vegan::metaMDS(var_imp %>% dplyr::select(-Taxon), distance = "bray")
var_nmds

vegan::stressplot(var_nmds)

# # Base plot of taxa
# ordiplot(var_nmds, type = "n")           # empty plot
# text(var_nmds, display = "sites",        # taxa positions
#      labels = rownames(var_imp),
#      cex = 0.8, col = "steelblue")
# 
# # Optional: add vectors for environmental variables
# fit <- envfit(var_nmds, var_imp, permutations = 999)
# plot(fit, col = "black")

nmds_scores <- as.data.frame(scores(var_nmds, display = "sites"))
nmds_scores <- nmds_scores %>%
  rownames_to_column("Species") %>%
  left_join(var_imp %>% 
              mutate(Species = rownames(var_imp)) %>%
              dplyr::select(Taxon, Species), by = "Species")
# Group median (alternative)
group_medians <- nmds_scores %>%
  group_by(Taxon) %>%
  summarise(NMDS1 = median(NMDS1),
            NMDS2 = median(NMDS2))

# vectors of predictors
fit <- envfit(var_nmds, var_imp %>% dplyr::select(-Taxon), permutations = 999)
env_vectors <- as.data.frame(scores(fit, display = "vectors"))
env_vectors$Variable <- rownames(env_vectors)

ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Taxon)) +
  xlim(c(-1.6, 1))+ ylim(c(-1.5, 1.3))+
  stat_ellipse(data = nmds_scores, aes(group = Taxon, color = Taxon), #%>% filter(Taxon != "Crassiclitellata")
               type = "t", level = 0.95, linetype = 2, linewidth = 1) +
  geom_hline(aes(yintercept = 0), color = "grey80")+
  geom_vline(aes(xintercept = 0), color = "grey80")+
  geom_text(aes(label = Species), cex = 1) +
  geom_segment(data = env_vectors,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "black",
               inherit.aes = FALSE) +
  geom_text(data = env_vectors,
            aes(x = NMDS1, y = NMDS2, label = Variable),
            color = "black", vjust = -0.5, inherit.aes = FALSE)+
  geom_point(data = group_medians, aes(x = NMDS1, y = NMDS2, fill = Taxon, shape = Taxon), size = 10, color = "black", alpha = 0.5) +
  scale_shape_manual(values=c(21:23, 21:25))+
  theme_minimal() +
  labs(title = "NMDS of Taxa based on Environmental Variable Importance",
       x = "NMDS1",
       y = "NMDS2")+
  theme(panel.grid = element_blank())
ggsave("_figures/NMDS_varImp.png")


#- - - - - - - - - - - - - - - - - - - - - -
## Load UNCERTAINTY (relevant for all plots) ####
#- - - - - - - - - - - - - - - - - - - - - -

# load uncertainty extent for all maps
extent_10 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_10.RData"))) #extent_df
extent_100 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_100.RData"))) #extent_df

# load uncertainty
uncertain_tif <- terra::rast(paste0("_results/SDM_Uncertainty_", Taxon_name, "_10.tif")) 
uncertain_10 <- terra::as.data.frame(uncertain_tif, xy=TRUE)

uncertain_tif <- terra::rast(paste0("_results/SDM_Uncertainty_", Taxon_name, "_100.tif")) 
uncertain_100 <- terra::as.data.frame(uncertain_tif, xy=TRUE)

rm(uncertain_tif)

#- - - - - - - - - - - - - - - - - - - - - -
## Uncertainty ####
#- - - - - - - - - - - - - - - - - - - - - -
# ESMs
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  geom_tile(data=uncertain_10 %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
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
ggsave(file=paste0("_figures/Uncertainty_ESM_", Taxon_name, ".pdf"), height = 5, width = 8)

if(Taxon_name != "Crassiclitellata"){
  # SDMs
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    geom_tile(data=uncertain_100 %>% filter(Mean!=0), aes(x=x, y=y, fill=Mean))+
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
  ggsave(file=paste0("_figures/Uncertainty_SDM_", Taxon_name, ".pdf"), height = 5, width = 8)
}

# uncertainty threshold
# if(Taxon_name == "Earthworms") temp_thresh <- 452.75
# if(Taxon_name == "Nematodes") temp_thresh10 <- 485.5714; temp_thresh100 <- 57.45907
temp_thresh10 <- stats::quantile(uncertain_10$Mean, 0.9, na.rm=TRUE)
temp_thresh100 <- stats::quantile(uncertain_100$Mean, 0.9, na.rm=TRUE)

# ESMs
ggplot()+
  geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
  coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                  ylim = c(extent_portugal[3], extent_portugal[4]))+
  
  geom_tile(data=uncertain_10 %>% filter(Mean<temp_thresh10), aes(x=x, y=y, fill=Mean))+
  ggtitle("Coefficient of variation averaged across SDMs")+
  scale_fill_viridis_c(option="E")+
  geom_tile(data=uncertain_10 %>% filter(Mean>=temp_thresh10), aes(x=x, y=y), fill="linen")+
  theme_bw()+
  theme(axis.title = element_blank(), legend.title = element_blank(),
        legend.position = "right",legend.direction = "vertical",
        axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
ggsave(file=paste0("_figures/Uncertainty_ESM_", Taxon_name, "_", temp_thresh10, ".pdf"), height = 5, width = 8)

if(Taxon_name != "Crassiclitellata"){
  # SDMs
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                    ylim = c(extent_portugal[3], extent_portugal[4]))+
    
    geom_tile(data=uncertain_100 %>% filter(Mean<temp_thresh100), aes(x=x, y=y, fill=Mean))+
    ggtitle("Coefficient of variation averaged across SDMs")+
    scale_fill_viridis_c(option="E")+
    geom_tile(data=uncertain_100 %>% filter(Mean>=temp_thresh100), aes(x=x, y=y), fill="linen")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "right",legend.direction = "vertical",
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  ggsave(file=paste0("_figures/Uncertainty_SDM_", Taxon_name, "_", temp_thresh100, ".pdf"), height = 5, width = 8)
}

#- - - - - - - - - - - - - - - - - - - - - -
## Map species uncertainty maps ####
# to plot uncertainty of species with temp_thresh, add commented text and
# change limits of scale_fill to c(0, temp_thresh) in both plots and legend
if(!(Taxon_name %in% c("Fungi", "Bacteria", "Protists", "Eukaryotes"))){
  uncertain_df <- uncertain_10; temp_thresh <- temp_thresh10
  
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
  
  pdf(file=paste0("_figures/Uncertainty_ESM_allSpecies_", Taxon_name, ".pdf"), height = 10, width = 10)
  do.call(grid.arrange, plots2)
  dev.off()
}

if(!(Taxon_name %in% c("Fungi", "Bacteria", "Protists", "Eukaryotes"))){
  # 100
  uncertain_df <- uncertain_100; temp_thresh <- temp_thresh100
  
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
  
  pdf(file=paste0("_figures/Uncertainty_SDM_allSpecies_", Taxon_name, ".pdf"), height = 10, width = 10)
  do.call(grid.arrange, plots2)
  dev.off()
}

#- - - - - - - - - - - - - - - - - - - - - -
## Species richness (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes", "Bacteria")){
  
  print(Taxon_name)
  # load uncertainty extent for all maps
  try(extent_10 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_10.RData")))) #extent_df
  try(extent_100 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_100.RData")))) #extent_df
  
  if(exists("extent_10") & exists("extent_100")){
    extent_df <- extent_10 %>% inner_join(extent_100)
  } else{
    if(exists("extent_10")) extent_df <- extent_10
    if(exists("extent_100")) extent_df <- extent_100
  }
  species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))
  species_stack <- terra::as.data.frame(species_stack, xy = TRUE)
  
  species_stack <- species_stack %>% mutate(x = round(x, 5), y = round(y,5))
  extent_df <- extent_df %>% mutate(x = round(x, 5), y = round(y,5))
  
  summary(species_stack$Richness)
  sd(species_stack$Richness)
  
  # species richness
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                    ylim = c(extent_portugal[3], extent_portugal[4]))+
  
    geom_tile(data=extent_df %>% left_join(species_stack %>% filter(Richness>0), by=c("x","y")),
              aes(x=x, y=y, fill=Richness))+
    ggtitle(paste0("Taxon richness (number of taxa) - ", Taxon_name))+
    labs(subtitle = paste0("Total number of modelled taxa: ", length(names(species_stack))-3))+
    scale_fill_viridis_c()+
    #geom_tile(data=extent_df %>% left_join(species_stack %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "right",legend.direction = "vertical",
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=30), legend.key.size = unit(1, 'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  ggsave(file=paste0("_figures/SpeciesRichness_all_", Taxon_name, ".pdf"),
         last_plot(),
         height = 5, width = 8)
}

# stack all groups
stack_all <- terra::rast()
for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes", "Bacteria")){try({
  species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))
  
  temp_richness <- species_stack$Richness
  names(temp_richness) <- Taxon_name
  stack_all <- c(stack_all, temp_richness)
})}
stack_all

stack_all <- terra::mask(stack_all, env_por[[1]])
terra::writeRaster(stack_all, "_results/_Maps/SDM_stack_richness_allTaxa.tif", overwrite=TRUE)

# individually for SDM and ESM
for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes", "Bacteria")){try({
  species100 <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))
  
  ## SDMs
  extent_100 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_100.RData"))) #extent_df
  
  species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))
  species_stack <- terra::as.data.frame(species_stack, xy = TRUE)
  
  species_stack <- species_stack %>% mutate(x = round(x, 5), y = round(y,5))
  extent_df <- extent_df %>% mutate(x = round(x, 5), y = round(y,5))
  
  species_stack <- species_stack %>% dplyr::select(x, y, any_of(pull(species100)), Richness)
  
  # species richness
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                    ylim = c(extent_portugal[3], extent_portugal[4]))+
    
    geom_tile(data=extent_df %>% left_join(species_stack %>% filter(Richness>0), by=c("x","y")),
              aes(x=x, y=y, fill=Richness))+
    ggtitle(paste0("Taxon richness (number of taxa) - ", Taxon_name))+
    labs(subtitle = paste0("Total number of modelled taxa: ", length(names(species_stack))-3))+
    scale_fill_viridis_c()+
    #geom_tile(data=extent_df %>% left_join(species_stack %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "right",legend.direction = "vertical",
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=30), legend.key.size = unit(1, 'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  ggsave(file=paste0("_figures/SpeciesRichness_SDM_", Taxon_name, ".pdf"),
         last_plot(),
         height = 5, width = 8)
  
})}

for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes", "Bacteria")){try({
  species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))
  
  ## ESMs
  extent_10 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_10.RData"))) #extent_df
  
  species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))
  species_stack <- terra::as.data.frame(species_stack, xy = TRUE)
  
  species_stack <- species_stack %>% dplyr::select(x, y, any_of(pull(species10)), Richness)
  
  species_stack <- species_stack %>% mutate(x = round(x, 5), y = round(y,5))
  extent_df <- extent_df %>% mutate(x = round(x, 5), y = round(y,5))
  
  summary(species_stack$Richness)
  sd(species_stack$Richness)
  
  # species richness
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80") +
    coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                    ylim = c(extent_portugal[3], extent_portugal[4]))+
    
    geom_tile(data=extent_df %>% left_join(species_stack %>% filter(Richness>0), by=c("x","y")),
              aes(x=x, y=y, fill=Richness))+
    ggtitle(paste0("Taxon richness (number of taxa) - ", Taxon_name))+
    labs(subtitle = paste0("Total number of modelled taxa: ", length(names(species_stack))-3))+
    scale_fill_viridis_c()+
    #geom_tile(data=extent_df %>% left_join(species_stack %>% filter(Richness==0), by=c("x","y")), aes(x=x, y=y), fill="grey60")+
    theme_bw()+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "right",legend.direction = "vertical",
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          legend.text = element_text(size=30), legend.key.size = unit(1, 'cm'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  ggsave(file=paste0("_figures/SpeciesRichness_ESM_", Taxon_name, ".pdf"),
         last_plot(),
         height = 5, width = 8)
})}

#- - - - - - - - - - - - - - - - - - - - - -
## Species distributions (current) ####
#- - - - - - - - - - - - - - - - - - - - - -

species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))
#species_stack <- terra::as.data.frame(species_stack, xy = TRUE)

occ_points <- read_csv(file=paste0("_intermediates/Occurrence_rasterized_1km_", Taxon_name, ".csv"))
occ_points <- occ_points %>% pivot_longer(cols = all_of(speciesSub), values_to = "occ")

extent_tif <- terra::rasterize(terra::vect(extent_df, geom = c("x", "y")), species_stack)

# Optionally, mask the cropped_layer to only keep values where layer1 == 1
species_stack <- mask(species_stack, extent_tif)
species_df <- terra::as.data.frame(species_stack, xy=TRUE)

# map binary species distributions
plots <- lapply(3:(ncol(species_df)-1), function(s) {try({
  print(s-2)
  temp_occ <- occ_points %>% 
    filter(name == colnames(species_df)[s]) %>% 
    dplyr::select(!(starts_with(c("PC", "Species"))))
  ggplot()+
    geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "white", color = "grey80") +
    coord_cartesian(xlim = c(extent_portugal[1], extent_portugal[2]),
                    ylim = c(extent_portugal[3], extent_portugal[4]))+
    
    geom_tile(data=species_df, 
              aes(x=x, y=y, fill=as.factor(species_df[,s])))+
    geom_point(data=temp_occ, 
               aes(x=x, y=y, color=as.factor(occ)), cex = 5, shape = 19)+
    geom_text(aes(x = -7.65, y = 42.2, label = "Spain"), size = 10, color = "grey60")+
    geom_text(aes(x = -8.05, y = 40.55, label = "Portugal"), size = 10, color = "grey60")+
    ggtitle(colnames(species_df)[s])+
    scale_fill_manual(values=c("1"="forestgreen", "0"="grey90","NA"="white"))+
    scale_color_manual(values=c("1"="black", "0"="grey60"))+
    theme_bw()+
    ggtitle(colnames(species_df)[s])+
    guides(fill = guide_legend(# title.hjust = 1, # adjust title if needed
      label.position = "bottom",
      label.hjust = 0.5))+
    theme(axis.title = element_blank(), legend.title = element_blank(),
          legend.position = "none", legend.direction = "horizontal",
          legend.key.size = unit(0.2, 'cm'),
          legend.text = element_text(size=10),
          title = element_text(size=30),
          axis.line = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
})
})

#pdf(file=paste0("_figures/DistributionMap_bestBinary_", Taxon_name, ".pdf"), )
png(file=paste0("_figures/DistributionMap_bestBinary_", Taxon_name, ".png"),width=3000, height=3300)
do.call(grid.arrange, plots)
dev.off()

#- - - - - - - - - - - - - - - - - - - - -
## Model performance ####
#- - - - - - - - - - - - - - - - - - - - -

data_eval <- lapply(taxaNames, function(x){
  y <- read_csv(paste0("_results/Model_evaluation_", x, ".csv"))
  
  # save top 10 predictors
  z <- y %>%
    mutate(Taxon = x,
           Model = ifelse(is.na(KAPPA), 10, 100))
  z
})
data_eval <- do.call(rbind, data_eval)
data_eval

data_eval %>% dplyr::select(-Species) %>% group_by(Model, Taxon) %>% summarize_all(median, na.rm=TRUE) #note: ESMs no kappa
data_eval %>% dplyr::select(-Species) %>% group_by(Model, Taxon) %>% summarize_all(sd, na.rm=TRUE) #note: ESMs no kappa

# point plot with lables, tss over roc
ggplot(data_eval, 
       aes(x=MaxTSS, y=AUC, color=Species))+
  #geom_point()+
  geom_text(aes(label = Species), nudge_x = 0, nudge_y = 0, check_overlap = F, cex=2)+
  facet_wrap(vars(Model))+
  xlim(0,1)+
  theme_bw()+
  theme(legend.position = "none")
ggsave(paste0("_figures/Model_performance_allTaxa.pdf"))

# boxplot, tss per algorithm
ggplot()+
  geom_boxplot(data = data_eval, aes(x=AUC, y = Taxon, fill = "AUC"), alpha = 0.2, outlier.alpha = 1, 
               outlier.shape = 18, outliers = TRUE, staplewidth = 0.2, lty = "longdash")+
  geom_jitter(data = data_eval, aes(x=AUC, y = Taxon, color = "AUC"), height = 0.2, alpha = 0.5)+
  geom_boxplot(data = data_eval, aes(x=MaxTSS, y = Taxon, fill = "TSS"), alpha = 0.2, outlier.alpha = 1, outlier.shape = 17, outliers = TRUE, staplewidth = 0.2)+
  geom_jitter(data = data_eval, aes(x=MaxTSS, y = Taxon, color = "TSS"), height = 0.2, alpha = 0.5)+
  facet_grid(vars(Model))+
  xlab("")+ylab("")+
  scale_x_continuous(expand=c(0,0), limits = c(0,1))+
  #coord_flip()+
  theme_bw()+
  theme(panel.grid.major.y = element_blank(),
        text = element_text(size = 15))
ggsave(paste0("_figures/Model_performance_allTaxa_boxplot.png"),
       height = 5, width = 8)

data_eval %>% group_by(Model) %>% count()
data_eval %>% filter(AUC>=0.7) %>% 
  group_by(Model, Taxon) %>% count()
data_eval %>% filter(AUC>=0.7) %>% 
  group_by(Model) %>% count()
data_eval %>% filter(AUC<0.7) %>% 
  group_by(Model) %>% arrange(AUC)

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


