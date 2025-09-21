#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#      Check output from Zonation           #
#          author: Romy Zeiss               #
#            date: 2025-09-17               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

zonation_dir <- "../SoilBioPrio_Zonation/POR"
list.files(zonation_dir, include.dirs = TRUE)

library(terra)
library(tidyverse)

approaches <- c("target", "complement", "prevent")
species_numbers <- c("", "_10", "_100")
group_weights <- c("", "_grW")
protected_areas <- c("", "_allPA")
protection_distances <- c("", "_dist")
degradation_weights <- c("", "_max", "_maxSubset") #"_sum", "_sumx10", "_maxdiv3"

# # move content from "subregion_1" (area of analysis) to scenario subfolder
# list_dirs <- list.dirs(zonation_dir, recursive = TRUE, full.names = TRUE)
# list_dirs <- list_dirs[str_detect(list_dirs, "subregion_1") & !(str_detect(list_dirs, "gui"))]
# 
# for(temp_list in list_dirs){
#   print(paste0("Move folder ", temp_list))
#   file.copy(temp_list, paste0(zonation_dir, "/", str_extract(temp_list, "(?<=/POR/)[^/]+")), recursive=TRUE)
#   
#   print(paste0("Done: ", str_extract(temp_list, "(?<=/POR/)[^/]+")))
# }

# load priority maps
priority_rast <- terra::rast()
for(approach in approaches){
  for(protection_distance in protection_distances){
  for(degradation_weight in degradation_weights){ 
  for(protected_area in protected_areas){
  for(species_number in species_numbers){
  for(group_weight in group_weights) try({
    # print(paste0(zonation_dir, "/", approach, 
    #              species_number, group_weight,
    #              protected_area, protection_distance, degradation_weight,
    #              "/subregion_1"))
    
    if(dir.exists(paste0(paste0(zonation_dir, "/", approach, 
                                species_number, group_weight,
                                protected_area, protection_distance, degradation_weight,
                                "/subregion_1")))){
      temp_rast <- terra::rast(paste0(zonation_dir, "/", approach, 
                                      species_number, group_weight,
                                      protected_area, protection_distance, degradation_weight, "/merged_rankmap.tif"))
    } else {
      temp_rast <- terra::rast(paste0(zonation_dir, "/", approach, 
                                      species_number, group_weight,
                                      protected_area, protection_distance, degradation_weight, "/rankmap.tif"))
    }
    names(temp_rast) <- paste0(approach,  
                               species_number, group_weight,
                               protected_area, protection_distance,
                               degradation_weight)
    
    priority_rast <- c(priority_rast, temp_rast)
    
    print(paste0("Priority map of ", paste0(approach,
                                            species_number, group_weight,
                                            protected_area, protection_distance, degradation_weight), " added."))
  }, silent = TRUE)}
}}}}
priority_rast

# mask to extent to get rid of strange orange background
env_por <- terra::rast( "_intermediates/EnvPredictor_1km_POR_normalized.tif")
#priority_rast <- terra::crop(priority_rast, env_por[[1]])
priority_rast <- terra::mask(priority_rast, env_por[[1]])

terra::plot(priority_rast, nly=30)
terra::plot(terra::diff(priority_rast))

#- - - - - - - - - - - - - - - - - - - - - -#
## Classify top areas ####

priority_rast_c <- priority_rast

# top 30, 10, 5, and 1% of area (if no protected areas)
m <- matrix(c(
  0,    0.7,  0,
  0.7,  0.9,  1,
  0.9,  0.95, 2,
  0.95, 0.99, 3,
  0.99, 1.0,  4
), ncol=3, byrow=TRUE)

# maps without protected areas
for (i in names(priority_rast_c)[names(priority_rast_c) %in% apply(expand.grid("target", species_numbers, group_weights), 1, paste, collapse = "")]) {
  
  priority_rast_c[[i]] <- classify(priority_rast_c[[i]], m)
  
  temp_name <- names(priority_rast_c[[i]])
  levels(priority_rast_c[[i]]) <- data.frame(
    ID    = 0:4,
    class = c(">30%", "30%","10%","5%","1%")
  )
  names(priority_rast_c[[i]]) <- temp_name
}

# top 30, 10, 5, and 1% of area WITH protected areas
protect_stack <- terra::rast(paste0("_intermediates/Protection_POR_coverage.tif"))

coverage_all <- 38569
coverage_pa <- unlist(global(protect_stack$PA==1, fun = "sum", na.rm = TRUE))
coverage_iucnpa <- unlist(global(protect_stack$IUCN_PA==1, fun = "sum", na.rm = TRUE))
print(paste0("Number of cells in total: ", coverage_all)); print(paste0("Number of cells under any protection: ", coverage_pa)); print(paste0("Number of cells under IUCN protection: ", coverage_iucnpa))

percentage_pa <- coverage_pa / coverage_all #Zonation: 0.75164
percentage_iucnpa <- coverage_iucnpa / coverage_all #Zonation: 0.883715
print(paste0("Percentage of area coverage by any PA: ", round(percentage_pa*100,2), "% and by IUCN PAs: ", round(percentage_iucnpa*100,2), "%."))

m_pa <- matrix(c(
  0, 1-percentage_pa-0.3,  0,
  1-percentage_pa-0.3,  1-percentage_pa-0.1,  1,
  1-percentage_pa-0.1,  1-percentage_pa-0.05, 2,
  1-percentage_pa-0.05, 1-percentage_pa-0.01, 3,
  1-percentage_pa-0.01, 1-percentage_pa+0.00001, 4,
  1-percentage_pa-0.00001, 1.0, 5
), ncol=3, byrow=TRUE)
m_pa; percentage_pa

m_iucnpa <- matrix(c(
  0, 1-percentage_iucnpa-0.3,  0,
  1-percentage_iucnpa-0.3,  1-percentage_iucnpa-0.1,  1,
  1-percentage_iucnpa-0.1,  1-percentage_iucnpa-0.05, 2,
  1-percentage_iucnpa-0.05, 1-percentage_iucnpa-0.01, 3,
  1-percentage_iucnpa-0.01, 1-percentage_iucnpa+0.00001, 4,
  1-percentage_iucnpa-0.00001, 1.0, 5
), ncol=3, byrow=TRUE)
m_iucnpa; percentage_iucnpa

# maps without protected areas
layer_names <- apply(expand.grid("complement", species_numbers, group_weights, protected_areas, protection_distances), 1, paste, collapse = "")
layer_names <- c(layer_names, apply(expand.grid("prevent", species_numbers, group_weights, protected_areas, protection_distances, degradation_weights), 1, paste, collapse = ""))
for (i in names(priority_rast_c)[names(priority_rast_c) %in% layer_names]) { 
  
  if(str_detect(i, "allPA")){
    priority_rast_c[[i]] <- classify(priority_rast_c[[i]], m_pa)
  } else {
    priority_rast_c[[i]] <- classify(priority_rast_c[[i]], m_iucnpa)
  }
  
  temp_name <- names(priority_rast_c[[i]])
  levels(priority_rast_c[[i]]) <- data.frame(
    ID    = 0:5,
    class = c(">30%", "30%","10%","5%","1%", "PA")
  )
  names(priority_rast_c[[i]]) <- temp_name
}

terra::plot(priority_rast_c)
terra::plot(terra::diff(priority_rast>=0.95)) #note: only diff to layer before

# save priority maps
terra::writeRaster(priority_rast, "_results/_Maps/Zonation_priorities_raw.tif", overwrite=TRUE)
terra::writeRaster(priority_rast_c, "_results/_Maps/Zonation_priorities.tif", overwrite=TRUE)

# plot classified priorities (top 30, 5 etc. %)
pdf("_figures/Zonation_priorities.pdf")
terra::plot(priority_rast_c, maxnl=50)
dev.off()

# # plot difference between top 5% areas
# m_diff <- matrix(c(
#   0,    0.7,  0,
#   0.7,  0.9,  1,
#   0.9,  0.95, 2,
#   0.95, 0.99, 3,
#   0.99, 1.0,  4
# ), ncol=3, byrow=TRUE)
# 
# priority_rast_diff <- terra::diff(priority_rast>0.95)
# 
# pdf("_figures/Zonation_priorities_diff_top5.pdf")
# terra::plot(priority_rast_diff)
# dev.off()

#- - - - - - - - - - - - - - - - - - - - - -
## Performance curves ####
approach <- "target"
group_curve <- read_delim(paste0(zonation_dir,  "/", approach, "/subregion_1/group_curves.csv"))
summary_curve <- read_delim(paste0(zonation_dir,  "/", approach, "/subregion_1/summary_curves.csv"))

names_list <- read_csv(paste0("_results/Zonation/TaxaNames_legend.csv"))

group_curve <- group_curve %>% 
  dplyr::select(-`...20`) %>%
  pivot_longer(
    cols = matches("\\.(\\d+)$"),      # all columns ending with .number
    names_to = c("metric", "group"),
    names_pattern = "(.*)\\.(\\d+)",   # split into metric and group number
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from  = metric,              # mean, weighted_mean, min
    values_from = value
  ) %>%
  mutate(group = as.integer(group)) 

group_curve <- group_curve %>% 
  full_join(names_list %>% dplyr::select(Taxon, wgrp) %>% unique(),
            by = c("group" = "wgrp"))

# feature_curve <- read_delim(paste0(zonation_dir, "/", approach,  "/subregion_1/feature_curves.csv"))
# density_curve <- read_delim(paste0(zonation_dir, "/", approach,  "/subregion_1/feature_density_curves.csv"))
# 
# density_curve$weighted_feature_density_pos <- 
#   as.numeric(gsub(density_curve$weighted_feature_density_pos, pattern = "\r", replacement = ""))
# 
# feature_curve <- feature_curve %>%
#   pivot_longer(
#     cols = -rank,
#     names_to = "feature",
#     values_to = "coverage"
#   )

ggplot(summary_curve)+
  geom_line(aes(x = rank, y = weighted_mean_positive))+
  geom_line(aes(x = rank, y = min, color = "red", alpha = 0.5))

ggplot(group_curve)+
  geom_line(aes(x = rank, y = weighted_mean, lty = Taxon))

# stats
top_legend <- data.frame(
  ID    = 0:5,
  class = c(">30%", "30%","10%","5%","1%", "PA"),
  m_pa_low = m_pa[,1],
  m_pa_upp = m_pa[,2],
  m_iucnpa_low = m_iucnpa[,1],
  m_iucnpa_upp = m_iucnpa[,2],
  m_low = c(m[,1], m_iucnpa[6,1]),
  m_upp = c(m[,2], m_iucnpa[6,2])
)
top_legend

coverage_df <- list()
for(approach in approaches){
  for(protection_distance in protection_distances){
    for(degradation_weight in degradation_weights){ 
      for(protected_area in protected_areas){
        for(species_number in species_numbers){
          for(group_weight in group_weights) try({
            
            temp_curve <- read_delim(paste0(zonation_dir, "/", approach, 
                                            species_number, group_weight,
                                            protected_area, protection_distance, degradation_weight, 
                                            "/subregion_1/summary_curves.csv"))
            
            for(temp_class in top_legend$class){
              if(approach == "target"){
                temp_thresh <- "m_low"
              } else {
                if(protected_area == "_allPA"){
                  temp_thresh <- "m_pa_low"
                } else {
                  temp_thresh <- "m_iucnpa_low"
                }
              }
              #temp_thresh
              temp_rank <- top_legend[top_legend$class==temp_class, temp_thresh]
              #temp_rank
              temp_cover <- max(temp_curve[temp_curve$rank>=temp_rank, "weighted_mean_positive"])
              #temp_cover
              
              coverage_df <- c(coverage_df, list( data.frame(
                "scenario" = paste0(approach, species_number, group_weight,
                                    protected_area, protection_distance, degradation_weight),
                "class" = temp_class,
                "rank" = temp_rank,
                "coverage" = temp_cover
              )))
            }
            
          }, silent = TRUE)}
      }}}}
coverage_df <- do.call(rbind, coverage_df)
coverage_df

write_csv(coverage_df, "_results/Zonation_performance_curve_coverage.csv")

#- - - - - - - - - - - - - - - - - - - - - -
## Comparison statistics ####
priority_rast_c <- terra::rast("_results/_Maps/Zonation_priorities.tif")

#priority_rast_diff <- priority_rast_c>=1 & priority_rast_c!=5 #30% and PAs
#priority_rast_diff <- priority_rast_c==4 #1% and PAs

#terra::plot(priority_rast_diff)

comparison_df_list <- list()
comparison_rast_list <- list()
for(temp_class in top_legend$class){
  
  temp_id <- top_legend[top_legend$class==temp_class, "ID"]
  
  if(temp_id != 5){
    priority_rast_diff <- priority_rast_c>=temp_id & priority_rast_c!=5 #class and PAs
  } else {
    priority_rast_diff <- priority_rast_c>=temp_id
  }
  
  # comparisons
  comparison_df <- data.frame(base_map = c("target_grW", "target_grW", "target_grW",
                                           "complement_grW", "complement_grW", 
                                           "complement_grW",
                                           "prevent_grW_maxSubset", "prevent_grW_maxSubset",
                                           "prevent_grW_maxSubset",
                                           "target_grW", "target_grW", "complement_grW"),
                              compare_map = c("target", "target_10_grW", "target_100_grW",
                                              "complement", "complement_grW_allPA",
                                              "complement_100_grW",
                                              "prevent_maxSubset", "prevent_grW_max",
                                              "prevent_grW_allPA_maxSubset",
                                              "complement_grW", "prevent_grW_maxSubset", "prevent_grW_maxSubset"),
                              name_diff = c("Group weights (t)", "ESM taxa (t)", "SDM taxa (t)",
                                            "Group weights (c)", "Protected areas (c)", 
                                            "SDM taxa (c)",
                                            "Group weights (p)", "Driver subset (p)",
                                            "Protected areas (p)",
                                            "c - t", "p - t", "p - c"),
                              class = temp_class,
                              n_cells = NA,
                              n_cell_diff = NA,
                              n_cell_diff_percent = NA,
                              n_cell_unchanged = NA,
                              n_cell_gain = NA,
                              n_cell_loss = NA)
  #comparison_df
  
  comparison_rast <- rast(ext(priority_rast_diff))
  for(i in 1:nrow(comparison_df)){
    temp_base <- priority_rast_diff[[comparison_df[i, "base_map"]]]
    temp_compare <- priority_rast_diff[[comparison_df[i, "compare_map"]]]
    temp_diff <- temp_compare - temp_base #diff. to compare, 1 = present in compare, -1 not present but in base
    names(temp_diff) <- comparison_df[i, "name_diff"]
    
    #terra::plot(temp_diff)
    temp_diff_df <- as.data.frame(table(values(temp_diff)))
    temp_diff_df
    
    comparison_df[i, "n_cells"] <- sum(values(temp_base), na.rm=TRUE)
    try(comparison_df[i, "n_cell_unchanged"] <- temp_diff_df[temp_diff_df$Var1==0, "Freq"] )
    try(comparison_df[i, "n_cell_gain"] <- temp_diff_df[temp_diff_df$Var1==1, "Freq"] )
    try(comparison_df[i, "n_cell_loss"] <- temp_diff_df[temp_diff_df$Var1==-1, "Freq"] )
    try(comparison_df[i, "n_cell_diff"] <- sum(c(unlist(comparison_df[i, "n_cell_loss"]), 
                                                        unlist(comparison_df[i, "n_cell_gain"])), 
                                               na.rm=TRUE))
    try(comparison_df[i, "n_cell_diff_percent"] <- comparison_df[i, "n_cell_diff"] / comparison_df[i, "n_cells"] )
    
    comparison_rast <- c(comparison_rast, temp_diff)
  }
  
  comparison_rast <- list(comparison_rast)
  names(comparison_rast) <- temp_class
  
  comparison_rast_list <- c(comparison_rast_list, comparison_rast)
  comparison_df_list <- c(comparison_df_list, list(comparison_df))
  #terra::plot(comparison_rast)
}
comparison_rast_list
comparison_df_list

comparison_df <- do.call(rbind, comparison_df_list) 

write_csv(comparison_df, "_results/Zonation_comparisons.csv")

## check out
comparison_df <- read_csv("_results/Zonation_comparisons.csv")

# coverage stats
coverage_df <- read_csv("_results/Zonation_performance_curve_coverage.csv")
coverage_df_sum <- coverage_df %>%
  dplyr::select(-rank) %>%
  mutate(class = paste0("coverage_", class)) %>%
  pivot_wider(names_from = class,
              values_from = coverage)
coverage_df_sum
write_csv(coverage_df_sum, "_results/Zonation_performance_curve_coverage_summary.csv")

## OLD: same coverage for scenarios except target against others
# comparison_df2 <- comparison_df %>% 
#   left_join(coverage_df_sum, by = c("base_map" = "scenario")) %>%
#   left_join(coverage_df_sum, by = c("compare_map" = "scenario"), suffix = c("_base", "_compare"))



## How well do PAs perform? ####
comparison_df %>% dplyr::filter(base_map == "target_grW" &
                                  compare_map == "complement_grW")
comparison_df %>% dplyr::filter(base_map == "complement_grW"  & 
                                  compare_map == "complement_grW_allPA")

coverage_df_sum %>% dplyr::filter(scenario == "target_grW" | scenario == "complement_grW" | 
                                    scenario == "complement_grW_allPA")
coverage_df_sum %>% dplyr::filter(scenario == "target_grW" | scenario == "complement_grW" | 
                                    scenario == "complement_grW_allPA") %>%
  mutate(across(where(is.numeric), function(x) x - coverage_PA))


# group weight to all same weight species at least targeted
comparison_df %>% dplyr::filter(base_map == "target_grW" &
                                  compare_map == "target")
coverage_df_sum %>% dplyr::filter(scenario == "target_grW" | scenario == "target")


# all species to 10 and 100 at least targeted
comparison_df %>% dplyr::filter(base_map == "target_grW" &
                                  compare_map == "target_100_grW" |
                                  base_map == "target_grW" &
                                  compare_map == "target_10_grW" )
coverage_df_sum %>% dplyr::filter(scenario == "target_grW" | 
                                    scenario == "target_10_grW" | 
                                    scenario == "target_100_grW")
#...maps

# - all PA to IUCN PA x group weight to same weight

# - target to complement: all and IUCN PA and group to same weight
# - prevent max and subset
# - prevent subset to complement grW

# read curves
curve_target <- read_delim(paste0(zonation_dir, "/", "target_grW",
                                "/subregion_1/summary_curves.csv"))
curve_complement <- read_delim(paste0(zonation_dir, "/", "complement_grW",
                                      "/subregion_1/summary_curves.csv"))
curve_complementAll <- read_delim(paste0(zonation_dir, "/", "complement_grW_allPA",
                                         "/subregion_1/summary_curves.csv"))
top_legend
temp_cover <- max(temp_curve[temp_curve$rank>=temp_rank, "weighted_mean_positive"])


# theme_curve <- theme_set(theme_bw()+
#                            theme(panel.border = element_blank(),
#                                  panel.grid.minor = element_blank(),
#                                  legend.position = "inside",
#                                  legend.position.inside = c(0.3, 0.3),
#                                  axis.title = element_text(size = 30),
#                                  axis.text = element_text(size = 20),
#                                  legend.text = element_text(size = 20),
#                                  legend.title = element_text(size = 30)))
# theme_set(theme_curve)
ggplot()+
  geom_line(data = curve_target, aes(x = rank, y = weighted_mean_positive, 
                                     color = "targeted", lty = "targeted"),
            lwd = 1.5)+
  geom_line(data = curve_complement, aes(x = rank, y = weighted_mean_positive,
                                  color = "complementing",
                                  lty = "complementing"),
            lwd = 1.5)+
  geom_line(data = curve_complementAll, aes(x = rank, y = weighted_mean_positive,
                                     color = "complementing (all PAs)",
                                     lty = "complementing (all PAs)"),
            lwd = 1.5)+
  xlab("Priority rank")+ 
  ylab("Weighted mean coverage of features")+
  scale_color_manual(values = c("targeted" = "#E18616",
                              "complementing" = "#059041",
                              "complementing (all PAs)" = "#059041"),
                     name = "Approach")+
  scale_linetype_manual(values = c("targeted" = "solid",
                                 "complementing" = "solid",
                                 "complementing (all PAs)" = "twodash"),
                        name = "Approach")+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.background = element_rect(color = "black"),
        legend.position.inside = c(0.4, 0.2),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 15))
ggsave("_figures/Zonation_performance.png", height = 5, width = 5)


# more performance curves
curve_target_noGrW <- read_delim(paste0(zonation_dir, "/", "target",
                                        "/subregion_1/summary_curves.csv"))
curve_target_10 <- read_delim(paste0(zonation_dir, "/", "target_10_grW",
                                        "/subregion_1/summary_curves.csv"))
curve_target_100 <- read_delim(paste0(zonation_dir, "/", "target_100_grW",
                                     "/subregion_1/summary_curves.csv"))

ggplot()+
  geom_line(data = curve_target_noGrW, aes(x = rank, y = weighted_mean_positive, 
                                           color = "targeted (no group weight)", 
                                           lty = "targeted (no group weight)",
                                           lwd = "targeted (no group weight)"))+
  geom_line(data = curve_target_100, aes(x = rank, y = weighted_mean_positive, 
                                           color = "targeted (SDMs)", 
                                           lty = "targeted (SDMs)",
                                          lwd = "targeted (SDMs)"))+
  geom_line(data = curve_target_10, aes(x = rank, y = weighted_mean_positive, 
                                           color = "targeted (ESMs)", 
                                           lty = "targeted (ESMs)",
                                           lwd = "targeted (ESMs)"))+
  geom_line(data = curve_target, aes(x = rank, y = weighted_mean_positive, 
                                     color = "targeted", lty = "targeted",
                                     lwd = "targeted"),)+
  geom_line(data = curve_complement, aes(x = rank, y = weighted_mean_positive,
                                         color = "complementing",
                                         lty = "complementing",
                                         lwd = "complementing"))+
  geom_line(data = curve_complementAll, aes(x = rank, y = weighted_mean_positive,
                                            color = "complementing (all PAs)",
                                            lty = "complementing (all PAs)",
                                            lwd = "complementing (all PAs)"))+
  xlab("Priority rank")+ 
  ylab("Weighted mean coverage of features")+
  scale_color_manual(values = c("targeted" = "#E18616",
                                "targeted (no group weight)" = "#E18616",
                                "targeted (SDMs)" = "#E18616",
                                "targeted (ESMs)" = "#E18616",
                                "complementing" = "#059041",
                                "complementing (all PAs)" = "#059041"),
                     name = "Approach")+
  scale_linetype_manual(values = c("targeted" = "solid",
                                   "targeted (no group weight)" = "twodash",
                                   "targeted (SDMs)" = "dotted",
                                   "targeted (ESMs)" = "dotdash",
                                   "complementing" = "solid",
                                   "complementing (all PAs)" = "twodash"),
                        name = "Approach")+
  scale_linewidth_manual(values = c("targeted" = 1.5,
                                   "targeted (no group weight)" = 1,
                                   "targeted (SDMs)" = 1,
                                   "targeted (ESMs)" = 1,
                                   "complementing" = 1.5,
                                   "complementing (all PAs)" = 1),
                        name = "Approach")+
  scale_x_continuous(limits = c(0,1), expand = c(0,0))+
  scale_y_continuous(limits = c(0,1), expand = c(0,0))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "inside",
        legend.background = element_rect(color = "black"),
        legend.position.inside = c(0.4, 0.2),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 15))
ggsave("_figures/Zonation_performance_Appendix.png", height = 5, width = 5)

# terra::plot(c(
#   terra::subset(priority_rast_c, c("target_grW", "complement_grW")),
#   comparison_rast$`c - t`), nc=3)


#- - - - - - - - - - - - - - - - - - - - - -
## Compare with richness maps ####
features <- read_delim(paste0(zonation_dir, "/features.txt"))
features_10 <- read_delim(paste0(zonation_dir, "/features_10.txt"))
features_100 <- read_delim(paste0(zonation_dir, "/features_100.txt"))

richness <- lapply(unique(features$group), function(x){
  temp_features <- features %>% filter(group == x)
  temp_rast <- terra::rast(paste0(zonation_dir, "/data/", temp_features$filename))
  temp_rast <- terra::app(temp_rast, sum)/1000
  names(temp_rast) <- str_extract(temp_features$filename[1], "(?<=\\.\\./data/)[^/]+")
  temp_rast
})
richness <- do.call(c, richness)

terra::plot(richness)
terra::writeRaster(richness, "_results/_Maps/Richness_allTaxa_features.tif", overwrite=TRUE)

pdf("_figures/Richness_allTaxa_features.pdf")
terra::plot(richness)
dev.off()
