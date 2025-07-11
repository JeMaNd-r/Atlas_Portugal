#- - - - - - - - - - - - - - - - - - - - - -#
#                                           #
#         Animated richness map             #
#          author: Romy Zeiss               #
#            date: 2025-06-26               #
#                                           #
#- - - - - - - - - - - - - - - - - - - - - -#

library(terra)
library(dplyr)
library(ggplot2)
library(gganimate)
library(viridis)
library(gifski)

for(Taxon_name in c("Eukaryotes")){ #"Crassiclitellata", "Nematodes", "Fungi", "Protists", , "Bacteria"
  
  print(Taxon_name)
   
  # Load the raster stack
  species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))
  
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
  
  try(if(!is.vector(species10)) species10 <- species10$species)
  try(if(!is.vector(species100)) species100 <- species100$species)
  
  try(species10 <- sort(species10[species10 %in% names(species_stack)]))
  try(species100 <- sort(species100[species100 %in% names(species_stack)]))
  #speciesSub

  if(Taxon_name == "Fungi"){
    # crop to right extent
    # use env. grid
    temp_extent <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_Protists.tif"))
    species_stack <- mask(species_stack, temp_extent)
    rm(temp_extent)
  }
  
  ## ESMs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(exists("species10") & length(species10)!=0){
    if(Taxon_name == "Crassiclitellata"){ # | Taxon_name == "Eukaryotes"
      # filter uncertainty
      extent_10 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_10.RData"))) #extent_df
      extent_tif <- terra::rasterize(terra::vect(extent_10, geom = c("x", "y")), species_stack)
      
      # Optionally, mask the cropped_layer to only keep values where layer1 == 1
      species_stack <- mask(species_stack, extent_tif)
    }
    
    # Pre-allocate raster for cumulative richness
    cumulative_raster <- species_stack[[1]] * 0
    
    # Prepare a list for storing minimal data per step
    cumulative_df_list <- list()
    
    # Loop efficiently without exploding memory
    for (i in seq_along(species10)) {
      # Add species layer
      cumulative_raster <- cumulative_raster + species_stack[[species10[i]]]
      
      # Sample cumulative raster at this step
      step_df <- as.data.frame(cumulative_raster, xy = TRUE)
      names(step_df)[3] <- "Richness"
      
      # Save only relevant data for this step
      step_df$step <- i
      step_df$species_added <- species10[i]
      
      cumulative_df_list[[i]] <- step_df
    }
    
    # Combine into one dataframe after the loop
    richness_animation_df <- dplyr::bind_rows(cumulative_df_list)
    
    if(Taxon_name != "Crassiclitellata" ){#& Taxon_name != "Eukaryotes"
      # filter uncertainty
      extent_10 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_10.RData"))) #extent_df
      richness_animation_df <- left_join(extent_10, richness_animation_df, by = c("x", "y"))
    }
    
    # Load world map
    world.inp <- map_data("world")
    
    # Plot limits
    plot_xlim <- c(min(richness_animation_df$x), max(richness_animation_df$x))
    plot_ylim <- c(min(richness_animation_df$y), max(richness_animation_df$y))
    
    # Build plot
    p <- ggplot() +
      geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80", color = NA) +
      coord_cartesian(xlim = plot_xlim, ylim = plot_ylim) +
      # Fixed color for richness = 0
      geom_tile(data = richness_animation_df %>% filter(Richness == 0),
                aes(x = x, y = y), fill = "grey60") +
      
      # Viridis color for richness > 0
      geom_tile(data = richness_animation_df %>% filter(Richness > 0),
                aes(x = x, y = y, fill = Richness)) +
      scale_fill_viridis_c() +
      theme_bw() +
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "right")
    
    # Animate
    anim_esm <- p +
      transition_manual(species_added, cumulative = TRUE) +
      view_follow(fixed_x = TRUE, fixed_y = TRUE) +
      labs(title = paste0("Number of ", Taxon_name," taxa (9<n<100)"),
           subtitle = paste0("Total number of taxa: ", length(species10)),
           fill = "Richness",
           caption = "Current species added: {current_frame}") 
    #anim_esm
  
    # Render animation
    gganimate::animate(anim_esm, fps = 10, width = 800, height = 600, 
            end_pause = 20,
            renderer = gifski_renderer(paste0("_figures/SpeciesRichness_ESM_", Taxon_name, ".gif")))
  }
  
  ## SDMs - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if(exists("species100") & length(species100)!=0){
    if(Taxon_name != "Crassiclitellata"){
      # Pre-allocate raster for cumulative richness
      cumulative_raster <- species_stack[[1]] * 0
      
      # Prepare a list for storing minimal data per step
      cumulative_df_list <- list()
      
      # Loop efficiently without exploding memory
      for (i in seq_along(species100)) {
        # Add species layer
        cumulative_raster <- cumulative_raster + species_stack[[species100[i]]]
        
        # Sample cumulative raster at this step
        step_df <- as.data.frame(cumulative_raster, xy = TRUE)
        names(step_df)[3] <- "Richness"
        
        # Save only relevant data for this step
        step_df$step <- i
        step_df$species_added <- species100[i]
        
        cumulative_df_list[[i]] <- step_df
      }
      
      # Combine into one dataframe after the loop
      richness_animation_df <- dplyr::bind_rows(cumulative_df_list)
      
      
      if(Taxon_name != "Crassiclitellata"){ #Taxon_name != "Eukaryotes" & 
        extent_100 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_100.RData"))) #extent_df
        richness_animation_df <- left_join(extent_100, richness_animation_df, by = c("x", "y"))
      }
    
      # Load world map
      world.inp <- map_data("world")
      
      # Plot limits
      plot_xlim <- c(min(richness_animation_df$x), max(richness_animation_df$x))
      plot_ylim <- c(min(richness_animation_df$y), max(richness_animation_df$y))
      
      # Build plot
      p <- ggplot() +
        geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80", color = NA) +
        coord_cartesian(xlim = plot_xlim, ylim = plot_ylim) +
        # Fixed color for richness = 0
        geom_tile(data = richness_animation_df %>% filter(Richness == 0),
                  aes(x = x, y = y), fill = "grey60") +
        
        # Viridis color for richness > 0
        geom_tile(data = richness_animation_df %>% filter(Richness > 0),
                  aes(x = x, y = y, fill = Richness)) +
        scale_fill_viridis_c() +
        theme_bw() +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.position = "right")
      
      # Animate
      anim_sdm <- p +
        transition_manual(species_added, cumulative = TRUE) +
        view_follow(fixed_x = TRUE, fixed_y = TRUE) +
        labs(title = paste0("Number of ", Taxon_name," taxa (n>99)"),
             subtitle = paste0("Total number of taxa: ", length(species100)), 
             fill = "Richness",
             caption = "Current species added: {current_frame}") 
      #anim_sdm
      
      # Render animation
      gganimate::animate(anim_sdm, fps = 10, width = 800, height = 600, 
              end_pause = 20,
              renderer = gifski_renderer(paste0("_figures/SpeciesRichness_SDM_", Taxon_name, ".gif")))
    }
  }
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## ALL TAXA TOEGTHER ####
# for(Taxon_name in c("Crassiclitellata", "Nematodes", "Fungi", "Protists", "Eukaryotes", "Bacteria")){
#   
#   print(Taxon_name)
#   
#   # load number of occurrences per species and focal species names
#   species100 <- read.csv(file=paste0("_intermediates/SDM_", Taxon_name, ".csv"))
#   if(nrow(species100) != 0){
#     species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv"))
#     speciesSub <-  species100 %>% 
#       pull(species) %>%
#       c(species10 %>% pull(species))
#   }else{
#     species10 <- read.csv(file=paste0("_intermediates/ESM_", Taxon_name, ".csv")) %>% pull(species)
#     speciesSub <- species10
#   }
#   #if(Taxon_name == "Fungi") species10 <- species10[species10 != "F02025"]
#   
#   species10
#   species100
#   speciesSub
#   
#   # Load the raster stack
#   species_stack <- terra::rast(paste0("_results/_Maps/SDM_stack_binary_", Taxon_name, ".tif"))
#   species_order <- names(species_stack)[-nlyr(species_stack)]
#   
#   if(Taxon_name == "Crassiclitellata"){
#     # filter uncertainty
#     extent_10 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_10.RData"))) #extent_df
#     extent_tif <- terra::rasterize(terra::vect(extent_10, geom = c("x", "y")), species_stack)
#     
#     # Optionally, mask the cropped_layer to only keep values where layer1 == 1
#     species_stack <- mask(species_stack, extent_tif)
#   }
#   
#   # Pre-allocate raster for cumulative richness
#   cumulative_raster <- species_stack[[1]] * 0
#   
#   # Prepare a list for storing minimal data per step
#   cumulative_df_list <- list()
#   
#   # Loop efficiently without exploding memory
#   for (i in seq_along(species_order)) {
#     # Add species layer
#     cumulative_raster <- cumulative_raster + species_stack[[species_order[i]]]
#     
#     # Sample cumulative raster at this step
#     step_df <- as.data.frame(cumulative_raster, xy = TRUE)
#     names(step_df)[3] <- "Richness"
#     
#     # Save only relevant data for this step
#     step_df$step <- i
#     step_df$species_added <- species_order[i]
#     
#     cumulative_df_list[[i]] <- step_df
#   }
#   
#   # Combine into one dataframe after the loop
#   richness_animation_df <- dplyr::bind_rows(cumulative_df_list)
#   
#   if(Taxon_name != "Crassiclitellata"){
#     # filter uncertainty
#     extent_10 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_10.RData"))) #extent_df
#     extent_100 <- get(load(file=paste0("_results/SDM_Uncertainty_extent_", Taxon_name, "_100.RData"))) #extent_df
#     
#     # Logical raster where both rasters are not NA
#     overlap_mask <- inner_join(extent_10, extent_100, by = c("x", "y"))
#     richness_animation_df <- left_join(overlap_mask, richness_animation_df, by = c("x", "y"))
#   }
#   
#   # Load world map
#   world.inp <- map_data("world")
#   
#   # Plot limits
#   plot_xlim <- c(min(richness_animation_df$x), max(richness_animation_df$x))
#   plot_ylim <- c(min(richness_animation_df$y), max(richness_animation_df$y))
#   
#   # Build plot
#   p <- ggplot() +
#     geom_map(data = world.inp, map = world.inp, aes(map_id = region), fill = "grey80", color = NA) +
#     coord_cartesian(xlim = plot_xlim, ylim = plot_ylim) +
#     # Fixed color for richness = 0
#     geom_tile(data = richness_animation_df %>% filter(Richness == 0),
#               aes(x = x, y = y), fill = "grey60") +
#     
#     # Viridis color for richness > 0
#     geom_tile(data = richness_animation_df %>% filter(Richness > 0),
#               aes(x = x, y = y, fill = Richness)) +
#     scale_fill_viridis_c() +
#     labs(title = paste0("Number of ", Taxon_name," taxa (9<n<100"), fill = "Richness") +
#     theme_bw() +
#     theme(axis.title = element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           legend.position = "right")
#   
#   # Animate
#   anim <- p +
#     transition_manual(species_added, cumulative = TRUE) +
#     view_follow(fixed_x = TRUE, fixed_y = TRUE)
#   #anim
#   
#   # Render animation
#   animate(anim, fps = 10, width = 800, height = 600, 
#           end_pause = 20,
#           renderer = gifski_renderer(paste0("_figures/SpeciesRichness_ESM_", Taxon_name, ".gif")))
# not all frames, compress
# frame_indices <- seq(1, length(species_order), length.out = 100)
# 
# for (i in frame_indices) {
#   # Add species as before
# }
