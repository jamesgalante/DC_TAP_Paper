# Script: distance_by_effect_size_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/distance_by_effect_size_plots.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggExtra)
})

message("Loading input files")
k562_discovery_results <- readRDS(snakemake@input$discovery_results[[1]])
wtc11_discovery_results <- readRDS(snakemake@input$discovery_results[[2]])

k562_final_sceptre_object <- readRDS(snakemake@input$final_sceptre_object[[1]])
wtc11_final_sceptre_object <- readRDS(snakemake@input$final_sceptre_object[[2]])

k562_distances <- read_tsv(snakemake@input$distances[[1]])
wtc11_distances <- read_tsv(snakemake@input$distances[[2]])


###  Doing Thing =========================================================

# Function to create and save plot for each reference
plot_distance_by_es <- function(discovery_results, distances, cell_type) {
  
  # Combine the discovery_results and distances file 
  merged <- distances %>% inner_join(discovery_results, by = c("response_id", "grna_group" = "grna_target"))
  
  # Convert the log_2_fold_change to log_10_fold_change
  merged <- merged %>%
    mutate(fold_change = 2^log_2_fold_change,
           percent_change = (fold_change - 1) * 100,  # Calculate percent change
           distance = abs(distance)) %>%             # Use absolute distance
    mutate(distance = abs(distance)) %>% # Make the distance the absolute value of the distance from gene
    arrange(significant) %>% # Order by significance so these points appear on top of the insignificant points
    filter(!is.na(significant)) %>% # Remove any NA significant points
    filter(target_type == "enh") %>% # Filter for only "enh"
    mutate(effect_color = case_when(
      significant == TRUE & percent_change > 0 ~ "Positive",
      significant == TRUE & percent_change < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    ))
  
  # Create the plot
  p <- ggplot(merged, aes(x = distance / 1e6, y = percent_change, color = effect_color)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Dotted line at y = 0
    scale_color_manual(values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "grey"),
                       labels = c("Positive" = "Pos.", "Negative" = "Neg.", "Non-significant" = "Non-sig")) +
    ylim(-100, 100) +
    labs(title = "Distance by Effect Size",
         subtitle = paste0(cell_type, " discovery pairs"),
         x = "Distance to TSS (Mb)", 
         y = "CRISPRi effect size \n(% change)", 
         color = "CRISPR") +
    theme_classic(base_size = 14) +  # Increase base font size
    theme(
      plot.title = element_text(vjust = -1, margin = margin(b = 10), size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.background = element_blank(),
      legend.position = "bottom"
    )
  
  p <- ggMarginal(p, type = "density", margins = "both", groupColour = TRUE, groupFill = FALSE)
  
  return(p)
}

# Create the plots
p1 <- plot_distance_by_es(k562_discovery_results, k562_distances, "K562")
p2 <- plot_distance_by_es(wtc11_discovery_results, wtc11_distances, "WTC11")


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = p1, filename = snakemake@output$k562_distance_by_es, device = "pdf", height = 5, width = 4.5)
ggsave(plot = p2, filename = snakemake@output$wtc11_distance_by_es, device = "pdf", height = 5, width = 4.5)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)