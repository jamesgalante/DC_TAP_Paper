# Script: distance_by_effect_size_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/distance_by_effect_size_plots_FDR_wo_pos_controls.rda"))
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
screen_results <- read_tsv(snakemake@input$screen_results)

# Filter for Random_DistalElement_Gene
results_with_element_gene_pair_categories_modified <- screen_results %>% 
  filter(significant_wo_pos_controls | power_at_effect_size_15_wo_pos_controls >= 0.8) %>%
  filter(Random_DistalElement_Gene)


###  Doing Thing =========================================================

# Function to create and save plot for each reference
plot_distance_by_es <- function(df, cell_type_to_plot) {
  
  # Filter the table for the cell type
  df <- df %>% 
    filter(cell_type == cell_type_to_plot)
  
  # Convert the log_2_fold_change to log_10_fold_change
  df <- df %>%
    mutate(percent_change = pct_change_effect_size,
           distance = abs(distance_to_gencode_gene_TSS)) %>%             # Use absolute distance
    arrange(significant) %>% # Order by significance so these points appear on top of the insignificant points
    mutate(effect_color = case_when(
      significant == TRUE & percent_change > 0 ~ "Positive",
      significant == TRUE & percent_change < 0 ~ "Negative",
      TRUE ~ "Non-significant"
    ))
  
  # Create the plot
  p <- ggplot(df, aes(x = distance / 1e6, y = percent_change, color = effect_color)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # Dotted line at y = 0
    scale_color_manual(values = c("Positive" = "blue", "Negative" = "red", "Non-significant" = "grey"),
                       labels = c("Positive" = "Pos.", "Negative" = "Neg.", "Non-significant" = "Non-sig")) +
    ylim(-100, 100) +
    labs(title = "Distance by Effect Size",
         subtitle = paste0(cell_type_to_plot, " discovery pairs"),
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
  
  # Suppress automatic PDF creation that ggMarginal can trigger
  pdf(NULL)
  p <- ggMarginal(p, type = "density", margins = "both", groupColour = TRUE, groupFill = FALSE)
  dev.off()
  
  return(p)
}

# Create the plots
p1 <- plot_distance_by_es(results_with_element_gene_pair_categories_modified, "K562")
p2 <- plot_distance_by_es(results_with_element_gene_pair_categories_modified, "WTC11")


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