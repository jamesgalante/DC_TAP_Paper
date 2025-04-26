# Script: plot_promoter_overlapping_elements.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/plot_promoter_overlapping_elements.rda"))
message("Saved Image")
stop("Manually Stopped Program after Saving Image")

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
  library(stats)
})

message("Loading input files")
Formatted_DC_TAP_Seq_Results_w_Categories <- read_tsv(snakemake@input$Formatted_DC_TAP_Seq_Results_w_Categories)


### BARPLOTS OF PROMOTER PAIRS ================================================

# Let's see how the selfPromoters worked
selfPromoters <- Formatted_DC_TAP_Seq_Results_w_Categories %>%
  filter(selfPromoter) %>%
  ggplot(aes(x = cell_type, y = pct_change_effect_size, fill = significant)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), color = "black", alpha = 0.8, size = 0.8) +
  labs(
    title = "selfPromoters",
    x = "Cell Type",
    y = "Effect Size (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )
ggsave(plot = selfPromoters, filename = snakemake@output$selfPromoter_boxplot,
       device = "pdf", height = 3, width = 3)

# Let's see the DistalPromoter_Gene pairs
DistalPromoter_Genes <- Formatted_DC_TAP_Seq_Results_w_Categories %>%
  filter(DistalPromoter_Gene) %>%
  filter(significant) %>%
  ggplot(aes(x = cell_type, y = pct_change_effect_size, fill = cell_type)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), color = "black", alpha = 0.8, size = 0.8) +
  labs(
    title = "DistalPromoter_Gene",
    x = "Cell Type",
    y = "Effect Size (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  )
ggsave(plot = DistalPromoter_Genes, filename = snakemake@output$DistalPromoter_Gene_boxplot,
       device = "pdf", height = 3, width = 3)


### DISTANCE BY EFFECT SIZE ===================================================

# Let's see these DistalPromoter_Gene pairs on a distance by effect size plot
# If these were trans effects, we wouldn't assume any distance preference
# We also wouldn't assume any effect size direction preference
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
         subtitle = paste0(cell_type_to_plot, " DistalPromoter_Gene pairs"),
         x = "Distance to TSS (Mb)", 
         y = "CRISPRi effect size \n(% change)", 
         color = "CRISPR") +
    theme_classic(base_size = 14) +  # Increase base font size
    theme(
      plot.title = element_text(vjust = -1, margin = margin(b = 10), size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      strip.background = element_blank(),
      legend.position = "bottom",
      aspect.ratio = 1
    )
  
  p <- ggMarginal(p, type = "density", margins = "both", groupColour = TRUE, groupFill = FALSE)
  
  return(p)
}

# Create the plots
p1 <- plot_distance_by_es(Formatted_DC_TAP_Seq_Results_w_Categories %>% filter(DistalPromoter_Gene), "K562")
p2 <- plot_distance_by_es(Formatted_DC_TAP_Seq_Results_w_Categories %>% filter(DistalPromoter_Gene), "WTC11")
# It looks like the density of pos, neg, and non-sig is different in WTC11 versus K562


### EXPLORING H3K27ac BY EFFECT SIZE ==========================================

# DistalPromoter_Gene pairs EffectSize correlates with h3k27ac
DistalPromoter_Gene_ES_by_Element_Category <- Formatted_DC_TAP_Seq_Results_w_Categories %>%
  filter(DistalPromoter_Gene) %>%
  filter(significant) %>%
  mutate(element_category = case_when(
    element_category == "H3K27ac low" ~ "low",
    element_category == "H3K27ac medium" ~ "medium",
    element_category == "H3K27ac high" ~ "high"
  )) %>%
  ggplot(aes(x = factor(element_category, level = c("low", "medium", "high")), y = abs(pct_change_effect_size), fill = cell_type)) + 
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5, outlier.size = 0) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.8, size = 0.8) +
  labs(
    title = "DistalPromoter_Genes by Element Category",
    x = "Element Category (H3K27ac Level)",
    y = "Absolute Effect Size (%)"
  ) +
  theme_classic()
ggsave(plot = DistalPromoter_Gene_ES_by_Element_Category, filename = snakemake@output$DistalPromoter_Gene_ES_by_Element_Category,
       device = "pdf", height = 3, width = 3.5)

# DistalPromoter_Gene pairs EffectSize correlates with h3k27ac
selfPromoter_ES_by_Element_Category <- Formatted_DC_TAP_Seq_Results_w_Categories %>%
  filter(selfPromoter) %>%
  filter(significant) %>%
  mutate(element_category = case_when(
    element_category == "H3K27ac low" ~ "low",
    element_category == "H3K27ac medium" ~ "medium",
    element_category == "H3K27ac high" ~ "high",
    element_category == "CTCF overlap" | element_category == "H3K27me3 overlap" ~ "other",
    TRUE ~ element_category
  )) %>%
  ggplot(aes(x = factor(element_category, level = c("low", "medium", "high", "other")), y = abs(pct_change_effect_size), fill = cell_type)) + 
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5, outlier.size = 0) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.8, size = 0.8) +
  labs(
    title = "DistalPromoter_Genes by Element Category",
    x = "Element Category (H3K27ac Level)",
    y = "Absolute Effect Size (%)"
  ) +
  theme_classic()
ggsave(plot = selfPromoter_ES_by_Element_Category, filename = snakemake@output$selfPromoter_ES_by_Element_Category,
       device = "pdf", height = 3, width = 3.5)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = p1, filename = snakemake@output$k562_distance_by_es_DistalPromoter_Gene, device = "pdf", height = 5, width = 4.5)
ggsave(plot = p2, filename = snakemake@output$wtc11_distance_by_es_DistalPromoter_Gene, device = "pdf", height = 5, width = 4.5)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)
