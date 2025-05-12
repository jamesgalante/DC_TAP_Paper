# Script: effect_size_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/effect_size_plots.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages(
  library(tidyverse)
)

message("Loading input files")
results_with_element_gene_pair_categories_modified <- read_tsv(snakemake@input$results_with_element_gene_pair_categories_modified)


### MAKE PLOT =================================================================

# Make a plot for all Random Distal Element Gene Pairs that are are significant and have a negative effect size
effect_size_boxplot <- results_with_element_gene_pair_categories_modified %>%
  filter(Random_DistalElement_Gene, significant, pct_change_effect_size < 0) %>%
  ggplot(aes(x = cell_type, y = pct_change_effect_size, fill = cell_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.9, alpha = 0.8) +
  labs(
    x = "Cell Type",
    y = "Effect Size (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none"
  )

# Let's look at all pairs labelled as tss_pos in the screen design
positive_control_self_promoters_effect_size_boxplot <- results_with_element_gene_pair_categories_modified %>%
  filter(Positive_Control_selfPromoter) %>%
  filter(significant | power_at_effect_size_15 >= 0.8) %>%
  ggplot(aes(x = cell_type, y = pct_change_effect_size)) +
  geom_boxplot(aes(fill = cell_type), outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(aes(color = significant), width = 0.2, size = 0.9, alpha = 0.8) +
  scale_color_manual(values = c("FALSE" = "darkgrey", "TRUE" = "black")) +
  labs(
    x = "Cell Type",
    y = "Effect Size (%)"
  ) +
  theme_classic()


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = effect_size_boxplot, filename = snakemake@output$effect_size_boxplot,
       device = "pdf", height = 2.5, width = 2.3)
ggsave(plot = positive_control_self_promoters_effect_size_boxplot, filename = snakemake@output$positive_control_self_promoters_effect_size_boxplot,
       device = "pdf", height = 2.5, width = 3)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)