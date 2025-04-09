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
experiment_summary_table <- read_tsv(snakemake@input$experiment_summary_table)


### MAKE PLOT =================================================================

effect_size_boxplot <- experiment_summary_table %>%
  filter(Random_DistalElement_Gene, Significant, EffectSize < 0) %>%
  ggplot(aes(x = ExperimentCellType, y = EffectSize*100, fill = ExperimentCellType)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.9, alpha = 0.8) +
  labs(
    x = "Cell Type",
    y = "Effect Size (%)"
  ) +
  ylim(-25, 0) +
  theme_classic() +
  theme(
    legend.position = "none"
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = effect_size_boxplot, filename = snakemake@output$effect_size_boxplot,
       device = "pdf", height = 2.5, width = 2.3)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)