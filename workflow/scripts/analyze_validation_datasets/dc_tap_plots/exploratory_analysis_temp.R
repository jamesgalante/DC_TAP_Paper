# Script: exploratory_analysis_temp.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/exploratory_analysis_temp.rda"))
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
Formatted_DC_TAP_Seq_Results_w_Categories <- read_tsv(snakemake@input$Formatted_DC_TAP_Seq_Results_w_Categories)


### X =========================================================

# Distribution of effect sizes across different elements
es_by_element <- Formatted_DC_TAP_Seq_Results_w_Categories %>%
  filter(significant) %>%
  ggplot(aes(x = element_category, y = pct_change_effect_size, fill = cell_type)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), color = "black", alpha = 0.8, size = 0.7) +
  labs(
    title = "",
    x = "Element Category",
    y = "Effect Size (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust=1)
  )
ggsave(plot = es_by_element, filename = "es_element.pdf",
       device = "pdf", height = 4, width = 4)

# Distribution of distance across different elements
dist_by_element <- Formatted_DC_TAP_Seq_Results_w_Categories %>%
  filter(significant) %>%
  ggplot(aes(x = element_category, y = abs(distance_to_gencode_gene_TSS), fill = cell_type)) +
  geom_boxplot(position = position_dodge(width = 0.75), alpha = 0.5, outlier.shape = NA) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), color = "black", alpha = 0.8, size = 0.7) +
  labs(
    title = "",
    x = "Element Category",
    y = "Distance"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust=1)
  )
ggsave(plot = dist_by_element, filename = "dist_element.pdf",
       device = "pdf", height = 4, width = 4)

# Identification of genes with multiple enhancers
Formatted_DC_TAP_Seq_Results_w_Categories %>%
  group_by(cell_type) %>%
  filter(gene_symbol)

# Identification of enhancers affecting multiple genes



### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)