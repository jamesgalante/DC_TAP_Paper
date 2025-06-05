# Script: slc2a3_effect_size_plot.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/slc2a3_effect_size_plot.rda"))
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
wtc11_singleton_results <- readRDS(snakemake@input$wtc11_singleton_results)
wtc11_calibration_check_results <- readRDS(snakemake@input$wtc11_calibration_check_results)
negative_controls_together <- readRDS(snakemake@input$negative_controls_together)


### FORMAT RESULTS ============================================================

# Filter for the pair of interest
pair_of_interest <- results_with_element_gene_pair_categories_modified %>% 
  filter(cell_type == "WTC11", element_gene_pair_identifier_hg38 == "SLC2A3|chr12:7960676-7960977")

# Filter for negative controls tested against gene of interest
neg_of_interest <- negative_controls_together %>%
  filter(startsWith(grna_target, "pseudo")) %>% 
  filter(response_id == "ENSG00000059804") %>%
  mutate(
    pct_change_effect_size = (fold_change - 1) * 100,
    standard_error_pct_change = se_fold_change * 100,
    lower_CI_95_pct_change = pct_change_effect_size - 1.96 * standard_error_pct_change,
    upper_CI_95_pct_change = pct_change_effect_size + 1.96 * standard_error_pct_change
  )

# Filter for guides of interest
guide_ids <- unlist(strsplit(pair_of_interest$guide_ids, ","))
guides_of_interest <- wtc11_singleton_results %>%
  filter(response_id == "ENSG00000059804") %>%
  filter(grna_id %in% guide_ids) %>%
  mutate(
    pct_change_effect_size = 100*(2^log_2_fold_change - 1),
    group = "chr12:8113272-8113573"  # Add group column
  )

# Filter for negative control guides tested against gene of interest
neg_guides <- wtc11_calibration_check_results %>%
  filter(response_id == "ENSG00000059804") %>%
  mutate(
    pct_change_effect_size = 100*(2^log_2_fold_change - 1),
    group = "non-targeting"  # Add group column
  )

### PREPARE PLOTTING DATA =====================================================

message("Preparing data for plotting")
# Create a data frame for the bar plot showing means and CIs
bar_data <- tibble(
  group = c("chr12:8113272-8113573", "non-targeting"),
  mean_effect = c(
    pair_of_interest$pct_change_effect_size, 
    neg_of_interest$pct_change_effect_size
  ),
  ci_lower = c(
    pair_of_interest$lower_CI_95_pct_change,
    neg_of_interest$lower_CI_95_pct_change
  ),
  ci_upper = c(
    pair_of_interest$upper_CI_95_pct_change,
    neg_of_interest$upper_CI_95_pct_change
  )
)

# Combine the guides data for plotting individual points
all_guides <- bind_rows(guides_of_interest, neg_guides)


### CREATE PLOT ===============================================================

message("Creating barplot with significance coloring")
# Set custom colors for groups and significance
group_colors <- c("chr12:8113272-8113573" = "#F8766D", "non-targeting" = "#00BFC4")
sig_colors <- c("TRUE" = "black", "FALSE" = "gray70")

# Create the barplot
slc2a3_effect_plot <- ggplot() +
  geom_bar(data = bar_data, aes(x = group, y = mean_effect, fill = group), stat = "identity",  width = 0.6, alpha = 0.7) +
  geom_jitter(data = all_guides, aes(x = group, y = pct_change_effect_size, color = significant), width = 0.2, size = 0.5, alpha = 0.7) +
  geom_errorbar(data = bar_data, aes(x = group, ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = sig_colors) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "",
    y = "Effect on SLC2A3 (%)",
    x = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  ylim(-50, 50)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(filename = snakemake@output$slc2a3_effect_size_barplot, plot = slc2a3_effect_plot, width = 2.5, height = 4, device = "pdf")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)