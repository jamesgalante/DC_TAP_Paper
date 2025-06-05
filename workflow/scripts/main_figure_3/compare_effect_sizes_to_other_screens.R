# Script: compare_effect_sizes_to_other_screens.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/compare_effect_sizes_to_other_screens.rda"))
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
gasperini_results <- read_tsv(snakemake@input$gasperini_results)
dc_tap_results <- read_tsv(snakemake@input$dc_tap_results)
klann <- read_tsv(snakemake@input$klann)
morrisv1 <- read_tsv(snakemake@input$morrisv1)
morrisv2 <- read_tsv(snakemake@input$morrisv2)
xie <- read_tsv(snakemake@input$xie)


### FORMAT RESUTLS ============================================================

# Filter for significant valid pairs
klann_filt <- klann %>% filter(ValidConnection == T, Significant, EffectSize < 0)
morrisv1_filt <- morrisv1 %>% filter(ValidConnection == T, Significant, EffectSize < 0)
morrisv2_filt <- morrisv2 %>% filter(ValidConnection == T, Significant, EffectSize < 0)
xie_filt <- xie %>% filter(ValidConnection == T, Significant, EffectSize < 0)
gasperini_filt <- gasperini_results %>% filter(DistalElement_Gene, significant, pct_change_effect_size < 0)
k562_dc_tap_filt <- dc_tap_results %>% filter(Random_DistalElement_Gene, significant, pct_change_effect_size < 0, cell_type == "K562")
wtc11_dc_tap_filt <- dc_tap_results %>% filter(Random_DistalElement_Gene, significant, pct_change_effect_size < 0, cell_type == "WTC11")

# Convert effect sizes to percentage for datasets where needed
klann_pct <- klann_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Klann et al. 2021", group = "K562 Perturb-Seq")
morrisv1_pct <- morrisv1_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Morris et al. 2023", group = "K562 Perturb-Seq")
morrisv2_pct <- morrisv2_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Morris et al. 2023", group = "K562 Perturb-Seq")
xie_pct <- xie_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Xie et al. 2019", group = "K562 Perturb-Seq")

# Prepare Gasperini data
gasperini_pct <- gasperini_filt %>% mutate(pct_change = pct_change_effect_size, dataset = "Gasperini et al. 2019", group = "K562 Perturb-Seq")

# Prepare DC-TAP data by cell type
k562_dc_tap_pct <- k562_dc_tap_filt %>% mutate(pct_change = pct_change_effect_size, dataset = "K562 DC TAP Seq", group = "DC TAP")
wtc11_dc_tap_pct <- wtc11_dc_tap_filt %>% mutate(pct_change = pct_change_effect_size, dataset = "WTC11 DC TAP Seq", group = "DC TAP")

# Combine all datasets
all_datasets <- bind_rows(
  klann_pct,
  morrisv1_pct,
  morrisv2_pct,
  xie_pct,
  gasperini_pct,
  k562_dc_tap_pct,
  wtc11_dc_tap_pct
) %>%
  # Select only needed columns
  select(dataset, pct_change, group)


### PLOT THE EFFECT SIZES ====================================================

# Set the order of datasets for plotting
dataset_order <- c("WTC11 DC TAP Seq", "K562 DC TAP Seq", "Gasperini et al. 2019", "Morris et al. 2023", "Xie et al. 2019", "Klann et al. 2021")
all_datasets$dataset <- factor(all_datasets$dataset, levels = dataset_order)

# Set the colors for the comparison
dataset_colors <- c(
  "K562 DC TAP Seq"       = "#993333",
  "WTC11 DC TAP Seq"      = "#000099",
  "Gasperini et al. 2019" = "#336633",
  "Morris et al. 2023"    = "#E69F00",
  "Xie et al. 2019"       = "#56B4E9",
  "Klann et al. 2021"     = "#F0E442"
)

# Create a boxplot with consistent widths across facets
effect_size_plot <- ggplot(all_datasets, aes(x = pct_change, y = dataset, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
  geom_jitter(height = 0.25, alpha = 0.5, size = 1, colour = "black") +
  facet_grid(group ~ ., scales = "free", space = "free", switch = "y") +
  theme_classic() +
  labs(
    title = "",
    x = "Effect Size (%)",
    y = NULL
  ) +
  theme(
    axis.text.y        = element_text(hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.border       = element_rect(fill = NA, colour = "black"),
    legend.position    = "none",
    strip.text.y       = element_blank(),
    strip.background   = element_blank()
  ) +
  # reverse x-axis: 0 % → –100 %
  scale_x_reverse(limits = c(0, -100), labels = scales::percent_format(scale = 1)) +
  # apply custom fill colours
  scale_fill_manual(values = dataset_colors)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(filename = snakemake@output$effect_size_comparison, plot = effect_size_plot, height = 4.5, width = 4, device = "pdf")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)