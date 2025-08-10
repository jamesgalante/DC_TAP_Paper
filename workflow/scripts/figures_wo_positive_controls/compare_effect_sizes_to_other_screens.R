# Script: compare_effect_sizes_to_other_screens.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/compare_effect_sizes_to_other_screens_FDR_wo_pos_controls.rda"))
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
  library(data.table)
})

message("Loading input files")
gasperini_results <- fread(snakemake@input$gasperini_results, sep = "\t")
dc_tap_results <- fread(snakemake@input$dc_tap_results, sep = "\t")
klann <- fread(snakemake@input$klann, sep = "\t")
morrisv1 <- fread(snakemake@input$morrisv1, sep = "\t")
morrisv2 <- fread(snakemake@input$morrisv2, sep = "\t")
xie <- fread(snakemake@input$xie, sep = "\t")


### FORMAT RESUTLS ============================================================

# Filter for significant valid pairs
klann_filt <- klann %>% filter(ValidConnection == T, Significant, EffectSize < 0)
morrisv1_filt <- morrisv1 %>% filter(ValidConnection == T, Significant, EffectSize < 0)
morrisv2_filt <- morrisv2 %>% filter(ValidConnection == T, Significant, EffectSize < 0)
xie_filt <- xie %>% filter(ValidConnection == T, Significant, EffectSize < 0)
gasperini_filt <- gasperini_results %>% filter(DistalElement_Gene, significant, pct_change_effect_size < 0)
k562_dc_tap_filt <- dc_tap_results %>% filter(Random_DistalElement_Gene, significant_wo_pos_controls, pct_change_effect_size < 0, cell_type == "K562")
wtc11_dc_tap_filt <- dc_tap_results %>% filter(Random_DistalElement_Gene, significant_wo_pos_controls, pct_change_effect_size < 0, cell_type == "WTC11")

# Convert effect sizes to percentage for datasets where needed
klann_pct <- klann_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Klann et al. (2021)", group = "K562 Perturb-Seq")
morrisv1_pct <- morrisv1_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Morris et al. (2023)", group = "K562 Perturb-Seq")
morrisv2_pct <- morrisv2_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Morris et al. (2023)", group = "K562 Perturb-Seq")
xie_pct <- xie_filt %>% mutate(pct_change = (2^EffectSize - 1)*100, dataset = "Xie et al. (2019)", group = "K562 Perturb-Seq")

# Prepare Gasperini data
gasperini_pct <- gasperini_filt %>% mutate(pct_change = pct_change_effect_size, dataset = "Gasperini et al. (2019)", group = "K562 Perturb-Seq")

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

# Panel labels like the second figure
all_datasets <- all_datasets %>%
  mutate(
    panel = ifelse(group == "DC TAP", "DC-TAP-seq", "Previous datasets")
  )

# Within-panel order (left panel first, then right panel)
order_left  <- c("K562 DC TAP Seq", "WTC11 DC TAP Seq")
order_right <- c("Gasperini et al. (2019)", "Morris et al. (2023)", "Klann et al. (2021)", "Xie et al. (2019)")

all_datasets <- all_datasets %>%
  mutate(
    dataset = factor(
      dataset,
      levels = c(order_left, order_right)
    ),
    # a compact x label (legend will carry the long names)
    dataset_short = forcats::fct_relabel(dataset, function(x) dplyr::case_when(
      x == "K562 DC TAP Seq"       ~ "K562",
      x == "WTC11 DC TAP Seq"      ~ "WTC11",
      x == "Gasperini et al. (2019)" ~ "Gasperini",
      x == "Morris et al. (2023)"    ~ "Morris",
      x == "Klann et al. (2021)"     ~ "Klann",
      x == "Xie et al. (2019)"       ~ "Xie",
      TRUE ~ x
    ))
  )

# Colors: DC-TAP colored; others greys (like the second panel)
dataset_colors <- c(
  "K562 DC TAP Seq"       = "#7A1E24",  # deep red
  "WTC11 DC TAP Seq"      = "#C00078",  # magenta
  "Gasperini et al. (2019)" = "#333333",
  "Morris et al. (2023)"    = "#666666",
  "Klann et al. (2021)"     = "#999999",
  "Xie et al. (2019)"       = "#BBBBBB"
)

# Build plot: vertical layout, two facet columns, compact x spacing
effect_size_plot <- ggplot(
  all_datasets,
  aes(x = dataset_short, y = pct_change, fill = dataset)
) +
  geom_boxplot(outlier.shape = NA, width = 0.75, alpha = 0.9, colour = "black", linewidth = 0.35) +
  geom_point(
    aes(x = dataset_short, y = pct_change),
    position = position_jitter(width = 0.25, height = 0),
    size = 0.7, alpha = 0.4, inherit.aes = FALSE, colour = "black"
  ) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x", switch = "y") +
  scale_y_continuous(
    limits = c(-100, 0),
    breaks = seq(0, -100, by = -25),
    labels = scales::percent_format(scale = 1)
  ) +
  scale_fill_manual(values = dataset_colors, guide = guide_legend(title = NULL, ncol = 1)) +
  scale_fill_manual(
    values = dataset_colors,
    breaks = c(order_left, order_right),
    drop = FALSE,
    guide = guide_legend(
      title = NULL,
      ncol = 1,
      override.aes = list(
        linetype = 0,  # no outline lines
        shape = 22,    # solid square
        size = 5,      # legend key size
        colour = NA    # no border
      )
    )
  ) +
  labs(x = NULL, y = "Effect Size (%)") +
  theme_classic(base_size = 11) +
  theme(
    strip.text.x.top   = element_text(face = "bold", size = 11),
    panel.border       = element_rect(fill = NA, colour = "black", linewidth = 0.4),
    axis.title.y       = element_text(size = 13),
    axis.text.x        = element_blank(),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold"),
    legend.position    = "right",
    legend.text        = element_text(size = 10),     # legend text size
    legend.key.size    = unit(1, "cm"),             # legend box size
    axis.ticks.x       = element_blank(),
    plot.margin        = margin(6, 6, 6, 6)
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(filename = snakemake@output$effect_size_comparison, plot = effect_size_plot, height = 2.8, width = 6, device = "pdf")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)