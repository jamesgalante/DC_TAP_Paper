# Script: power_to_detect_change_for_all_datasets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/power_to_detect_change_for_all_datasets.rda"))
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


### HARMONIZE THE DISTANCES ===================================================

# The klann, morrisv1, morrisv2, and xie datasets only test pairs up to 1Mb
gasperini_results <- gasperini_results %>% filter(distance_to_gencode_gene_TSS < 1e6)
dc_tap_results <- dc_tap_results %>% filter(distance_to_gencode_gene_TSS < 1e6)


### COMBINE & PREP DATASETS ===================================================

# Merge Morris datasets
morris <- bind_rows(morrisv1, morrisv2)

# Split DC-TAP by cell type
k562_dc_tap <- dc_tap_results %>% filter(cell_type == "K562")
wtc11_dc_tap <- dc_tap_results %>% filter(cell_type == "WTC11")

# Order the datasets
datasets <- list(
  "K562" = k562_dc_tap,
  "WTC11" = wtc11_dc_tap,
  "Gasperini" = gasperini_results,
  "Morris" = morris,
  "Klann" = klann,
  "Xie" = xie
)


### POWER CALCULATION =========================================================

get_power_col <- function(data, es) {
  regex <- paste0("(?i)(PowerAtEffectSize|power_at_effect_size_)", es, "$")
  col <- names(data)[grepl(regex, names(data), perl = TRUE)]
  if (length(col)) col[1] else NA_character_
}

calc_power_props <- function(data, dataset_name) {
  map_dfr(c(10, 15, 20, 25, 50), function(es) {
    col <- get_power_col(data, es)
    if (is.na(col)) {
      warning(dataset_name, ": no column for effect-size ", es)
      return(tibble(effect_size = es,
                    well_powered = NA_integer_,
                    total_pairs = NA_integer_,
                    proportion = NA_real_,
                    dataset = dataset_name))
    }
    well <- sum(data[[col]] >= 0.8, na.rm = TRUE)
    total <- sum(!is.na(data[[col]]))
    tibble(effect_size = es,
           well_powered = well,
           total_pairs = total,
           proportion = well / total,
           dataset = dataset_name)
  })
}

message("Calculating power proportions for all effect sizes (2%, 3%, 5%, 10%, 15%, 20%, 25%, 50%)")
all_power_props <- imap_dfr(datasets, calc_power_props)


### RENAME DATASETS FOR PLOTTING =============================================

all_power_props <- all_power_props %>%
  mutate(dataset = case_when(
    dataset == "K562" ~ "K562 DC TAP Seq",
    dataset == "WTC11" ~ "WTC11 DC TAP Seq",
    dataset == "Gasperini" ~ "Gasperini et al. 2019",
    dataset == "Morris" ~ "Morris et al. 2023",
    dataset == "Klann" ~ "Klann et al. 2021",
    dataset == "Xie" ~ "Xie et al. 2019",
    TRUE ~ dataset
  ))

dataset_order <- c(
  "K562 DC TAP Seq",
  "WTC11 DC TAP Seq",
  "Gasperini et al. 2019",
  "Morris et al. 2023",
  "Klann et al. 2021",
  "Xie et al. 2019"
)
all_power_props$dataset <- factor(all_power_props$dataset, levels = dataset_order)


### VISUALISATION =============================================================

dataset_colors <- c(
  "K562 DC TAP Seq" = "#993333",
  "WTC11 DC TAP Seq" = "#000099",
  "Gasperini et al. 2019" = "#336633",
  "Morris et al. 2023" = "#E69F00",
  "Klann et al. 2021" = "#F0E442",
  "Xie et al. 2019" = "#56B4E9"
)

# 1. Create percentage/proportion plot
message("Creating proportion plot")
percentage_plot <- ggplot(all_power_props, aes(x = factor(effect_size), y = proportion * 100, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = dataset_colors) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20)) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  labs(title = "Proportion of Well-Powered Pairs",
       x = "Effect Size",
       y = "Proportion of Well-Powered Pairs (%)") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

# 2. Create count plot
message("Creating count plot")
count_plot <- ggplot(all_power_props, aes(x = factor(effect_size), y = well_powered, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = dataset_colors) +
  scale_x_discrete(labels = function(x) paste0(x, "%")) +
  labs(title = "Number of Well-Powered Pairs",
       x = "Effect Size",
       y = "Number of Well-Powered Pairs") +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))



### SAVE OUTPUT ===============================================================

# Save output files
message("Saving proportion plot")
ggsave(filename = snakemake@output$power_for_all_datasets_by_proportion, 
       plot = percentage_plot, width = 10, height = 6, device = "pdf")

message("Saving count plot")
ggsave(filename = snakemake@output$power_for_all_datasets_by_count, 
       plot = count_plot, width = 10, height = 6, device = "pdf")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)