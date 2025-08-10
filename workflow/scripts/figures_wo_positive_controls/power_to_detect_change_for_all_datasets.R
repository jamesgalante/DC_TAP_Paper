# Script: power_to_detect_change_for_all_datasets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/power_to_detect_change_for_all_datasets_FDR_wo_pos_controls.rda"))
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

# Prefer *_wo_pos_controls for DC-TAP; otherwise use regular names (snake or camel)
get_power_col <- function(nms, es, prefer_wo = FALSE) {
  base_snake <- paste0("power_at_effect_size_", es)
  base_camel <- paste0("PowerAtEffectSize", es)
  wo_snake   <- paste0(base_snake, "_wo_pos_controls")
  wo_camel   <- paste0(base_camel, "_wo_pos_controls")
  
  # Try wo_pos_controls first if requested
  if (prefer_wo) {
    if (wo_snake %in% nms) return(wo_snake)
    if (wo_camel %in% nms) return(wo_camel)
  }
  # Otherwise, or fallback
  if (base_snake %in% nms) return(base_snake)
  if (base_camel %in% nms) return(base_camel)
  
  NA_character_
}

calc_power_props <- function(df, dataset_name, prefer_wo = FALSE) {
  es_list <- c(10, 15, 20, 25, 50)
  
  purrr::map_dfr(es_list, function(es) {
    col <- get_power_col(names(df), es, prefer_wo = prefer_wo)
    if (is.na(col)) {
      warning(dataset_name, ": no column for effect-size ", es)
      return(tibble(effect_size = es,
                    well_powered = NA_integer_,
                    total_pairs  = NA_integer_,
                    proportion   = NA_real_,
                    dataset      = dataset_name))
    }
    well  <- sum(df[[col]] >= 0.8, na.rm = TRUE)
    total <- sum(!is.na(df[[col]]))
    tibble(effect_size = es,
           well_powered = well,
           total_pairs  = total,
           proportion   = well / total,
           dataset      = dataset_name)
  })
}

message("Calculating power proportions (wo_pos_controls for DC-TAP only)")
all_power_props <- dplyr::bind_rows(
  calc_power_props(k562_dc_tap,    "K562",      prefer_wo = TRUE),
  calc_power_props(wtc11_dc_tap,   "WTC11",     prefer_wo = TRUE),
  calc_power_props(gasperini_results, "Gasperini", prefer_wo = FALSE),
  calc_power_props(morris,         "Morris",    prefer_wo = FALSE),
  calc_power_props(klann,          "Klann",     prefer_wo = FALSE),
  calc_power_props(xie,            "Xie",       prefer_wo = FALSE)
)


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

fill_order <- c(
  "K562 DC TAP Seq", "WTC11 DC TAP Seq",
  "Gasperini et al. 2019", "Morris et al. 2023",
  "Klann et al. 2021", "Xie et al. 2019"
)

dataset_colors <- c(
  "K562 DC TAP Seq"       = "#7A1E24",  # deep red
  "WTC11 DC TAP Seq"      = "#C00078",  # magenta
  "Gasperini et al. 2019" = "#333333",
  "Morris et al. 2023"    = "#666666",
  "Klann et al. 2021"     = "#999999",
  "Xie et al. 2019"       = "#BBBBBB"
)

# --- ensure clean, explicit x labels (no more "…") ---
# make effect_size a factor with fixed levels
all_power_props <- all_power_props %>%
  dplyr::mutate(effect_size = factor(effect_size, levels = c(10, 15, 20, 25, 50)))

# map levels -> labels explicitly (use plain hyphen to avoid font issues)
x_labels <- c("10" = "-10%", "15" = "-15%", "20" = "-20%", "25" = "-25%", "50" = "-50%")

# axis text size knobs
axis_text_size  <- 14
axis_title_size <- 15

message("Creating proportion plot (styled like panel b)")
percentage_plot <- ggplot(
  all_power_props,
  aes(x = effect_size, y = proportion * 100, fill = dataset)
) +
  geom_col(
    # wider bars + tighter dodge compresses space between groups
    position = position_dodge2(width = 0.9, padding = 0.03, preserve = "single"),
    width = 0.8,
    colour = NA
  ) +
  scale_fill_manual(values = dataset_colors, breaks = fill_order) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_discrete(labels = x_labels, drop = FALSE) +
  labs(x = "Effect Size", y = "Well-Powered DE-G Pairs (%)") +
  theme_classic(base_size = 11) +
  theme(
    legend.position  = "none",
    plot.title       = element_blank(),
    axis.text.x      = element_text(size = axis_text_size),
    axis.text.y      = element_text(size = axis_text_size),
    axis.title.x     = element_text(size = axis_title_size),
    axis.title.y     = element_text(size = axis_title_size),
    panel.grid       = element_blank(),
    axis.line.x      = element_line(linewidth = 0.6),
    axis.line.y      = element_line(linewidth = 0.6)
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving proportion plot")
ggsave(filename = snakemake@output$power_for_all_datasets_by_proportion, 
       plot = percentage_plot, width = 6, height = 3.5, device = "pdf")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)