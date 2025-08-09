# Script: add_results_wo_positive_controls.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/add_results_wo_positive_controls.rda"))
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
element_gene_pairs <- read_tsv(snakemake@input$results_with_element_gene_pair_categories_modified)

# Read in power simulation data for ALL effect sizes
effect_sizes <- c(2, 3, 5, 10, 15, 20, 25, 50)

# Load K562 power simulations for all effect sizes using map_dfr
k562_power_simulation <- map_dfr(seq_along(effect_sizes), function(i) {
  fread(snakemake@input$combined_power_analysis_output_K562[[i]]) %>%
    as_tibble() %>%
    mutate(cell_type = "K562", effect_size = effect_sizes[i])
})

# Load WTC11 power simulations for all effect sizes using map_dfr
wtc11_power_simulation <- map_dfr(seq_along(effect_sizes), function(i) {
  fread(snakemake@input$combined_power_analysis_output_WTC11[[i]]) %>%
    as_tibble() %>%
    mutate(cell_type = "WTC11", effect_size = effect_sizes[i])
})

# Combine all power simulations
combined_power_simulation <- rbind(k562_power_simulation, wtc11_power_simulation)


### CUSTOM FDR CORRECTION WITHOUT POSITIVE CONTROLS ==========================

message("Applying custom FDR correction excluding positive controls")

# Apply FDR correction by cell type, excluding positive controls
final_pairs <- element_gene_pairs %>%
  group_by(cell_type) %>%
  mutate(
    # Create a flag for pairs to include in FDR correction
    include_in_fdr = (selfPromoter == FALSE & Positive_Control_DistalElement_Gene == FALSE & Positive_Control_selfPromoter == FALSE),
    
    # Apply Benjamini-Hochberg correction only to the subset
    sceptre_adj_p_value_wo_pos_controls = ifelse(
      include_in_fdr & !is.na(sceptre_p_value),
      p.adjust(ifelse(include_in_fdr & !is.na(sceptre_p_value), sceptre_p_value, NA), method = "BH"),
      NA_real_
    ),
    
    # Call significance at FDR 0.1 (equivalent to alpha 0.05)
    significant_wo_pos_controls = case_when(
      include_in_fdr & !is.na(sceptre_adj_p_value_wo_pos_controls) ~ sceptre_adj_p_value_wo_pos_controls <= 0.1,
      TRUE ~ NA
    )
  ) %>%
  ungroup()


### RECALCULATE POWER FOR WO_POS_CONTROLS FOR ALL EFFECT SIZES ===============

message("Recalculating power for wo_pos_controls significance threshold for all effect sizes")

# Get max nominal p-value from wo_pos_controls significant pairs for each cell type
max_nom_p_val_K562 <- final_pairs %>% 
  filter(significant_wo_pos_controls == TRUE, cell_type == "K562") %>% 
  pull(sceptre_p_value) %>% 
  max(na.rm = TRUE)

max_nom_p_val_WTC11 <- final_pairs %>% 
  filter(significant_wo_pos_controls == TRUE, cell_type == "WTC11") %>% 
  pull(sceptre_p_value) %>% 
  max(na.rm = TRUE)

message(paste("Max nominal p-value for K562 wo_pos_controls:", max_nom_p_val_K562))
message(paste("Max nominal p-value for WTC11 wo_pos_controls:", max_nom_p_val_WTC11))

# Function to calculate power for a given effect size
calculate_power_for_effect_size <- function(df, effect_size_val, max_p_val, cell_type_val) {
  df %>%
    filter(cell_type == cell_type_val, effect_size == effect_size_val, !is.na(log_2_fold_change)) %>%
    group_by(grna_target, response_id) %>%
    summarize(
      !!paste0("power_at_effect_size_", effect_size_val, "_wo_pos_controls") := 
        mean(p_value < max_p_val & log_2_fold_change < 0),
      .groups = "drop"
    )
}

# Calculate power for all effect sizes for K562
power_wo_pos_controls_K562 <- map(effect_sizes, function(es) {
  calculate_power_for_effect_size(combined_power_simulation, es, max_nom_p_val_K562, "K562")
}) %>%
  reduce(left_join, by = c("grna_target", "response_id")) %>%
  mutate(cell_type = "K562")

# Calculate power for all effect sizes for WTC11
power_wo_pos_controls_WTC11 <- map(effect_sizes, function(es) {
  calculate_power_for_effect_size(combined_power_simulation, es, max_nom_p_val_WTC11, "WTC11")
}) %>%
  reduce(left_join, by = c("grna_target", "response_id")) %>%
  mutate(cell_type = "WTC11")

# Combine power results
power_wo_pos_controls <- rbind(power_wo_pos_controls_K562, power_wo_pos_controls_WTC11)

# Calculate mean perturbation cells across all effect sizes for both cell types
mean_pert_cells_summary <- combined_power_simulation %>%
  group_by(grna_target, response_id, cell_type) %>%
  summarize(
    mean_sim_pert_cells = mean(num_pert_cells, na.rm = TRUE),
    .groups = "drop"
  )

# Merge both power calculations and mean_pert_cells back to final_pairs
final_pairs <- final_pairs %>%
  left_join(
    power_wo_pos_controls,
    by = c("design_file_target_name" = "grna_target", "gene_id" = "response_id", "cell_type")
  ) %>%
  left_join(
    mean_pert_cells_summary,
    by = c("design_file_target_name" = "grna_target", "gene_id" = "response_id", "cell_type")
  )


### MAKE POWER_WO_POS_CONTROLS NA FOR ALL EFFECT SIZES =======================

# Power at all effect sizes without positive controls should be NA where significance wasn't calculated
# Because Positive controls weren't included in the FDR correction
for(es in effect_sizes) {
  power_col <- paste0("power_at_effect_size_", es, "_wo_pos_controls")
  final_pairs[[power_col]] <- ifelse(is.na(final_pairs$significant_wo_pos_controls), 
                                     NA, 
                                     final_pairs[[power_col]])
}

# Set mean_sim_pert_cells of pairs which power isn't calculated for to NA
final_pairs <- final_pairs %>% 
  mutate(mean_sim_pert_cells = ifelse(is.na(significant_wo_pos_controls), NA, mean_sim_pert_cells))


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(final_pairs, snakemake@output$results_wo_pos_controls)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)