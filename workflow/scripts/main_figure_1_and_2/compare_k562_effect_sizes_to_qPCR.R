# Script: compare_k562_effect_sizes_to_qPCR.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/compare_k562_effect_sizes_to_qPCR.rda"))
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
  library(Seurat)
  library(sceptre)
  library(ggrepel)
})

message("Loading input files")
# Pilot DC TAP cell ranger outputs for main guide comparison
cell_ranger_outputs <- Read10X(snakemake@input$cell_ranger_outputs)
gene_matrix <- cell_ranger_outputs[[1]]
guide_matrix <- cell_ranger_outputs[[2]]

# qPCR results
qPCR_results <- read_tsv(snakemake@input$qPCR_results)
extra_qPCRs <- read_tsv(snakemake@input$extra_qPCRs)

# Full singleton DC TAP results for extra guide comparison
k562_singleton_diffex_results <- readRDS(snakemake@input$k562_singleton_diffex_results)


### SCEPTRE SETUP =============================================================

# Creating the grna_target_data_frame - (grna_id, grna_target) 
guide_names <- rownames(guide_matrix)
grna_target_data_frame <- data.frame(
  grna_id = guide_names,
  grna_target = ifelse(1:length(guide_names) <= 98, "non-targeting", guide_names)
)

# These guides are not present in the fastq files
guide_names_not_in_screen <- c("safe_1", "safe_2", "safe_6", "safe_17", "safe_19", "safe_22", "safe_26", "safe_41", "nontarget_4", "nontarget_32", "nontarget_33", "nontarget_34", "nontarget_35")
grna_target_data_frame <- grna_target_data_frame %>% filter(!grna_id %in% guide_names_not_in_screen)

# Creating the discovery_pairs dataframe - all by all - (grna_target, response_id)
testing_target_names <- grna_target_data_frame %>% filter(grna_target != "non-targeting") %>% pull(grna_id)
response_id_names <- rownames(gene_matrix)
discovery_pairs <- expand.grid(
  grna_target = testing_target_names,
  response_id = response_id_names,
  stringsAsFactors = FALSE
)


### CREATE SCEPTRE OBJECTS ====================================================

# Create sceptre_object
sceptre_object <- import_data(
  response_matrix = gene_matrix,
  grna_matrix = guide_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "low"
)

# Set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  side = "both",
  grna_integration_strategy = "union",
)

print(sceptre_object)


### RUN SCEPTRE DIFFEX ========================================================

sceptre_object <- sceptre_object %>% 
  assign_grnas(method = "mixture") %>%
  run_qc() %>%
  run_calibration_check() %>%
  run_discovery_analysis(parallel = FALSE)


### FORMAT SCEPTRE WITH qPCR RESULTS ==========================================

# Let's get the discovery results to compare to the qPCR
discovery_results <- get_result(
  sceptre_object = sceptre_object,
  analysis = "run_discovery_analysis"
)

# Here we're comparing the results of the guides from a pilot screen with qPCR results
# Let's look at all guides that we tested that are in discovery pairs
filt_qPCR_results <- qPCR_results %>% 
  mutate(type = "TSS") %>%
  mutate(qPCR_EffectSize = expr_remaining - 1) %>%
  left_join(discovery_results, by = c("gene" = "response_id", "guide" = "grna_target")) %>%
  mutate(pctChange = 2^log_2_fold_change - 1) %>%
  mutate(Source = "Original") %>%
  select(-expr_remaining) %>%
  mutate(plotting_name = guide)


### FORMAT EXTRA RESULTS ======================================================

response_id_to_symbol_conversion <- c(
  "ENSG00000106367" = "AP1S1",
  "ENSG00000102007" = "PLP2",
  "ENSG00000158578" = "ALAS2",
  "ENSG00000235169" = "SMIM1",
  "ENSG00000118640" = "VAMP8",
  "ENSG00000147689" = "FAM83A",
  "ENSG00000127445" = "PIN1",
  "ENSG00000160932" = "LY6E"
)

# Here we're comparing the results of some guides in the full screen that were additionally tested with qPCR
joined_extra_qPCRs <- extra_qPCRs %>%
  left_join(k562_singleton_diffex_results %>% group_by(grna_id) %>% dplyr::slice(1), by = c("Name" = "grna_id")) %>% # I checked each response_id, and they are the expected gene_symbol
  mutate(pctChange = 2^log_2_fold_change - 1) %>%
  dplyr::rename(guide = "Name", sem = `SE_qPCR`) %>%
  mutate(type = case_when(
    grepl("tss", `Gene-Target`) ~ "TSS",
    grepl("enh", `Gene-Target`) ~ "Enh"
  )) %>%
  mutate(gene = response_id_to_symbol_conversion[response_id]) %>%
  mutate(Source = "Extra") %>%
  mutate(qPCR_EffectSize = -qPCR) %>%
  mutate(plotting_name = gene) %>%
  select(-c("response_id", "grna_target")) %>%
  select(colnames(filt_qPCR_results))


### INITIAL PLOTTING RESULTS ==================================================

# Merge the two datasets together
merged_data <- rbind(
  filt_qPCR_results,
  joined_extra_qPCRs
) 

comparison_of_all_by_type <- merged_data %>%
  ggplot(aes(x = qPCR_EffectSize*100, y = pctChange*100, label = plotting_name, color = type)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3.8, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(
    x = "qPCR Effect Size (%)",
    y = "SCEPTRE Effect Size (%)",
    color = "Target Type"
  ) +
  theme_minimal(base_size = 15) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        aspect.ratio = 1) +
  coord_cartesian(xlim = c(-100, 0), ylim = c(-100, 0))

# Save the plot
ggsave(plot = comparison_of_all_by_type,
       filename = snakemake@output$comparison_of_all_by_type,
       device = "pdf", height = 6, width = 6)


### FIGURES FOR PAPER =========================================================

# Subsetting the data for two plots
positive_controls <- merged_data %>%
  filter(gene %in% c("HDAC6", "GATA1", "PLP2", "FAM83A"))
other_tests <- merged_data %>%
  filter(!gene %in% positive_controls$gene)

# Plotting the positive controls
positive_controls[positive_controls$gene == "PLP2", "plotting_name"] <- "PLP2-TSS"
positive_controls[positive_controls$gene == "FAM83A", "plotting_name"] <- "FMA83A-TSS"

# Plot the positive controls comparison
comparison_positive_controls <- positive_controls %>%
  ggplot(aes(x = qPCR_EffectSize*100, y = pctChange*100, label = plotting_name)) +
  geom_point(size = 3, color = "red") +
  geom_text_repel(size = 4.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(
    x = "qPCR Effect Size (%)",
    y = "SCEPTRE Effect Size (%)"
  ) +
  theme_minimal(base_size = 15) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        aspect.ratio = 1) +
  coord_cartesian(xlim = c(-100, 0), ylim = c(-100, 0))

# Save the plot
ggsave(plot = comparison_positive_controls,
       filename = snakemake@output$comparison_positive_controls,
       device = "pdf", height = 5, width = 5)

# Plot the other comparisons
comparison_others <- other_tests %>%
  ggplot(aes(x = qPCR_EffectSize*100, y = pctChange*100, label = plotting_name, color = type)) +
  geom_point(size = 3) +
  geom_text_repel(size = 4.5, color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  labs(
    x = "qPCR Effect Size (%)",
    y = "SCEPTRE Effect Size (%)",
    color = "Target Type"
  ) +
  theme_minimal(base_size = 15) + 
  theme(panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        aspect.ratio = 1) +
  coord_cartesian(xlim = c(-100, 0), ylim = c(-100, 0))

# Save the plot
ggsave(plot = comparison_others,
       filename = snakemake@output$comparison_others,
       device = "pdf", height = 6, width = 6)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(grna_target_data_frame, snakemake@output$grna_target_data_frame_filtered_guides)
write_tsv(merged_data, snakemake@output$formatted_qPCR_Sceptre_table)
write_outputs_to_directory(sceptre_object = sceptre_object, directory = dirname(snakemake@output$results_run_discovery_analysis))


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)