# Script: compare_k562_effect_sizes_to_qPCR.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/compare_k562_effect_sizes_to_qPCR.rda"))
message("Saved Image")
stop("Manually Stopped Program after Saving Image")

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
  library(openxlsx)
})

message("Loading input files")
# Pilot DC TAP cell ranger outputs for main guide comparison
cell_ranger_outputs <- Read10X(snakemake@input$cell_ranger_outputs)
gene_matrix <- cell_ranger_outputs[[1]]
guide_matrix <- cell_ranger_outputs[[2]]

# qPCR results
qPCR_results <- read.xlsx(snakemake@input$qPCR_results)
extra_qPCRs <- read_tsv(snakemake@input$extra_qPCRs)

# Full singleton DC TAP results for extra guide comparison
k562_singleton_diffex_results <- readRDS(snakemake@input$k562_singleton_diffex_results)


### SCEPTRE SETUP =============================================================

# Creating the grna_target_data_frame - (grna_id, grna_target) 
guide_names <- rownames(guide_matrix)
grna_target_data_frame <- data.frame(
  grna_id = guide_names,
  grna_target = ifelse(1:length(guide_names) <= 94, "non-targeting", guide_names)
)

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
  filter(guides %in% testing_target_names) %>%
  filter(grepl("NC", guides) | sapply(seq_along(guides), function(i) grepl(genes[i], guides[i]))) %>% # Only keep `genes` - `guides` pairs where the `guides` contains "NC" or the `guides` contains the `genes` gene name
  mutate(type = case_when(
    grepl("TSS", guides) ~ "TSS",
    grepl("NC", guides) ~ "NC",
    grepl("e-", guides) | grepl("-e", guides) ~ "Enh"
  )) %>%
  left_join(discovery_results, by = c("genes" = "response_id", "guides" = "grna_target")) %>%
  mutate(pctChange = 2^log_2_fold_change - 1)


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
  dplyr::rename(guides = "Name", means = "qPCR", sds = `Stdev qPCR` ) %>%
  mutate(ics = NA, p_values = NA, n_nc_cells = NA) %>%
  mutate(type = case_when(
    grepl("tss", `Gene-Target`) ~ "TSS",
    grepl("enh", `Gene-Target`) ~ "enh"
  )) %>%
  mutate(genes = response_id_to_symbol_conversion[response_id]) %>%
  select(colnames(filt_qPCR_results))


### PLOT RESULTS ==============================================================

# Merge the two datasets together
merged_data <- rbind(
  filt_qPCR_results,
  joined_extra_qPCRs
) %>%
mutate(means = -(1 - means)) # Convert the "means" to "pctChange"

# Plot the correlation between "means" and "pctChange"
merged_data %>%
  ggplot(aes(x = means*100, y = pctChange*100, color = type)) +
  geom_point() +
  labs(
    title = "qPCR v. Sceptre Effect Sizes",
    y = "Sceptre Effect Size (%)",
    x = "qPCR Effect Size (%)"
  ) +
  theme_classic() +
  theme(
    aspect.ratio = 1
  )

### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_outputs_to_directory(sceptre_object = sceptre_object, directory = dirname(snakemake@output$results_run_discovery_analysis))


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)