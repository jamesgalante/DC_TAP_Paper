# Script: compare_replicates.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/compare_replicates_", snakemake@wildcards$sample, ".rda"))
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
  library(sceptre)
})

message("Loading input files")
# Import the main counts data matrix
dge <- readRDS(snakemake@input$dge)
# Import the perturbation status dataframe
perturb_status <- readRDS(snakemake@input$perturb_status)
# Import the discovery pairs
gene_gRNA_group_pairs <- readRDS(snakemake@input$gene_gRNA_group_pairs)
# Import the gRNA-target reference file
gRNA_groups_table <- readRDS(snakemake@input$gRNA_groups_table)
# Import the metadata file
metadata <- readRDS(snakemake@input$metadata)


### SPLIT THE BATCHES =========================================================

# first, get the unique batch ids as characters to avoid factor level issues
unique_batches <- as.character(unique(metadata$batch))
num_batches <- length(unique_batches)

# Randomly assign batches to two groups
set.seed(123)  # for reproducibility
batch_group1 <- sample(unique_batches, size = floor(num_batches/2))
batch_group2 <- setdiff(unique_batches, batch_group1)

# Print which batches are in each group
print("Batch Group 1:")
print(batch_group1)
print("Batch Group 2:")
print(batch_group2)

# Subset the metadata, dge, and perturb_status into the different batches
metadata1 <- metadata %>% filter(batch %in% batch_group1) %>% mutate(batch = as.factor(as.character(batch)))
metadata2 <- metadata %>% filter(batch %in% batch_group2) %>% mutate(batch = as.factor(as.character(batch)))

perturb_status1 <- perturb_status[,colnames(perturb_status) %in% rownames(metadata1)]
perturb_status2 <- perturb_status[,colnames(perturb_status) %in% rownames(metadata2)]

dge1 <- dge[,colnames(dge) %in% rownames(metadata1)]
dge2 <- dge[,colnames(dge) %in% rownames(metadata2)]


### SETUP DIFFEX ==============================================================

# Create sceptre_object1
sceptre_object1 <- import_data(
  response_matrix = dge1,
  grna_matrix = perturb_status1,
  grna_target_data_frame = gRNA_groups_table,
  moi = "high",
  extra_covariates = metadata1
)

# Create sceptre_object2
sceptre_object2 <- import_data(
  response_matrix = dge2,
  grna_matrix = perturb_status2,
  grna_target_data_frame = gRNA_groups_table,
  moi = "high",
  extra_covariates = metadata2
)

### ========================================

# Set analysis parameters1
sceptre_object1 <- set_analysis_parameters(
  sceptre_object = sceptre_object1,
  discovery_pairs = gene_gRNA_group_pairs,
  side = "both",
  grna_integration_strategy = "union",
)

# Set analysis parameters2
sceptre_object2 <- set_analysis_parameters(
  sceptre_object = sceptre_object2,
  discovery_pairs = gene_gRNA_group_pairs,
  side = "both",
  grna_integration_strategy = "union",
)


### RUN DIFFEX ================================================================

# Run differential expression1
sceptre_object1 <- sceptre_object1 %>% 
  assign_grnas(method = "thresholding", threshold = 5) %>%
  run_qc() %>%
  run_calibration_check() %>%
  run_discovery_analysis(parallel = FALSE)

# Run differential expression2
sceptre_object2 <- sceptre_object2 %>% 
  assign_grnas(method = "thresholding", threshold = 5) %>%
  run_qc() %>%
  run_calibration_check() %>%
  run_discovery_analysis(parallel = FALSE)


### SAVE OUTPUT ===============================================================

# Create plots directories inside each sample's output directory
plots_dir1 <- file.path(dirname(snakemake@output$discovery_results1), "plots")
plots_dir2 <- file.path(dirname(snakemake@output$discovery_results2), "plots")

# Check if directory exists before creating it
if (!dir.exists(plots_dir1)) {dir.create(plots_dir1, recursive = TRUE)}
if (!dir.exists(plots_dir2)) {dir.create(plots_dir2, recursive = TRUE)}

# Define ggsave wrapper that takes directory as parameter
save_plot <- function(plot, file_name, plots_dir) {
  file_path <- file.path(plots_dir, paste0(file_name, ".pdf"))
  ggsave(
    filename = file_path,
    plot = plot,
    device = "pdf",
    width = 5,
    height = 5
  )
}

# Save plots for replicate 1
message("Saving plots for batch group 1")
save_plot(plot_grna_count_distributions(sceptre_object1), "plot_grna_count_distributions", plots_dir1)
save_plot(plot_assign_grnas(sceptre_object1), "plot_assign_grnas", plots_dir1)
save_plot(plot_run_calibration_check(sceptre_object1), "plot_run_calibration_check", plots_dir1)
save_plot(plot_run_discovery_analysis(sceptre_object1), "plot_run_discovery_analysis", plots_dir1)
save_plot(plot_covariates(sceptre_object1), "plot_covariates", plots_dir1)
save_plot(plot_run_qc(sceptre_object1), "plot_run_qc", plots_dir1)

# Save plots for replicate 2
message("Saving plots for batch group 2")
save_plot(plot_grna_count_distributions(sceptre_object2), "plot_grna_count_distributions", plots_dir2)
save_plot(plot_assign_grnas(sceptre_object2), "plot_assign_grnas", plots_dir2)
save_plot(plot_run_calibration_check(sceptre_object2), "plot_run_calibration_check", plots_dir2)
save_plot(plot_run_discovery_analysis(sceptre_object2), "plot_run_discovery_analysis", plots_dir2)
save_plot(plot_covariates(sceptre_object2), "plot_covariates", plots_dir2)
save_plot(plot_run_qc(sceptre_object2), "plot_run_qc", plots_dir2)

# Use directories from the Snakemake output paths
batch1_dir <- dirname(snakemake@output$discovery_results1)
batch2_dir <- dirname(snakemake@output$discovery_results2)

# Write all outputs to directories
message("Saving output files for batch group 1")
write_outputs_to_directory(sceptre_object1, batch1_dir)

message("Saving output files for batch group 2")
write_outputs_to_directory(sceptre_object2, batch2_dir)

# Save the final sceptre objects and discovery results
message("Saving final sceptre objects and discovery results")
saveRDS(sceptre_object1, snakemake@output$final_sceptre_object1)
saveRDS(sceptre_object2, snakemake@output$final_sceptre_object2)
saveRDS(get_result(sceptre_object1, analysis = "run_discovery_analysis"), snakemake@output$discovery_results1)
saveRDS(get_result(sceptre_object2, analysis = "run_discovery_analysis"), snakemake@output$discovery_results2)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)