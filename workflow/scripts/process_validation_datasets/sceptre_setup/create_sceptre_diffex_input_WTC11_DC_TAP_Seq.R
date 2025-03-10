# Script: # Script: create_sceptre_diffex_input_WTC11_DC_TAP_Seq.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_sceptre_diffex_input_WTC11_DC_TAP_Seq.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

# Download Sceptre
# library(devtools)
# message("Installing Sceptre")
# devtools::install_github("katsevich-lab/sceptre")
# message("Sceptre Installation Complete")

message("Loading in packages")
suppressPackageStartupMessages({
  library(sceptre)
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  source(file.path(snakemake@scriptdir, "gene_target_pairing_functions.R"))
})

message("Loading input files")
# Import the main counts data matrix
dge <- readRDS(snakemake@input$dge)
# Import the perturbation status dataframe
perturb_status <- readRDS(snakemake@input$perturb_status)
# Import the annotation file
annot <- import(snakemake@input$annot)
# Import the guide_targets file
guide_targets <- read_tsv(snakemake@input$guide_targets)


### CREATE GENE_GRNA_GROUP_PAIRS ==============================================

# Create the file
target_search_results <- find_genes_near_targets(guide_targets %>% filter(!is.na(target_chr)), annot, rownames(dge))
gene_gRNA_group_pairs <- target_search_results[[1]]
errors <- target_search_results[[2]]

# Modify for Sceptre Input specifications
gene_gRNA_group_pairs <- gene_gRNA_group_pairs %>%
  select(grna_group, response_id) %>%
  dplyr::rename(grna_target = grna_group)


### CHECK THE ERRORED GENE NAMES ==============================================

# Some of the genes were negative controls, so we don't expect them to be within 2Mb of any target
# Check which of the genes that don't get paired up with any target are/aren't a negative control

negative_control_genes <- as.vector(snakemake@params$negative_control_genes)
genes_paired_with_no_targets <- errors[[2]]$gene_name

message("Genes, which aren't negative control genes, that weren't paired up with any target within 2Mb")
print(setdiff(genes_paired_with_no_targets, negative_control_genes))


### CREATE GRNA_GROUPS_TABLE ================================================== 

# Create the grna_groups_table which is just the guide name matched to the target name
gRNA_groups_table <- guide_targets %>% 
  select(name, target_name) %>%
  dplyr::rename(grna_id = name, grna_target = target_name) %>%
  mutate(
    grna_target = case_when(
      grna_target == "safe_targeting" ~ "non-targeting",
      grna_target == "negative_control" ~ "non-targeting",
      TRUE ~ grna_target
    )
  )

# Subset for guides that are also in the perturb_status guide counts matrix - as not all are
message("Guides in the guide_targets file that aren't in perturb_status matrix: ")
print(gRNA_groups_table %>% filter(!grna_id %in% rownames(perturb_status)))
gRNA_groups_table <- gRNA_groups_table %>% 
  filter(grna_id %in% rownames(perturb_status))


### CREATE THE METADATA FILE ==================================================

# Create the metadata indicating the batch for each cell barcode
cell_barcode <- colnames(dge)
cell_batches <- as.factor(str_extract(cell_barcode, "\\d+$"))

# Create the metadata dataframe
metadata <- data.frame(row.names = cell_barcode, batch = cell_batches)


### CREATE SCEPTRE OBJECT =====================================================

# Create sceptre_object
sceptre_object <- import_data(
  response_matrix = dge,
  grna_matrix = perturb_status,
  grna_target_data_frame = gRNA_groups_table,
  moi = "high",
  extra_covariates = metadata
)

# Set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = gene_gRNA_group_pairs,
  side = "both",
  grna_integration_strategy = "union",
)

print(sceptre_object)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(gene_gRNA_group_pairs, snakemake@output$gene_gRNA_group_pairs)
saveRDS(gRNA_groups_table, snakemake@output$gRNA_groups_table)
saveRDS(metadata, snakemake@output$metadata)
saveRDS(sceptre_object, snakemake@output$sceptre_diffex_input)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)