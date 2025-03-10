# Script: process_cell_ranger_outputs_K562_DC_TAP_Seq.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/process_cell_ranger_outputs_K562_DC_TAP_Seq.rda"))
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
  library(Seurat)
  library(data.table)
})

message("Loading input files")
cell_ranger_output <- Read10X(snakemake@input$cell_ranger_directory, gene.column = 1)


### GET INDIVIDUAL FILES ======================================================

gene_matrix <- cell_ranger_output$`Gene Expression`
guide_matrix <- cell_ranger_output$`CRISPR Guide Capture`


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(gene_matrix, snakemake@output$gene_matrix)
saveRDS(guide_matrix, snakemake@output$guide_matrix)
message("Files successfully saved")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)