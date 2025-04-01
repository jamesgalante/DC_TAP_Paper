# Script: understand_pair_drop_out.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/understand_pair_drop_out.rda"))
message("Saved Image")
stop("Manually Stopped Program after Saving Image")

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
# Read in the files that are split
discovery_pairs <- lapply(snakemake@input$discovery_pairs, readRDS)
diffex_output <- lapply(snakemake@input$diffex_output, readRDS)
power_analysis_output <- lapply(snakemake@input$power_analysis_output, read_tsv)
create_encode_dataset_output <- lapply(snakemake@input$create_encode_dataset_output, read_tsv)
liftover_output <-  lapply(snakemake@input$liftover_output, read_tsv)
filter_crispr_dataset_output <- lapply(snakemake@input$filter_crispr_dataset_output, read_tsv)

# Read in the files that are combined
create_ensemble_encode_output <- read_tsv(snakemake@input$create_ensemble_encode_output)
create_ensemble_epbenchmarking_output <- read_tsv(snakemake@input$create_ensemble_epbenchmarking_output)


### GET DROP OUT NUMS =========================================================

# Make sure to understand type of pair that dropped out

dim(discovery_pairs[[1]])
dim(diffex_output[[1]])
dim(power_analysis_output[[1]])
dim(create_encode_dataset_output[[1]])
dim(liftover_output[[1]])
dim(filter_crispr_dataset_output[[1]])
dim(create_ensemble_encode_output %>% filter(Dataset == "K562_DC_TAP_Seq"))
dim(create_ensemble_epbenchmarking_output %>% filter(Dataset == "K562_DC_TAP_Seq"))

# > dim(discovery_pairs[[1]])
# [1] 7556    2
# > dim(diffex_output[[1]])
# [1] 7556    8
# > dim(power_analysis_output[[1]]) 
# [1] 7516   24
# > dim(create_encode_dataset_output[[1]])
# [1] 7516   31
# > dim(liftover_output[[1]])
# [1] 7498   31
# > dim(filter_crispr_dataset_output[[1]])
# [1] 7498   31
# > dim(create_ensemble_encode_output %>% filter(Dataset == "K562_DC_TAP_Seq"))
# [1] 7498   30
# > dim(create_ensemble_epbenchmarking_output %>% filter(Dataset == "K562_DC_TAP_Seq"))
# [1] 7498   21

dim(discovery_pairs[[2]])
dim(diffex_output[[2]])
dim(power_analysis_output[[2]])
dim(create_encode_dataset_output[[2]])
dim(liftover_output[[2]])
dim(filter_crispr_dataset_output[[2]])
dim(create_ensemble_encode_output %>% filter(Dataset == "WTC11_DC_TAP_Seq"))
dim(create_ensemble_epbenchmarking_output %>% filter(Dataset == "WTC11_DC_TAP_Seq"))

# > dim(discovery_pairs[[2]])
# [1] 6843    2
# > dim(diffex_output[[2]])
# [1] 6843    8
# > dim(power_analysis_output[[2]])
# [1] 6581   24
# > dim(create_encode_dataset_output[[2]])
# [1] 6581   31
# > dim(liftover_output[[2]])
# [1] 6581   31
# > dim(filter_crispr_dataset_output[[2]])
# [1] 6581   31
# > dim(create_ensemble_encode_output %>% filter(Dataset == "WTC11_DC_TAP_Seq"))
# [1] 6581   30
# > dim(create_ensemble_epbenchmarking_output %>% filter(Dataset == "WTC11_DC_TAP_Seq"))
# [1] 6581   21


############# QUANTIFY EXACTLY HOW MANY PAIRS OF EACH TYPE DROP OUT AT EACH STEP - K562

# 40 drop out after power analysis - 
missing_pairs <- anti_join(
  diffex_output[[1]] %>% dplyr::select(grna_target, response_id), 
  power_analysis_output[[1]] %>% dplyr::select(perturbation, gene),
  by = c("grna_target" = "perturbation", "response_id" = "gene")
)
diffex_output[[1]] %>% filter(grna_target %in% missing_pairs$grna_target, response_id %in% missing_pairs$response_id)
# None of these pairs pass QC by sceptre - in fact these are the only pairs that don't pass qc

# After liftover, 18 drop out - all pairs with `chr8:145537346-145537700` liftOver error message: "Partially deleted in new: Sequence insufficiently intersects one chain"


############# QUANTIFY EXACTLY HOW MANY PAIRS OF EACH TYPE DROP OUT AT EACH STEP - WTC11

# For WTC11, pairs only dropped out after sceptre, so these are from not passing Sceptre QCs
sum(!diffex_output[[2]]$pass_qc) ==  nrow(diffex_output[[2]]) - nrow(power_analysis_output[[2]]) 






### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)