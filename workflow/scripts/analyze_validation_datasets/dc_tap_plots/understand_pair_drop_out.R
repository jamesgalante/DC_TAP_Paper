# Script: understand_pair_drop_out.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/understand_pair_drop_out.rda"))
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
# Read in the files that are split
discovery_pairs <- lapply(snakemake@input$discovery_pairs, readRDS)
diffex_output <- lapply(snakemake@input$diffex_output, readRDS)
power_analysis_output <- lapply(snakemake@input$power_analysis_output, read_tsv)
create_encode_dataset_output <- lapply(snakemake@input$create_encode_dataset_output, read_tsv)
liftover_output <-  lapply(snakemake@input$liftover_output, readRDS)
filter_crispr_dataset_output <- lapply(snakemake@input$filter_crispr_dataset_output, read_tsv)

# Read in the files that are combined
create_ensemble_encode_output <- read_tsv(snakemake@input$create_ensemble_encode_output)
create_ensemble_epbenchmarking_output <- read_tsv(snakemake@input$create_ensemble_epbenchmarking_output)


### GET DROP OUT NUMS =========================================================

# Make sure to understand type of pair that dropped out

dim(discovery_pairs)
dim(diffex_output)
dim(power_analysis_output)
dim(create_encode_dataset_output)
dim(liftover_output)
dim(filter_crispr_dataset_output)
dim(create_ensemble_encode_output)
dim(create_ensemble_epbenchmarking_output)



# > dim(discovery_pairs)
# [1] 7964    2
# > dim(diffex_output)
# [1] 7964    8
# > dim(power_analysis_output)
# [1] 7919   24
# > dim(create_encode_dataset_output)
# [1] 7919   31
# > dim(liftover_output)
# [1] 7901   31
# > dim(filter_crispr_dataset_output)
# [1] 7901   31
# > dim(create_ensemble_encode_output)
# [1] 7525   30
# > dim(create_ensemble_epbenchmarking_output)
# [1] 7525   21


############# QUANTIFY EXACTLY HOW MANY PAIRS OF EACH TYPE DROP OUT AT EACH STEP


### Why do 45 pairs drop out in power analysis?
# ???
# Num tss v enh
num_tss_before <- sum(power_analysis_output$target_type == "TSSCtrl")
num_tss_after <- create_encode_dataset_output %>% filter(str_detect(ValidConnection, "TSS target")) %>% nrow()
print(num_tss_before - num_tss_after)
# None of these are TSS


### Why do 18 pairs drop out in liftover_output
# ???
# Num tss v enh
num_tss_before <- create_encode_dataset_output %>% filter(str_detect(ValidConnection, "TSS target")) %>% nrow()
num_tss_after <- liftover_output %>% filter(str_detect(ValidConnection, "TSS target")) %>% nrow()
print(num_tss_before - num_tss_after)
# None of these are TSS


### Why do 282 pairs drop out in filter_crispr_dataset_output* `filter_crispr_dataset.R`  !!!!! line 45 distance filter (*fixed*)
# distance filter is a hard cutoff at 1000bp so none of the TSS Ctrls are going to go through if they're targeting their own gene
# How many of these were TSS targeting versus enh
num_tss_before <- liftover_output %>% filter(str_detect(ValidConnection, "TSS target")) %>% nrow()
num_tss_after <- filter_crispr_dataset_output %>% filter(str_detect(ValidConnection, "TSS target")) %>% nrow()
print(num_tss_before - num_tss_after)
# So it seems like 150/282 pairs that dropped out at this step were TSS targets


### Why do 355 pairs drop out in create_ensemble_encode_output* `create_ensemble_dataset.R` !!!!! line 115 reduce gets rid of everything
# This is targets with overlapping regions that are paired with the same gene
# WHY are there overlapping targets -> did these exist before liftOver??
# How many tss duplicates are removed in the create_ensemble_encode_output step
num_tss_before <- filter_crispr_dataset_output %>% filter(str_detect(ValidConnection, "TSS target")) %>% nrow()
num_tss_after <- create_ensemble_encode_output %>% filter(str_detect(ValidConnection, "TSS target")) %>% nrow()
print(num_tss_before - num_tss_after)
# Apparently 149 / 355 pairs removed were genes tested against a TSS target
# Might have to check what's happening here -> like if fixing the filter crispr dataset thing doesn't help then it could be that enhancers overlapped with tss targets and they had more power or significance or something





### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)