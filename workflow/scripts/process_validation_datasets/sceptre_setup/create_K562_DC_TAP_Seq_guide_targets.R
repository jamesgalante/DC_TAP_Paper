# Script: create_K562_DC_TAP_Seq_guide_targets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_K562_DC_TAP_Seq_guide_targets.rda"))
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
  library(GenomicRanges)
})

message("Loading input files")
guide_design_file <- read_tsv(snakemake@input$guide_design_file)



################################# TEMP UNTIL ADD SCREEN DESIGN TO PIPELINE
positive_controls <- read_tsv(snakemake@input$positive_controls, col_names = FALSE)
old_guide_targets <- read_tsv(snakemake@input$old_guide_targets)
old_pre_merged_guide_targets <- read_tsv(snakemake@input$old_pre_merged_guide_targets)
################################# TEMP UNTIL ADD SCREEN DESIGN TO PIPELINE



### GET THE TARGET TYPES ======================================================






################################# TEMP UNTIL ADD SCREEN DESIGN TO PIPELINE

# Parse the positive_controls for distal enhancers
distal_element_controls <- positive_controls %>%
  filter(startsWith(X4, "chr")) %>% # Select all the DEs
  filter(X4 %in% guide_design_file$guideSet) %>% # Remove the one DE that is extra and not in the screen
  pull(X4)

# Which guides target "tss_pos"
tss_pos <- old_guide_targets %>%
  filter(target_type == "tss_pos") %>%
  pull(name)

tss_random <- old_guide_targets %>%
  filter(target_type == "tss_random") %>%
  pull(name)

# Add target_type to the guide_design file
guide_design_file <- guide_design_file %>%
  mutate(
    target_type = case_when(
      name %in% tss_pos ~ "tss_pos",
      name %in% tss_random ~ "tss_random",
      guideSet %in% distal_element_controls ~ "DE",
      startsWith(guideSet, "chr") ~ "enh",
      TRUE ~ guideSet # Should only be safe and non targeting
    )
  )
################################# TEMP UNTIL ADD SCREEN DESIGN TO PIPELINE





# Summarize screen stats by target_type
design_summary <- guide_design_file %>%
  group_by(target_type) %>%
  summarize(
    n_gRNAs = n_distinct(name),
    n_guideSets = n_distinct(guideSet)
  )
message("Summary stats for the screen design by target type: ")
print(design_summary)


### CREATE THE TARGETS RANGES ================================================

# For all "enh" this should just be the guideSet name, split up into target_chr, target_start, target_end
# For all non-targeting guides, this should be "NA"
# For all TSS of the same type, get the range of all guides targeting the same guideSet

safe_min <- function(x) {
  if(all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  if(all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
}

processed_guide_design_file <- guide_design_file %>%
  mutate(
    # Initialize columns with NA (explicitly as proper types)
    target_chr = NA_character_,
    target_start = NA_character_,  # temporarily as character for extraction
    target_end = NA_character_
  ) %>%
  # Extract chr, start, and end for "enh" guides
  mutate(
    target_chr = case_when(
      target_type %in% c("DE", "enh") ~ sub(":.*", "", guideSet),
      target_type %in% c("negative_control", "safe_targeting") ~ NA_character_,
      TRUE ~ target_chr
    ),
    target_start = case_when(
      target_type %in% c("DE", "enh") ~ sub(".*:(\\d+)-.*", "\\1", guideSet),
      target_type %in% c("negative_control", "safe_targeting") ~ NA_character_,
      TRUE ~ target_start
    ),
    target_end = case_when(
      target_type %in% c("DE", "enh") ~ sub(".*-(\\d+)", "\\1", guideSet),
      target_type %in% c("negative_control", "safe_targeting") ~ NA_character_,
      TRUE ~ target_end
    )
  ) %>%

  # Convert extracted start and end values to numeric
  mutate(
    target_start = as.numeric(target_start),
    target_end = as.numeric(target_end)
  ) %>%

  # Compute genomic ranges for DE and TSS guides using safe_min/max
  group_by(guideSet) %>%
  mutate(
    target_start = case_when(
      target_type %in% c("tss_random", "tss_pos") ~ safe_min(start),
      TRUE ~ target_start
    ),
    target_end = case_when(
      target_type %in% c("tss_random", "tss_pos") ~ safe_max(end),
      TRUE ~ target_end
    ),
    target_chr = case_when(
      target_type %in% c("tss_random", "tss_pos") ~ if(length(unique(chr[!is.na(chr)])) == 0) NA_character_ else unique(chr[!is.na(chr)]),
      TRUE ~ target_chr
    )
  ) %>%
  ungroup()


### CREATE THE GUIDE TARGETS FILE ============================================

# Create the guide_targets file
guide_targets <- data.frame(
  chr = processed_guide_design_file$chr,
  start = processed_guide_design_file$start,
  end = processed_guide_design_file$end,
  name = processed_guide_design_file$name,
  strand = processed_guide_design_file$strand,
  spacer = processed_guide_design_file$GuideSequenceWithPAM,
  target_chr = processed_guide_design_file$target_chr,
  target_start = processed_guide_design_file$target_start,
  target_end = processed_guide_design_file$target_end,
  target_name = processed_guide_design_file$guideSet,
  target_strand = ".",
  target_type = processed_guide_design_file$target_type
)


### MERGING ANY OVERLAPPING TARGETS ===========================================

# There are some targets which are overlapping, so we should merge these targets before running differential expression
# Let's first make sure that target_region start and end are correct for everything (specifically for TSS)
  # Because currently I define a TSS region as something that is limited to where the guides are
  # The TSS region should be 500 bp centered on the summit of the dnase peak representing the TSS


# In the following, we're getting the TSS target region from the old guide targets file - making sure this is 500 bp
################################# TEMP UNTIL ADD SCREEN DESIGN TO PIPELINE
# Extract the desired name order
name_order <- guide_targets$name
# Reorder the rows of old_pre_merged_guide_targets
pre_merged_ordered <- old_pre_merged_guide_targets[match(name_order, old_pre_merged_guide_targets$name), ]

# Add the target_start and target_end from the old_pre_merged_guide_targets file to the official file for the positive controls only
pre_merged_ordered %>% filter(target_type %in% c("tss_pos", "tss_random"))
guide_targets %>% filter(target_type %in% c("tss_pos", "tss_random"))

# For these pairs, modify the guide_targets file target_start and target_end to match that of the pre merged ordered file
# Filter down to just the positive/random TSS entries in both dataframes
to_update <- pre_merged_ordered %>%
  filter(target_type %in% c("tss_pos", "tss_random")) %>%
  select(name, target_start, target_end)

# Join updated values into guide_targets and overwrite target_start/end
guide_targets_updated <- guide_targets %>%
  left_join(to_update, by = "name", suffix = c("", ".new")) %>%
  mutate(
    target_start = if_else(target_type %in% c("tss_pos", "tss_random"), target_start.new, target_start),
    target_end   = if_else(target_type %in% c("tss_pos", "tss_random"), target_end.new, target_end)
  ) %>%
  select(-target_start.new, -target_end.new)
################################# TEMP UNTIL ADD SCREEN DESIGN TO PIPELINE


# Now that we have the width of each target confirmed, we can overlap the "target_chr", "target_start", and "target_end" of each target
# We can then get a list of all target_names that need to be merged - along with the target types

# First get a list of all target names, their regions, and the target type, and remove any non targeting
all_targets <- guide_targets_updated %>% 
  select(target_name, target_chr, target_start, target_end, target_type) %>%
  filter(!duplicated(target_name)) %>%
  filter(target_type %in% c("enh", "tss_pos", "tss_random", "DE"))


# Step 1: Create GRanges from all_targets
gr <- GRanges(
  seqnames = all_targets$target_chr,
  ranges = IRanges(start = all_targets$target_start, end = all_targets$target_end),
  target_name = all_targets$target_name,
  target_type = all_targets$target_type
)

# Step 2: Reduce overlapping regions into merged "blocks"
reduced_gr <- reduce(gr)

# Step 3: Map original targets to their reduced blocks
hits <- findOverlaps(gr, reduced_gr)

# Step 4: Build mapping of each merged block to all targets that fall into it
merged_map <- as.data.frame(hits) %>%
  mutate(
    original_target_name = mcols(gr)$target_name[queryHits],
    original_target_type = mcols(gr)$target_type[queryHits],
    original_chr = as.character(seqnames(gr)[queryHits]),
    original_start = start(gr)[queryHits],
    original_end = end(gr)[queryHits],
    merged_chr = as.character(seqnames(reduced_gr)[subjectHits]),
    merged_start = start(reduced_gr)[subjectHits],
    merged_end = end(reduced_gr)[subjectHits]
  )


# Step 5: Assign merged name/type logic
merged_targets <- merged_map %>%
  group_by(subjectHits) %>%
  summarise(
    target_chr = dplyr::first(merged_chr),
    target_start = dplyr::first(merged_start),
    target_end = dplyr::first(merged_end),
    original_names = paste(unique(original_target_name), collapse = "_AND_"),
    types = paste(unique(original_target_type), collapse = ","),
    target_type = if (length(unique(original_target_type)) == 1) {
      unique(original_target_type)
    } else {
      "tss_random"
    },
    target_name = if (all(unique(original_target_type) == "enh")) {
      paste0(dplyr::first(merged_chr), ":", dplyr::first(merged_start), "-", dplyr::first(merged_end))
    } else {
      original_names
    },
    .groups = "drop"
  ) %>%
  select(target_name, target_chr, target_start, target_end, target_type)

# Verify that the number of rows in merged_targets matches the number of unique subjectHits
# This ensures our row_number() approach is valid
message("Validating merged target mapping...")
if(length(unique(merged_map$subjectHits)) != nrow(merged_targets)) {
  warning("Number of unique subject hits doesn't match number of merged targets. Check the mapping.")
}

# Create target mapping 
target_mapping <- merged_map %>%
  select(original_target_name, subjectHits) %>%
  # Join with merged_targets to get the merged target information
  left_join(
    merged_targets %>% mutate(subjectHits = row_number()),
    by = "subjectHits"
  ) %>%
  select(original_target_name, 
         merged_name = target_name,
         merged_chr = target_chr,
         merged_start = target_start,
         merged_end = target_end,
         merged_type = target_type)

# Apply the mapping to update guide_targets_updated
guide_targets_final <- guide_targets_updated %>%
  left_join(
    target_mapping,
    by = c("target_name" = "original_target_name")
  ) %>%
  mutate(
    # Update target information with merged information
    target_name = ifelse(!is.na(merged_name), merged_name, target_name),
    target_chr = ifelse(!is.na(merged_chr), merged_chr, target_chr),
    target_start = ifelse(!is.na(merged_start), merged_start, target_start),
    target_end = ifelse(!is.na(merged_end), merged_end, target_end),
    target_type = ifelse(!is.na(merged_type), merged_type, target_type)
  ) %>%
  # Select only the original columns
  select(chr, start, end, name, strand, spacer,
         target_chr, target_start, target_end,
         target_name, target_strand, target_type)

# Summarize screen stats by target_type
design_summary <- guide_targets_final %>%
  group_by(target_type) %>%
  summarize(
    n_gRNAs = n_distinct(name),
    n_guideSets = n_distinct(target_name)
  )
message("Summary stats for the screen design by target type after merging: ")
print(design_summary)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(guide_targets_final, snakemake@output$guide_targets_file)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)