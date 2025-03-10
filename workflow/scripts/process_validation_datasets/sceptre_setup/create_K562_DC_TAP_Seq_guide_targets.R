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
suppressPackageStartupMessages(
  library(tidyverse)
)

message("Loading input files")
guide_design_file <- read_tsv(snakemake@input$guide_design_file)



################################# TEMP UNTIL ADD SCREEN DESIGN TO PIPELINE
positive_controls <- read_tsv(snakemake@input$positive_controls, col_names = FALSE)
old_guide_targets <- read_tsv(snakemake@input$old_guide_targets)
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


### CHECKING TARGET SIZES ====================================================

# Let's visualize the target sizes
data <- data.frame(value = unique(guide_targets$target_end) - unique(guide_targets$target_start))

p <- ggplot(data, aes(x = value)) +
  geom_histogram(bins = 40, fill = "steelblue", color = "black") +
  scale_y_log10() +
  labs(title = "Histogram with Log10 Y-Axis", x = "Value", y = "Count (log10 scale)") +
  theme_classic()
# print(p)

# Filter and process target sizes
filtered_targets <- guide_targets %>%
  filter(target_end - target_start != 301) %>%
  mutate(size_diff = target_end - target_start) %>%
  group_by(target_name, target_chr, target_start, target_end) %>%  # Include chr, start, and end
  summarize(
    num_guides = n(),
    size_diff_values = paste(unique(size_diff), collapse = ", "),
    min_size = min(size_diff),  # Extract the smallest value for sorting
    .groups = "drop"
  ) %>%
  arrange(min_size) %>%
  select(target_name, num_guides, target_chr, target_start, target_end, size_diff_values)

# Print results
print(filtered_targets)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
write_tsv(guide_targets, snakemake@output$guide_targets_file)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)