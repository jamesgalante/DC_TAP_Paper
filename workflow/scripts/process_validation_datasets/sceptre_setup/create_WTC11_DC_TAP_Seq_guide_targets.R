# Script: create_WTC11_DC_TAP_Seq_guide_targets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/create_WTC11_DC_TAP_Seq_guide_targets.rda"))
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


### =========================================================






############################## TEMP UNTIL SCREEN DESIGN ADDED TO PIPELINE

guide_design_file <- guide_design_file %>%
  mutate(tss_gene = case_when(
    str_detect(guideSet, "_TSS") ~ str_replace(guideSet, "^(.*?)(?:_(?:\\d+))?_TSS(?:_(?:\\d+))?$", "\\1"),
    TRUE ~ NA
    )
  )
           
# Set these genes as "tss_pos" and the rest as "tss_random": 
# See slides for why these genes are tss_pos: https://docs.google.com/presentation/d/1VqcZunNG8Gk_XWgI4YgDISJdOzLP29swu-FjGVnmbXc/edit#slide=id.g332dac9883a_1_8
# Evvie's file referenced in the slides is here: "resources/process_validation_datasets/sceptre_setup/WTC11_DC_TAP_Seq/wtc11_full_gene_set_use_real_w_ensemble_id.txt"
tss_pos_genes <- c("KCNH2", "KCNH2", "MESP1", "MYBPC3", "BAG3", "GSK3B", "POU5F1", "POU5F1", "ROCK1", "SOX2")


# Add target_type to the guide_design file
guide_design_file <- guide_design_file %>%
  mutate(
    target_type = case_when(
      tss_gene %in% tss_pos_genes ~ "tss_pos",
      str_detect(guideSet, "_TSS") ~ "tss_random",
      str_detect(guideSet, "_DE")  ~ "DE",
      startsWith(guideSet, "chr") ~ "enh",
      TRUE ~ guideSet # Should only be safe and non targeting
    )
  ) %>%
  select(-tss_gene)

############################## TEMP UNTIL SCREEN DESIGN ADDED TO PIPELINE






# Summarize screen stats by target_type
design_summary <- guide_design_file %>%
  group_by(target_type) %>%
  summarize(
    n_gRNAs = n_distinct(name),
    n_guideSets = n_distinct(guideSet)
  )
message("Summary stats for the screen design by target type: ")
print(design_summary)


### RENAME GUIDES TO MATCH MATRIX ============================================

# Name format in perturb_status.rds is WTC11_Random_Screen_Crop_1 - WTC11_Random_Screen_Crop_16197
# The current format in guide_design_file is WTC11_Random_Screen_sgOpti-DC_1 - WTC11_Random_Screen_sgOpti-DC_16197
# Need to convert the current format to match up with the format in perturb_status
guide_design_file$name <- gsub("sgOpti-DC", "Crop", guide_design_file$name)


### PROCESS EXTRA GUIDES =====================================================

# Add in 10 guides from the Jones group - obtained from previous TF-Perturb-Seq experiment
jones_group_guides <- data.frame(
  name = c(
    "WTC11_Random_Screen_Crop_SOX2-1553", "WTC11_Random_Screen_Crop_TCF4", "WTC11_Random_Screen_Crop_YBX1-821",
    "WTC11_Random_Screen_Crop_YBX1-823", "WTC11_Random_Screen_Crop_OCT4", "WTC11_Random_Screen_Crop_NANOG",
    "WTC11_Random_Screen_Crop_SOX2", "WTC11_Random_Screen_Crop_BAG3", "WTC11_Random_Screen_Crop_ROCK1",
    "WTC11_Random_Screen_Crop_GSK3B"
  ),
  GuideSequence = c(
    "GCGGCAGGATCGGCCAGAGG", "GCATCACCATGGACTCCCCCG", "GCCATAGAGTACCCGAGAGGA", # For GCATCACCATGGACTCCCCCG, I had GCATCACCATGGACTCCCCC before (no G)
    "GCCGCCGGCCGCCATAGAGAC", "GGGTGAAATGAGGGCTTGCGA", "GCCAGCAGAACGTTAAAATCC",
    "GCCCTGACAGCCCCCGTCACA", "GTTCATAAAGGTGCCCGGCGC", "GCGGGGCGCGGACGCTCGGAA",
    "GGGGTGGCTCGGAGATGCGAC"
  )
)


# Function to add "NGG" and trim to max 23 characters
adjust_guide_sequence <- function(seq) {
  if (is.na(seq)) return(NA)  # Handle NAs safely
  seq_pam <- paste0(seq, "NGG")  # Add NGG
  len <- nchar(seq_pam)
  
  # If longer than 23, trim from the start
  if (len > 23) {
    seq_pam <- substr(seq_pam, len - 22, len)  # Keep last 23 characters
  }
  
  return(seq_pam)
}

# Function to create "seq" by removing last three letters ("NGG") from PAM column
remove_pam_ngg <- function(seq_pam) {
  if (is.na(seq_pam) || nchar(seq_pam) < 3) return(NA)  # Handle short sequences
  return(substr(seq_pam, 1, nchar(seq_pam) - 3))  # Remove last 3 characters ("NGG")
}

# Process `jones_group_guides`
jones_group_guides <- jones_group_guides %>%
  mutate(
    GuideSequenceWithPAM = sapply(GuideSequence, adjust_guide_sequence),  # Add NGG & trim
    seq = sapply(GuideSequenceWithPAM, remove_pam_ngg)  # Remove "NGG" to get seq
  )

# BLAT results - performed on each "seq" individually through https://genome.ucsc.edu/cgi-bin/hgBlat with hg19
# GCGGCAGGATCGGCCAGAGG chr3   +   181429913 181429932
# CATCACCATGGACTCCCCCG chr18  +    53256940  53256959
# CCATAGAGTACCCGAGAGGA chr1   -    43147833  43147852
# CCGCCGGCCGCCATAGAGAC chr1   -    43147886  43147905
# GGTGAAATGAGGGCTTGCGA chr6   +    31138429  31138448
# CCAGCAGAACGTTAAAATCC chr12  -   7942016   7942035 (chr12:7789420-7789439) & chr15  +  35377487  35377506 (chr15:35085286-35085305)
# CCCTGACAGCCCCCGTCACA chr3   +   181429691 181429710
# TTCATAAAGGTGCCCGGCGC chr10  +   121411051 121411070
# CGGGGCGCGGACGCTCGGAA chr18  +    18691782  18691801
# GGGTGGCTCGGAGATGCGAC chr3   +   119812768 119812787

# For now, we'll keep the guides that didn't map or multi-mapped and we can remove these in downstream files (create_sceptre_diffex_input)

# The chr start-end is then the "seq" BLAT region where start = (BLAT start - 1) and end = (BLAT end)
# The guideSet was added by overlapping each added guide with all other guides in guide_design_file and using the overlapped guideSet - some didn't overlap any other - these are mentioned below
jones_group_guides <- jones_group_guides %>%
  mutate(
    chr = c("chr3", "chr18", "chr1", "chr1", "chr6", "chr12", "chr3", "chr10", "chr18", "chr3"),
    strand = c("+", "+", "-", "-", "+", "-", "+", "+", "+", "+"),
    start = c(181429912, 53256939, 43147832, 43147885, 31138428, 7942015, 181429690, 121411050, 18691781, 119812767),  # BLAT start - 1
    end = c(181429932, 53256959, 43147852, 43147905, 31138448, 7942035, 181429710, 121411070, 18691801, 119812787),  # BLAT end
    guideSet = c("SOX2_TSS", "TCF4_TSS", "YBX1_TSS", "YBX1_TSS", "POU5F1_TSS_2", "NANOG_TSS", "SOX2_TSS", "BAG3_TSS", "ROCK1_TSS", "GSK3B_TSS"),
    target_type = c("tss_pos", "tss_pos", "tss_pos", "tss_pos", "tss_pos", "tss_random", "tss_pos", "tss_pos", "tss_pos", "tss_pos")
  )

# Note that chr1:43147832-43147852 and chr1:43147885-43147905 didn't overlap any other guides in guide_design_file. In IGV these are on YBX1 promoter and the guide name indicates it's targeting YBX1 so I added the name here
# For chr12:7942015-7942035, this overlaps the NANOG promoter, so I added "NANOG_TSS" - note that this also overlaps a NANOG paralog
# For chr3:119812767-119812787, this correctly overlaps the GSK3B promoter when checked in IGV, but it's a bit further from other guides targeting that promoter
# The one NA value, which didn't map through BLAT to hg19 is labelled in the guide name as targeting TCF4, however no other guide in the experiment targets that gene (guide_design_file %>% filter(str_detect(guideSet, "TCF4")) %>% pull(guideSet))


### ADD IN EXTRA GUIDES ======================================================

# Ensure `jones_group_guides` has the same columns as `guide_design_file`
missing_cols <- setdiff(names(guide_design_file), names(jones_group_guides))

# Add missing columns as NA in `jones_group_guides`
for (col in missing_cols) {
  jones_group_guides[[col]] <- NA
}

# Ensure column order matches `guide_design_file`
jones_group_guides <- jones_group_guides %>% select(names(guide_design_file))

# Append `jones_group_guides` to `guide_design_file`
guide_design_file <- bind_rows(guide_design_file, jones_group_guides)


### CREATE THE TARGETS RANGES ================================================

# For all "enh" this should just be the guideSet name, split up into target_chr, target_start, target_end
# For all non-targeting guides, this should be "NA"
# For all DE of the same type, get the range of all guides targeting the same guideSet
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
      target_type == "enh" ~ sub(":.*", "", guideSet),
      target_type %in% c("negative_control", "safe_targeting") ~ NA_character_,
      TRUE ~ target_chr
    ),
    target_start = case_when(
      target_type == "enh" ~ sub(".*:(\\d+)-.*", "\\1", guideSet),
      target_type %in% c("negative_control", "safe_targeting") ~ NA_character_,
      TRUE ~ target_start
    ),
    target_end = case_when(
      target_type == "enh" ~ sub(".*-(\\d+)", "\\1", guideSet),
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
      target_type %in% c("DE", "tss_pos", "tss_random") ~ safe_min(start),
      TRUE ~ target_start
    ),
    target_end = case_when(
      target_type %in% c("DE", "tss_pos", "tss_random") ~ safe_max(end),
      TRUE ~ target_end
    ),
    target_chr = case_when(
      target_type %in% c("DE", "tss_pos", "tss_random") ~ if(length(unique(chr[!is.na(chr)])) == 0) NA_character_ else unique(chr[!is.na(chr)]),
      TRUE ~ target_chr
    )
  ) %>%
  ungroup()


### CREATE THE GUIDE TARGETS FILE ============================================

# Create the guide_targets file
guide_targets <- data.frame(
  chr = processed_guide_design_file$chr, # 
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