# Script: calculate_summary_statistics_for_screen.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/calculate_summary_statistics_for_screen.rda"))
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
  library(cowplot)
})

message("Loading input files")
combined_validation_unfiltered <- read_tsv(snakemake@input$combined_validation_unfiltered)
combined_training_unfiltered <- read_tsv(snakemake@input$combined_training_unfiltered)

# Load in the negative control genes from the config file
k562_negative_control_genes <- as.vector(snakemake@params$k562_negative_control_genes)
wtc11_negative_control_genes <- as.vector(snakemake@params$wtc11_negative_control_genes)

# Load in the guide_targets files
k562_guide_targets <- read_tsv(snakemake@input$guide_targets[[1]])
wtc11_guide_targets <- read_tsv(snakemake@input$guide_targets[[2]])


### ARE ANY NEGATIVE CONTROL PAIRS REPRESENTED ================================

# Separate the wtc11 and k562 pairs
wtc11 <- combined_validation_unfiltered %>% filter(category == "WTC11 DC TAP Seq")
k562 <- combined_validation_unfiltered %>% filter(category == "K562 DC TAP Seq")

# Check how many (if any) negative control genes are represented in the pairs 
# nrow(wtc11 %>% filter(ValidConnection == TRUE) %>% filter(Regulated == TRUE) %>% filter(measuredGeneSymbol %in% wtc11_negative_control_genes))
# nrow(k562 %>% filter(ValidConnection == TRUE) %>% filter(Regulated == TRUE) %>% filter(measuredGeneSymbol %in% k562_negative_control_genes))

# Both of these have 0 rows meaning that there were no enhancer-gene pairs (Regulated == TRUE & ValidConnection == TRUE) discovered in the screen between the targets and negative control genes


### PROMOTER CASES ============================================================

# Are there any instances where ValidConnection contains "TSS targeting guide(s)" but not "overlaps potential promoter"?
combined_validation_unfiltered %>% 
  filter(str_detect(ValidConnection, "TSS targeting guide") & !str_detect(ValidConnection, "overlaps potential promoter"))

# There are three cases where we have a TSS targeting guide that doesn't overlap a potential promoter - all are K562
  # chr2 55048528 55048911 ERLEC1|chr2:55275664-55276047:. 0.02978613 chr2 53786680 53787180 ERLEC1 FALSE
  # chr2 55137588 55138076 ERLEC1|chr2:55364724-55365212:. 0.04190005 chr2 53786680 53787180 ERLEC1 FALSE
  # chr2 55137588 55138076 RTN4|chr2:55364724-55365212:. 0.01525616 chr2 55050348 55050848 RTN4 FALSE

# These only involve two purported TSS targets (don't be distracted by ERLEC1 being represented twice - this doesn't have to do with the target being a TSS control)
  # chr2 55048528 55048911 (hg38)
  # chr2 55137588 55138076 (hg38)

# Both of these instances are with RTN4 - this locus is a bit strange given how many annotated TSSs there are. There also only seems to be one active TSS
# I think, because these TSSs are so close together, that I will keep these as being labelled as TSSs, but it's good to keep this in mind when looking at results


### CREATE COMBINED TABLE W/ TSS NAMES ========================================

guide_targets_all <- bind_rows(
  k562_guide_targets,
  wtc11_guide_targets
) %>%
  group_by(target_name) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  mutate(
    # For rows where target_type is "TSS" (or the target_name contains "TSS"),
    # compute TSS_control_gene; otherwise, leave it as NA.
    TSS_control_gene = if_else(target_type == "TSS" | str_detect(target_name, "TSS"),
                               str_replace(target_name, "^(.*?)(?:_(?:\\d+))?_TSS(?:_(?:\\d+))?$", "\\1"),
                               NA_character_)
  )

# Parse combined_validation_unfiltered: extract coordinate information
combined_parsed <- combined_validation_unfiltered %>%
  # Split 'name' at the "|" and discard the gene part (we already have measuredGeneSymbol)
  separate(name, into = c(NA, "hg19_target_coords"), sep = "\\|", remove = FALSE) %>%
  # Remove trailing ":." from the coordinates and extract components
  mutate(
    hg19_target_coords = str_remove(hg19_target_coords, ":\\.$"),
    hg19_target_chr   = str_extract(hg19_target_coords, "chr[^:]+"),
    hg19_target_start = as.numeric(str_extract(hg19_target_coords, "(?<=:)[0-9]+")),
    hg19_target_end   = as.numeric(str_extract(hg19_target_coords, "(?<=-)[0-9]+"))
  )

# Join the parsed combined data with the combined TSS controls by matching coordinates.
combined_joined <- combined_parsed %>%
  left_join(
    guide_targets_all %>% select(target_name, target_type, target_chr, target_start, target_end, TSS_control_gene),
    by = c("hg19_target_chr" = "target_chr",
           "hg19_target_start" = "target_start",
           "hg19_target_end" = "target_end")
  )


### CALCULATE STATISTICS FOR CATEGORY 1 =======================================

message("~~~~~ 'DE-G' pairs (w/ ENCODE filtering, so filter for ValidConnection == TRUE)")
message("\tfor all pairs after encode filtering, what are the stats of sig and not")
cat1_significant_pairs <- combined_validation_unfiltered %>%
  filter(ValidConnection == TRUE, Significant == TRUE) %>%
  dplyr::rename(Downregulated = Regulated) %>%
  select(category, Downregulated) %>%
  table()
print(cat1_significant_pairs)

cat1_nonsignificant_pairs <- combined_validation_unfiltered %>%
  filter(Significant == FALSE) %>%
  filter(str_detect(ValidConnection, "TRUE")) %>%
  mutate("Underpowered" = str_detect(ValidConnection, "< 80% power at 15% effect size")) %>%
  select(category, Underpowered) %>%
  table()
print(cat1_nonsignificant_pairs)


### CALCULATE STATISTICS FOR CATEGORY 2a ======================================

message('~~~~~ "DP-G" pairs (distal promoters effects on another nearby gene)')
message("\tfor all targets that overlap a promoter in the ENCODE pipeline AND which aren't paired with the gene of the promoter they're overlapping, what are the stats of sig and not")
pairs_overlapping_nonself_promoter <- combined_validation_unfiltered %>%
  filter(str_detect(ValidConnection, "overlaps potential promoter")) %>%
  filter(!str_detect(ValidConnection, "distance <1000")) %>% # Filter out any pairs that are <1kb in distance (this will remove a lot of the TSS pos controls)
  rowwise() %>%
  mutate(promoters_overlapping_target = list({
    parts <- str_split(ValidConnection, ";")[[1]] %>% str_trim()
    pp <- parts[str_detect(parts, "overlaps potential promoter")]
    genes <- str_remove(pp, "overlaps potential promoter: ") %>% str_split("\\|") %>% unlist() %>% str_trim()
    genes
    })) %>%
  ungroup() %>%
  filter(!(measuredGeneSymbol %in% promoters_overlapping_target))

cat2_prom_significant_pairs <- pairs_overlapping_nonself_promoter %>%
  filter(Significant == TRUE) %>%
  dplyr::rename(Downregulated = Regulated) %>%
  select(category, Downregulated) %>%
  table()
print(cat2_prom_significant_pairs)

cat2_prom_nonsignificant_pairs <- pairs_overlapping_nonself_promoter %>%
  filter(Significant == FALSE) %>%
  mutate("Underpowered" = str_detect(ValidConnection, "< 80% power at 15% effect size")) %>%
  select(category, Underpowered) %>%
  table()
print(cat2_prom_nonsignificant_pairs)


### CALCULATE STATISTICS FOR CATEGORY 2b ======================================

# We want to now just filter out all instances where the measuredGeneSymbol and the TSS_control_gene are the same
message("\tfor all targets that are a TSS control in the guide design AND which aren't paired with the gene of the TSS target, what are the stats of sig and not")
combined_joined_tss_only <- combined_joined %>%
  filter(!str_detect(ValidConnection, "distance <1000")) %>% # Filter out any pairs that are <1kb in distance (this will remove a lot of the TSS pos controls)
  filter(!is.na(TSS_control_gene)) %>%
  filter(measuredGeneSymbol != TSS_control_gene)

# Now we calculate the summary statistics on this table : it's all the TSSs targeted and their gene pairs (not including the gene of the TSS)
cat2_tss_significant_pairs <- combined_joined_tss_only %>%
  filter(Significant == TRUE) %>%
  dplyr::rename(Downregulated = Regulated) %>%
  select(category, Downregulated) %>%
  table()
print(cat2_tss_significant_pairs)

cat2_tss_nonsignificant_pairs <- combined_joined_tss_only %>%
  filter(Significant == FALSE) %>%
  mutate(Underpowered = str_detect(ValidConnection, "< 80% power at 15% effect size")) %>%
  select(category, Underpowered) %>%
  table()
print(cat2_tss_nonsignificant_pairs)


### CALCULATE STATISTICS FOR CATEGORY 3a ======================================

message("~~~~~ X# Promoter effects on target genes")
message("\tfor all targets that overlap a promoter in the ENCODE pipeline AND ARE paired with the gene of the promoter they're overlapping, what are the stats of sig and not")
pairs_overlapping_self_promoter <- combined_validation_unfiltered %>%
  filter(str_detect(ValidConnection, "overlaps potential promoter")) %>%
  filter(!str_detect(ValidConnection, "distance <1000")) %>% # Filter out any pairs that are <1kb in distance (this will remove a lot of the TSS pos controls)
  rowwise() %>%
  mutate(promoters_overlapping_target = list({
    parts <- str_split(ValidConnection, ";")[[1]] %>% str_trim()
    pp <- parts[str_detect(parts, "overlaps potential promoter")]
    genes <- str_remove(pp, "overlaps potential promoter: ") %>% str_split("\\|") %>% unlist() %>% str_trim()
    genes
  })) %>%
  ungroup() %>%
  filter(measuredGeneSymbol %in% promoters_overlapping_target)

cat3_prom_significant_pairs <- pairs_overlapping_self_promoter %>%
  filter(Significant == TRUE) %>%
  dplyr::rename(Downregulated = Regulated) %>%
  select(category, Downregulated) %>%
  table()
print(cat3_prom_significant_pairs)

cat3_prom_nonsignificant_pairs <- pairs_overlapping_self_promoter %>%
  filter(Significant == FALSE) %>%
  mutate("Underpowered" = str_detect(ValidConnection, "< 80% power at 15% effect size")) %>%
  select(category, Underpowered) %>%
  table()
print(cat3_prom_nonsignificant_pairs)


### CALCULATE STATISTICS FOR CATEGORY 3b ======================================

message("\tfor all targets that are a TSS control in the guide design AND ARE paired with the gene of the promoter they're overlapping, what are the stats of sig and not")
combined_joined_tss_only <- combined_joined %>%
  filter(!is.na(TSS_control_gene)) %>%
  filter(measuredGeneSymbol == TSS_control_gene)

cat3_tss_significant_pairs <- combined_joined_tss_only %>%
  filter(Significant == TRUE) %>%
  dplyr::rename(Downregulated = Regulated) %>%
  select(category, Downregulated) %>%
  table()
print(cat3_tss_significant_pairs)

cat3_tss_nonsignificant_pairs <- combined_joined_tss_only %>%
  filter(Significant == FALSE) %>%
  mutate(Underpowered = str_detect(ValidConnection, "< 80% power at 15% effect size")) %>%
  select(category, Underpowered) %>%
  table()
print(cat3_tss_nonsignificant_pairs)


### CALCULATE STATISTICS FOR CATEGORY 4 =======================================

message("~~~~~ X# Positive control distal enhancer-gene pairs")

# Use combined_joined (this has the target_type) for each target
# Have to figure out the target_type for K562 - guide design file doesn't have it -> probably have to get elsewhere or something
# Would like to have labelled - positive control (TSS), positive control (DE), Locus gene TSS (TSS), Enh -> ask Judhajeet


### CREATE GRAPHS & TABLES ====================================================

# Function to create a bar plot for significant pairs
# Assumes the table has two dimensions: category and a binary variable (e.g. "Downregulated")
plot_sig <- function(table_obj, subtitle = "") {
  # Convert the table to a data frame and rename columns
  df <- as.data.frame(table_obj)
  # Rename columns: assuming dimnames are in order: category and Downregulated
  colnames(df) <- c("Category", "Group", "Count")
  
  # Set factor order for Group (adjust levels if needed)
  df$Group <- factor(df$Group, levels = c("FALSE", "TRUE"))
  
  p <- ggplot(df, aes(x = Category, y = Count, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Significant Pairs",
         subtitle = subtitle,
         x = "Category",
         y = "Count",
         fill = "Downregulated") +
    theme_classic()
  
  return(p)
}

# Function to create a bar plot for nonsignificant pairs
# Assumes the table has two dimensions: category and a binary variable (e.g. "Underpowered")
plot_non <- function(table_obj, subtitle = "") {
  # Convert the table to a data frame and rename columns
  df <- as.data.frame(table_obj)
  # Rename columns: assuming dimnames are in order: category and Underpowered
  colnames(df) <- c("Category", "Group", "Count")
  
  # Set factor order for Group (adjust levels if needed)
  df$Group <- factor(df$Group, levels = c("FALSE", "TRUE"))
  
  p <- ggplot(df, aes(x = Category, y = Count, fill = Group)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    labs(title = "Nonsignificant Pairs",
         subtitle = subtitle,
         x = "Category",
         y = "Count",
         fill = "Underpowered") +
    theme_classic()
  
  return(p)
}


# category 1
sig_cat1 <- plot_sig(cat1_significant_pairs, subtitle = "cat1")
non_cat1 <- plot_non(cat1_nonsignificant_pairs, subtitle = "cat1")

# category 2 prom
sig_cat2a <- plot_sig(cat1_significant_pairs, subtitle = "cat2a")
non_cat2a <- plot_non(cat1_nonsignificant_pairs, subtitle = "cat2a")
# category 2 tss
sig_cat2b <- plot_sig(cat1_significant_pairs, subtitle = "cat2b")
non_cat2b <- plot_non(cat1_nonsignificant_pairs, subtitle = "cat2b")

# category 3 prom
sig_cat3a <- plot_sig(cat1_significant_pairs, subtitle = "cat3a")
non_cat3a <- plot_non(cat1_nonsignificant_pairs, subtitle = "cat3a")
# category 3 tss
sig_cat3b <- plot_sig(cat1_significant_pairs, subtitle = "cat3b")
non_cat3b <- plot_non(cat1_nonsignificant_pairs, subtitle = "cat3b")

combined_grid <- plot_grid(
  sig_cat1, non_cat1,
  sig_cat2a, non_cat2a,
  sig_cat2b, non_cat2b,
  sig_cat3a, non_cat3a,
  sig_cat3b, non_cat3b,
  ncol = 2,
  align = "v",
  labels = "AUTO"
)

# Display the grid
# print(combined_grid)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
saveRDS(combined_joined, snakemake@output$unfiltered_validation_with_guide_targets)

# Save a text file with all of these summary statistics
# Save a bunch of plots which display all of these numbers


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)