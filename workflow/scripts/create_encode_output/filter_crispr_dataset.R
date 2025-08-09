## Filter ENCODE CRISPR dataset based on minimum statistical power for specified effect size. Adds
## a low power flag to the ValidConnection column if a pair is a negative below a specified power
## threshold (positives pass filter regardless of power). If specified removes all invalid
## connections based on ValidConnection column before writing to output. Also removes any pairs not
## within specified distance to TSS range
save.image(paste0("RDA_objects/filter_crispr_dataset", snakemake@wildcards$sample, ".rda"))

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# column types in ENCODE CRISPR data
crispr_data_cols <- cols(
  .default = col_double(),  # for power columns, which can vary in names
  chrom = col_character(),
  chromStart = col_integer(),
  chromEnd = col_integer(),
  name = col_character(),
  EffectSize = col_double(),
  strandPerturbationTarget = col_character(),
  PerturbationTargetID = col_character(),
  chrTSS = col_character(),
  startTSS = col_integer(),
  endTSS = col_integer(),
  strandGene = col_character(),
  EffectSize95ConfidenceIntervalLow = col_double(),
  EffectSize95ConfidenceIntervalHigh = col_double(),
  measuredGeneSymbol = col_character(),
  measuredEnsemblID = col_character(),
  guideSpacerSeq = col_character(),
  guideSeq = col_character(),
  Significant = col_logical(),
  pValue = col_double(),
  pValueAdjusted = col_double(),
  ValidConnection = col_character(),
  Notes = col_character(),
  Reference = col_character()
)

# load unfiltered dataset
dat <- read_tsv(snakemake@input[[1]], col_types = crispr_data_cols, progress = FALSE)

# # hard filter based on specified distance to TSS range
# dist_filter <- as.integer(snakemake@params$tss_to_dist)
# dat <- filter(dat, abs(distToTSS) > dist_filter[[1]], abs(distToTSS) <= dist_filter[[2]])

# This line was filtering out the positive controls
# Instead, let's just include a line in the ValidConnection

# hard "annotation" based on specified distance to TSS range
dist_filter <- as.integer(snakemake@params$tss_to_dist)
min_dist <- dist_filter[[1]]
max_dist <- dist_filter[[2]]
dat <- dat %>% 
  mutate(ValidConnection = case_when(
    abs(distToTSS) < min_dist ~ paste(ValidConnection, paste0("distance <", min_dist), sep = "; "),
    abs(distToTSS) > max_dist ~ paste(ValidConnection, paste0("distance >", max_dist), sep = "; "),
    TRUE ~ ValidConnection
  ))

# power filter
power_filt_col <- paste0("PowerAtEffectSize", snakemake@wildcards$es)
min_power <- as.numeric(snakemake@wildcards$pwr)

# ValidConnection flag to add to low power pairs
power_flag <- paste0("< ", min_power * 100, "% power at ", snakemake@wildcards$es, "% effect size")

# add filter data based on power
output <- dat %>% 
  mutate(ValidConnection = case_when(
    (!!sym(power_filt_col) < 0.8 | is.na(!!sym(power_filt_col))) & Significant != TRUE ~ 
      paste(ValidConnection, power_flag, sep = "; "),
    TRUE ~ ValidConnection
  ))

# write output to file
write_tsv(output, file = snakemake@output[[1]])
