# Script: filter_loci.R
### SETUP =====================================================================
# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/filter_loci.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")

### LOADING FILES =============================================================
# Required packages and functions
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
})

# Load input files
message("Loading input files...")
locus_stats_df <- readRDS(snakemake@input$locus_stats)

# Get parameters
dhs_threshold <- snakemake@params$dhs_threshold
kb100_threshold <- snakemake@params$kb100_threshold
q75_threshold <- snakemake@params$q75_threshold

### ADD QUANTILE INFORMATION =================================================
message("Adding quantile information to loci...")
# Add peak quantile information to loci
locus_stats_df$q75counts <- NA
locus_stats_df$bottom75 <- NA
locus_stats_df$q75counts_p <- NA
locus_stats_df$bottom75_p <- NA

for (i in 1:nrow(locus_stats_df)) {
  quant_list <- locus_stats_df$dhs_quants[[i]]
  
  # Count peaks in top 75% quantile
  q75_count <- sum(quant_list == "75_1", na.rm = TRUE)
  locus_stats_df$q75counts[i] <- q75_count
  
  # Total number of peaks
  total_peaks <- length(quant_list)
  
  # Calculate proportion of peaks in top 75% quantile
  if (total_peaks > 0) {
    locus_stats_df$q75counts_p[i] <- q75_count / total_peaks
  } else {
    locus_stats_df$q75counts_p[i] <- 0
  }
  
  # Count peaks in bottom 75% quantiles
  locus_stats_df$bottom75[i] <- total_peaks - q75_count
  locus_stats_df$bottom75_p[i] <- (total_peaks - q75_count) / total_peaks
}

### ADD TSS DISTANCE INFORMATION =============================================
message("Adding TSS distance information to loci...")
# Add TSS distance information to loci
locus_stats_df$g100kb <- NA
locus_stats_df$l100kb <- NA
locus_stats_df$g100kb_p <- NA
locus_stats_df$l100kb_p <- NA

for (i in 1:nrow(locus_stats_df)) {
  distances <- locus_stats_df$dhs_nearest_tss_distance[[i]]
  
  if (length(distances) == 0) {
    locus_stats_df$g100kb[i] <- 0
    locus_stats_df$l100kb[i] <- 0
    locus_stats_df$g100kb_p[i] <- 0
    locus_stats_df$l100kb_p[i] <- 0
  } else {
    # Count peaks more than 100kb from TSS
    g100kb_count <- sum(distances > 100000, na.rm = TRUE)
    locus_stats_df$g100kb[i] <- g100kb_count
    
    # Count peaks less than 100kb but more than 1kb from TSS
    l100kb_count <- sum(distances <= 100000 & distances > 1000, na.rm = TRUE)
    locus_stats_df$l100kb[i] <- l100kb_count
    
    # Calculate proportions
    locus_stats_df$g100kb_p[i] <- g100kb_count / length(distances)
    locus_stats_df$l100kb_p[i] <- l100kb_count / length(distances)
  }
}

### FILTER LOCI ==============================================================
message("Filtering loci based on criteria...")
# Filter loci based on criteria
usable_loci <- locus_stats_df %>%
  filter(
    q75counts >= q75_threshold & 
      bottom75 >= q75_threshold & 
      g100kb >= kb100_threshold & 
      l100kb >= 30
  )

message(paste("Filtered down to", nrow(usable_loci), "usable loci out of", 
              nrow(locus_stats_df), "total loci"))
message(paste("Percentage of usable loci:", round(nrow(usable_loci) / nrow(locus_stats_df) * 100, 1), "%"))

### VISUALIZATION ============================================================
message("Creating filtering visualization plots...")

# Plot quantile distribution
q75_plot <- ggplot(locus_stats_df, aes(x = q75counts)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = q75_threshold, color = "red") +
  labs(title = "Distribution of peaks in top 75% quantile",
       x = "Number of peaks in top 75% quantile",
       y = "Count") +
  theme_bw()

# Plot >100kb distance distribution
g100kb_plot <- ggplot(locus_stats_df, aes(x = g100kb)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = kb100_threshold, color = "red") +
  labs(title = "Distribution of peaks >100kb from TSS",
       x = "Number of peaks >100kb from TSS",
       y = "Count") +
  theme_bw()

# Compare filtered vs unfiltered loci
locus_stats_df$filtered <- ifelse(
  locus_stats_df$locus %in% usable_loci$locus, 
  "Usable", 
  "Filtered out"
)

comparison_plot <- ggplot(locus_stats_df, aes(x = q75counts, y = g100kb, color = filtered)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = q75_threshold, linetype = "dashed") +
  geom_hline(yintercept = kb100_threshold, linetype = "dashed") +
  scale_color_manual(values = c("Usable" = "blue", "Filtered out" = "gray70")) +
  labs(title = "Loci filtering criteria",
       x = "Peaks in top 75% quantile",
       y = "Peaks >100kb from TSS",
       color = "Status") +
  theme_bw()

# Combine plots
combined_plot <- plot_grid(
  q75_plot, g100kb_plot, comparison_plot,
  ncol = 1,
  labels = c("A", "B", "C")
)

### SAVE OUTPUT ==============================================================
message("Saving output files")
saveRDS(usable_loci, file = snakemake@output$filtered_loci)
ggsave(snakemake@output$filtering_summary, plot = combined_plot, width = 10, height = 12)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)