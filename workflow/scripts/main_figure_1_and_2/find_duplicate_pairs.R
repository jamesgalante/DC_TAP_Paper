# Script: find_duplicate_pairs.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/find_duplicate_pairs.rda"))
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
  library(Seurat)
  library(SingleCellExperiment)
})

message("Loading input files")
gasperini_results <- read_tsv(snakemake@input$gasperini_results)
dc_tap_results <- read_tsv(snakemake@input$dc_tap_results) %>% filter(cell_type == "K562")


### REMOVE UNDERPOWERED NONSIGNIFICANT PAIRS ==================================

gasperini_results <- gasperini_results %>%
  filter(significant | power_at_effect_size_15 >= 0.8)

dc_tap_results <- dc_tap_results %>%
  filter(significant | power_at_effect_size_15 >= 0.8)


### SET UP GRANGES ============================================================

# Create the GRanges object
gr_gasperini <- GRanges(seqnames = gasperini_results$targeting_chr_hg38,  
                        ranges = IRanges(start = gasperini_results$targeting_start_hg38,  
                                         end = gasperini_results$targeting_end_hg38),  
                        measuredGeneSymbol = gasperini_results$gene_symbol,  
                        training = gasperini_results)

gr_dc_tap <- GRanges(seqnames = dc_tap_results$targeting_chr_hg38, 
                     ranges = IRanges(start = dc_tap_results$targeting_start_hg38, 
                                      end = dc_tap_results$targeting_end_hg38), 
                     measuredGeneSymbol = dc_tap_results$gene_symbol, 
                     dc_tap = dc_tap_results)


### OVERLAP THE GR OBJECTS ====================================================

# Find the overlaps and subset
overlaps <- findOverlaps(gr_dc_tap, gr_gasperini)
dc_tap_overlaps <- gr_dc_tap[queryHits(overlaps)]
gasperini_overlaps <- gr_gasperini[subjectHits(overlaps)]

# Only keep the pairs overlaps that have the same measuredGeneSymbol
matches <- dc_tap_overlaps$measuredGeneSymbol == gasperini_overlaps$measuredGeneSymbol
dc_tap_duplicates <- dc_tap_overlaps[matches] %>% as.data.frame() %>% mutate(id = seq_along(1:sum(matches)))
gasperini_duplicates <- gasperini_overlaps[matches] %>% as.data.frame() %>% mutate(id = seq_along(1:sum(matches)))

# Merge the two dataframes by ID (added in lines above)
duplicates <- dc_tap_duplicates %>%
  left_join(gasperini_duplicates, by = "id", suffix = c(".dc_tap", ".training"))


### CREATE DUPLICATES PLOT ====================================================

# Both dataframes are now filtered for valid pairs that pass filtering - lets plot their effect sizes and color by regulated status
comparing_all_duplicate_pairs <- duplicates %>%
  filter(dc_tap.DistalElement_Gene & training.DistalElement_Gene) %>%
  mutate(color_code = case_when(
    (dc_tap.significant & dc_tap.pct_change_effect_size < 0) & (training.significant & training.pct_change_effect_size < 0) ~ "darkblue",
    TRUE ~ "grey"
  )) %>%
  {
    data <- .
    
    # Calculate r value only for significant points (purple points where both are significant with negative effect)
    sig_points <- data %>% filter(color_code == "darkblue")
    
    if(nrow(sig_points) > 1) {
      r_value <- cor(sig_points$dc_tap.pct_change_effect_size, 
                     sig_points$training.pct_change_effect_size, 
                     method = "pearson", 
                     use = "complete.obs")
    } else {
      r_value <- NA
    }
    
    # Create the plot with the correlation annotation
    ggplot(data, aes(y = dc_tap.pct_change_effect_size, x = training.pct_change_effect_size, color = color_code)) +
      geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
      geom_point(size = 1.5) +
      # Add R value directly on the plot
      annotate("text", x = -45, y = 15, 
               label = paste0("r = ", round(r_value, 3)), 
               hjust = 0, size = 3.5, color = "red") +
      labs(
        title = "Comparing Duplicate Pairs",
        y = "K562 DC TAP Seq Effect Size (%)",
        x = "Gasperini et al. 2019 Effect Size (%)",
        color = "Regulated Status"
      ) +
      theme_bw() +
      ylim(-50, 20) +
      xlim(-50, 20) +
      scale_color_manual(
        values = c("darkblue" = "darkblue"),
        labels = c(
          "darkblue" = "DC TAP & Gasperini")
      ) +
      theme(
        aspect.ratio = 1,
        panel.grid = element_blank(),
        legend.position = "none"
      )
  }


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = comparing_all_duplicate_pairs, filename = snakemake@output$duplicate_pairs_plot, device = "pdf", height = 3, width = 3.5)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)