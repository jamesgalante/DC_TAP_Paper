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
gasperini_MAST_and_Sceptre <- readRDS(snakemake@input$gasperini_MAST_and_Sceptre)
combined_validation <- read_tsv(snakemake@input$combined_validation) %>% filter(Reference == "K562_DC_TAP_Seq") %>% filter(pred_id == "ENCODE_rE2G")
gasperini_sceptre_results_w_symbol <- readRDS(snakemake@input$gasperini_sceptre_results_w_symbol)

# Load for getting guide information
perturb_sce <- readRDS(snakemake@input$perturb_sce)
dctap_guide_targets <- read_tsv(snakemake@input$dc_tap_guide_targets)

# Load the singleton differential expression results to get guide effect sizes
singleton_dctap <- readRDS(snakemake@input$singleton_dctap)
singleton_gasperini <- readRDS(snakemake@input$singleton_gasperini)


### UNDERSTAND GASPERINI SCEPTRE V MAST ======================================

# Create a color column based on the classification of significance b/t MAST and Sceptre
gasperini_MAST_and_Sceptre <- gasperini_MAST_and_Sceptre %>%
  mutate(
    color_category = case_when(
      significant == TRUE & Sceptre_pctChange < 0 & Regulated == TRUE ~ "purple",   # Sceptre and MAST
      significant == TRUE & Sceptre_pctChange < 0 & Regulated == FALSE ~ "blue",  # Sceptre Only
      Regulated == TRUE & significant == FALSE ~ "red",                            # MAST Only
      TRUE ~ "grey"                                                                 # Neither
    )
  )

# Separate grey points and colored points to control point z axis
grey_points <- gasperini_MAST_and_Sceptre %>% filter(color_category == "grey")
colored_points <- gasperini_MAST_and_Sceptre %>% filter(color_category != "grey")

gasperini_MAST_v_Sceptre_plot <- gasperini_MAST_and_Sceptre %>%
  ggplot() +
  geom_point(data = grey_points, aes(x = MAST_pctChange, y = Sceptre_pctChange), color = "grey", alpha = 1, size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 0.7, aes(linetype = "y ~ x")) +
  geom_point(data = colored_points, aes(x = MAST_pctChange, y = Sceptre_pctChange, color = color_category),
             alpha = 1, size = 0.5) +
  geom_smooth(data = gasperini_MAST_and_Sceptre, aes(x = MAST_pctChange, y = Sceptre_pctChange), method = "lm",
              se = TRUE, color = "black", linetype = "solid") +
  geom_smooth(data = colored_points %>% filter(color_category == "purple"), 
              aes(x = MAST_pctChange, y = Sceptre_pctChange), method = "lm",
              se = TRUE, color = "purple", linetype = "solid") +
  labs(
    title = "Correlation b/t MAST and Sceptre Results for Gasperini",
    x = "MAST Effect Size (%)",
    y = "Sceptre Effect Size (%)",
    color = "Significance"
  ) +
  # Define manual colors and labels
  scale_color_manual(
    values = c("purple" = "purple", "blue" = "blue", "red" = "red"),
    labels = c(
      "purple" = "Sceptre and MAST",
      "blue" = "Sceptre Only",
      "red" = "MAST Only"
    )
  ) +
  # Set axis limits
  coord_cartesian(xlim = c(-80, 45), ylim = c(-80, 45)) +
  # Apply theme
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold"),
    aspect.ratio = 1
  )

# Save the plot
ggsave(plot = gasperini_MAST_v_Sceptre_plot, 
       filename = snakemake@output$gasperini_sceptre_v_mast_comparison, 
       device = "pdf", height = 5, width = 6)


### SET UP GRANGES ============================================================

# First define "Regulated" in gasperini_sceptre_results_w_symbol for downstream purpose
gasperini_sceptre_results_w_symbol <- gasperini_sceptre_results_w_symbol %>%
  mutate(Regulated = Sceptre_pctChange < 0 & significant == TRUE)

# Create the GRanges object
gr_gasperini <- GRanges(seqnames = gasperini_sceptre_results_w_symbol$chrom,  
                       ranges = IRanges(start = gasperini_sceptre_results_w_symbol$chromStart,  
                                        end = gasperini_sceptre_results_w_symbol$chromEnd),  
                       measuredGeneSymbol = gasperini_sceptre_results_w_symbol$gene_name,  
                       training = gasperini_sceptre_results_w_symbol)

gr_dc_tap <- GRanges(seqnames = combined_validation$chrom, 
                     ranges = IRanges(start = combined_validation$chromStart, 
                                      end = combined_validation$chromEnd), 
                     measuredGeneSymbol = combined_validation$measuredGeneSymbol, 
                     dc_tap = combined_validation)


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


### PLOT THE DUPLICATES =======================================================

# Without filtering for ValidConnection == TRUE, there are a few points that don't seem to line up between Gasperini and DC TAP
comparing_all_duplicate_pairs <- duplicates %>%
  ggplot(aes(y = dc_tap.EffectSize * 100, x = (2^training.log_2_fold_change - 1) * 100)) +
  geom_point(size = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Comparing Duplicate Pairs",
    y = "K562 DC TAP Seq Effect Size (%)",
    x = "Gasperini et al. 2019 Effect Size (%)"
  ) +
  theme_bw() +
  ylim(-50, 20) +
  xlim(-50, 20) +
  scale_color_brewer(palette = "Set1") +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank()
  )
ggsave(plot = comparing_all_duplicate_pairs,
       filename = snakemake@output$comparing_all_duplicate_pairs,
       device = "pdf", height = 3, width = 3)

# If we look further into some of the points that drop out after filtering for ValidConnection == TRUE, we find disparities between DC TAP and Gapserini
comparing_all_duplicate_pairs_with_color <- duplicates %>%
  mutate(color_code = case_when(
    dc_tap.Regulated == TRUE & training.Regulated == TRUE ~ "purple",
    dc_tap.Regulated == TRUE & training.Regulated == FALSE ~ "blue",
    dc_tap.Regulated == FALSE & training.Regulated == TRUE ~ "red",
    TRUE ~ "grey"
  )) %>%
  ggplot(aes(y = dc_tap.EffectSize * 100, x = training.Sceptre_pctChange, color = color_code)) +
  geom_point(size = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
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
    values = c("purple" = "purple", "blue" = "blue", "red" = "red"),
    labels = c(
      "purple" = "DC TAP & Gasperini",
      "blue" = "DC TAP Only",
      "red" = "Gasperini Only")
  ) +
  theme(
    aspect.ratio = 1,
    panel.grid.minor = element_blank()
  )
ggsave(plot = comparing_all_duplicate_pairs_with_color,
       filename = snakemake@output$comparing_all_duplicate_pairs_with_color,
       device = "pdf", height = 4, width = 4.5)


### SUBSETTING WEIRD PAIRS ====================================================

# All the pairs where the "color_code" is "red" are significant in Gasperini but not in DC TAP - let's look at these on a guide level
targets <- duplicates %>% filter(dc_tap.Regulated != training.Regulated) %>% pull(training.grna_target)

# Now we can get the coordinates of each of these guides for gasperini
gasp_grna_perts <- rowData(altExp(perturb_sce, "grna_perts"))
gasp_guide_information <- as.tibble(gasp_grna_perts) %>% filter(target_name %in% targets) %>% select(chr, start, end, name, target_name)

# Save bed file
write_tsv(gasp_guide_information %>% select(chr, start, end), snakemake@output$gasperini_guides_bed_file, col_names = FALSE)

# Now let's create a bed file so we can liftOver then visualize these guides in IGV
# The resulting liftedOver coordinates are here:
# chr8	129689617	129689640
# chr11	5280667	5280690
# chr11	5288289	5288312
# chr11	5280569	5280592
# chr2	55140483	55140506
# chr2	55140626	55140649
# chr11	5280877	5280900
# chr11	5280666	5280689
# chr11	5280802	5280825
# chr11	5284680	5284703
# chr8	129689534	129689557
# chr2	55140795	55140818
# chr2	55140665	55140688

# We can also get all guide information on the dc tap guides
dctap_targets <- duplicates %>% 
  filter(dc_tap.Regulated != training.Regulated) %>%
  mutate(dc_tap.grna_target = gsub("^[^|]*\\|([^:]*:[^:]*-[^:]*)[:].*", "\\1", dc_tap.name)) %>%
  pull(dc_tap.grna_target)
dctap_guide_information <- dctap_guide_targets %>% filter(target_name %in% dctap_targets) %>% select(chr, start, end, name, target_name)


### FOCUSING IN ON HBD ========================================================

plot_guide_effects <- function(gasp_target, dc_target, singleton_gasperini, singleton_dctap,
                               gasp_guide_info, dctap_guide_info, gene_id, gene_name) {
  
  # Get guides for these regions
  g_guides <- gasp_guide_info %>% filter(target_name == gasp_target) %>% pull(name)
  d_guides <- dctap_guide_info %>% filter(target_name == dc_target) %>% pull(name)
  
  # Get data for these guides targeting our gene
  g_data <- singleton_gasperini %>% 
    filter(grna_id %in% g_guides, response_id == gene_id) %>%
    mutate(pctChange = 2^log_2_fold_change - 1) %>%
    left_join(gasp_guide_info, by = c("grna_id" = "name")) %>%
    mutate(midpoint = (start + end)/2, dataset = "Gasperini")
  
  d_data <- singleton_dctap %>% 
    filter(grna_id %in% d_guides, response_id == gene_id) %>%
    mutate(pctChange = 2^log_2_fold_change - 1) %>%
    left_join(dctap_guide_info, by = c("grna_id" = "name")) %>%
    mutate(midpoint = (start + end)/2, dataset = "DCTAP")
  
  # Combine and plot
  combined_data <- bind_rows(g_data, d_data) %>% filter(!is.na(pctChange))
  
  ggplot(combined_data, aes(x = midpoint, y = pctChange, color = dataset)) +
    geom_errorbarh(aes(xmin = start, xmax = end), height = 0.02, alpha = 0.6) +
    geom_point(aes(shape = significant), size = 2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    scale_color_manual(values = c("Gasperini" = "red", "DCTAP" = "blue")) +
    scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
    labs(
      title = paste0(gene_name, " gRNA Effect Sizes"),
      x = "Genomic Position",
      y = "Percent Change"
    ) +
    theme_minimal()
}

# First pair (a lot of gasp guides)
hbd_plot1 <- plot_guide_effects(
  gasp_target = "chr11:5301716-5302216", dc_target = "chr11:5301767-5302068", 
  singleton_gasperini = singleton_gasperini, singleton_dctap = singleton_dctap,
  gasp_guide_info = gasp_guide_information, dctap_guide_info = dctap_guide_information,
  gene_id = "ENSG00000223609", gene_name = "HBD"
)

# Second pair (gasp guide off center)
hbd_plot2 <- plot_guide_effects(
  gasp_target = "chr11:5305780-5306280", dc_target = "chr11:5305872-5306173", 
  singleton_gasperini = singleton_gasperini, singleton_dctap = singleton_dctap,
  gasp_guide_info = gasp_guide_information, dctap_guide_info = dctap_guide_information,
  gene_id = "ENSG00000223609", gene_name = "HBD"
)

# Third pair (one gasp guide on center)
hbd_plot3 <- plot_guide_effects(
  gasp_target = "chr11:5309248-5309748", dc_target = "chr11:5309369-5309670", 
  singleton_gasperini = singleton_gasperini, singleton_dctap = singleton_dctap,
  gasp_guide_info = gasp_guide_information, dctap_guide_info = dctap_guide_information,
  gene_id = "ENSG00000223609", gene_name = "HBD"
)

# Save the plots if needed
ggsave(filename = file.path(dirname(snakemake@output$comparing_all_duplicate_pairs_with_color), "hbd_plot1.pdf"), plot = hbd_plot1, width = 6, height = 5)
ggsave(filename = file.path(dirname(snakemake@output$comparing_all_duplicate_pairs_with_color), "hbd_plot2.pdf"), plot = hbd_plot2, width = 6, height = 5)
ggsave(filename = file.path(dirname(snakemake@output$comparing_all_duplicate_pairs_with_color), "hbd_plot3.pdf"), plot = hbd_plot3, width = 6, height = 5)


### FOCUSING IN ON CHR8 LOCUS ========================================================

# Filter guides for the new locus
gasperini_guides <- gasp_guide_information %>% 
  filter(target_name == "chr8:130701583-130702083") %>% 
  pull(name)

dctap_guides <- dctap_guide_information %>% 
  filter(target_name == "chr8:130701786-130702087") %>% 
  pull(name)

target_gene_id <- "ENSG00000229140"

# Get singleton differential expression results for the guides targeting this gene
gasp_chr8_guides <- singleton_gasperini %>% 
  filter(grna_id %in% gasperini_guides) %>% 
  filter(response_id == target_gene_id)

dctap_chr8_guides <- singleton_dctap %>% 
  filter(grna_id %in% dctap_guides) %>% 
  filter(response_id == target_gene_id)

# It seems like dctap guides were ineffective for this locus - very few n_nonzero_trt cells


### MAKE FINAL PLOT ===========================================================

# Let's make the final plot with all pairs that are significant between both gasperini and dc tap
# Let's plot the correlation R value
# Let's color gray if non significant and darkblue if significant in both
# Create the final correlation plot
final_duplicate_pairs_plot <- duplicates %>%
  # Apply ValidConnection filter
  # filter(dc_tap.ValidConnection == TRUE) %>%
  # Only keep points where both are significant or both not significant
  
  filter((dc_tap.Regulated == TRUE & training.Regulated == TRUE) |
           (dc_tap.Regulated == FALSE & training.Regulated == FALSE)) %>%
  
  # Create color coding
  mutate(significance = if_else(dc_tap.Regulated == TRUE & training.Regulated == TRUE, 
                                "Both", "Neither")) %>%
  # Calculate correlation for significant points only
  {
    data <- .
    
    # Calculate r value only for significant points
    
    sig_points <- data %>% filter(significance == "Both")
    
    if(nrow(sig_points) > 1) {
      r_value <- cor(sig_points$dc_tap.EffectSize * 100, 
                     sig_points$training.Sceptre_pctChange, 
                     method = "pearson", 
                     use = "complete.obs")
    } else {
      r_value <- NA
    }
    
    # Create the plot
    ggplot(data, aes(y = dc_tap.EffectSize * 100, x = training.Sceptre_pctChange, 
                     color = significance)) +
      geom_point(size = 1.4) +
      geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
      # Add R value directly on the plot
      annotate("text", x = -45, y = 15, 
               label = paste0("r = ", round(r_value, 3)), 
               hjust = 0, size = 3.5) +
      labs(
        title = "",
        y = "K562 DC TAP Seq Effect Size (%)",
        x = "Gasperini et al. 2019 Effect Size (%)",
        color = "Regulated Status"
      ) +
      scale_color_manual(
        values = c("Both" = "darkblue", "Neither" = "grey70")
      ) +
      theme_bw() +
      ylim(-50, 20) +
      xlim(-50, 20) +
      theme(
        aspect.ratio = 1,
        panel.grid.minor = element_blank(),
        legend.position = "right",
        panel.grid = element_blank()
      )
  }

# Save the plot
ggsave(plot = final_duplicate_pairs_plot,
       filename = file.path(dirname(snakemake@output$comparing_all_duplicate_pairs_with_color), "final_correlation_plot.pdf"),
       device = "pdf", height = 4, width = 4.5)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)