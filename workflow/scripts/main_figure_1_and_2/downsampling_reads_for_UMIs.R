# Script: downsampling_reads_for_UMIs.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/downsampling_reads_for_UMIs.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(rhdf5)
  library(Seurat)
})

message("Loading input files")
# Get dataset paths from Snakemake input object
datasets <- list(
  list(
    name = "DC-TAP-seq (MOI6)",
    mol_info = snakemake@input$moi6_mol_info,
    filtered_h5 = snakemake@input$moi6_filtered_h5
  ),
  list(
    name = "CRISPRi Direct Capture",
    mol_info = snakemake@input$crisprdi_mol_info,
    filtered_h5 = snakemake@input$crisprdi_filtered_h5
  ),
  list(
    name = "DC-TAP-seq (MOI3)",
    mol_info = snakemake@input$moi3_mol_info,
    filtered_h5 = snakemake@input$moi3_filtered_h5
  )
)

message(paste("Will process", length(datasets), "datasets:"))
for(dataset in datasets) {
  message(paste("  -", dataset$name))
}


### ANALYSIS FUNCTIONS ========================================================

get_complexity_df <- function(read_df, dataset_name, fracs = NULL) {
  message(paste("Processing reads data for", dataset_name))
  
  # Get the raw reads
  read_df$reads <- as.integer(read_df$count)
  read_df <- read_df %>% rename(cell = barcode_idx)
  
  # Create a dataframe with each read represented by a row
  message(paste("Expanding reads for", dataset_name, "..."))
  df_repeat <- read_df[rep(row.names(read_df), read_df$reads), c("umi", "cell", "gene_name")]
  
  # Define downsampling fractions (10% to 100%)
  if(is.null(fracs)) {
    fracs <- seq(0.1, 1.0, 0.1)
  }
  
  # Create empty results dataframe
  results <- data.frame(
    proportion = numeric(length(fracs)),
    mean_reads_per_cell = numeric(length(fracs)),
    mean_umis_per_cell_per_gene = numeric(length(fracs)),
    dataset = character(length(fracs)),
    stringsAsFactors = FALSE
  )
  
  # Downsample at each fraction
  for(i in seq_along(fracs)) {
    frac <- fracs[i]
    message(paste("Processing fraction", frac, "for", dataset_name))
    
    # Take a sampling of our reads at that fraction
    set.seed(42)
    sample_indices <- sample(nrow(df_repeat), size = floor(nrow(df_repeat) * frac))
    sample_df <- df_repeat[sample_indices, ]
    
    # Total reads leftover after sampling
    total_reads <- nrow(sample_df)
    
    # Drop the duplicate reads now to count unique molecular identifiers per cell per gene
    umis_cell_gene <- sample_df %>%
      distinct() %>%
      group_by(cell, gene_name) %>%
      summarise(umi_count = n(), .groups = "drop")
    
    # Remake a cell counts matrix (cell x gene with umi counts per gene in that cell)
    cell_counts_matrix <- umis_cell_gene %>%
      pivot_wider(names_from = gene_name, values_from = umi_count, values_fill = 0)
    
    # Get metrics from that matrix
    mean_umis_per_cell_per_gene <- cell_counts_matrix %>%
      select(-cell) %>%
      summarise_all(mean) %>%
      unlist() %>%
      mean()
    
    # Total reads per cell considers all in the sampling
    mean_reads_per_cell <- total_reads / length(unique(sample_df$cell))
    
    # Store results
    results$proportion[i] <- frac
    results$mean_reads_per_cell[i] <- mean_reads_per_cell
    results$mean_umis_per_cell_per_gene[i] <- mean_umis_per_cell_per_gene
    results$dataset[i] <- dataset_name
  }
  
  return(results)
}

process_dataset <- function(dataset) {
  name <- dataset$name
  mol_info_path <- dataset$mol_info
  filtered_h5_path <- dataset$filtered_h5
  
  message(paste("Processing", name, "..."))
  
  # Read in final output of counts using Seurat
  counts <- Read10X_h5(filtered_h5_path)
  if(is.list(counts)) {
    counts <- counts$`Gene Expression`  # Take gene expression if multiple assays
  }
  
  # Get gene info - handle both list and matrix cases
  counts_no_names <- Read10X_h5(filtered_h5_path, use.names = FALSE)
  if(is.list(counts_no_names)) {
    gene_ids <- rownames(counts_no_names$`Gene Expression`)
  } else {
    gene_ids <- rownames(counts_no_names)
  }
  
  gene_info <- data.frame(
    gene_ids = gene_ids,
    gene_names = rownames(counts),
    stringsAsFactors = FALSE
  )
  
  # Read in molecule_info.h5
  message("Reading molecule info file...")
  
  # Read molecular info tables
  barcode_idx <- h5read(mol_info_path, "barcode_idx")
  umi <- h5read(mol_info_path, "umi")
  count <- h5read(mol_info_path, "count")
  feature_idx <- h5read(mol_info_path, "feature_idx")
  
  # Read barcodes and gene IDs
  barcodes <- h5read(mol_info_path, "barcodes")
  gene_ids <- h5read(mol_info_path, "features/id")
  
  h5closeAll()
  
  # Create molecule info dataframe
  molecule_info <- data.frame(
    barcode_idx = barcode_idx,
    umi = umi,
    count = count,
    feature_idx = feature_idx,
    stringsAsFactors = FALSE
  )
  
  # Map barcode indices to actual barcodes (vectorized - much faster!)
  barcode_lookup <- as.character(barcodes)
  idx_1based <- molecule_info$barcode_idx + 1  # Convert to 1-based indexing
  valid_mask <- idx_1based >= 1 & idx_1based <= length(barcode_lookup)
  molecule_info$barcode_idx <- NA_character_
  molecule_info$barcode_idx[valid_mask] <- barcode_lookup[idx_1based[valid_mask]]
  
  # Map feature indices to gene IDs (vectorized - much faster!)
  gene_lookup <- as.character(gene_ids)
  gene_idx_1based <- molecule_info$feature_idx + 1  # Convert to 1-based indexing
  gene_valid_mask <- gene_idx_1based >= 1 & gene_idx_1based <= length(gene_lookup)
  molecule_info$feature_idx <- NA_character_
  molecule_info$feature_idx[gene_valid_mask] <- gene_lookup[gene_idx_1based[gene_valid_mask]]
  
  # Map to gene names from the counts matrix
  gene_map <- setNames(gene_info$gene_names, gene_info$gene_ids)
  molecule_info$gene_name <- gene_map[molecule_info$feature_idx]
  
  # Filter to only the 93 amplified genes that are in your counts matrix
  # This removes the other ~35,000 genes you don't care about
  molecule_info <- molecule_info[!is.na(molecule_info$gene_name), ]
  
  # Filter to cells in the counts matrix
  cell_barcodes <- gsub("-1$", "", colnames(counts))
  reads_df <- molecule_info[molecule_info$barcode_idx %in% cell_barcodes, ]
  
  message(paste("Final dataset has", nrow(reads_df), "molecular observations for analysis"))
  
  # Run the downsampling analysis
  return(get_complexity_df(reads_df, name))
  
}


### MAIN ANALYSIS =============================================================

message("Processing all datasets")

# Process all datasets
all_results <- list()

for(i in seq_along(datasets)) {
  results <- process_dataset(datasets[[i]])
  if(!is.null(results)) {
    all_results[[i]] <- results
  }
}

# Remove NULL results
all_results <- all_results[!sapply(all_results, is.null)]


### CREATE PLOTS ==============================================================

message("Creating plots")

# Combine all results
combined_results <- do.call(rbind, all_results)

# Create plot
saturation_plot <- ggplot(combined_results, aes(x = mean_reads_per_cell, y = mean_umis_per_cell_per_gene, color = dataset)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("#4477AA", "#EE6677", "#228833")) +
  labs(
    title = "DC-TAP-seq Effect Size Analysis",
    x = "Mean Total Reads per Cell",
    y = "Mean UMIs per Cell per Gene",
    color = "Dataset"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

message("Analysis complete. Saving outputs...")
  

### SAVE OUTPUT ===============================================================

message("Saving output files")

# Save the raw results
write_csv(combined_results, "combined_results.tsv")
write_csv(combined_results, snakemake@output$combined_results)

# Save plots
ggsave(filename = snakemake@output$plot_pdf, plot = saturation_plot, width = 10, height = 6, device = "pdf")
saveRDS(saturation_plot, file = "saturation_plot.rds")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)
