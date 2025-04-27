# Script: QC_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/QC_plots.rda"))
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
# Load in the K562 files
k_gene_matrix <- fread(snakemake@input$k_gene_matrix)
k_guide_matrix <- fread(snakemake@input$k_guide_matrix)

# Load in the WTC11 files
w_gene_matrix <- fread(snakemake@input$w_gene_matrix)
w_guide_matrix <- fread(snakemake@input$w_guide_matrix)

# Create output plots directory if needed
output_directory <- dirname(snakemake@output$filler_plot)
# Check if the directory exists, and create it only if it doesn't
if (!file.exists(output_directory)) { dir.create(output_directory, recursive = TRUE) }


### UMI PLOTS ================================================================

## Distribution of UMIs per cell
# - genes

genes_umis_hist <- ggplot(data.frame(value = c(colSums(k_gene_matrix), colSums(w_gene_matrix)), 
                                     group = rep(c("K562", "WTC11"), c(ncol(k_gene_matrix), ncol(w_gene_matrix)))), 
                          aes(x = value, fill = group, color = group)) +  
  geom_histogram(alpha = 0.3, linewidth = 0.8, bins = 30, position = "identity") +  # Histogram with overlapping bars
  scale_fill_manual(values = c("K562" = "green", "WTC11" = "blue")) +  
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +  
  labs(title = "Histogram of Gene Expression", 
       x = "UMIs per Cell", y = "Number of Cells") +  # Update y-axis label
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +  # Fix log x-axis labels
  scale_y_continuous(expand = c(0, 0)) +  # Ensure y-axis starts at 0
  theme_classic() +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank()
  )  # Correct legend positioning
ggsave(plot = genes_umis_hist, filename = file.path(output_directory, "genes_umis_hist.pdf"), device = "pdf", width = 4.5, height = 4)

# - guides

guides_umis_hist <- ggplot(data.frame(value = c(colSums(k_guide_matrix), colSums(w_guide_matrix)), 
                                      group = rep(c("K562", "WTC11"), c(ncol(k_guide_matrix), ncol(w_guide_matrix)))), 
                           aes(x = value, fill = group, color = group)) +  
  geom_histogram(alpha = 0.3, linewidth = 0.8, bins = 30, position = "identity") +  # Histogram with overlapping bars
  scale_fill_manual(values = c("K562" = "green", "WTC11" = "blue")) +  
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +  
  labs(title = "Histogram of Guide Expression", 
       x = "UMIs per Cell", y = "Number of Cells") +  # Update y-axis label
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +  # Fix log x-axis labels
  scale_y_continuous(expand = c(0, 0)) +  # Ensure y-axis starts at 0
  theme_classic() +
  theme(
    legend.position = "bottom", 
    legend.title = element_blank()
  )  # Correct legend positioning
ggsave(plot = guides_umis_hist, filename = file.path(output_directory, "guides_umis_hist.pdf"), device = "pdf", width = 4.5, height = 4)

## Distribution of UMIs per cell per batch
# - genes

genes_umis_distr_w_batch <- ggplot(data.frame(
  value = c(colSums(k_gene_matrix), colSums(w_gene_matrix)),  # Sum UMIs per cell
  group = rep(c("K562", "WTC11"), c(ncol(k_gene_matrix), ncol(w_gene_matrix))),  # Label cell type
  batch = c(as.numeric(sub(".*-", "", colnames(k_gene_matrix))), 
            as.numeric(sub(".*-", "", colnames(w_gene_matrix))))  # Extract batch number
), 
aes(x = value, color = group, group = batch)) +  # Fill based on batch gradient
  geom_density(alpha = 0.8, linewidth = 1) +  # Separate density lines for each batch
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +  # Color by cell type
  labs(title = "Density Plot of Gene Expression by Batch", 
       x = "UMIs per Cell", y = "Density", color = "Cell Type") +  
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = c(0, 0)) +  
  theme_classic() +
  facet_grid(rows = vars(group)) +  # Split into two panels (one for each cell type)
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )  
ggsave(plot = genes_umis_distr_w_batch, filename = file.path(output_directory, "genes_umis_distr_w_batch.pdf"), device = "pdf", width = 4.5, height = 4)

# - guides

guides_umis_distr_w_batch <- ggplot(data.frame(
  value = c(colSums(k_guide_matrix), colSums(w_guide_matrix)),  # Sum UMIs per cell
  group = rep(c("K562", "WTC11"), c(ncol(k_guide_matrix), ncol(w_guide_matrix))),  # Label cell type
  batch = c(as.numeric(sub(".*-", "", colnames(k_guide_matrix))), 
            as.numeric(sub(".*-", "", colnames(w_guide_matrix))))  # Extract batch number
), 
aes(x = value, color = group, group = batch)) +  # Fill based on batch gradient
  geom_density(alpha = 0.8, linewidth = 1) +  # Separate density lines for each batch
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +  # Color by cell type
  labs(title = "Density Plot of Gene Expression by Batch", 
       x = "UMIs per Cell", y = "Density", color = "Cell Type") +  
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = c(0, 0)) +  
  theme_classic() +
  facet_grid(rows = vars(group)) +  # Split into two panels (one for each cell type)
  theme(
    legend.position = "bottom", 
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank()
  )  
ggsave(plot = guides_umis_distr_w_batch, filename = file.path(output_directory, "guides_umis_distr_w_batch.pdf"), device = "pdf", width = 4.5, height = 4)


### GENE PLOTS ===============================================================

## Number of unique genes per cell
# - genes

genes_hist <- ggplot(data.frame(value = c(colSums(k_gene_matrix > 0), colSums(w_gene_matrix > 0)), 
                                group = rep(c("K562", "WTC11"), c(ncol(k_gene_matrix), ncol(w_gene_matrix)))), 
                     aes(x = value, fill = group, color = group)) +  
  geom_histogram(alpha = 0.3, linewidth = 0.8, bins = 30, position = "identity") +
  scale_fill_manual(values = c("K562" = "green", "WTC11" = "blue")) +
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +
  labs(title = "Histogram of Unique Genes per Cell", 
       x = "Unique Genes per Cell", y = "Number of Cells") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(plot = genes_hist, filename = file.path(output_directory, "genes_hist.pdf"), device = "pdf", width = 4.5, height = 4)

# - guides

guides_hist <- ggplot(data.frame(value = c(colSums(k_guide_matrix > 0), colSums(w_guide_matrix > 0)), 
                                 group = rep(c("K562", "WTC11"), c(ncol(k_guide_matrix), ncol(w_guide_matrix)))), 
                      aes(x = value, fill = group, color = group)) +  
  geom_histogram(alpha = 0.3, linewidth = 0.8, bins = 30, position = "identity") +
  scale_fill_manual(values = c("K562" = "green", "WTC11" = "blue")) +
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +
  labs(title = "Histogram of Unique Guides per Cell", 
       x = "Unique Guides per Cell", y = "Number of Cells") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(plot = guides_hist, filename = file.path(output_directory, "guides_hist.pdf"), device = "pdf", width = 4.5, height = 4)


## Number of unique genes per cell per batch
# - genes

genes_distr_w_batch <- ggplot(data.frame(
  value = c(colSums(k_gene_matrix > 0), colSums(w_gene_matrix > 0)),
  group = rep(c("K562", "WTC11"), c(ncol(k_gene_matrix), ncol(w_gene_matrix))),
  batch = c(as.numeric(sub(".*-", "", colnames(k_gene_matrix))), 
            as.numeric(sub(".*-", "", colnames(w_gene_matrix))))
), 
aes(x = value, color = group, group = batch)) +  
  geom_density(alpha = 0.8, linewidth = 1) +
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +
  labs(title = "Density Plot of Unique Genes per Cell by Batch", 
       x = "Unique Genes per Cell", y = "Density", color = "Cell Type") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  facet_grid(rows = vars(group)) +  
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_blank())

ggsave(plot = genes_distr_w_batch, filename = file.path(output_directory, "genes_distr_w_batch.pdf"), device = "pdf", width = 4.5, height = 4)

# - guides

guides_distr_w_batch <- ggplot(data.frame(
  value = c(colSums(k_guide_matrix > 0), colSums(w_guide_matrix > 0)),
  group = rep(c("K562", "WTC11"), c(ncol(k_guide_matrix), ncol(w_guide_matrix))),
  batch = c(as.numeric(sub(".*-", "", colnames(k_guide_matrix))), 
            as.numeric(sub(".*-", "", colnames(w_guide_matrix))))
), 
aes(x = value, color = group, group = batch)) +  
  geom_density(alpha = 0.8, linewidth = 1) +
  scale_color_manual(values = c("K562" = "darkgreen", "WTC11" = "darkblue")) +
  labs(title = "Density Plot of Unique Guides per Cell by Batch", 
       x = "Unique Guides per Cell", y = "Density", color = "Cell Type") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  facet_grid(rows = vars(group)) +  
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.background = element_blank(), strip.text = element_blank())

ggsave(plot = guides_distr_w_batch, filename = file.path(output_directory, "guides_distr_w_batch.pdf"), device = "pdf", width = 4.5, height = 4)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)