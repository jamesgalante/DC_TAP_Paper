# Script: loci_qc.R
### SETUP =====================================================================
# Saving image for debugging
if (!file.exists("RDA_objects/k562_creating_targets")) { dir.create("RDA_objects/k562_creating_targets", recursive = TRUE) }
save.image(paste0("RDA_objects/k562_creating_targets/loci_qc.rda"))
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
  library(rmarkdown)
})

# Load input files
message("Loading input files...")
sampled_loci <- readRDS(snakemake@input$sampled_loci)
locus_stats_df <- readRDS(snakemake@input$locus_stats)

### COMPARE LOCI STATISTICS ==================================================
message("Comparing sampled loci with overall locus statistics...")
# Mark sampled loci in overall locus statistics
locus_stats_df$is_sampled <- locus_stats_df$locus %in% sampled_loci$locus

### CREATE DENSITY COMPARISON PLOTS ==========================================
message("Creating comparison plots...")

# Plot quantile distribution comparison
q75_plot <- ggplot(locus_stats_df, aes(x = q75counts_p, color = is_sampled, fill = is_sampled)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Distribution of peaks in top 75% quantile",
       x = "Proportion of peaks in top 75% quantile",
       y = "Density",
       color = "Sampled",
       fill = "Sampled") +
  theme_bw()

# Plot TSS distance distribution comparison
g100kb_plot <- ggplot(locus_stats_df, aes(x = g100kb_p, color = is_sampled, fill = is_sampled)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Distribution of peaks >100kb from TSS",
       x = "Proportion of peaks >100kb from TSS",
       y = "Density",
       color = "Sampled",
       fill = "Sampled") +
  theme_bw()

# Plot DHS count distribution comparison
dhs_plot <- ggplot(locus_stats_df, aes(x = dhs, color = is_sampled, fill = is_sampled)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Distribution of DHS count per locus",
       x = "Number of DHS peaks",
       y = "Density",
       color = "Sampled",
       fill = "Sampled") +
  theme_bw() +
  scale_x_log10()

# Plot gene count distribution comparison
genes_plot <- ggplot(locus_stats_df, aes(x = genes, color = is_sampled, fill = is_sampled)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Distribution of gene count per locus",
       x = "Number of genes",
       y = "Density",
       color = "Sampled",
       fill = "Sampled") +
  theme_bw() +
  scale_x_log10()

# Plot genes above TPM threshold distribution comparison
genes_tpm_plot <- ggplot(locus_stats_df, aes(x = genes_above_tpms, color = is_sampled, fill = is_sampled)) +
  geom_density(alpha = 0.3) +
  scale_color_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Distribution of genes above TPM threshold per locus",
       x = "Number of genes above TPM threshold",
       y = "Density",
       color = "Sampled",
       fill = "Sampled") +
  theme_bw()

### CREATE BOXPLOT COMPARISON PLOTS ==========================================
# Compare with boxplots
q75_boxplot <- ggplot(locus_stats_df, aes(x = is_sampled, y = q75counts_p, fill = is_sampled)) +
  geom_boxplot() +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Proportion of peaks in top 75% quantile",
       x = "Sampled",
       y = "Proportion") +
  theme_bw()

g100kb_boxplot <- ggplot(locus_stats_df, aes(x = is_sampled, y = g100kb_p, fill = is_sampled)) +
  geom_boxplot() +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Proportion of peaks >100kb from TSS",
       x = "Sampled",
       y = "Proportion") +
  theme_bw()

dhs_boxplot <- ggplot(locus_stats_df, aes(x = is_sampled, y = dhs, fill = is_sampled)) +
  geom_boxplot() +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Number of DHS peaks",
       x = "Sampled",
       y = "Count") +
  theme_bw()

genes_boxplot <- ggplot(locus_stats_df, aes(x = is_sampled, y = genes, fill = is_sampled)) +
  geom_boxplot() +
  scale_fill_manual(values = c(`TRUE` = "blue", `FALSE` = "gray50")) +
  labs(title = "Number of genes",
       x = "Sampled",
       y = "Count") +
  theme_bw()

### COMBINE PLOTS ===========================================================
# Combine plots
combined_plot_density <- plot_grid(
  q75_plot, g100kb_plot, dhs_plot, genes_plot, genes_tpm_plot,
  ncol = 2,
  labels = c("A", "B", "C", "D", "E")
)

combined_plot_boxplot <- plot_grid(
  q75_boxplot, g100kb_boxplot, dhs_boxplot, genes_boxplot,
  ncol = 2,
  labels = c("F", "G", "H", "I")
)

combined_plot <- plot_grid(
  combined_plot_density, combined_plot_boxplot,
  ncol = 1,
  rel_heights = c(3, 2)
)

### GENERATE SUMMARY STATISTICS ==============================================
message("Generating summary statistics...")
summary_stats <- data.frame(
  Metric = c(
    "Total loci", 
    "Sampled loci",
    "Mean DHS per locus (sampled)", 
    "Mean genes per locus (sampled)", 
    "Mean genes above TPM threshold per locus (sampled)",
    "Mean proportion of peaks in top 75% quantile (sampled)",
    "Mean proportion of peaks >100kb from TSS (sampled)"
  ),
  Value = c(
    nrow(locus_stats_df),
    nrow(sampled_loci),
    mean(sampled_loci$dhs),
    mean(sampled_loci$genes),
    mean(sampled_loci$genes_above_tpms),
    mean(sampled_loci$q75counts_p, na.rm = TRUE),
    mean(sampled_loci$g100kb_p, na.rm = TRUE)
  )
)

### CREATE QC REPORT ========================================================
message("Creating QC report...")
report_rmd <- tempfile(fileext = ".Rmd")
cat(
  '---
title: "CRISPRi Screen Design Loci QC Report"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
library(knitr)
library(ggplot2)
library(cowplot)
```

## Summary Statistics

```{r summary}
kable(summary_stats)
```

## Loci Distribution Comparisons

```{r plots, fig.width=12, fig.height=10}
combined_plot
```

## Sampled Loci Details

```{r loci_details}
sampled_loci_summary <- sampled_loci %>%
  select(locus, chr, dhs, genes, genes_above_tpms) %>%
  arrange(chr, locus)

kable(sampled_loci_summary)
```

## Chromosome Distribution

```{r chr_dist}
ggplot(sampled_loci, aes(x = chr)) +
  geom_bar() +
  labs(title = "Chromosome distribution of sampled loci",
       x = "Chromosome", y = "Count") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```
',
  file = report_rmd
)

# Render RMarkdown to HTML
rmarkdown::render(
  report_rmd,
  output_file = snakemake@output$qc_report,
  envir = new.env()
)

### SAVE OUTPUT ==============================================================
message("Saving output files")
ggsave(snakemake@output$qc_plots, plot = combined_plot, width = 12, height = 15)

### CLEAN UP =================================================================
message("Closing log file")
sink()
sink(type = "message")
close(log)