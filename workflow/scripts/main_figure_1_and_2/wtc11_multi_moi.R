# Script: wtc11_multi_moi.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/wtc11_multi_moi.rda"))
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
  library(data.table)
})

message("Loading input files")
wtc11_multi_moi <- read_tsv(snakemake@input$wtc11_multi_moi)


### CREATE PLOT ===============================================================

# Data reshaping - all together
tab_tss <- wtc11_multi_moi[grepl(wtc11_multi_moi$perturbation, pattern = "TSS"), ]

# Keep only rows where gene name appears in the perturbation
keep_rows <- c()
for (i in 1:nrow(tab_tss)){
  gene <- tab_tss[i, "gene"]
  if (grepl(tab_tss$perturbation[i], pattern = gene)){
    keep_rows <- c(keep_rows, i)
  }
}
tss_tab1 <- tab_tss[keep_rows, ]

# Filter for significant results in MOI1
tss_tab1_use <- tss_tab1[tss_tab1$pval_adj_moi1 < 0.01, ]

# Convert to proper data.table (this is the key step)
tss_tab1_use_dt <- as.data.table(tss_tab1_use)

# Reshape data from wide to long format using data.table's melt
wtc11_tab <- data.table::melt(data = tss_tab1_use_dt,
                              id.vars = c("perturbation", "gene", "logFC_moi1"),
                              measure.vars = c("logFC_moi1", "logFC_moi2", "logFC_moi3",
                                               "logFC_moi6", "logFC_moi10"),
                              variable.name = "variable",
                              value.name = "value")

# Normalize all values to MOI1 - showing as positive KD values
wtc11_tab[, normalized_to_moi1 := value / logFC_moi1]

# Create and print the plot with modified y-axis
wtc11_multi_moi <- ggplot(wtc11_tab, aes(x = variable, 
                                         y = as.numeric(normalized_to_moi1),
                                         fill = variable)) +
  geom_boxplot(lwd=1) +
  scale_x_discrete(labels = c("moi1", "moi2", "moi3", "moi6", "moi10")) +
  theme_classic() + 
  theme(
    axis.text = element_text(size = 20), 
    axis.title = element_text(size = 20),
    legend.position = "none"
  ) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(labels = c("moi1", "moi2", "moi3", "moi6", "moi10"),
                    values = c("#FFFFB2", "#FFEDA0", "#FED976", "#FECC5C", "#FC4E2A")) +
  labs(fill = "MOI", 
       y = "KD Normalized to mean MOI1 KD", 
       x = "") +
  ylim(0, 2)

### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = wtc11_multi_moi, filename = snakemake@output$wtc11_multi_moi_plot, height = 5, width = 5, device = "pdf")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)