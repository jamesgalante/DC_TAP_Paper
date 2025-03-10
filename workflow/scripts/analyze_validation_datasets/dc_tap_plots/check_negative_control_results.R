# Script: check_negative_control_results.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/check_negative_control_results.rda"))
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
k562_discovery_results <- readRDS(snakemake@input$discovery_results[[1]])
wtc11_discovery_results <- readRDS(snakemake@input$discovery_results[[2]])

k562_final_sceptre_object <- readRDS(snakemake@input$final_sceptre_object[[1]])
wtc11_final_sceptre_object <- readRDS(snakemake@input$final_sceptre_object[[2]])

k562_negative_control_genes <- as.vector(snakemake@params$k562_negative_control_genes)
wtc11_negative_control_genes <- as.vector(snakemake@params$wtc11_negative_control_genes)

k562_features <- read_tsv(snakemake@input$features[[1]], col_names = FALSE)
wtc11_features <- read_tsv(snakemake@input$features[[2]], col_names = FALSE)


### VISUALIZE K562 RESULTS ====================================================

# Get the ENSEMBLE ID for each response_id in discovery_results
message("Processing K562")
k562_negative_control_results <- k562_discovery_results %>%
  left_join(k562_features, by = c("response_id" = "X1")) %>%
  filter(X2 %in% k562_negative_control_genes)

# How many significant hits are there for the negative control genes
message(paste0("There are ", sum(k562_negative_control_results$significant, na.rm = TRUE), " significant hits for the negative control genes"))

# Print the Gene-Target pairs that are significant
k562_negative_control_results %>%
  filter(significant == TRUE) %>%
  select(response_id, grna_target, X2)

# Volcano plot of hits
dashed_line <- k562_negative_control_results %>%
  filter(significant == TRUE) %>%
  pull(p_value) %>%
  max()

k562 <- k562_negative_control_results %>%
  ggplot(aes(x = log_2_fold_change, y = -log10(p_value), color = significant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = rgb(128/255, 51/255, 218/255), "FALSE" = rgb(69/255, 142/255, 247/255))) +  # Color by significance
  geom_hline(yintercept = -log10(dashed_line), linetype = "dashed", color = "black") +  # Dashed significance threshold line
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.2, 1.2)) +  # Explicit x-axis ticks
  labs(
    x = "Log fold change",
    y = "-log10(P-value)",
    title = "Discovery volcano plot K562",
    subtitle = "Negative Controls Only"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Remove legend if not needed
  )

# The points with the lowest p values are
# ENSG00000198668 (CALM1) - GATA1_TSS
# ENSG00000179115 (FARSA) - GATA1_TSS
# ENSG00000179115 (FARSA) - chr8:128972548-128972849


### VISUALIZE WTC11 RESULTS ===================================================

# Get the ENSEMBLE ID for each response_id in discovery_results
message("Processing WTC11")
wtc11_negative_control_results <- wtc11_discovery_results %>%
  left_join(wtc11_features, by = c("response_id" = "X1")) %>%
  filter(X2 %in% wtc11_negative_control_genes)

# How many significant hits are there for the negative control genes
message(paste0("There are ", sum(wtc11_negative_control_results$significant, na.rm = TRUE), " significant hits for the negative control genes"))

# Print the Gene-Target pairs that are significant
wtc11_negative_control_results %>%
  filter(significant == TRUE) %>%
  select(response_id, grna_target, X2)

# Volcano plot of hits
dashed_line <- wtc11_negative_control_results %>%
  filter(significant == TRUE) %>%
  pull(p_value) %>%
  max()

wtc11 <- wtc11_negative_control_results %>%
  ggplot(aes(x = log_2_fold_change, y = -log10(p_value), color = significant)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = rgb(128/255, 51/255, 218/255), "FALSE" = rgb(69/255, 142/255, 247/255))) +  # Color by significance
  geom_hline(yintercept = -log10(dashed_line), linetype = "dashed", color = "black") +  # Dashed significance threshold line
  scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.2, 1.2)) +  # Explicit x-axis ticks
  labs(
    x = "Log fold change",
    y = "-log10(P-value)",
    title = "Discovery volcano plot WTC11",
    subtitle = "Negative Controls Only"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.position = "none"  # Remove legend if not needed
  )

# The points with the lowest p values are
# ENSG00000141682 (PMAIP1) - chr8:145701156-145701457 
# ENSG00000141682 (PMAIP1) -                OCT4_DE_1  
# ENSG00000141682 (PMAIP1) - chr8:145700492-145700793  
# ENSG00000141682 (PMAIP1) -             POU5F1_TSS_2  
# ENSG00000141682 (PMAIP1) -                OCT4_DE_2  
# ENSG00000141682 (PMAIP1) - chr8:145699133-145699434  


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = k562, file = snakemake@output$k562_volcano, device = "pdf", height = 4, width = 4)
ggsave(plot = wtc11, file = snakemake@output$wtc11_volcano, device = "pdf", height = 4, width = 4)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)