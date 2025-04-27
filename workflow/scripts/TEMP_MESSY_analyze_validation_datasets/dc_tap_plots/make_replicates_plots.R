# Script: make_replicates_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/make_replicates_plots.rda"))
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
rep1 <- lapply(snakemake@input$discovery_results1, readRDS)
rep2 <- lapply(snakemake@input$discovery_results2, readRDS)


### CREATE REP FILES -=========================================================

# Bind the cell types together
k562 <- rep1[[1]] %>% left_join(rep2[[1]], by = c("response_id", "grna_target"))
wtc11 <- rep1[[2]] %>% left_join(rep2[[2]], by = c("response_id", "grna_target"))

# Filter for only pairs that pass qc
k562 <- k562 %>% filter(pass_qc.x == TRUE & pass_qc.y == TRUE)
wtc11 <- wtc11 %>% filter(pass_qc.x == TRUE & pass_qc.y == TRUE)

# Convert log_2_fold_change into pctChange
k562 <- k562 %>% 
  mutate(pctChange.x = 2^log_2_fold_change.x - 1) %>%
  mutate(pctChange.y = 2^log_2_fold_change.y - 1)

wtc11 <- wtc11 %>%
  mutate(pctChange.x = 2^log_2_fold_change.x - 1) %>%
  mutate(pctChange.y = 2^log_2_fold_change.y - 1)
  

### CREATE REP CORRELATION PLOTS ==============================================

# Calculate Pearson correlation coefficient
k562_filtered <- k562 %>% filter(significant.x == TRUE & significant.y == TRUE)
k562_cor <- cor(k562_filtered$pctChange.x, k562_filtered$pctChange.y, method = "pearson")
k562_cor_text <- paste("r =", round(k562_cor, 2))

# Create the plot for K562
k562_plot <- k562_filtered %>%
  ggplot(aes(x = pctChange.x * 100, y = pctChange.y * 100)) +
  geom_point(color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", 
           x = -Inf, y = Inf, 
           label = k562_cor_text, 
           hjust = -0.1, vjust = 1.5, 
           size = 4, color = "red") +
  labs(
    subtitle = "K562 Significant Pairs",
    title = "Replicate Comparison",
    x = "Replicate 1 Effect Size (%)",
    y = "Replicate 2 Effect Size (%)"
  ) +
  scale_x_continuous(limits = c(-60, 70), breaks = c(-60, -30, 0, 30, 60)) +
  scale_y_continuous(limits = c(-60, 70), breaks = c(-60, -30, 0, 30, 60)) +
  theme_classic() +
  theme(
    aspect.ratio = 1
  )

# Calculate Pearson correlation coefficient
wtc11_filtered <- wtc11 %>% filter(significant.x == TRUE & significant.y == TRUE)
wtc11_cor <- cor(wtc11_filtered$pctChange.x, wtc11_filtered$pctChange.y, method = "pearson")
wtc11_cor_text <- paste("r =", round(wtc11_cor, 2))

# Create the plot for WTC11
wtc11_plot <- wtc11_filtered %>%
  ggplot(aes(x = pctChange.x * 100, y = pctChange.y * 100)) +
  geom_point(color = "darkblue") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", 
           x = -Inf, y = Inf, 
           label = wtc11_cor_text, 
           hjust = -0.1, vjust = 1.5, 
           size = 4, color = "red") +
  labs(
    subtitle = "WTC11 Significant Pairs",
    title = "Replicate Comparison",
    x = "Replicate 1 Effect Size (%)",
    y = "Replicate 2 Effect Size (%)"
  ) +
  scale_x_continuous(limits = c(-60, 70), breaks = c(-60, -30, 0, 30, 60)) +
  scale_y_continuous(limits = c(-60, 70), breaks = c(-60, -30, 0, 30, 60)) +
  theme_classic() +
  theme(
    aspect.ratio = 1
  )


### CREATE REP NEG ONLY CORRELATION PLOTS =====================================

# Calculate Pearson correlation coefficient
k562_filtered <- k562 %>% filter(significant.x == TRUE & significant.y == TRUE) %>% filter(pctChange.x < 0 & pctChange.y < 0)
k562_cor <- cor(k562_filtered$pctChange.x, k562_filtered$pctChange.y, method = "pearson")
k562_cor_text <- paste("r =", round(k562_cor, 2))

# Create the plot for K562
k562_neg_plot <- k562_filtered %>%
  ggplot(aes(x = pctChange.x * 100, y = pctChange.y * 100)) +
  geom_point(color = "darkblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", 
           x = -Inf, y = Inf, 
           label = k562_cor_text, 
           hjust = -0.1, vjust = 1.5, 
           size = 4, color = "red") +
  labs(
    subtitle = "K562 Significant Pairs (<0)",
    title = "Replicate Comparison",
    x = "Replicate 1 Effect Size (%)",
    y = "Replicate 2 Effect Size (%)"
  ) +
  scale_x_continuous(limits = c(-60, 0), breaks = c(-60, -30, 0)) +
  scale_y_continuous(limits = c(-60, 0), breaks = c(-60, -30, 0)) +
  theme_classic() +
  theme(
    aspect.ratio = 1
  )

# Calculate Pearson correlation coefficient
wtc11_filtered <- wtc11 %>% filter(significant.x == TRUE & significant.y == TRUE) %>% filter(pctChange.x < 0 & pctChange.y < 0)
wtc11_cor <- cor(wtc11_filtered$pctChange.x, wtc11_filtered$pctChange.y, method = "pearson")
wtc11_cor_text <- paste("r =", round(wtc11_cor, 2))

# Create the plot for WTC11
wtc11_neg_plot <- wtc11_filtered %>%
  filter(significant.x == TRUE & significant.y == TRUE) %>%
  ggplot(aes(x = pctChange.x * 100, y = pctChange.y * 100)) +
  geom_point(color = "darkblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  annotate("text", 
           x = -Inf, y = Inf, 
           label = wtc11_cor_text, 
           hjust = -0.1, vjust = 1.5, 
           size = 4, color = "red") +
  labs(
    subtitle = "WTC11 Significant Pairs (<0)",
    title = "Replicate Comparison",
    x = "Replicate 1 Effect Size (%)",
    y = "Replicate 2 Effect Size (%)"
  ) +
  scale_x_continuous(limits = c(-60, 0), breaks = c(-60, -30, 0)) +
  scale_y_continuous(limits = c(-60, 0), breaks = c(-60, -30, 0)) +
  theme_classic() +
  theme(
    aspect.ratio = 1
  )


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
ggsave(plot = k562_plot, filename = snakemake@output$k562_plot, device = "pdf", height = 3, width = 3)
ggsave(plot = wtc11_plot, filename = snakemake@output$wtc11_plot, device = "pdf", height = 3, width = 3)
ggsave(plot = k562_neg_plot, filename = snakemake@output$k562_neg_plot, device = "pdf", height = 3, width = 3)
ggsave(plot = wtc11_neg_plot, filename = snakemake@output$wtc11_neg_plot, device = "pdf", height = 3, width = 3)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)