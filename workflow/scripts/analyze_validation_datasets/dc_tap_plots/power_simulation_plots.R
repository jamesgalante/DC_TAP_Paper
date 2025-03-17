# Script: power_simulation_plots.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/power_simulation_plots.rda"))
message("Saved Image")
stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(rtracklayer)
  library(cowplot)
})

message("Loading input files")
k562_power_output <- read_tsv(snakemake@input$k562_power_output)
wtc11_power_output <- read_tsv(snakemake@input$wtc11_power_output)
gasperini_power_output <- read_tsv(snakemake@input$gasperini_power_output)
k562_tpm_file <- read_tsv(snakemake@input$k562_tpm_file)
annot <- import(snakemake@input$annot)

# Process annot file to get gene name conversions
annot <- annot[annot$type == "gene"]
annot <- annot[!grepl("_PAR_Y", annot$gene_id),]
gene_name_conversion <- annot %>% 
  as_tibble() %>% 
  mutate(ensemble_id = sub("\\..*", "", gene_id)) %>% 
  dplyr::select(gene_name, ensemble_id)

# Convert gene names in tpm file to ensemble ids to match with gasperini power output
k562_tpm_file_w_ens <- k562_tpm_file %>%
  left_join(gene_name_conversion, by = c("gene" = "gene_name"))

# There are some gene symbols that have multiple ensemble_ids
# Let's remove these ensemble_ids to avoid confusion with the mapping
k562_tpm_file_w_ens_filt <- k562_tpm_file_w_ens %>% 
  filter(!ensemble_id %in% (k562_tpm_file_w_ens[which(duplicated(k562_tpm_file_w_ens$gene)), ] %>% pull(ensemble_id)))


### POWER SIM RESULTS =========================================================

# Plot the number of pairs tested that have 80% power to detect each effect size [10, 15, 20, 25, 50]
# Function to calculate power summary for each dataset (enhancers only)
calculate_power_summary <- function(data, cell_type) {
  data %>%
    filter(target_type == "enh") %>%
    summarize(
      `10` = sum(PowerAtEffectSize10 >= 0.8) / n() * 100,
      `15` = sum(PowerAtEffectSize15 >= 0.8) / n() * 100,
      `20` = sum(PowerAtEffectSize20 >= 0.8) / n() * 100,
      `25` = sum(PowerAtEffectSize25 >= 0.8) / n() * 100,
      `50` = sum(PowerAtEffectSize50 >= 0.8) / n() * 100
    ) %>%
    pivot_longer(
      cols = everything(),
      names_to = "effect_size",
      values_to = "percentage"
    ) %>%
    mutate(
      cell_type = cell_type,
      effect_size_label = paste0(effect_size, "%"),
      effect_size_num = as.numeric(effect_size),
      effect_size = factor(effect_size, levels = c("10", "15", "20", "25", "50"))
    )
}

# Calculate power summary for both cell types
k562_summary <- calculate_power_summary(k562_power_output, "K562")
wtc11_summary <- calculate_power_summary(wtc11_power_output, "WTC11")

# Combine data for plotting
combined_summary <- bind_rows(k562_summary, wtc11_summary) %>%
  mutate(cell_type = factor(cell_type, levels = c("K562", "WTC11")))

# Define common theme elements
common_theme <- theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text = element_text(size = 9),
    aspect.ratio = 1
  )

# Define common color scheme
cell_colors <- c("K562" = "firebrick", "WTC11" = "darkblue")

# Create the bar plot with grouped bars
bar_plot <- ggplot(combined_summary, 
                   aes(x = effect_size, y = percentage, fill = cell_type)) +
  # Grouped bars side by side
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # Labels and title
  labs(
    title = "Tested pairs with High Power",
    x = "Effect Size",
    y = "Percentage of Pairs (%)",
    fill = "Cell Type"
  ) +
  # Set y-axis limits and breaks
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  # Set x-axis labels and spacing
  scale_x_discrete(labels = function(x) paste0(x, "%"), 
                   expand = c(0.1, 0.1)) +
  # Set colors for cell types
  scale_fill_manual(values = cell_colors) +
  common_theme

# Create the line plot
line_plot <- ggplot(combined_summary, 
                    aes(x = effect_size, y = percentage, color = cell_type, group = cell_type)) +
  # Add points
  geom_point(size = 3) +
  # Add connecting lines
  geom_line(linewidth = 1) +
  # Labels and title
  labs(
    title = "Tested pairs with High Power",
    x = "Effect Size",
    y = "Percentage of Pairs (%)",
    color = "Cell Type"
  ) +
  # Set y-axis limits and breaks
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  # Set x-axis labels and spacing
  scale_x_discrete(labels = function(x) paste0(x, "%"), 
                   expand = c(0.1, 0.1)) +
  # Set colors for cell types
  scale_color_manual(values = cell_colors) +
  common_theme

# Let's save both plots
ggsave(plot = line_plot, 
       filename = snakemake@output$power_sim_results_line_plot,
       device = "pdf",
       height = 4,
       width = 4)
ggsave(plot = bar_plot,
       filename = snakemake@output$power_sim_results_bar_plot,
       device = "pdf",
       height = 4,
       width = 4)


### GASPERINI V K562 DC TAP ===================================================

# We want to compare, for the same genes, the statistical power to detect a 15% effect size between Gasperini and K562 DC TAP Seq

# Function to process Gasperini data - with column renaming
process_gasperini <- function(power_data, tpm_data) {
  # Get the key columns needed for comparison
  processed <- power_data %>%
    # Ensure we're working with enhancer-gene pairs
    filter(target_type == "enh") %>%
    # Add TPM data by matching on ensemble_id
    left_join(tpm_data %>% dplyr::select(ensemble_id, TPM), by = c("gene" = "ensemble_id")) %>%
    # Rename power columns to match K562 format (camel case)
    dplyr::rename(
      PowerAtEffectSize10 = "power_effect_size_10",
      PowerAtEffectSize15 = "power_effect_size_15",
      PowerAtEffectSize20 = "power_effect_size_20",
      PowerAtEffectSize25 = "power_effect_size_25",
      PowerAtEffectSize50 = "power_effect_size_50"
    )
  
  # Report percentage of genes with TPM values
  tpm_count <- sum(!is.na(processed$TPM))
  total_count <- nrow(processed)
  message(paste0("Gasperini: Found TPM values for ", tpm_count, " out of ", total_count, 
                 " genes (", round(tpm_count/total_count*100, 1), "%)"))
  
  # Return processed data
  return(processed)
}

# Function to process K562 data
process_k562 <- function(power_data, tpm_data) {
  # K562 has gene column with ensemble IDs already
  processed <- power_data %>%
    # Ensure we're working with enhancer-gene pairs
    filter(target_type == "enh") %>%
    # Add TPM data by matching ensemble_id with TPM data
    left_join(k562_tpm_file_w_ens_filt %>% dplyr::select(ensemble_id, TPM), 
              by = c("gene" = "ensemble_id"))
  
  # Report percentage of genes with TPM values
  tpm_count <- sum(!is.na(processed$TPM))
  total_count <- nrow(processed)
  message(paste0("K562: Found TPM values for ", tpm_count, " out of ", total_count, 
                 " genes (", round(tpm_count/total_count*100, 1), "%)"))
  
  # Return processed data
  return(processed)
}

# Process Gasperini and K562 data
gasperini_processed <- process_gasperini(gasperini_power_output, k562_tpm_file_w_ens_filt)
k562_processed <- process_k562(k562_power_output, k562_tpm_file_w_ens_filt)

# Define TPM ranges for comparison
tpm_ranges <- list(
  c(5, 15),
  c(15, 25),
  c(25, 35),
  c(35, 45),
  c(45, 55)
)

# Store plots in a list but don't save them yet
tpm_plots <- list()

# Define the effect size to analyze (10, 15, 20, 25, or 50)
effect_size <- 20  # Change this value to analyze different effect sizes

# For each TPM range, create a plot
for (i in seq_along(tpm_ranges)) {
  low <- tpm_ranges[[i]][1]
  high <- tpm_ranges[[i]][2]
  
  # Construct the power column name based on effect size
  power_col <- paste0("PowerAtEffectSize", effect_size)
  
  # Filter data for current TPM range
  k562_filtered <- k562_processed %>% 
    filter(TPM > low & TPM < high) %>%
    mutate(source = "K562 DC TAP")
  
  gasperini_filtered <- gasperini_processed %>% 
    filter(TPM > low & TPM < high) %>%
    mutate(source = "Gasperini et al. 2019")
  
  # Skip this range if no data in either dataset
  if (nrow(k562_filtered) == 0 && nrow(gasperini_filtered) == 0) {
    message(paste0("No data for TPM range ", low, "-", high))
    next
  }
  
  # Bin the data to reduce noise - especially important for high cell counts
  k562_binned <- k562_filtered %>%
    mutate(
      # More granular binning at lower cell counts, wider bins at higher counts
      bin_width = case_when(
        cells < 200 ~ 100,  # 20-cell bins for < 200 cells
        cells < 500 ~ 200,  # 50-cell bins for 200-500 cells
        cells < 1000 ~ 300, # 100-cell bins for 500-1000 cells
        TRUE ~ 400         # 200-cell bins for >= 1000 cells
      ),
      cells_binned = bin_width * round(cells / bin_width)
    ) %>%
    group_by(source, cells_binned) %>%
    summarise(
      power_value = mean(!!sym(power_col), na.rm = TRUE),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    filter(n_samples >= 3)  # Minimum sample threshold
  
  gasperini_binned <- gasperini_filtered %>%
    mutate(
      bin_width = case_when(
        cells < 200 ~ 100,
        cells < 500 ~ 200,
        cells < 1000 ~ 300,
        cells < 2000 ~ 600,
        TRUE ~ 3000         # Wider bins for the very high cell counts
      ),
      cells_binned = bin_width * round(cells / bin_width)
    ) %>%
    group_by(source, cells_binned) %>%
    summarise(
      power_value = mean(!!sym(power_col), na.rm = TRUE),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    filter(n_samples >= 3)  # Minimum sample threshold
  
  # Combine binned data
  combined_binned <- bind_rows(k562_binned, gasperini_binned)
  
  # Find first point where each source reaches 80% power - avoid using slice()
  k562_threshold <- k562_binned %>%
    filter(power_value * 100 >= 80) %>%
    arrange(cells_binned)
  
  gasperini_threshold <- gasperini_binned %>%
    filter(power_value * 100 >= 80) %>%
    arrange(cells_binned)
  
  # Get the first row if available
  k562_threshold_cells <- if(nrow(k562_threshold) > 0) k562_threshold$cells_binned[1] else NULL
  gasperini_threshold_cells <- if(nrow(gasperini_threshold) > 0) gasperini_threshold$cells_binned[1] else NULL
  
  # Create the plot
  p <- ggplot(combined_binned, aes(x = cells_binned, y = power_value * 100, color = source)) +
    geom_line(size = 1.3) +
    geom_point() +  # Standard size points, no scaling by sample count
    geom_hline(yintercept = 80, linetype = "dashed") +
    labs(
      title = paste0("Power at ", effect_size, "% Effect Size (TPM ", low, "-", high, ")"),
      #subtitle = paste0("Genes: K562 n=", nrow(k562_filtered), ", Gasperini n=", nrow(gasperini_filtered)),
      x = "Number of Perturbed Cells",
      y = "Power (%)",
      color = "Source"
    ) +
    scale_x_log10() +
    scale_color_manual(values = c("K562 DC TAP" = "firebrick", "Gasperini et al. 2019" = "darkblue")) +
    theme_classic() +
    theme(
      aspect.ratio = 1,
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
  # # Add vertical lines for 80% power threshold if available
  # if (!is.null(k562_threshold_cells)) {
  #   p <- p + geom_vline(xintercept = k562_threshold_cells,
  #                       linetype = "dashed",
  #                       color = "firebrick",
  #                       alpha = 0.7)
  # 
  #   # Add label for K562 threshold
  #   p <- p + annotate("text",
  #                     x = k562_threshold_cells,
  #                     y = 95,  # Position near top of plot
  #                     label = paste0("K562: ", k562_threshold_cells, " cells"),
  #                     color = "firebrick",
  #                     hjust = -0.1,  # Adjust horizontal position
  #                     size = 3.5)
  # }
  # 
  # if (!is.null(gasperini_threshold_cells)) {
  #   p <- p + geom_vline(xintercept = gasperini_threshold_cells,
  #                       linetype = "dashed",
  #                       color = "darkblue",
  #                       alpha = 0.7)
  # 
  #   # Add label for Gasperini threshold
  #   p <- p + annotate("text",
  #                     x = gasperini_threshold_cells,
  #                     y = 90,  # Position near top of plot but below K562 label
  #                     label = paste0("Gasperini: ", gasperini_threshold_cells, " cells"),
  #                     color = "darkblue",
  #                     hjust = -0.1,  # Adjust horizontal position
  #                     size = 3.5)
  # }
  
  # Store plot in the list
  tpm_plots[[i]] <- p
  
  # Print the plot
  # print(p)
  
  # Log the number of pairs in each range and the 80% threshold points
  message(paste0("TPM range ", low, "-", high, ": ",
                 "K562 pairs = ", nrow(k562_filtered),
                 ", Gasperini pairs = ", nrow(gasperini_filtered)))
  
  if (!is.null(k562_threshold_cells)) {
    message(paste0("K562 reaches 80% power at ", k562_threshold_cells, " cells"))
  } else {
    message("K562 does not reach 80% power in this TPM range")
  }
  
  if (!is.null(gasperini_threshold_cells)) {
    message(paste0("Gasperini reaches 80% power at ", gasperini_threshold_cells, " cells"))
  } else {
    message("Gasperini does not reach 80% power in this TPM range")
  }
}


### PLOT ALL TPMS TOGETHER ===================================================

# Create empty dataframe to store all binned data
all_binned_data <- data.frame()
effect_size <- 20  # Change this value to analyze different effect sizes

# Process each TPM range and add to the combined dataset
for (i in seq_along(tpm_ranges)) {
  low <- tpm_ranges[[i]][1]
  high <- tpm_ranges[[i]][2]
  
  # Create a label for this TPM range
  tpm_label <- paste0(low, "-", high)
  
  # Construct the power column name based on effect size
  power_col <- paste0("PowerAtEffectSize", effect_size)
  
  # Filter data for current TPM range
  k562_filtered <- k562_processed %>% 
    filter(TPM > low & TPM < high) %>%
    mutate(source = "K562 DC TAP")
  
  gasperini_filtered <- gasperini_processed %>% 
    filter(TPM > low & TPM < high) %>%
    mutate(source = "Gasperini et al. 2019")
  
  # Skip this range if no data in either dataset
  if (nrow(k562_filtered) == 0 && nrow(gasperini_filtered) == 0) {
    message(paste0("No data for TPM range ", low, "-", high))
    next
  }
  
  # Bin the data to reduce noise - especially important for high cell counts
  k562_binned <- k562_filtered %>%
    mutate(
      bin_width = case_when(
        cells < 200 ~ 100,
        cells < 500 ~ 200,
        cells < 1000 ~ 300,
        TRUE ~ 400
      ),
      cells_binned = bin_width * round(cells / bin_width)
    ) %>%
    group_by(source, cells_binned) %>%
    summarise(
      power_value = mean(!!sym(power_col), na.rm = TRUE),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    filter(n_samples >= 3) %>%
    mutate(tpm_range = tpm_label,
           tpm_numeric = low)  # Store numeric value for color gradient
  
  gasperini_binned <- gasperini_filtered %>%
    mutate(
      bin_width = case_when(
        cells < 200 ~ 100,
        cells < 500 ~ 200,
        cells < 1000 ~ 300,
        cells < 2000 ~ 600,
        TRUE ~ 3000
      ),
      cells_binned = bin_width * round(cells / bin_width)
    ) %>%
    group_by(source, cells_binned) %>%
    summarise(
      power_value = mean(!!sym(power_col), na.rm = TRUE),
      n_samples = n(),
      .groups = 'drop'
    ) %>%
    filter(n_samples >= 3) %>%
    mutate(tpm_range = tpm_label,
           tpm_numeric = low)  # Store numeric value for color gradient
  
  # Add to combined dataset
  all_binned_data <- bind_rows(all_binned_data, k562_binned, gasperini_binned)
  
  # Log the number of pairs in each range
  message(paste0("TPM range ", low, "-", high, ": ",
                 "K562 pairs = ", nrow(k562_filtered),
                 ", Gasperini pairs = ", nrow(gasperini_filtered)))
}

# Make sure tpm_range is a factor with levels in the right order
all_binned_data <- all_binned_data %>%
  mutate(tpm_range = factor(tpm_range, 
                            levels = c("5-15", "15-25", "25-35", "35-45", "45-55")))

# Create separate dataframes for each source to control colors independently
k562_data <- all_binned_data %>% filter(source == "K562 DC TAP")
gasperini_data <- all_binned_data %>% filter(source == "Gasperini et al. 2019")


# Create an alternative version with a color scale for each source
p_alt_colors <- ggplot(all_binned_data, 
                aes(x = cells_binned, y = power_value * 100, 
                    color = interaction(source, tpm_range),
                    group = interaction(source, tpm_range))) +
  geom_line(size = 1.3) +
  geom_point() +
  geom_hline(yintercept = 80, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(
    name = "Dataset & TPM Range",
    values = c(
      # K562 red gradient (light to dark)
      "K562 DC TAP.5-15" = "#FFCCCC",
      "K562 DC TAP.15-25" = "#FF9999",
      "K562 DC TAP.25-35" = "#FF6666",
      "K562 DC TAP.35-45" = "#FF3333",
      "K562 DC TAP.45-55" = "#CC0000",
      
      # Gasperini blue gradient (light to dark)
      "Gasperini et al. 2019.5-15" = "#CCCCFF",
      "Gasperini et al. 2019.15-25" = "#9999FF",
      "Gasperini et al. 2019.25-35" = "#6666FF",
      "Gasperini et al. 2019.35-45" = "#3333FF",
      "Gasperini et al. 2019.45-55" = "#0000CC"
    )
  ) +
  labs(
    title = paste0("Power at ", effect_size, "% Effect Size by TPM Range"),
    x = "Number of Perturbed Cells",
    y = "Power (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    aspect.ratio = 1
  )

# Print the alternative plot
# print(p_alt_colors)

# Save the plots


### P_ALT W/ LINES INSTEAD OF VALUES ==========================================

# Create an improved version with separate color for source and linetype for TPM
p_alt_lines <- ggplot(all_binned_data, 
                aes(x = cells_binned, y = power_value * 100, 
                    color = source,         # Use source for main color distinction
                    linetype = tpm_range,   # Use TPM ranges for linetype (grayscale effect)
                    group = interaction(source, tpm_range))) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 80, linetype = "dashed", color = "gray30") +
  scale_x_log10() +
  # Just two colors for the datasets
  scale_color_manual(
    name = "Dataset",
    values = c(
      "K562 DC TAP" = "firebrick",
      "Gasperini et al. 2019" = "darkblue"
    )
  ) +
  # Grayscale effect through line patterns for TPM ranges - from dotted (low TPM) to solid (high TPM)
  scale_linetype_manual(
    name = "TPM Range",
    values = c(
      "5-15" = "dotted",
      "15-25" = "dotdash",
      "25-35" = "dashed", 
      "35-45" = "longdash",
      "45-55" = "solid"
    )
  ) +
  labs(
    title = paste0("Power at ", effect_size, "% Effect Size by TPM Range"),
    x = "Number of Perturbed Cells",
    y = "Power (%)"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    aspect.ratio = 1
  ) +
  # Order legends and make them more compact
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    linetype = guide_legend(order = 2)
  )

# Save the plots


### PLOT THE TPM RANGES BY EFFECT SIZE ========================================

# Set the base directory for output plots
base_dir <- dirname(snakemake@output$power_sim_results_line_plot)
message(paste0("Using base directory: ", base_dir))

# Define the effect sizes to analyze
effect_sizes <- c(10, 15, 20, 25, 50)

# Loop through each effect size
for (effect_size in effect_sizes) {
  message(paste0("Processing effect size: ", effect_size, "%"))
  
  # Create empty dataframe to store all binned data for this effect size
  all_binned_data <- data.frame()
  
  # Process each TPM range and add to the combined dataset for this effect size
  for (i in seq_along(tpm_ranges)) {
    low <- tpm_ranges[[i]][1]
    high <- tpm_ranges[[i]][2]
    
    # Create a label for this TPM range
    tpm_label <- paste0(low, "-", high)
    
    # Construct the power column name based on current effect size
    power_col <- paste0("PowerAtEffectSize", effect_size)
    
    # Filter data for current TPM range
    k562_filtered <- k562_processed %>% 
      filter(TPM > low & TPM < high) %>%
      mutate(source = "K562 DC TAP")
    
    gasperini_filtered <- gasperini_processed %>% 
      filter(TPM > low & TPM < high) %>%
      mutate(source = "Gasperini et al. 2019")
    
    # Skip this range if no data in either dataset
    if (nrow(k562_filtered) == 0 && nrow(gasperini_filtered) == 0) {
      message(paste0("No data for TPM range ", low, "-", high))
      next
    }
    
    # Bin the data to reduce noise - especially important for high cell counts
    k562_binned <- k562_filtered %>%
      mutate(
        bin_width = case_when(
          cells < 200 ~ 100,
          cells < 500 ~ 200,
          cells < 1000 ~ 300,
          TRUE ~ 400
        ),
        cells_binned = bin_width * round(cells / bin_width)
      ) %>%
      group_by(source, cells_binned) %>%
      summarise(
        power_value = mean(!!sym(power_col), na.rm = TRUE),
        n_samples = n(),
        .groups = 'drop'
      ) %>%
      filter(n_samples >= 3) %>%
      mutate(tpm_range = tpm_label,
             tpm_numeric = low)  # Store numeric value for color gradient
    
    gasperini_binned <- gasperini_filtered %>%
      mutate(
        bin_width = case_when(
          cells < 200 ~ 100,
          cells < 500 ~ 200,
          cells < 1000 ~ 300,
          cells < 2000 ~ 600,
          TRUE ~ 3000
        ),
        cells_binned = bin_width * round(cells / bin_width)
      ) %>%
      group_by(source, cells_binned) %>%
      summarise(
        power_value = mean(!!sym(power_col), na.rm = TRUE),
        n_samples = n(),
        .groups = 'drop'
      ) %>%
      filter(n_samples >= 3) %>%
      mutate(tpm_range = tpm_label,
             tpm_numeric = low)  # Store numeric value for color gradient
    
    # Add to combined dataset
    all_binned_data <- bind_rows(all_binned_data, k562_binned, gasperini_binned)
    
    # Log the number of pairs in each range
    message(paste0("TPM range ", low, "-", high, ": ",
                   "K562 pairs = ", nrow(k562_filtered),
                   ", Gasperini pairs = ", nrow(gasperini_filtered)))
  }
  
  # Make sure tpm_range is a factor with levels in the right order
  all_binned_data <- all_binned_data %>%
    mutate(tpm_range = factor(tpm_range, 
                              levels = c("5-15", "15-25", "25-35", "35-45", "45-55")))
  
  # Create separate dataframes for each source to control colors independently
  k562_data <- all_binned_data %>% filter(source == "K562 DC TAP")
  gasperini_data <- all_binned_data %>% filter(source == "Gasperini et al. 2019")
  
  # Create p_alt_colors with a color scale for each source
  p_alt_colors <- ggplot(all_binned_data, 
                         aes(x = cells_binned, y = power_value * 100, 
                             color = interaction(source, tpm_range),
                             group = interaction(source, tpm_range))) +
    geom_line(size = 0.8) +
    geom_point(size = 1) +
    geom_hline(yintercept = 80, linetype = "dashed") +
    scale_x_log10() +
    scale_color_manual(
      name = "Dataset & TPM Range",
      values = c(
        # K562 red gradient (light to dark)
        "K562 DC TAP.5-15" = "#FFCCCC",
        "K562 DC TAP.15-25" = "#FF9999",
        "K562 DC TAP.25-35" = "#FF6666",
        "K562 DC TAP.35-45" = "#FF3333",
        "K562 DC TAP.45-55" = "#CC0000",
        
        # Gasperini blue gradient (light to dark)
        "Gasperini et al. 2019.5-15" = "#CCCCFF",
        "Gasperini et al. 2019.15-25" = "#9999FF",
        "Gasperini et al. 2019.25-35" = "#6666FF",
        "Gasperini et al. 2019.35-45" = "#3333FF",
        "Gasperini et al. 2019.45-55" = "#0000CC"
      )
    ) +
    labs(
      title = paste0("Power at ", effect_size, "% Effect Size"),
      x = "Number of Perturbed Cells",
      y = "Power (%)"
    ) +
    theme_classic() +
    ylim(0, 100) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      aspect.ratio = 1
    )
  
  # Create p_alt_lines with separate color for source and linetype for TPM
  p_alt_lines <- ggplot(all_binned_data, 
                        aes(x = cells_binned, y = power_value * 100, 
                            color = source,
                            linetype = tpm_range,
                            group = interaction(source, tpm_range))) +
    geom_line(size = 0.8) +
    geom_point(size = 1) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "gray30") +
    scale_x_log10() +
    # Just two colors for the datasets
    scale_color_manual(
      name = "Dataset",
      values = c(
        "K562 DC TAP" = "firebrick",
        "Gasperini et al. 2019" = "darkblue"
      )
    ) +
    # Grayscale effect through line patterns for TPM ranges - from dotted (low TPM) to solid (high TPM)
    scale_linetype_manual(
      name = "TPM Range",
      values = c(
        "5-15" = "dotted",
        "15-25" = "dotdash",
        "25-35" = "dashed", 
        "35-45" = "longdash",
        "45-55" = "solid"
      )
    ) +
    labs(
      title = paste0("Power at ", effect_size, "% Effect Size"),
      x = "Number of Perturbed Cells",
      y = "Power (%)"
    ) +
    ylim(0, 100) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      aspect.ratio = 1
    ) +
    # Order legends and make them more compact
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 3)),
      linetype = guide_legend(order = 2)
    )
  
  # Save the plots for this effect size
  # Create filenames with effect size included
  colors_filename <- file.path(base_dir, paste0("colors_effect_size_", effect_size, ".pdf"))
  lines_filename <- file.path(base_dir, paste0("lines_effect_size_", effect_size, ".pdf"))
  
  # Save the plots
  message(paste0("Saving plots for effect size ", effect_size, "%"))
  ggsave(plot = p_alt_colors, 
         filename = colors_filename,
         device = "pdf",
         height = 4,
         width = 4.5)
  ggsave(plot = p_alt_lines,
         filename = lines_filename,
         device = "pdf",
         height = 4,
         width = 4.5)
}


### POWER RESULTS W GAPSERINI (BAR) ===========================================

# Make same plot as first plot but for gasperini and k562 dctap


### COMPARISONS ===============================================================

# For each TPM range, how many perturbed cells are needed to match the power of DC TAP
# What's the increase in efficiency for detecting small effect sizes from these 


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)