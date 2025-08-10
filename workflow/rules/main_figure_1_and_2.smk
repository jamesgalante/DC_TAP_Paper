
# Run sceptre differential expression with "singleton"
rule sceptre_singleton_differential_expression:
  input:
    sceptre_diffex_input = "results/process_validation_datasets/{sample}/differential_expression/sceptre_diffex_input.rds",
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds"
  output:
    discovery_results = "results/process_validation_datasets/{sample}/singleton_differential_expression/results_run_discovery_analysis.rds",
    final_sceptre_object = "results/process_validation_datasets/{sample}/singleton_differential_expression/final_sceptre_object.rds"
  log: "results/process_validation_datasets/{sample}/singleton_differential_expression/singleton_discovery_results.log"
  conda:
    "../envs/sceptre_dev_for_CIs.yml"
  resources:
    mem = "32G",
    time = "12:00:00"
  script:
    "../scripts/process_validation_datasets/sceptre_power_analysis/sceptre_singleton_differential_expression.R"

# Create distance by effect size plots
rule distance_by_effect_size_plots:
  input:
    results_with_element_gene_pair_categories_modified = "results/formatted_dc_tap_results/results_with_element_gene_pair_categories_modified.tsv"
  output:
    k562_distance_by_es = "results/main_figure_1_and_2/k562_distance_by_es.pdf",
    wtc11_distance_by_es = "results/main_figure_1_and_2/wtc11_distance_by_es.pdf"
  log: "results/main_figure_1_and_2/logs/distance_by_effect_size_plots.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/distance_by_effect_size_plots.R"

# Effect Size Plots
rule effect_size_plots:
  input:
    results_with_element_gene_pair_categories_modified = "results/formatted_dc_tap_results/results_with_element_gene_pair_categories_modified.tsv"
  output:
    effect_size_boxplot = "results/main_figure_1_and_2/effect_size_boxplot.pdf",
    positive_control_self_promoters_effect_size_boxplot = "results/main_figure_1_and_2/positive_control_self_promoters_effect_size_boxplot.pdf"
  log: "results/main_figure_1_and_2/logs/effect_size_plots.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/effect_size_plots.R" 
    
# Group the non-targeting guides under one grna_group (not as a negative control) to get the se_fold_change of all negative controls together
rule negative_targeting_confidence_interval:
  input:
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/{sample}/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/{sample}/metadata.rds",
    dge = "results/process_validation_datasets/{sample}/raw_counts/dge.rds",
    perturb_status = "results/process_validation_datasets/{sample}/raw_counts/perturb_status.rds"
  output:
    discovery_results = "results/main_figure_1_and_2/negative_controls_differential_expression_w_confidence_intervals_{sample}/results_run_discovery_analysis.rds",
    final_sceptre_object = "results/main_figure_1_and_2/negative_controls_differential_expression_w_confidence_intervals_{sample}/final_sceptre_object.rds"
  log: "results/formatted_dc_tap_results/logs/negative_controls_differential_expression_w_confidence_intervals_{sample}.log"
  conda:
    "../envs/sceptre_dev_for_CIs.yml"
  resources:
    mem = "32G",
    time = "12:00:00"
  script:
    "../scripts/main_figure_1_and_2/negative_controls_differential_expression_w_confidence_intervals.R"
    
# Plot for SLC2A3|chr12:7960676-7960977 Effect Sizes compared to Negative Controls
rule slc2a3_effect_size_plot:
  input:
    results_with_element_gene_pair_categories_modified = "results/formatted_dc_tap_results/results_with_element_gene_pair_categories_modified.tsv",
    wtc11_calibration_check_results = "results/process_validation_datasets/WTC11_DC_TAP_Seq/singleton_differential_expression/results_run_calibration_check.rds",
    wtc11_singleton_results = "results/process_validation_datasets/WTC11_DC_TAP_Seq/singleton_differential_expression/results_run_discovery_analysis.rds",
    negative_controls_together = "results/main_figure_1_and_2/negative_controls_differential_expression_w_confidence_intervals_WTC11_DC_TAP_Seq/results_run_discovery_analysis.rds"
  output:
    slc2a3_effect_size_barplot = "results/main_figure_1_and_2/slc2a3_effect_size_barplot.pdf"
  log: "results/main_figure_1_and_2/logs/slc2a3_effect_size_plot.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/slc2a3_effect_size_plot.R" 
    
# Replicates Analysis
rule compare_replicates:
  input:
    dge = "results/process_validation_datasets/{sample}/raw_counts/dge.rds",
    perturb_status = "results/process_validation_datasets/{sample}/raw_counts/perturb_status.rds",
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/{sample}/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/{sample}/metadata.rds",
  output:
    discovery_results1 = "results/main_figure_1_and_2/compare_replicates/{sample}_1/results_run_discovery_analysis.rds",
    final_sceptre_object1 = "results/main_figure_1_and_2/compare_replicates/{sample}_1/final_sceptre_object.rds",
    discovery_results2 = "results/main_figure_1_and_2/compare_replicates/{sample}_2/results_run_discovery_analysis.rds",
    final_sceptre_object2 = "results/main_figure_1_and_2/compare_replicates/{sample}_2/final_sceptre_object.rds"
  log: "results/main_figure_1_and_2/logs/compare_replicates_{sample}.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/compare_replicates.R"

# Make plots for the replicates analysis
rule make_replicates_plots:
  input:
    discovery_results1 = expand("results/main_figure_1_and_2/compare_replicates/{sample}_1/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    discovery_results2 = expand("results/main_figure_1_and_2/compare_replicates/{sample}_2/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    k562_plot = "results/main_figure_1_and_2/compare_replicates/k562_plot.pdf",
    wtc11_plot = "results/main_figure_1_and_2/compare_replicates/wtc11_plot.pdf",
    k562_neg_plot = "results/main_figure_1_and_2/compare_replicates/k562_neg_plot.pdf",
    wtc11_neg_plot = "results/main_figure_1_and_2/compare_replicates/wtc11_neg_plot.pdf"
  log: "results/main_figure_1_and_2/logs/make_replicates_plots.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/make_replicates_plots.R"

# DC TAP Seq x qPCR Effect Size Comparison
rule compare_k562_effect_sizes_to_qPCR:
  input:
    cell_ranger_outputs = "resources/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR",
    qPCR_results = "resources/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/qPCR_results.tsv",
    extra_qPCRs = "resources/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/extra_qPCR.tsv",
    k562_singleton_diffex_results = "results/process_validation_datasets/K562_DC_TAP_Seq/singleton_differential_expression/results_run_discovery_analysis.rds"
  output:
    results_run_discovery_analysis = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/results_run_discovery_analysis.rds",
    comparison_of_all_by_type = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/comparison_of_all_by_type.pdf",
    comparison_positive_controls = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/comparison_positive_controls.pdf",
    comparison_others = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/comparison_others.pdf",
    formatted_qPCR_Sceptre_table = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/formatted_qPCR_Sceptre_table.tsv",
    grna_target_data_frame_filtered_guides = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/grna_target_data_frame_filtered_guides.tsv"
  log: "results/main_figure_1_and_2/logs/compare_k562_effect_sizes_to_qPCR.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR.R"   
    
# Downsampling Analysis Rule
rule downsampling_reads_for_UMIs:
  input:
    # CRISPRi Direct Capture files
    crisprdi_mol_info = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/crisprdi_direct_capture/molecule_info.h5",
    crisprdi_filtered_h5 = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/crisprdi_direct_capture/filtered_feature_bc_matrix.h5",
    # DC-TAP-seq (MOI3) files
    moi3_mol_info = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/moi3_fixed/molecule_info.h5",
    moi3_filtered_h5 = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/moi3_fixed/filtered_feature_bc_matrix.h5"
  output:
    plot_pdf = "results/main_figure_1_and_2/downsampling_reads_for_UMIs/read_umi_saturation_comparison.pdf",
    combined_results = "results/main_figure_1_and_2/downsampling_reads_for_UMIs/combined_results.csv"
  log: "results/main_figure_1_and_2/logs/downsampling_reads_for_UMIs.log"
  conda:
    "../envs/downsampling.yml"
  resources:
    mem = "32G",
    time = "4:00:00"
  script:
    "../scripts/main_figure_1_and_2/downsampling_reads_for_UMIs.R"

# Duplicates Analysis
# Because this analysis requires many rules, most of the analysis was moved to duplicate_pairs_analysis.smk
rule find_duplicate_pairs:
  input:
    gasperini_results = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv",
    dc_tap_results = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv"
  output:
    duplicate_pairs_plot = "results/main_figure_1_and_2/duplicate_pairs_analysis/duplicate_pairs.pdf"
  log: "results/main_figure_1_and_2/logs/find_duplicate_pairs.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/find_duplicate_pairs.R"

rule wtc11_multi_moi:
  input:
    wtc11_multi_moi = "resources/main_figure_1_and_2/wtc11_multi_moi/wtc11_all_moi_mast_tab.txt"
  output:
    wtc11_multi_moi_plot = "results/main_figure_1_and_2/wtc11_multi_moi.pdf"
  log: "results/main_figure_1_and_2/logs/wtc11_multi_moi.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/wtc11_multi_moi.R"
