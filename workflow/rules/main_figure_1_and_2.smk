
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
    effect_size_boxplot = "results/main_figure_1_and_2/effect_size_boxplot.pdf"
  log: "results/main_figure_1_and_2/logs/effect_size_plots.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../scripts/main_figure_1_and_2/effect_size_plots.R" 

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
    # DC-TAP-seq (MOI6) files
    moi6_mol_info = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/moi6_fixed/molecule_info.h5",
    moi6_filtered_h5 = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/moi6_fixed/filtered_feature_bc_matrix.h5",
    # CRISPRi Direct Capture files
    crisprdi_mol_info = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/crisprdi_direct_capture/molecule_info.h5",
    crisprdi_filtered_h5 = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/crisprdi_direct_capture/filtered_feature_bc_matrix.h5",
    # DC-TAP-seq (MOI3) files
    moi3_mol_info = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/moi3_fixed/molecule_info.h5",
    moi3_filtered_h5 = "resources/main_figure_1_and_2/downsampling_reads_for_UMIs/moi3_fixed/filtered_feature_bc_matrix.h5"
  output:
    plot_pdf = "results/main_figure_1_and_2/downsampling_reads_for_UMIs/read_umi_saturation_comparison.pdf",
    plot_png = "results/main_figure_1_and_2/downsampling_reads_for_UMIs/read_umi_saturation_comparison.png",
    combined_results = "results/main_figure_1_and_2/downsampling_reads_for_UMIs/combined_results.csv"
  log: "results/main_figure_1_and_2/logs/downsampling_reads_for_UMIs.log"
  conda:
    "../envs/python_env.yml"
  resources:
    mem = "32G",
    time = "4:00:00"
  script:
    "../scripts/main_figure_1_and_2/downsampling_reads_for_UMIs.py"

# Duplicates Analysis
# Because this analysis requires many rules, it was moved to duplicate_pairs_analysis.smk
# The output of this analysis is included in the rule below

# Main Figure 2
rule create_main_figure_1_and_2:
  input:
    # Figure 1 plot
    dc_tap_x_qpcr_fig1 = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/comparison_positive_controls.pdf",
    # Figure 2 plots
    k562_distance_by_es_fig2 = "results/main_figure_1_and_2/k562_distance_by_es.pdf",
    wtc11_distance_by_es_fig2 = "results/main_figure_1_and_2/wtc11_distance_by_es.pdf",
    effect_size_plot_fig2 = "results/main_figure_1_and_2/effect_size_boxplot.pdf",
    k562_replicates_fig2 = "results/main_figure_1_and_2/compare_replicates/k562_plot.pdf", 
    wtc11_replicates_fig2 = "results/main_figure_1_and_2/compare_replicates/wtc11_plot.pdf", 
    dc_tap_x_qpcr_fig2 = "results/main_figure_1_and_2/compare_k562_effect_sizes_to_qPCR/comparison_others.pdf",
    duplicates_plot_fig2 = "results/main_figure_1_and_2/duplicate_pairs_analysis/final_correlation_plot.pdf"
  output:
    touch("results/main_figure_1_and_2/create_main_figure_1_and_2.done")
