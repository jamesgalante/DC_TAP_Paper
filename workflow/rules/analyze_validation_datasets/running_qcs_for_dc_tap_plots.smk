
# Rule to get numbers for the paper
rule calculate_summary_statistics_for_screen:
  input:
    combined_validation = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt",
    guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    summary_stats = "results/analyze_validation_datasets/dc_tap_plots/summary_stats.txt"
  params:
    k562_negative_control_genes = config["process_validation_datasets"]["sceptre_setup"]["negative_controls"]["K562_DC_TAP_Seq"],
    wtc11_negative_control_genes = config["process_validation_datasets"]["sceptre_setup"]["negative_controls"]["WTC11_DC_TAP_Seq"]
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/calculate_summary_statistics_for_screen.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/calculate_summary_statistics_for_screen.R"
    
# Rule to calculate pair drop out at different steps of the pipeline
rule understand_pair_drop_out:
  input:
    discovery_pairs = expand("results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    diffex_output = expand("results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    power_analysis_output = expand("results/process_validation_datasets/{sample}/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    create_encode_dataset_output = expand("results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_{sample}_0.13gStd_Sceptre_perCRE_hg19.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    liftover_output = expand("results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_{sample}_0.13gStd_Sceptre_perCRE_GRCh38.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    filter_crispr_dataset_output = expand("results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    create_ensemble_encode_output = "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz",
    create_ensemble_epbenchmarking_output = "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  output:
    drop_out_stats = "results/analyze_validation_datasets/dc_tap_plots/drop_out_stats.txt"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/understand_pair_drop_out.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/understand_pair_drop_out.R"
    
# Create volcano plots for the negative controls results
rule check_negative_control_results:
  input:
    discovery_results = expand("results/process_validation_datasets/{sample}/differential_expression_w_negative_controls/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    final_sceptre_object = expand("results/process_validation_datasets/{sample}/differential_expression_w_negative_controls/final_sceptre_object.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    features = expand("resources/process_validation_datasets/sceptre_setup/{sample}/cell_ranger_output/features.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    k562_volcano = "results/analyze_validation_datasets/dc_tap_plots/check_negative_control_results/k562_volcano.pdf",
    wtc11_volcano = "results/analyze_validation_datasets/dc_tap_plots/check_negative_control_results/wtc11_volcano.pdf"
  params:
    k562_negative_control_genes = config["process_validation_datasets"]["sceptre_setup"]["negative_controls"]["K562_DC_TAP_Seq"],
    wtc11_negative_control_genes = config["process_validation_datasets"]["sceptre_setup"]["negative_controls"]["WTC11_DC_TAP_Seq"]
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/check_negative_control_results.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/check_negative_control_results.R"
    
rule compare_k562_effect_sizes_to_qPCR:
  input:
    cell_ranger_outputs = "resources/analyze_validation_datasets/running_qcs",
    qPCR_results = "resources/analyze_validation_datasets/running_qcs/qPCR_results.tsv",
    extra_qPCRs = "resources/analyze_validation_datasets/running_qcs/extra_qPCR.tsv",
    k562_singleton_diffex_results = "results/process_validation_datasets/K562_DC_TAP_Seq/singleton_differential_expression/results_run_discovery_analysis.rds"
  output:
    results_run_discovery_analysis = "results/analyze_validation_datasets/dc_tap_plots/qPCR_qc/results_run_discovery_analysis.rds",
    comparison_of_all_by_type = "results/analyze_validation_datasets/dc_tap_plots/qPCR_qc/comparison_of_all_by_type.pdf",
    comparison_positive_controls = "results/analyze_validation_datasets/dc_tap_plots/qPCR_qc/comparison_positive_controls.pdf",
    comparison_others = "results/analyze_validation_datasets/dc_tap_plots/qPCR_qc/comparison_others.pdf",
    formatted_qPCR_Sceptre_table = "results/analyze_validation_datasets/dc_tap_plots/qPCR_qc/formatted_qPCR_Sceptre_table.tsv"
  params:
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/compare_k562_effect_sizes_to_qPCR.log"
  conda:
    "../../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/compare_k562_effect_sizes_to_qPCR.R"   
    
  
rule compare_replicates:
  input:
    dge = "results/process_validation_datasets/{sample}/raw_counts/dge.rds",
    perturb_status = "results/process_validation_datasets/{sample}/raw_counts/perturb_status.rds",
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/{sample}/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/{sample}/metadata.rds",
  output:
    discovery_results1 = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/{sample}_1/results_run_discovery_analysis.rds",
    final_sceptre_object1 = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/{sample}_1/final_sceptre_object.rds",
    discovery_results2 = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/{sample}_2/results_run_discovery_analysis.rds",
    final_sceptre_object2 = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/{sample}_2/final_sceptre_object.rds"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/compare_replicates_{sample}.log"
  conda:
    "../../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/compare_replicates.R"
  
rule make_replicates_plots:
  input:
    discovery_results1 = expand("results/analyze_validation_datasets/dc_tap_plots/compare_replicates/{sample}_1/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    discovery_results2 = expand("results/analyze_validation_datasets/dc_tap_plots/compare_replicates/{sample}_2/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    k562_plot = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/k562_plot.pdf",
    wtc11_plot = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/wtc11_plot.pdf",
    k562_neg_plot = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/k562_neg_plot.pdf",
    wtc11_neg_plot = "results/analyze_validation_datasets/dc_tap_plots/compare_replicates/wtc11_neg_plot.pdf"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/make_replicates_plots.log"
  conda:
    "../../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/make_replicates_plots.R"  
  
  


  
