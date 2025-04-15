
# Rule to get numbers for the paper
rule calculate_summary_statistics_for_screen:
  input:
    combined_validation = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/validation_expt_pred_merged_annot.txt.gz",
    guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    encode_training_dataset = "resources/analyze_validation_datasets/expt_pred_merged_annot/training_expt_pred_merged_annot.txt.gz",
    create_ensemble_encode_input = expand("results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    abc_canonical_tss = "resources/benchmark_validation_datasets/crispr_benchmarking/genome_annotations/CollapsedGeneBounds.hg38.TSS500bp.bed"
  output:
    combined_joined_w_categories = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/combined_joined_w_categories.tsv",
    summarized_categories = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summarized_categories.tsv",
    igvf_formatted_file = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/igvf_formatted_file.tsv",
    summary_K562 = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summary_K562.tsv",
    summary_WTC11 = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summary_WTC11.tsv"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/calculate_summary_statistics_for_screen.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/calculate_summary_statistics_for_screen.R"
    
# Rule to deal with specific pairs
rule modify_specific_pairs_in_final_file:
  input:
    summarized_categories = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/summarized_categories.tsv"
  output:
    Formatted_DC_TAP_Seq_Results = "results/analyze_validation_datasets/dc_tap_plots/summary_stats/Formatted_DC_TAP_Seq_Results.tsv"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/modify_specific_pairs_in_final_file.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/modify_specific_pairs_in_final_file.R"
    
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
  
rule igv_file_creation:
  input:
    combined_joined_w_categories = "results/analyze_validation_datasets/dc_tap_plots/combined_joined_w_categories.tsv",
    k562_guide_targets = "results/process_validation_datasets/K562_DC_TAP_Seq/guide_targets.tsv",
    wtc11_guide_targets = "results/process_validation_datasets/WTC11_DC_TAP_Seq/guide_targets.tsv",
    hg19ToHg38_chain_file = "resources/benchmark_validation_datasets/create_encode_output/hg19ToHg38.over.chain",
    annot = "resources/process_validation_datasets/sceptre_setup/genome_annotation_files/gencode.v32.annotation.gtf.gz",
    k562_gene_table = "resources/analyze_validation_datasets/igv_file_creation/k562_gene_table.txt",
    wtc11_gene_table = "resources/analyze_validation_datasets/igv_file_creation/wtc11_gene_table.txt"
  output:
    k562_guides = "results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_guide_targets_hg38.bed",
    k562_targets = "results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_target_regions_hg38.bed",
    k562_arcs = "results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_guide_target_arcs_hg38.bedpe",
    wtc11_guides = "results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_guide_targets_hg38.bed",
    wtc11_targets = "results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_target_regions_hg38.bed",
    wtc11_arcs = "results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_guide_target_arcs_hg38.bedpe",
    guide_color_map = "results/analyze_validation_datasets/dc_tap_plots/igv_files/guide_color_map.tsv",
    target_color_map = "results/analyze_validation_datasets/dc_tap_plots/igv_files/target_color_map.tsv",
    k562_TAP_seq_genes = "results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_TAP_seq_genes.bed",
    wtc11_TAP_seq_genes = "results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_TAP_seq_genes.bed",
    k562_signif_arcs = "results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_signif_arcs.bedpe",
    wtc11_signif_arcs = "results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_signif_arcs.bedpe"
  log: "results/analyze_validation_datasets/dc_tap_plots/logs/igv_file_creation.log"
  conda:
    "../../envs/all_packages.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/dc_tap_plots/igv_file_creation.R"  
    
# Lines to copy all file to mitra
# Copy K562 files
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_guide_targets_hg38.bed \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_K562/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_target_regions_hg38.bed \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_K562/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_guide_target_arcs_hg38.bedpe \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_K562/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_TAP_seq_genes.bed \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_K562/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/k562_signif_arcs.bedpe \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_K562/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/guide_color_map.tsv \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_K562/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/target_color_map.tsv \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_K562/

# Copy WTC11 files
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_guide_targets_hg38.bed \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_WTC11/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_target_regions_hg38.bed \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_WTC11/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_guide_target_arcs_hg38.bedpe \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_WTC11/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_TAP_seq_genes.bed \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_WTC11/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/wtc11_signif_arcs.bedpe \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_WTC11/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/guide_color_map.tsv \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_WTC11/
# cp results/analyze_validation_datasets/dc_tap_plots/igv_files/target_color_map.tsv \
#    /oak/stanford/groups/engreitz/public/RayJagoda2024/DC-TAPseq/new_WTC11/

  
# Mitra links:
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_K562/k562_guide_targets_hg38.bed
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_K562/k562_guide_target_arcs_hg38.bedpe
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_K562/k562_TAP_seq_genes.bed
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_K562/k562_signif_arcs.bedpe
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_K562/k562_target_regions_hg38.bed
# 
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_WTC11/wtc11_guide_targets_hg38.bed
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_WTC11/wtc11_guide_target_arcs_hg38.bedpe
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_WTC11/wtc11_TAP_seq_genes.bed
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_WTC11/wtc11_signif_arcs.bedpe
# https://mitra.stanford.edu/engreitz/oak/public/RayJagoda2024/DC-TAPseq/new_WTC11/wtc11_target_regions_hg38.bed
  
