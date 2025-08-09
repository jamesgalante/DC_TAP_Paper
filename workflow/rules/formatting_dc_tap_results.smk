

# Download the TSS bed file to add ABC promoter annotations
rule download_TSS_500bp_bed_file:
  output: "results/genome_annotation_files/CollapsedGeneBounds.hg38.TSS500bp.bed"
  params:
    url = config["benchmark_validation_datasets"]["download_urls"]["tss_bed"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# Run differential expression with sceptre's dev branch to get confidence intervals
rule sceptre_dev_differential_expression_w_confidence_intervals:
  input:
    sceptre_diffex_input = "results/process_validation_datasets/{sample}/differential_expression/sceptre_diffex_input.rds"
  output:
    discovery_results = "results/formatted_dc_tap_results/differential_expression_w_confidence_intervals_{sample}/results_run_discovery_analysis.rds",
    final_sceptre_object = "results/formatted_dc_tap_results/differential_expression_w_confidence_intervals_{sample}/final_sceptre_object.rds"
  log: "results/formatted_dc_tap_results/logs/sceptre_dev_differential_expression_w_confidence_intervals_{sample}.log"
  conda:
    "../envs/sceptre_dev_for_CIs.yml"
  resources:
    mem = "32G",
    time = "12:00:00"
  script:
    "../scripts/format_dc_tap_results/sceptre_dev_differential_expression_w_confidence_intervals.R"

# Rule to get numbers for the paper
rule adding_design_file_information:
  input:
    combined_validation = "results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz",
    guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    results_with_design_file_features = "results/formatted_dc_tap_results/results_with_design_file_features.tsv"
  log: "results/formatted_dc_tap_results/logs/adding_design_file_information.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/adding_design_file_information.R"

# Add more categories for understanding element overlap with different genomic features
rule adding_genomic_feature_overlaps:
  input:
    results_with_design_file_features = "results/formatted_dc_tap_results/results_with_design_file_features.tsv",
    abc_canonical_tss = "results/genome_annotation_files/CollapsedGeneBounds.hg38.TSS500bp.bed"
  output:
    results_with_design_file_and_genomic_feature_overlaps = "results/formatted_dc_tap_results/results_with_design_file_and_genomic_feature_overlaps.tsv"
  log: "results/formatted_dc_tap_results/logs/adding_genomic_feature_overlaps.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/adding_genomic_feature_overlaps.R"
    
# Creating categories to define the Random Set, Promoters, Valid Distal Element Gene pairs, etc.
rule adding_element_gene_pair_categories:
  input:
    results_with_design_file_and_genomic_feature_overlaps = "results/formatted_dc_tap_results/results_with_design_file_and_genomic_feature_overlaps.tsv",
    guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    create_ensemble_encode_input = expand("results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    results_with_element_gene_pair_categories = "results/formatted_dc_tap_results/results_with_element_gene_pair_categories.tsv"
  log: "results/formatted_dc_tap_results/logs/adding_element_gene_pair_categories.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/adding_element_gene_pair_categories.R"

rule modify_specific_pairs_in_final_file:
  input:
    results_with_element_gene_pair_categories = "results/formatted_dc_tap_results/results_with_element_gene_pair_categories.tsv",
    discovery_results_w_CIs = expand("results/formatted_dc_tap_results/differential_expression_w_confidence_intervals_{sample}/results_run_discovery_analysis.rds", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"])
  output:
    results_with_element_gene_pair_categories_modified = "results/formatted_dc_tap_results/results_with_element_gene_pair_categories_modified.tsv"
  log: "results/formatted_dc_tap_results/logs/modify_specific_pairs_in_final_file.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/modify_specific_pairs_in_final_file.R"

# Add the significance results with FDR correction on non-positive controls
rule add_results_wo_positive_controls:
  input:
    results_with_element_gene_pair_categories_modified = "results/formatted_dc_tap_results/results_with_element_gene_pair_categories_modified.tsv",
    combined_power_analysis_output_K562 = expand("results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/combined_power_analysis_output_es_{es}.tsv", es = [0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.50]),
    combined_power_analysis_output_WTC11 = expand("results/process_validation_datasets/WTC11_DC_TAP_Seq/power_analysis/combined_power_analysis_output_es_{es}.tsv", es = [0.02, 0.03, 0.05, 0.10, 0.15, 0.20, 0.25, 0.50])
  output:
    results_wo_pos_controls = "results/formatted_dc_tap_results/results_wo_pos_controls.tsv"
  log: "results/formatted_dc_tap_results/logs/add_results_wo_positive_controls.log"
  conda:
    "../envs/sceptre_dev_for_CIs.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/add_results_wo_positive_controls.R"

# Resize the elements to 500bp and merge overlapping elements for compatibility with chromatin category overlap pipeline
rule resize_and_merge_dc_tap_elements_for_chromatin_categories:
  input:
    results_wo_pos_controls = "results/formatted_dc_tap_results/results_wo_pos_controls.tsv",
  output:
    resized_and_merged_input_for_chromatin_categorization_pipeline = "results/formatted_dc_tap_results/resized_and_merged_input_for_chromatin_categorization_pipeline.tsv"
  log: "results/formatted_dc_tap_results/logs/resize_and_merge_dc_tap_elements.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/resize_and_merge_dc_tap_elements.R"
  
# Add Maya's epigenetic categories for each pair
rule add_element_chromatin_categories:
  input:
    resized_and_merged_input_for_chromatin_categorization_pipeline = "results/formatted_dc_tap_results/resized_and_merged_input_for_chromatin_categorization_pipeline.tsv",
    categorized_data = "resources/formatting_dc_tap_results/all_dc_tap.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.tsv.gz"
  output:
    Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements = "results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv"
  log: "results/formatted_dc_tap_results/logs/add_element_epigenetic_categories.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/add_element_epigenetic_categories.R"
