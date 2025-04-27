
# Download the TSS bed file to add ABC promoter annotations
rule download_TSS_500bp_bed_file:
  output: "results/genome_annotation_files/CollapsedGeneBounds.hg38.TSS500bp.bed"
  params:
    url = config["benchmark_validation_datasets"]["download_urls"]["tss_bed"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# Rule to get numbers for the paper
rule calculate_summary_statistics_for_screen:
  input:
    combined_validation = "results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz",
    guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    create_ensemble_encode_input = expand("results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    abc_canonical_tss = "results/genome_annotation_files/CollapsedGeneBounds.hg38.TSS500bp.bed"
  output:
    combined_joined_w_categories = "results/formatted_dc_tap_results/combined_joined_w_categories.tsv",
    summarized_categories = "results/formatted_dc_tap_results/summarized_categories.tsv",
    igvf_formatted_file = "results/formatted_dc_tap_results/igvf_formatted_file.tsv"
  log: "results/formatted_dc_tap_results/logs/calculate_summary_statistics_for_screen.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/calculate_summary_statistics_for_screen.R"
    
# Rule to deal with specific pairs
rule modify_specific_pairs_in_final_file:
  input:
    summarized_categories = "results/formatted_dc_tap_results/summarized_categories.tsv"
  output:
    Formatted_DC_TAP_Seq_Results = "results/formatted_dc_tap_results/Formatted_DC_TAP_Seq_Results.tsv",
    summary_K562 = "results/formatted_dc_tap_results/summary_K562.tsv",
    summary_WTC11 = "results/formatted_dc_tap_results/summary_WTC11.tsv"
  log: "results/formatted_dc_tap_results/logs/modify_specific_pairs_in_final_file.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/modify_specific_pairs_in_final_file.R"
    
# Resize DC TAP elements for epigenetic categories overlap
rule resize_dc_tap_elements:
  input:
    Formatted_DC_TAP_Seq_Results = "results/formatted_dc_tap_results/Formatted_DC_TAP_Seq_Results.tsv"
  output:
    resized_Formatted_DC_TAP_Seq_Results = "results/formatted_dc_tap_results/resized_Formatted_DC_TAP_Seq_Results.tsv"
  log: "results/formatted_dc_tap_results/logs/resize_dc_tap_elements.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/resize_dc_tap_elements.R"

# Add Maya's epigenetic categories for each pair
rule add_element_epigenetic_categories:
  input:
    resized_Formatted_DC_TAP_Seq_Results = "results/formatted_dc_tap_results/resized_Formatted_DC_TAP_Seq_Results.tsv",
    categorized_data = "resources/formatting_dc_tap_results/all_DC_TAP.h3k27me3_quantile_50.ratio_quantile_35.h3k27ac_quantiles_60_90.WTC11_as_WTC11.tsv"
  output:
    Formatted_DC_TAP_Seq_Results_w_Categories = "results/formatted_dc_tap_results/Formatted_DC_TAP_Seq_Results_w_Categories.tsv"
  log: "results/formatted_dc_tap_results/logs/add_element_epigenetic_categories.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/format_dc_tap_results/add_element_epigenetic_categories.R"
