# Here we're analyzing duplicate pairs in the datasets and subsequently processing gasperini with Sceptre

# This rule is for taking the CRISPRi pipeline gasperini files and creating inputs for the sceptre analysis
rule create_gasperini_sceptre_inputs:
  input:
    perturb_sce = "resources/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/perturb_sce.rds",
    crispr_pipeline_output = "resources/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/output_0.13gStd_MAST_perCRE.tsv.gz"
  output:
    raw_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/raw_counts.rds",
    binarized_guide_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/binarized_guide_counts.rds",
    response_id_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/response_id_target_pairs.tsv",
    guide_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/guide_target_pairs.tsv"
  log: "results/analyze_validation_datasets/duplicate_pairs_analysis/logs/create_gasperini_sceptre_inputs.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "64G",
    time = "3:00:00"
  script:
    "../../scripts/analyze_validation_datasets/duplicate_pairs_analysis/create_gasperini_sceptre_inputs.R"


# # This rule is to process the gasperini dataset with Sceptre to see if the MAST and Sceptre Effect sizes are correlated
# rule analyze_gasperini_with_sceptre:
#   input:
#     raw_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/raw_counts.rds",
#     binarized_guide_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/binarized_guide_counts.rds",
#     response_id_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/response_id_target_pairs.tsv",
#     guide_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/guide_target_pairs.tsv"
#   output:
#     discovery_results = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_analysis/results_run_discovery_analysis.rds"
#   log: "results/analyze_validation_datasets/duplicate_pairs_analysis/logs/analyze_gasperini_with_sceptre.log"
#   conda:
#     "../../envs/sceptre_power_simulations.yml"
#   resources:
#     mem = "72G",
#     time = "6:00:00"
#   script:
#     "../../scripts/analyze_validation_datasets/duplicate_pairs_analysis/analyze_gasperini_with_sceptre.R"
    
# # Get the per guide effect size of the gasperini guides
rule gasperini_per_guide_effect_sizes:
  input:
    raw_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/raw_counts.rds",
    binarized_guide_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/binarized_guide_counts.rds",
    response_id_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/response_id_target_pairs.tsv",
    guide_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/guide_target_pairs.tsv"
  output:
    discovery_results = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_analysis/per_guide_diffex/results_run_discovery_analysis.rds"
  log: "results/analyze_validation_datasets/duplicate_pairs_analysis/logs/gasperini_per_guide_effect_sizes.log"
  conda:
    "../../envs/sceptre_power_simulations.yml"
  resources:
    mem = "72G",
    time = "15:00:00"
  script:
    "../../scripts/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_per_guide_effect_sizes.R"
    
    
# Download annotation file for getting gene symbols for each ensemble id in gasperini
rule download_gasperini_annotation_file:
  output:
    downloaded_annot = "resources/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/gencode.v26lift37.annotation.gtf.gz" 
  params:
    annot = config["analyze_validation_datasets"]["duplicate_pairs_analysis"]["download_gasperini_annotation_file"]["annot"]
  resources:
    mem = "16G",
    time = "1:00:00"
  shell:
    """
    wget -O {output.downloaded_annot} {params.annot}
    """

# Compare gasperini effect sizes between Sceptre and MAST
rule compare_gasperini_Sceptre_and_MAST:
  input:
    gasperini_sceptre_results = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_analysis/results_run_discovery_analysis.rds",
    combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt",
    annotation_file = "resources/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/gencode.v26lift37.annotation.gtf.gz" 
  output:
    gasperini_MAST_and_Sceptre = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_analysis/gasperini_MAST_and_Sceptre.rds"
  log: "results/analyze_validation_datasets/duplicate_pairs_analysis/logs/compare_gasperini_Sceptre_and_MAST.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/duplicate_pairs_analysis/compare_gasperini_Sceptre_and_MAST.R"

# This rule is to visualize the DC TAP dataset per guide effect sizes and compare with duplicate pairs in the training datasets
rule analyze_duplicate_k562_dc_tap_pairs:
  input:
    combined_unfiltered_k562_dc_tap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_unfiltered_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt",
    gasperini_MAST_and_Sceptre = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_analysis/gasperini_MAST_and_Sceptre.rds",
    per_guide_effect_sizes_unfiltered_k562_dc_tap = "resources/analyze_validation_datasets/duplicate_pairs_analysis/dc_tap_data/k562_dc_tap_per_guide_effect_sizes.txt",
    grna_target_table = "resources/analyze_validation_datasets/duplicate_pairs_analysis/dc_tap_data/gRNA_groups_table.txt",
    bonferroni_corrected_k562_dc_tap = "resources/analyze_validation_datasets/duplicate_pairs_analysis/dc_tap_data/k562_dc_tap_bonferroni_integration.rds",
    k562_dc_tap_discovery_results = "resources/analyze_validation_datasets/duplicate_pairs_analysis/dc_tap_data/k562_dc_tap_discovery_results.txt",
    k562_dc_tap_gene_mapping = "resources/analyze_validation_datasets/duplicate_pairs_analysis/dc_tap_data/k562_dc_tap_gene_mapping.tsv"
  output:
    "results/analyze_validation_datasets/duplicate_pairs_analysis/analyze_duplicate_k562_dc_tap_pairs.html"
  conda: 
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "48G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/duplicate_pairs_analysis/analyze_duplicate_k562_dc_tap_pairs.Rmd"


