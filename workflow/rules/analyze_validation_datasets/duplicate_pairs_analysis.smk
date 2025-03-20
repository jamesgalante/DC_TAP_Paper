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

# This rule is to process the gasperini dataset with Sceptre to see if the MAST and Sceptre Effect sizes are correlated
rule analyze_gasperini_with_sceptre:
  input:
    raw_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/raw_counts.rds",
    binarized_guide_counts = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/binarized_guide_counts.rds",
    response_id_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/response_id_target_pairs.tsv",
    guide_target_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/guide_target_pairs.tsv"
  output:
    discovery_results = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_analysis/results_run_discovery_analysis.rds"
  log: "results/analyze_validation_datasets/duplicate_pairs_analysis/logs/analyze_gasperini_with_sceptre.log"
  conda:
    "../../envs/sceptre_power_simulations.yml"
  resources:
    mem = "72G",
    time = "6:00:00"
  script:
    "../../scripts/analyze_validation_datasets/duplicate_pairs_analysis/analyze_gasperini_with_sceptre.R"
    
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
    downloaded_annot = "resources/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_data/gencode.v26lift37.annotation.gtf.gz" # Just for getting the gene names
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

# Analyze the duplicate pairs between DC TAP Seq and Gasperini
rule find_duplicate_pairs:
  input:
    gasperini_MAST_and_Sceptre = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_analysis/gasperini_MAST_and_Sceptre.rds",
    combined_validation = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_expt_pred_merged_annot.txt"
  output:
    gasperini_sceptre_v_mast_comparison = "results/analyze_validation_datasets/duplicate_pairs_analysis/gasperini_sceptre_v_mast_comparison.pdf",
    comparing_all_duplicate_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/comparing_all_duplicate_pairs.pdf",
    comparing_valid_duplicate_pairs = "results/analyze_validation_datasets/duplicate_pairs_analysis/comparing_valid_duplicate_pairs.pdf",
    comparing_all_duplicate_pairs_with_color = "results/analyze_validation_datasets/duplicate_pairs_analysis/comparing_all_duplicate_pairs_with_color.pdf"
  conda: 
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "48G",
    time = "2:00:00"
  script:
    "../../scripts/analyze_validation_datasets/duplicate_pairs_analysis/find_duplicate_pairs.R" 

    














