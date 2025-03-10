# This rule file processes the expt_pred_merged_annot outputs for the subset_upsampling_analysis

# Rule to download BED files
rule download_bed_files:
  output:
    "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}.bed"
  params:
    lambda wildcards: config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]['download_bed_files'][wildcards.bam_type][wildcards.cell_type]
  resources:
    mem = "24G",
    time = "4:00:00"
  shell:
    """
    wget -O {output}.gz {params}
    gunzip -f {output}.gz
    """
    
    
# An intermediate file used in the ABC pipeline is used in downstream analysis for overlapping h3k27ac and dnase values with enhancers
# This file already existed for wtc11, k562
  # K562: /oak/stanford/groups/engreitz/Projects/Benchmarking/Revisions/Predictors/ABC/dnase_h3k27ac_intactHiC/K562/Neighborhoods/EnhancerList.txt
  # WTC11: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_temp/ENCODE_rE2G/results/2024_0612_WTC11_for_EJ/wtc11/Neighborhoods/EnhancerList.txt
# I then copied these intermediate files from the ABC pipeline into the resources folder in resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/ABC_outputs


### EQTL Questions
# CD14+ monocyte:
  # H3K27me3 experiment: https://www.encodeproject.org/experiments/ENCSR000ASK/
  # CTCF experiment: https://www.encodeproject.org/experiments/ENCSR000ATN/
  # EnhancerList: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_v1.0.0/ENCODE_rE2G/results/2024_0110_monocyte_T/CD14_pos_monocyte/Neighborhoods/EnhancerList.txt
# CD4+ T-cell:
  # H3K27me3 experiment: https://www.encodeproject.org/experiments/ENCSR043SBG/
  # CTCF experiment: https://www.encodeproject.org/experiments/ENCSR470KCE/
  # EnhancerList: /oak/stanford/groups/engreitz/Users/sheth/ENCODE_rE2G_v1.0.0/ENCODE_rE2G/results/2024_0110_monocyte_T/CD4_pos_T/Neighborhoods/EnhancerList.txt


rule overlap_h3k27me3_and_ctcf_with_enhancers:
  input:
    training_expt_pred_merged_annot = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/training_expt_pred_merged_annot.txt.gz",
    validation_expt_pred_merged_annot = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/validation_expt_pred_merged_annot.txt.gz",
    peak_bed_file = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/{bam_type}_peaks/{cell_type}.bed", bam_type = ["ctcf", "h3k27me3"], cell_type = ["k562", "wtc11"])
  output:
    training_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/h3k27me3_ctcf_overlap/training_expt_pred_merged_annot.txt",
    validation_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/h3k27me3_ctcf_overlap/validation_expt_pred_merged_annot.txt"
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/overlap_h3k27me3_and_ctcf_with_enhancers.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/overlap_h3k27me3_and_ctcf_with_enhancers.R"


rule overlap_dnase_and_h3k27ac_with_enhancers:
  input:
    training_expt_pred_merged_annot = "resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/training_expt_pred_merged_annot.txt.gz",
    validation_expt_pred_merged_annot = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/validation_expt_pred_merged_annot.txt.gz",
    enhancer_list = expand("resources/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/ABC_outputs/{cell_type}_EnhancerList.txt", cell_type = ["k562", "wtc11"])
  output:
    training_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/dnase_h3k27ac_overlap/training_expt_pred_merged_annot.txt",
    validation_overlap = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/dnase_h3k27ac_overlap/validation_expt_pred_merged_annot.txt"
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/overlap_dnase_and_h3k27ac_with_enhancers.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "4:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/overlap_dnase_and_h3k27ac_with_enhancers.R"


bam_types = ["h3k27me3_ctcf", "dnase_h3k27ac"]
rule combine_all_tracks:
  input:
    validation_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/validation_expt_pred_merged_annot.txt", bam_type = bam_types),
    training_overlaps = expand("results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/{bam_type}_overlap/training_expt_pred_merged_annot.txt", bam_type = bam_types)
  output:
    combined_validation = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_validation_expt_pred_merged_annot.txt",
    combined_training = "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/expt_pred_merged_annot/combined_training_expt_pred_merged_annot.txt"
  params:
    model_threshold = config["analyze_validation_datasets"]["process_bam_files_and_combined_w_EG_results"]["combine_all_tracks"]["model_threshold"]
  log: "results/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/logs/combine_all_tracks.log"
  conda:
    "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "16G",
    time = "1:00:00"
  script:
    "../../scripts/analyze_validation_datasets/process_bam_files_and_combined_w_EG_results/combine_all_tracks.R"  





