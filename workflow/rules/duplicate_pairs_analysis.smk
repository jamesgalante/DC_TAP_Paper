# This script is for processing Gasperini et al. 2019 results with SCEPTRE and a SCEPTRE-based power simulation

# Get the raw counts and binarized grna counts from the perturb_sce object
rule process_perturb_sce_Gasperini:
  input:
    perturb_sce = "resources/main_figure_1_and_2/duplicate_pairs_analysis/perturb_sce.rds",
  output:
    raw_counts = "results/main_figure_1_and_2/duplicate_pairs_analysis/raw_counts.rds",
    binarized_guide_counts = "results/main_figure_1_and_2/duplicate_pairs_analysis/binarized_guide_counts.rds",
    guide_targets = "results/main_figure_1_and_2/duplicate_pairs_analysis/guide_targets.tsv"
  log: "results/main_figure_1_and_2/logs/process_perturb_sce_Gasperini.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "72G",
    time = "2:00:00"
  script:
    "../scripts/main_figure_1_and_2/process_perturb_sce_Gasperini.R"

# Create the SCEPTRE differential expression input object
rule create_sceptre_diffex_input_Gasperini:
  input:
    raw_counts = "results/main_figure_1_and_2/duplicate_pairs_analysis/raw_counts.rds",
    binarized_guide_counts = "results/main_figure_1_and_2/duplicate_pairs_analysis/binarized_guide_counts.rds",
    guide_targets = "results/main_figure_1_and_2/duplicate_pairs_analysis/guide_targets.tsv",
    annot = "results/genome_annotation_files/gencode.v26lift37.annotation.gtf.gz"
  output:
    gene_gRNA_group_pairs = "results/main_figure_1_and_2/{sample}/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/main_figure_1_and_2/{sample}/gRNA_groups_table.rds",
    metadata = "results/main_figure_1_and_2/{sample}/metadata.rds",
    sceptre_diffex_input = "results/main_figure_1_and_2/{sample}/differential_expression/sceptre_diffex_input.rds",
    distances = "results/main_figure_1_and_2/{sample}/distances.tsv"
  log: "results/main_figure_1_and_2/logs/create_sceptre_diffex_input_Gasperini_{sample}.log"
  conda:
    "../envs/all_packages.yml"
  resources:
    mem = "72G",
    time = "2:00:00"
  script:
    "../scripts/main_figure_1_and_2/create_sceptre_diffex_input_Gasperini.R"
    
# Run differential expression with sceptre's dev branch to get confidence intervals
rule sceptre_differential_expression_Gasperini:
  input:
    sceptre_diffex_input = "results/main_figure_1_and_2/{sample}/differential_expression/sceptre_diffex_input.rds"
  output:
    discovery_results = "results/main_figure_1_and_2/{sample}/differential_expression/results_run_discovery_analysis.rds",
    final_sceptre_object = "results/main_figure_1_and_2/{sample}/differential_expression/final_sceptre_object.rds"
  log: "results/main_figure_1_and_2/logs/sceptre_differential_expression_Gasperini_{sample}.log"
  conda:
    "../envs/sceptre_dev_for_CIs.yml"
  resources:
    mem = "64G",
    time = "12:00:00"
  script:
    "../scripts/main_figure_1_and_2/sceptre_differential_expression_Gasperini.R"

# Create the sce object from SCEPTRE object for simulations
rule create_sce_Gasperini:
  input:
    final_sceptre_object = "results/main_figure_1_and_2/{sample}/differential_expression/final_sceptre_object.rds"
  output:
    perturb_sce = "results/main_figure_1_and_2/{sample}/perturb_sce.rds"
  log: "results/main_figure_1_and_2/logs/create_sce_Gasperini_{sample}.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "128G",
    time = "4:00:00"
  script:
     "../scripts/process_validation_datasets/sceptre_power_analysis/create_sce_object.R"
    
# Create and split the discovery pairs file
G_BATCHES = 300 # There are around 5000 unique grna_groups, so this processes ~15 in each batch
rule split_target_response_pairs_Gasperini:
  input:
    gene_gRNA_group_pairs = "results/main_figure_1_and_2/{sample}/gene_gRNA_group_pairs.rds"
  output:
    splits = expand("results/main_figure_1_and_2/{{sample}}/pair_splits/gene_gRNA_group_pairs_{split}.txt", split = range(1, G_BATCHES + 1))
  params:
    batches = G_BATCHES
  log: "results/main_figure_1_and_2/logs/split_target_response_pairs_Gasperini_{sample}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "24G",
    time = "2:00:00"
  script:
    "../scripts/process_validation_datasets/sceptre_power_analysis/split_target_response_pairs.R"

# Run the power simulation with sceptre for each split
rule sceptre_power_analysis_Gasperini:
  input:
    gene_gRNA_group_pairs_split = "results/main_figure_1_and_2/{sample}/pair_splits/gene_gRNA_group_pairs_{split}.txt",
    final_sceptre_object = "results/main_figure_1_and_2/{sample}/differential_expression/final_sceptre_object.rds",
    gRNA_groups_table = "results/main_figure_1_and_2/{sample}/gRNA_groups_table.rds",
    perturb_sce = "results/main_figure_1_and_2/{sample}/perturb_sce.rds"
  output:
    power_analysis_output = "results/main_figure_1_and_2/{sample}/power_analysis/effect_size_{effect_size}/power_analysis_output_{split}.tsv"
  params:
    reps = config["process_validation_datasets"]["power_analysis"]["n_reps"]
  log: "results/main_figure_1_and_2/logs/sceptre_power_analysis_Gasperini_{sample}_es{effect_size}_split{split}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../scripts/process_validation_datasets/sceptre_power_analysis/sceptre_power_analysis.R"

# Combine the split outputs of the power analysis
rule combine_sceptre_power_analysis_Gasperini:
  input:
    splits = expand("results/main_figure_1_and_2/{{sample}}/power_analysis/effect_size_{{effect_size}}/power_analysis_output_{split}.tsv", split = range(1, G_BATCHES + 1))
  output:
    combined_power_analysis_output = "results/main_figure_1_and_2/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv"
  log: "results/main_figure_1_and_2/logs/combine_sceptre_power_analysis_Gasperini_{sample}_es{effect_size}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/process_validation_datasets/sceptre_power_analysis/combine_sceptre_power_analysis.R"

# Compute the power from the power simulations
rule compute_power_from_simulations_Gasperini:
  input:
    combined_power_analysis_output = "results/main_figure_1_and_2/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv",
    discovery_results = "results/main_figure_1_and_2/{sample}/differential_expression/results_run_discovery_analysis.rds"
  output:
    power_analysis_results = "results/main_figure_1_and_2/{sample}/power_analysis/power_analysis_results_es_{effect_size}.tsv"
  log: "results/main_figure_1_and_2/logs/compute_power_from_simulations_Gasperini_{sample}_es{effect_size}.log"
  conda:
    "../envs/sceptre_power_simulations.yml"
  resources:
    mem = "24G",
    time = "1:00:00"
  script:
    "../scripts/process_validation_datasets/sceptre_power_analysis/compute_power_from_simulations.R"

# Format sceptre output for compatibility with ENCODE pipelines
rule format_sceptre_output_Gasperini:
  input:
    power_analysis_results = expand("results/main_figure_1_and_2/{{sample}}/power_analysis/power_analysis_results_es_{effect_size}.tsv", effect_size = [0.10, 0.15, 0.20, 0.25, 0.50]),
    discovery_results = "results/main_figure_1_and_2/{sample}/differential_expression/results_run_discovery_analysis.rds",
    gene_gRNA_group_pairs = "results/main_figure_1_and_2/{sample}/gene_gRNA_group_pairs.rds",
    distances = "results/main_figure_1_and_2/{sample}/distances.tsv"
  output:
    final_output = "results/main_figure_1_and_2/{sample}/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv"
  log: "results/main_figure_1_and_2/logs/format_sceptre_output_Gasperini_{sample}.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/process_validation_datasets/sceptre_power_analysis/format_sceptre_output.R"

# Overall rule to run everything for the Gasperini dataset
rule run_gasperini_analysis:
  input:
    "results/main_figure_1_and_2/duplicate_pairs_analysis/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv"
  output:
    "results/main_figure_1_and_2/gasperini_analysis_complete.txt"
  shell:
    "touch {output}"
