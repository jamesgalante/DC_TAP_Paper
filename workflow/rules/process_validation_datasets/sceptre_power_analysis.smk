# Snakemake rules to run power analysis using Sceptre

# Run sceptre differential expression with "union"
rule sceptre_differential_expression:
  input:
    sceptre_diffex_input = "results/process_validation_datasets/{sample}/differential_expression/sceptre_diffex_input.rds"
  output:
    discovery_results = "results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds",
    final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression/final_sceptre_object.rds"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/sceptre_differential_expression_{sample}.log"
  conda:
    "../../envs/all_packages.yml"
  resources:
    mem = "32G",
    time = "12:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_power_analysis/sceptre_differential_expression.R"
    
# # Run sceptre differential expression with "mixture"
# rule sceptre_differential_expression_w_mixture:
#     input:
#       sceptre_diffex_input = "results/process_validation_datasets/{sample}/differential_expression/sceptre_diffex_input.rds"
#     output:
#       discovery_results = "results/process_validation_datasets/{sample}/differential_expression_w_mixture/results_run_discovery_analysis.rds",
#       final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression_w_mixture/final_sceptre_object.rds"
#     log: "results/process_validation_datasets/sceptre_power_analysis/logs/sceptre_differential_expression_w_mixture_{sample}.log"
#     conda:
#       "../../envs/all_packages.yml"
#     resources:
#       mem = "32G",
#       time = "12:00:00"
#     script:
#       "../../scripts/process_validation_datasets/sceptre_power_analysis/sceptre_differential_expression_w_mixture.R"
# 
# # Run sceptre differential expression include all the negative controls tested against all elements
# rule sceptre_differential_expression_w_negative_controls:
#   input:
#     sceptre_diffex_input = "results/process_validation_datasets/{sample}/differential_expression/sceptre_diffex_input.rds",
#     gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds",
#     gRNA_groups_table = "results/process_validation_datasets/{sample}/gRNA_groups_table.rds",
#     features = "resources/process_validation_datasets/sceptre_setup/{sample}/cell_ranger_output/features.tsv.gz"
#   output:
#     discovery_results = "results/process_validation_datasets/{sample}/differential_expression_w_negative_controls/results_run_discovery_analysis.rds",
#     final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression_w_negative_controls/final_sceptre_object.rds"
#   params:
#     negative_control_genes = lambda wildcards: config["process_validation_datasets"]["sceptre_setup"]["negative_controls"][wildcards.sample]
#   log: "results/process_validation_datasets/{sample}/differential_expression_w_negative_controls/sceptre_differential_expression_w_negative_controls.log"
#   conda:
#     "../../envs/all_packages.yml"
#   resources:
#     mem = "32G",
#     time = "12:00:00"
#   script:
#     "../../scripts/process_validation_datasets/sceptre_power_analysis/sceptre_diffex_with_negative_controls.R"
# 
# # Run sceptre differential expression with "singleton"
# rule sceptre_singleton_differential_expression:
#   input:
#     sceptre_diffex_input = "results/process_validation_datasets/{sample}/differential_expression/sceptre_diffex_input.rds",
#     gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds"
#   output:
#     discovery_results = "results/process_validation_datasets/{sample}/singleton_differential_expression/results_run_discovery_analysis.rds",
#     final_sceptre_object = "results/process_validation_datasets/{sample}/singleton_differential_expression/final_sceptre_object.rds"
#   log: "results/process_validation_datasets/{sample}/singleton_differential_expression/singleton_discovery_results.log"
#   conda:
#     "../../envs/all_packages.yml"
#   resources:
#     mem = "32G",
#     time = "12:00:00"
#   script:
#     "../../scripts/process_validation_datasets/sceptre_power_analysis/sceptre_singleton_differential_expression.R"

rule create_sce:
  input:
    final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression/final_sceptre_object.rds"
  output:
    perturb_sce = "results/process_validation_datasets/{sample}/perturb_sce.rds"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/create_sce_{sample}.log"
  conda: "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "64G",
    time = "4:00:00"
  script:
     "../../scripts/process_validation_datasets/sceptre_power_analysis/create_sce_object.R"
    
# Create and split the discovery pairs file
N_BATCHES = 100
rule split_target_response_pairs:
  input:
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds"
  output:
    splits = expand("results/process_validation_datasets/{{sample}}/pair_splits/gene_gRNA_group_pairs_{split}.txt", split = range(1, N_BATCHES + 1)) # Default to 100 splits - this only works if every dataset has more than 100 genes, which is true (error will be thrown if not)
  params:
    batches = N_BATCHES
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/split_target_response_pairs_{sample}.log"
  conda:
    "../../envs/sceptre_power_simulations.yml"
  resources:
    mem = "24G",
    time = "2:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_power_analysis/split_target_response_pairs.R"
    

# Run the power simulation with sceptre for each split
N_REPS = 100
rule sceptre_power_analysis:
  input:
    gene_gRNA_group_pairs_split = "results/process_validation_datasets/{sample}/pair_splits/gene_gRNA_group_pairs_{split}.txt",
    final_sceptre_object = "results/process_validation_datasets/{sample}/differential_expression/final_sceptre_object.rds",
    gRNA_groups_table = "results/process_validation_datasets/{sample}/gRNA_groups_table.rds",
    perturb_sce = "results/process_validation_datasets/{sample}/perturb_sce.rds"
  output:
    power_analysis_output = "results/process_validation_datasets/{sample}/power_analysis/effect_size_{effect_size}/power_analysis_output_{split}.tsv"
  params:
    reps = N_REPS
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/sceptre_power_analysis_es{effect_size}_split{split}_{sample}.log"
  conda:
    "../../envs/sceptre_power_simulations.yml"
  resources:
    mem = "8G",
    time = "1:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_power_analysis/sceptre_power_analysis.R"

# Combine the split outputs of the power analysis
rule combine_sceptre_power_analysis:
 input:
   splits = expand("results/process_validation_datasets/{{sample}}/power_analysis/effect_size_{{effect_size}}/power_analysis_output_{split}.tsv", split = range(1, N_BATCHES + 1)),
 output:
   combined_power_analysis_output = "results/process_validation_datasets/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv"
 log: "results/process_validation_datasets/sceptre_power_analysis/logs/combine_sceptre_power_analysis_es{effect_size}_{sample}.log"
 conda:
   "../../envs/sceptre_power_simulations.yml"
 resources:
   mem = "32G",
   time = "2:00:00"
 script:
   "../../scripts/process_validation_datasets/sceptre_power_analysis/combine_sceptre_power_analysis.R"

# Compute the power from the power simulations
rule compute_power_from_simulations:
  input:
    combined_power_analysis_output = "results/process_validation_datasets/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv",
    discovery_results = "results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds"
  output:
    power_analysis_results = "results/process_validation_datasets/{sample}/power_analysis/power_analysis_results_es_{effect_size}.tsv"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/compute_power_from_simulations_es{effect_size}_{sample}.log"
  conda:
    "../../envs/sceptre_power_simulations.yml"
  resources:
    mem = "24G",
    time = "1:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_power_analysis/compute_power_from_simulations.R"

# Format sceptre output for compatibility with ENCODE pipelines
rule format_sceptre_output:
  input:
    power_analysis_results = expand("results/process_validation_datasets/{{sample}}/power_analysis/power_analysis_results_es_{effect_size}.tsv", effect_size = [0.10, 0.15, 0.20, 0.25, 0.50]),
    discovery_results = "results/process_validation_datasets/{sample}/differential_expression/results_run_discovery_analysis.rds",
    gene_gRNA_group_pairs = "results/process_validation_datasets/{sample}/gene_gRNA_group_pairs.rds",
    distances = "results/process_validation_datasets/{sample}/distances.tsv",
    guide_targets = "results/process_validation_datasets/{sample}/guide_targets.tsv"
  output:
    final_output = "results/process_validation_datasets/{sample}/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv"
  log: "results/process_validation_datasets/sceptre_power_analysis/logs/format_sceptre_output_{sample}.log"
  conda: "../../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_power_analysis/format_sceptre_output.R"
