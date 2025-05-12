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
    "../envs/all_packages.yml"
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

# # Compute the power from the power simulations
# rule compute_power_from_simulations_Gasperini:
#   input:
#     combined_power_analysis_output = "results/main_figure_1_and_2/{sample}/power_analysis/combined_power_analysis_output_es_{effect_size}.tsv",
#     discovery_results = "results/main_figure_1_and_2/{sample}/differential_expression/results_run_discovery_analysis.rds"
#   output:
#     power_analysis_results = "results/main_figure_1_and_2/{sample}/power_analysis/power_analysis_results_es_{effect_size}.tsv"
#   log: "results/main_figure_1_and_2/logs/compute_power_from_simulations_Gasperini_{sample}_es{effect_size}.log"
#   conda:
#     "../envs/sceptre_power_simulations.yml"
#   resources:
#     mem = "24G",
#     time = "1:00:00"
#   script:
#     "../scripts/process_validation_datasets/sceptre_power_analysis/compute_power_from_simulations.R"

# Format sceptre output for compatibility with ENCODE pipelines
rule format_sceptre_output_Gasperini:
  input:
    power_analysis_results = expand("results/main_figure_1_and_2/{{sample}}/power_analysis/power_analysis_results_es_{effect_size}.tsv", effect_size = [0.10, 0.15, 0.20, 0.25, 0.50]),
    discovery_results = "results/main_figure_1_and_2/{sample}/differential_expression/results_run_discovery_analysis.rds",
    gene_gRNA_group_pairs = "results/main_figure_1_and_2/{sample}/gene_gRNA_group_pairs.rds",
    guide_targets = "resources/main_figure_1_and_2/duplicate_pairs_analysis/guide_targets.tsv",
    distances = "results/main_figure_1_and_2/{sample}/distances.tsv"
  output:
    final_output = "results/main_figure_1_and_2/{sample}/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv"
  log: "results/main_figure_1_and_2/logs/format_sceptre_output_Gasperini_{sample}.log"
  conda: "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "5:00:00"
  script:
    "../scripts/main_figure_1_and_2/format_sceptre_output_Gasperini.R"
    
# compile output files in ENCODE format
rule create_encode_dataset_Gasperini:
  input:
    results = "results/main_figure_1_and_2/{sample}/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv",
    annot = "results/genome_annotation_files/gencode.v26lift37.annotation.gtf.gz",
    guide_targets = "resources/main_figure_1_and_2/duplicate_pairs_analysis/guide_targets.tsv"
  output: "results/main_figure_1_and_2/{sample}/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz"
  params:
    tss_min_dist = config["benchmark_validation_datasets"]["encode_datasets"]["dist_to_TSS"][0],
    gene_ids = "ensembl",
    tss_ctrl_tag = "TSSCtrl",
    padj_threshold = config["process_validation_datasets"]["differential_expression"]["padj_threshold"],
    reference = lambda wildcards: wildcards.sample
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    time = "2:00:00",
    mem = "32G"
  script:
    "../scripts/create_encode_output/create_encode_dataset.R"
    
# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_enhancers_Gasperini:
  input:
    results = "results/main_figure_1_and_2/{sample}/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz", 
    chain = "results/genome_annotation_files/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/main_figure_1_and_2/{sample}/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg19.bed",
    hg38 = "results/main_figure_1_and_2/{sample}/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    unlifted = "results/main_figure_1_and_2/{sample}/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_unlifted.bed"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    time = "1:00:00",
    mem = "64G"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"
    
# liftover EP benchmarking dataset from hg19 to hg38
rule liftover_crispr_dataset_Gasperini:
  input:
    results = "results/main_figure_1_and_2/{sample}/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz",
    enh_hg38 = "results/main_figure_1_and_2/{sample}/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    annot_hg38 = "results/genome_annotation_files/gencode.v26.annotation.gtf.gz"
  output: "results/main_figure_1_and_2/{sample}/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    time = "1:00:00",
    mem = "64G"
  script:
    "../scripts/create_encode_output/liftover_crispr_dataset.R"    
    
# filter EP benchmarking datasets for distance to TSS and minimum power
rule filter_crispr_dataset_Gasperini:
  input: "results/main_figure_1_and_2/{sample}/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
  output: temp("results/main_figure_1_and_2/{sample}/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{pwr}pwrAt{es}effect_{genome}.tsv.gz")
  params:
    tss_to_dist = config["benchmark_validation_datasets"]["encode_datasets"]["dist_to_TSS"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_encode_output/filter_crispr_dataset.R"
    
# create ensembl CRISPR dataset in both ENCODE format    
rule create_ensemble_encode_Gasperini:
  input:
    duplicate_pairs_analysis = "results/main_figure_1_and_2/{sample}/ENCODE/EPCrisprBenchmark/ENCODE_duplicate_pairs_analysis_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz"
  output: "results/main_figure_1_and_2/{sample}/ENCODE/ENCODE_duplicate_pairs_analysis_GRCh38.tsv.gz"
  params:
    effect_size = {"duplicate_pairs_analysis":"log2FC"},
    cell_types = {"duplicate_pairs_analysis": "K562"}
  conda: 
    "../envs/r_process_crispr_data.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/create_encode_output/create_ensemble_dataset.R"

# convert ensembl CRISPR dataset from ENCODE to EPBenchmarking format file  
rule create_ensemble_epbenchmarking_Gasperini:
  input: "results/main_figure_1_and_2/{sample}/ENCODE/ENCODE_duplicate_pairs_analysis_GRCh38.tsv.gz"
  output: "results/main_figure_1_and_2/{sample}/ENCODE/EPCrisprBenchmark/ENCODE_duplicate_pairs_analysis_GRCh38.tsv.gz"
  params:
    effect_size = "pctChange",
    min_pct_change = None
  conda: 
    "../envs/r_process_crispr_data.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/create_encode_output/create_ep_benchmarking_dataset.R"
    
    
    
    
    
    
    
    
    
    # Rule to get numbers for the paper
rule adding_design_file_information_Gasperini:
  input:
    gasperini = "results/main_figure_1_and_2/duplicate_pairs_analysis/ENCODE/EPCrisprBenchmark/ENCODE_duplicate_pairs_analysis_GRCh38.tsv.gz",
    guide_targets = "resources/main_figure_1_and_2/duplicate_pairs_analysis/guide_targets.tsv"
  output:
    results_with_design_file_features = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_design_file_features.tsv"
  log: "results/main_figure_1_and_2/logs/adding_design_file_information_Gasperini.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/main_figure_1_and_2/adding_design_file_information_Gasperini.R"

# Add more categories for understanding element overlap with different genomic features
rule adding_genomic_feature_overlaps_Gasperini:
  input:
    results_with_design_file_features = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_design_file_features.tsv",
    abc_canonical_tss = "results/genome_annotation_files/CollapsedGeneBounds.hg38.TSS500bp.bed"
  output:
    results_with_design_file_and_genomic_feature_overlaps = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_design_file_and_genomic_feature_overlaps.tsv"
  log: "results/main_figure_1_and_2/logs/adding_genomic_feature_overlaps_Gasperini.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/main_figure_1_and_2/adding_genomic_feature_overlaps_Gasperini.R"
    
# Creating categories to define the Random Set, Promoters, Valid Distal Element Gene pairs, etc.
rule adding_element_gene_pair_categories_Gasperini:
  input:
    results_with_design_file_and_genomic_feature_overlaps = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_design_file_and_genomic_feature_overlaps.tsv",
    guide_targets = "resources/main_figure_1_and_2/duplicate_pairs_analysis/guide_targets.tsv",
    create_ensemble_encode_input = "results/main_figure_1_and_2/duplicate_pairs_analysis/ENCODE/EPCrisprBenchmark/ENCODE_duplicate_pairs_analysis_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz",
    discovery_results = "results/main_figure_1_and_2/duplicate_pairs_analysis/differential_expression/results_run_discovery_analysis.rds"
  output:
    results_with_element_gene_pair_categories = "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv"
  log: "results/main_figure_1_and_2/logs/adding_element_gene_pair_categories_Gasperini.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/main_figure_1_and_2/adding_element_gene_pair_categories_Gasperini.R"

# Overall rule to run everything for the Gasperini dataset
rule run_gasperini_analysis:
  input:
    "results/main_figure_1_and_2/duplicate_pairs_analysis/results_with_element_gene_pair_categories.tsv"
  output:
    "results/main_figure_1_and_2/gasperini_analysis_complete.txt"
  shell:
    "touch {output}"
