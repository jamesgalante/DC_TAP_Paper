# Script: crispr_benchmarking.smk
    
# merge predictions with experimental data
rule mergePredictionsWithExperiment_noGeneRemoval:
  input:
    predictions = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["pred"].values(),
    experiment = "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz",
    tss_universe = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["tss_universe"],
    gene_universe = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["gene_universe"],
    pred_config = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["pred_config"],
    cell_type_mapping = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["cell_type_mapping"].values(),
    expressed_genes = lambda wildcards: config["benchmark_validation_datasets"]["crispr_benchmarking"]["comparisons"][wildcards.dataset]["expressed_genes"]
  output:
    merged = "results/benchmark_validation_datasets/crispr_benchmarking/expt_pred_merged_annot/{dataset}_expt_pred_merged_annot.txt.gz"
  params:
    pos_col = "Regulated",
    include_col = "include",
    filter_include_col = False
  log: 
    "results/benchmark_validation_datasets/crispr_benchmarking/logs/mergePredictionsWithExperiment_{dataset}.log"
  conda: 
    "../../envs/r_crispr_comparison.yml"
  resources:
    mem_mb = 72000
  script:
    "../../scripts/benchmark_validation_datasets/crispr_benchmarking/mergePredictionsWithExperiment.R"
