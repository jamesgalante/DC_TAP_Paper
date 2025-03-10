## Create input files for ENCODE and distal regulation CRISPR benchmarking pipeline

ruleorder: liftover_crispr_dataset > create_encode_dataset

# function to get samples that require liftover from hg19 to GRCh38
def liftover_samples(config):
  genome_builds = config["benchmark_validation_datasets"]["create_encode_output"]["encode_datasets"]["genome_build"].items()
  liftover_samples = list(dict(filter(lambda x: x[1] == "hg19", genome_builds)).keys())
  return(liftover_samples)

# download UCSC hg19 to hg38 liftover chain file
rule download_chain_file:
  output: "resources/benchmark_validation_datasets/create_encode_output/hg19ToHg38.over.chain.gz"
  params:
    url = config["benchmark_validation_datasets"]["create_encode_output"]["download_urls"]["liftover_chain"]
  conda: "../../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# download gencode annotations
rule download_gencode_annotations:
  output: "resources/benchmark_validation_datasets/create_encode_output/{annot}.annotation.gtf.gz"
  params:
    url = lambda wildcards: config["benchmark_validation_datasets"]["create_encode_output"]["download_urls"][wildcards.annot]
  conda: "../../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# Need to add a rule here to copy the guide targets over from the previous outputs into the benchmarking resources directory

# compile output files in ENCODE format
rule create_encode_dataset:
  input:
    results = "results/process_validation_datasets/{sample}/power_analysis/output_{sd}gStd_{method}_{strategy}.tsv",
    annot = lambda wildcards: config["benchmark_validation_datasets"]["create_encode_output"]["encode_datasets"]["annot"][wildcards.sample],
    guide_targets = "results/process_validation_datasets/{sample}/guide_targets.tsv"
  output: "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
  params:
    ignore_txs = lambda wildcards: config["benchmark_validation_datasets"]["create_encode_output"]["encode_datasets"]["ignore_transcripts"][wildcards.sample],
    tss_min_dist = config["benchmark_validation_datasets"]["create_encode_output"]["encode_datasets"]["dist_to_TSS"][0],
    gene_ids = lambda wildcards: config["benchmark_validation_datasets"]["create_encode_output"]["metadata"][wildcards.sample]["gene_ids"],
    tss_ctrl_tag = lambda wildcards: config["benchmark_validation_datasets"]["create_encode_output"]["metadata"][wildcards.sample]["tss_ctrl_tag"],
    padj_threshold = config["process_validation_datasets"]["differential_expression"]["padj_threshold"],
    reference = lambda wildcards: config["benchmark_validation_datasets"]["create_encode_output"]["metadata"][wildcards.sample]["reference"]
  conda: "../../envs/r_process_crispr_data.yml"
  resources:
    time = "2:00:00",
    mem = "32G"
  script:
    "../../scripts/benchmark_validation_datasets/encode_datasets/create_encode_dataset.R"
    
## Liftover CRISPRi datasets -----------------------------------------------------------------------

# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_enhancers:
  input:
    results = "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz", 
    chain = "resources/benchmark_validation_datasets/create_encode_output/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/benchmark_validation_datasets/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg19.bed",
    hg38 = "results/benchmark_validation_datasets/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    unlifted = "results/benchmark_validation_datasets/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_unlifted.bed"
  conda: "../../envs/r_process_crispr_data.yml"
  resources:
    time = "1:00:00",
    mem = "64G"
  shell:
    "zcat {input.results} | "
    """awk 'BEGIN {{OFS = "\\t"}} (NR>1) {{print $1, $2, $3, $1":"$2"-"$3, 0, "."}}' | """
    "sort | uniq > {output.hg19} ; "
    "liftOver {output.hg19} {input.chain} {output.hg38} {output.unlifted}"
    
# liftover EP benchmarking dataset from hg19 to hg38
rule liftover_crispr_dataset:
  input:
    results = "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz",
    enh_hg38 = "results/benchmark_validation_datasets/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    annot_hg38 = "resources/benchmark_validation_datasets/create_encode_output/gencode.v32.annotation.gtf.gz"
  output: "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_GRCh38.tsv.gz"
  wildcard_constraints:
    sample = "|".join(liftover_samples(config))
  conda: "../../envs/r_process_crispr_data.yml"
  resources:
    time = "1:00:00",
    mem = "64G"
  script:
    "../../scripts/benchmark_validation_datasets/encode_datasets/liftover_crispr_dataset.R"    

## Create EPBenchmarking CRISPR data files ---------------------------------------------------------

# filter EP benchmarking datasets for distance to TSS and minimum power
rule filter_crispr_dataset:
  input: "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
  output: temp("results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{pwr}pwrAt{es}effect_{genome}.tsv.gz")
  params:
    tss_to_dist = config["benchmark_validation_datasets"]["create_encode_output"]["encode_datasets"]["dist_to_TSS"],
    remove_filtered_pairs = False
  conda: "../../envs/r_process_crispr_data.yml"
  script:
    "../../scripts/benchmark_validation_datasets/encode_datasets/filter_crispr_dataset.R"

# convert ENCODE format files to EPBenchmarking format files
rule create_ep_benchmarking_dataset:
  input: "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_{sd}gStd_Sceptre_perCRE_{pwr}pwrAt{es}effect_{genome}.tsv.gz"
  output: "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/EPCrisprBenchmark_{sample}_{sd}gStd_{pwr}pwrAt{es}effect_{genome}.tsv.gz"
  params:
    effect_size = "log2FC", # Sceptre specific
    min_pct_change = None,
    cell_type = lambda wildcards: config["benchmark_validation_datasets"]["create_encode_output"]["metadata"][wildcards.sample]["cell_type"]
  conda: "../../envs/r_process_crispr_data.yml"
  script:
    "../../scripts/benchmark_validation_datasets/encode_datasets/create_ep_benchmarking_dataset.R"

## Create ensemble dataset -------------------------------------------------------------------------

# create ensembl CRISPR dataset in both ENCODE format    
rule create_ensemble_encode:
  input:
    K562_DC_TAP_Seq = "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_K562_DC_TAP_Seq_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz",
    WTC11_DC_TAP_Seq = "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_WTC11_DC_TAP_Seq_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz",
  output: "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  params:
    effect_size = {"K562_DC_TAP_Seq":"log2FC", "WTC11_DC_TAP_Seq":"log2FC"},
    cell_types = {"K562_DC_TAP_Seq": "K562", "WTC11_DC_TAP_Seq": "WTC11"}
  conda: "../../envs/r_process_crispr_data.yml"
  script:
    "../../scripts/benchmark_validation_datasets/encode_datasets/create_ensemble_dataset.R"
    
# convert ensembl CRISPR dataset from ENCODE to EPBenchmarking format file  
rule create_ensemble_epbenchmarking:
  input: "results/benchmark_validation_datasets/create_encode_output/ENCODE/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  output: "results/benchmark_validation_datasets/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  params:
    effect_size = "pctChange",
    min_pct_change = None
  conda: "../../envs/r_process_crispr_data.yml"
  script:
    "../../scripts/benchmark_validation_datasets/encode_datasets/create_ep_benchmarking_dataset.R"   
