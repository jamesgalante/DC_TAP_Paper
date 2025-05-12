## Create input files for ENCODE and distal regulation CRISPR benchmarking pipeline

# First call create_encode_dataset for the datasets in hg19 then liftover the datasets to hg38
ruleorder: liftover_crispr_dataset > create_encode_dataset

# function to get samples that require liftover from hg19 to GRCh38
def liftover_samples(config):
  genome_builds = config["benchmark_validation_datasets"]["encode_datasets"]["genome_build"].items()
  liftover_samples = list(dict(filter(lambda x: x[1] == "hg19", genome_builds)).keys())
  return(liftover_samples)

# download UCSC hg19 to hg38 liftover chain file
rule download_chain_file:
  output: "results/genome_annotation_files/hg19ToHg38.over.chain.gz"
  params:
    url = config["benchmark_validation_datasets"]["download_urls"]["liftover_chain"]
  conda: "../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"

# compile output files in ENCODE format
rule create_encode_dataset:
  input:
    results = "results/process_validation_datasets/{sample}/power_analysis/output_{sd}gStd_{method}_{strategy}.tsv",
    annot = "results/genome_annotation_files/gencode.v32lift37.annotation.gtf.gz",
    guide_targets = "results/process_validation_datasets/{sample}/guide_targets.tsv"
  output: "results/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
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
    
## Liftover CRISPRi datasets -----------------------------------------------------------------------

# lift enhancer coordinates from hg19 to hg38 using UCSC's liftOver software    
rule liftover_enhancers:
  input:
    results = "results/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz", 
    chain = "results/genome_annotation_files/hg19ToHg38.over.chain.gz"
  output:
    hg19 = "results/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg19.bed",
    hg38 = "results/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    unlifted = "results/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_unlifted.bed"
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
rule liftover_crispr_dataset:
  input:
    results = "results/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_hg19.tsv.gz",
    enh_hg38 = "results/create_encode_output/ENCODE/liftover/{sample}_{sd}gStd_{method}_{strategy}/enh_hg38.bed",
    annot_hg38 = "results/genome_annotation_files/gencode.v26.annotation.gtf.gz"
  output: "results/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_GRCh38.tsv.gz"
  wildcard_constraints:
    sample = "|".join(liftover_samples(config))
  conda: "../envs/r_process_crispr_data.yml"
  resources:
    time = "1:00:00",
    mem = "64G"
  script:
    "../scripts/create_encode_output/liftover_crispr_dataset.R"    

## Create EPBenchmarking CRISPR data files ---------------------------------------------------------

# filter EP benchmarking datasets for distance to TSS and minimum power
rule filter_crispr_dataset:
  input: "results/create_encode_output/ENCODE/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{genome}.tsv.gz"
  output: temp("results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_{sample}_{sd}gStd_{method}_{strategy}_{pwr}pwrAt{es}effect_{genome}.tsv.gz")
  params:
    tss_to_dist = config["benchmark_validation_datasets"]["encode_datasets"]["dist_to_TSS"]
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_encode_output/filter_crispr_dataset.R"

## Create ensemble dataset -------------------------------------------------------------------------

# create ensembl CRISPR dataset in both ENCODE format    
rule create_ensemble_encode:
  input:
    K562_DC_TAP_Seq = "results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_K562_DC_TAP_Seq_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz",
    WTC11_DC_TAP_Seq = "results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_WTC11_DC_TAP_Seq_0.13gStd_Sceptre_perCRE_0.8pwrAt15effect_GRCh38.tsv.gz",
  output: "results/create_encode_output/ENCODE/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  params:
    effect_size = {"K562_DC_TAP_Seq":"log2FC", "WTC11_DC_TAP_Seq":"log2FC"},
    cell_types = {"K562_DC_TAP_Seq": "K562", "WTC11_DC_TAP_Seq": "WTC11"}
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_encode_output/create_ensemble_dataset.R"
    
# convert ensembl CRISPR dataset from ENCODE to EPBenchmarking format file  
rule create_ensemble_epbenchmarking:
  input: "results/create_encode_output/ENCODE/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  output: "results/create_encode_output/ENCODE/EPCrisprBenchmark/ENCODE_Combined_Validation_Datasets_GRCh38.tsv.gz"
  params:
    effect_size = "pctChange",
    min_pct_change = None
  conda: "../envs/r_process_crispr_data.yml"
  script:
    "../scripts/create_encode_output/create_ep_benchmarking_dataset.R"
