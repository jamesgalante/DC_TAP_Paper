# This snakemake file is for running sceptre on shraivogel to compare to DC TAP Seq

# combine guide targets file from both experiments
rule combine_tapseq_guide_targets:
  input:
    chr8 = "resources/TAPseqChr8/guide_targets_hg38.tsv",
    chr11 = "resources/TAPseqChr11/guide_targets_hg38.tsv"
  output:
    guide_targets = "resources/TAPseq/guide_targets.tsv",
    guides_bed = "resources/TAPseq/guide_positions.bed",
    targets_bed = "resources/TAPseq/guide_targets.bed"
  conda: "../envs/r_create_tapseq_input.yml"
  script:
    "../scripts/tapseq_dataset/combine_tapseq_guide_targets.R"
    
# add 10x lanes to SCE colData and filter out cells from bad lanes
rule filter_tapseq_10x_lanes:
  input: "resources/{sample}/perturb_sce_unfilt.rds"
  output: "resources/{sample}/perturb_sce.rds"
  wildcard_constraints:
    sample = "TAPseqChr.+"
  params:
    remove_lanes = lambda wildcards: config["process_tapseq"]["remove_lanes"][wildcards.sample]
  log: "resources/{sample}/logs/filter_tapseq_10x_lanes.log"
  conda: "../envs/r_create_tapseq_input.yml"
  script:
    "../scripts/tapseq_dataset/filter_tapseq_10x_lanes.R"
