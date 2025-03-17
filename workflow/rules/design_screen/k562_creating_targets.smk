# Snakemake workflow for ENCODE CRISPRi screen design (K562 creating targets)

# Hardcoded configuration parameters
TPM_DATA = "resources/design_screen/k562_creating_targets/tpm.csv"
GENOME_ANNOTATION = "results/design_screen/annotation_files/gencode.v29.annotation.gtf.gz"
DNASE_PEAKS = "resources/design_screen/k562_creating_targets/ENCFF185XRG_w_ENCFF325RTP_q30_sorted.txt"
TPM_THRESHOLD = 50
LOCUS_WIDTH = 2000000
DHS_THRESHOLD = 20
KB100_THRESHOLD = 10
Q75_THRESHOLD = 20
NUM_LOCI = 25
RANDOM_SEED = 69167

# Define the final outputs that should be created
# rule all:
#     input:
#         "results/design_screen/k562_creating_targets/final_screen/selected_loci.tsv",
#         "results/design_screen/k562_creating_targets/final_screen/candidate_peaks.bed",
#         "results/design_screen/k562_creating_targets/final_screen/loci_gene_info.tsv",
#         "results/design_screen/k562_creating_targets/qc/loci_distributions.pdf"

# # download gencode annotations
# rule download_gencode_annotations:
#   output:"results/design_screen/annotation_files/gencode.v29.annotation.gtf.gz"
#   params:
#     url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz"
#   conda: "../../envs/r_process_crispr_data.yml"
#   shell:
#     "wget -O {output} {params.url}"

# Process TPM data to identify genes above threshold
rule process_tpm_data:
    input:
        tpm_file = TPM_DATA
    output:
        tpm_filtered = "results/design_screen/k562_creating_targets/gene_selection/tpm_filtered.rds",
        tpm_plot = "results/design_screen/k562_creating_targets/gene_selection/tpm_distribution.pdf"
    params:
        tpm_threshold = TPM_THRESHOLD
    log:
        "results/design_screen/k562_creating_targets/logs/process_tpm_data.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/process_tpm_data.R"

# Process genome annotations
rule process_annotations:
    input:
        annotation_file = GENOME_ANNOTATION
    output:
        processed_annotation = "results/design_screen/k562_creating_targets/annotations/processed_annotation.rds"
    log:
        "results/design_screen/k562_creating_targets/logs/process_annotations.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/process_annotations.R"

# Process DNase-seq peaks
rule process_dnase_peaks:
    input:
        dnase_peaks = DNASE_PEAKS
    output:
        processed_peaks = "results/design_screen/k562_creating_targets/enhancers/processed_peaks.rds",
        peak_distribution = "results/design_screen/k562_creating_targets/enhancers/peak_distribution.pdf"
    log:
        "results/design_screen/k562_creating_targets/logs/process_dnase_peaks.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/process_dnase_peaks.R"

# Calculate locus statistics
rule calculate_locus_stats:
    input:
        tpm_filtered = "results/design_screen/k562_creating_targets/gene_selection/tpm_filtered.rds",
        processed_annotation = "results/design_screen/k562_creating_targets/annotations/processed_annotation.rds",
        processed_peaks = "results/design_screen/k562_creating_targets/enhancers/processed_peaks.rds"
    output:
        locus_stats = "results/design_screen/k562_creating_targets/locus_stats/locus_statistics.rds",
        locus_plots = "results/design_screen/k562_creating_targets/locus_stats/locus_statistics_plots.pdf"
    params:
        locus_width = LOCUS_WIDTH
    log:
        "results/design_screen/k562_creating_targets/logs/calculate_locus_stats.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/calculate_locus_stats.R"

# Filter loci based on multiple criteria
rule filter_loci:
    input:
        locus_stats = "results/design_screen/k562_creating_targets/locus_stats/locus_statistics.rds"
    output:
        filtered_loci = "results/design_screen/k562_creating_targets/locus_filtering/filtered_loci.rds",
        filtering_summary = "results/design_screen/k562_creating_targets/locus_filtering/filtering_summary.pdf"
    params:
        dhs_threshold = DHS_THRESHOLD,
        kb100_threshold = KB100_THRESHOLD,
        q75_threshold = Q75_THRESHOLD
    log:
        "results/design_screen/k562_creating_targets/logs/filter_loci.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/filter_loci.R"

# Sample loci for screen design
rule sample_loci:
    input:
        filtered_loci = "results/design_screen/k562_creating_targets/locus_filtering/filtered_loci.rds"
    output:
        sampled_loci = "results/design_screen/k562_creating_targets/loci_sampling/sampled_loci.rds",
        sampling_plots = "results/design_screen/k562_creating_targets/loci_sampling/sampling_plots.pdf"
    params:
        num_loci = NUM_LOCI,
        seed = RANDOM_SEED
    log:
        "results/design_screen/k562_creating_targets/logs/sample_loci.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/sample_loci.R"

# Quality control of selected loci
rule loci_qc:
    input:
        sampled_loci = "results/design_screen/k562_creating_targets/loci_sampling/sampled_loci.rds",
        locus_stats = "results/design_screen/k562_creating_targets/locus_stats/locus_statistics.rds"
    output:
        qc_report = "results/design_screen/k562_creating_targets/qc/loci_qc_report.html",
        qc_plots = "results/design_screen/k562_creating_targets/qc/loci_distributions.pdf"
    log:
        "results/design_screen/k562_creating_targets/logs/loci_qc.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/loci_qc.R"

# Generate final output files for screen design
rule generate_outputs:
    input:
        sampled_loci = "results/design_screen/k562_creating_targets/loci_sampling/sampled_loci.rds",
        processed_annotation = "results/design_screen/k562_creating_targets/annotations/processed_annotation.rds"
    output:
        selected_loci = "results/design_screen/k562_creating_targets/final_screen/selected_loci.tsv",
        candidate_peaks = "results/design_screen/k562_creating_targets/final_screen/candidate_peaks.bed",
        loci_gene_info = "results/design_screen/k562_creating_targets/final_screen/loci_gene_info.tsv"
    log:
        "results/design_screen/k562_creating_targets/logs/generate_outputs.log"
    conda:
        "../../envs/design_screen.yml"
    script:
        "../../scripts/design_screen/k562_creating_targets/generate_outputs.R"
