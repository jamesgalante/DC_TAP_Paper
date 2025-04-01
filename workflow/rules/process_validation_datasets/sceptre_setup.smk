
# download gencode annotations
rule download_gencode_v32lift37_K562_DC_TAP:
  output: "resources/process_validation_datasets/sceptre_setup/K562_DC_TAP_Seq/gencode.v32lift37.annotation.gtf.gz"
  params:
    url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/gencode.v32lift37.annotation.gtf.gz"
  conda: "../../envs/r_process_crispr_data.yml"
  shell:
    "wget -O {output} {params.url}"
    
    
# The process_cell_ranger_outputs scripts are essentially the same, but for qc'ing the immediate results, I feel I can separate these steps for now and easily combine them later
rule process_cell_ranger_outputs_WTC11_DC_TAP_Seq:
  input:
    cell_ranger_directory = "resources/process_validation_datasets/sceptre_setup/WTC11_DC_TAP_Seq/cell_ranger_output/"
  output:
    guide_matrix = "results/process_validation_datasets/WTC11_DC_TAP_Seq/raw_counts/perturb_status.rds",
    gene_matrix = "results/process_validation_datasets/WTC11_DC_TAP_Seq/raw_counts/dge.rds"
  log: "results/process_validation_datasets/WTC11_DC_TAP_Seq/logs/process_cell_ranger_outputs_WTC11_DC_TAP_Seq.log"
  conda:
    "../../envs/seurat_for_cell_ranger_outputs.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_setup/process_cell_ranger_outputs_WTC11_DC_TAP_Seq.R"
  
rule process_cell_ranger_outputs_K562_DC_TAP_Seq:
  input:
    cell_ranger_directory = "resources/process_validation_datasets/sceptre_setup/K562_DC_TAP_Seq/cell_ranger_output/"
  output:
    guide_matrix = "results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/perturb_status.rds",
    gene_matrix = "results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/dge.rds" 
  log: "results/process_validation_datasets/K562_DC_TAP_Seq/logs/process_cell_ranger_outputs_K562_DC_TAP_Seq.log"
  conda:
    "../../envs/seurat_for_cell_ranger_outputs.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_setup/process_cell_ranger_outputs_K562_DC_TAP_Seq.R"

# The guide design files were copied from [,] into the resources directory 
  # /oak/stanford/groups/engreitz/Users/dulguun/CRISPRDesign/220308_Random_Screen_FinalDesignFiles/220308_K562_Random_Screen_Crop.design.txt
  # /oak/stanford/groups/engreitz/Users/dulguun/CRISPRDesign/220308_Random_Screen_FinalDesignFiles/220308_WTC11_Random_Screen_Crop.design.txt
  
  
# I think i can drastically shorten this process with this file: "https://github.com/EngreitzLab/DC_TAP_paper/blob/main/inputs_for_Table_SX_RandomScreenDesign_guides_3_4/new_merged_k562_grna_groups_table_fixed.txt"
# Everything except the merging
rule create_K562_DC_TAP_Seq_guide_targets:
  input:
    guide_design_file = "resources/process_validation_datasets/sceptre_setup/K562_DC_TAP_Seq/220308_K562_Random_Screen_Crop.design.txt",
    
    ### TEMP UNTIL SCREEN DESIGN IN PIPELINE
    positive_controls = "resources/process_validation_datasets/sceptre_setup/K562_DC_TAP_Seq/ChosenGenes.AllRegions.bed",
    old_guide_targets = "resources/process_validation_datasets/sceptre_setup/K562_DC_TAP_Seq/old_guide_targets.tsv",
    old_pre_merged_guide_targets = "resources/process_validation_datasets/sceptre_setup/K562_DC_TAP_Seq/even_older_guide_targets.tsv"
    ### TEMP UNTIL SCREEN DESIGN IN PIPELINE
    
  output:
    guide_targets_file = "results/process_validation_datasets/K562_DC_TAP_Seq/guide_targets.tsv"
  log: "results/process_validation_datasets/K562_DC_TAP_Seq/logs/create_K562_DC_TAP_Seq_guide_targets.log"
  conda:
    "../../envs/r_process_crispr_data.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_setup/create_K562_DC_TAP_Seq_guide_targets.R"

rule create_WTC11_DC_TAP_Seq_guide_targets:
  input:
    guide_design_file = "resources/process_validation_datasets/sceptre_setup/WTC11_DC_TAP_Seq/220308_WTC11_Random_Screen_sgOpti-DC.design.txt",
    
    ### TEMP UNTIL SCREEN DESIGN IN PIPELINE
    old_pre_merged_guide_targets = "resources/process_validation_datasets/sceptre_setup/WTC11_DC_TAP_Seq/even_older_guide_targets.tsv",
    chosen_genes_all_regions = "resources/process_validation_datasets/sceptre_setup/WTC11_DC_TAP_Seq/ChosenGenes.AllRegions.bed",
    chosen_genes_distal_elements = "resources/process_validation_datasets/sceptre_setup/WTC11_DC_TAP_Seq/ChosenGenes.DistalElements.bed"
    ### TEMP UNTIL SCREEN DESIGN IN PIPELINE
    
  output:
    guide_targets_file = "results/process_validation_datasets/WTC11_DC_TAP_Seq/guide_targets.tsv"
  log: "results/process_validation_datasets/WTC11_DC_TAP_Seq/logs/create_WTC11_DC_TAP_Seq_guide_targets.log"
  conda:
    "../../envs/r_process_crispr_data.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_setup/create_WTC11_DC_TAP_Seq_guide_targets.R"



rule create_sceptre_diffex_input_K562_DC_TAP_Seq:
  input:
    dge = "results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/dge.rds",
    perturb_status = "results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/perturb_status.rds",
    guide_targets = "results/process_validation_datasets/K562_DC_TAP_Seq/guide_targets.tsv",
    annot = "resources/process_validation_datasets/sceptre_setup/genome_annotation_files/gencode.v32lift37.annotation.gtf.gz"
  output:
    gene_gRNA_group_pairs = "results/process_validation_datasets/K562_DC_TAP_Seq/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/K562_DC_TAP_Seq/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/K562_DC_TAP_Seq/metadata.rds",
    sceptre_diffex_input = "results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/sceptre_diffex_input.rds",
    distances = "results/process_validation_datasets/K562_DC_TAP_Seq/distances.tsv"
  params:
    negative_control_genes = config["process_validation_datasets"]["sceptre_setup"]["negative_controls"]["K562_DC_TAP_Seq"]
  log: "results/process_validation_datasets/K562_DC_TAP_Seq/logs/create_sceptre_diffex_input_K562_DC_TAP_Seq.log"
  conda:
    "../../envs/sceptre_power_simulations.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_setup/create_sceptre_diffex_input_K562_DC_TAP_Seq.R"
    
### SOMEHOW THERE ARE DUPLICATES IN DISTANCES - FIGURE OUT WHY AND FIX AND FIX IN OTHER SCRIPTS TOO

rule create_sceptre_diffex_input_WTC11_DC_TAP_Seq:
  input:
    dge = "results/process_validation_datasets/WTC11_DC_TAP_Seq/raw_counts/dge.rds",
    perturb_status = "results/process_validation_datasets/WTC11_DC_TAP_Seq/raw_counts/perturb_status.rds",
    guide_targets = "results/process_validation_datasets/WTC11_DC_TAP_Seq/guide_targets.tsv",
    annot = "resources/process_validation_datasets/sceptre_setup/genome_annotation_files/gencode.v32lift37.annotation.gtf.gz"
  output:
    gene_gRNA_group_pairs = "results/process_validation_datasets/WTC11_DC_TAP_Seq/gene_gRNA_group_pairs.rds",
    gRNA_groups_table = "results/process_validation_datasets/WTC11_DC_TAP_Seq/gRNA_groups_table.rds",
    metadata = "results/process_validation_datasets/WTC11_DC_TAP_Seq/metadata.rds",
    sceptre_diffex_input = "results/process_validation_datasets/WTC11_DC_TAP_Seq/differential_expression/sceptre_diffex_input.rds",
    distances = "results/process_validation_datasets/WTC11_DC_TAP_Seq/distances.tsv"
  params:
    negative_control_genes = config["process_validation_datasets"]["sceptre_setup"]["negative_controls"]["WTC11_DC_TAP_Seq"]
  log: "results/process_validation_datasets/WTC11_DC_TAP_Seq/logs/create_sceptre_diffex_input_WTC11_DC_TAP_Seq.log"
  conda:
    "../../envs/sceptre_power_simulations.yml"
  resources:
    mem = "64G",
    time = "2:00:00"
  script:
    "../../scripts/process_validation_datasets/sceptre_setup/create_sceptre_diffex_input_WTC11_DC_TAP_Seq.R"
    
    



