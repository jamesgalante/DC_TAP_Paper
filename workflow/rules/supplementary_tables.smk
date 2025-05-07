
# All rules here create a supplementary table or file in the paper

rule summary_of_element_gene_categories:
  input:
    wtc11_summary = "results/formatted_dc_tap_results/summary_WTC11.tsv",
    k562_summary = "results/formatted_dc_tap_results/summary_K562.tsv"
  output:
    summary_of_element_gene_categories_supplementary_table = "results/supplementary_tables/summary_of_element_gene_categories.tsv"
  log: "results/supplementary_tables/logs/summary_of_element_gene_categories.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "32G",
    time = "2:00:00"
  script:
    "../scripts/supplementary_tables/summary_of_element_gene_categories.R"
    
# Modify the guide targets files with the Guide Sequence without PAM or G
rule modify_guide_targets:
  input:
    guide_targets = expand("results/process_validation_datasets/{sample}/guide_targets.tsv", sample = ["K562_DC_TAP_Seq", "WTC11_DC_TAP_Seq"]),
    k562_design_file = "resources/process_validation_datasets/K562_DC_TAP_Seq/220308_K562_Random_Screen_Crop.design.txt",
    wtc11_design_file = "resources/process_validation_datasets/WTC11_DC_TAP_Seq/220308_WTC11_Random_Screen_sgOpti-DC.design.txt"
  output:
    k562_guide_targets_supp_table = "results/supplementary_tables/k562_guide_targets_supp_table.tsv",
    wtc11_guide_targets_supp_table = "results/supplementary_tables/wtc11_guide_targets_supp_table.tsv"
  log: "results/supplementary_tables/logs/modify_guide_targets.log"
  conda:
    "../envs/analyze_crispr_screen.yml"
  resources:
    mem = "8G",
    time = "2:00:00"
  script:
    "../scripts/supplementary_tables/modify_guide_targets.R"    
  
  
rule create_all_supplementary_tables:
  input:
    summary_of_element_gene_categories_supplementary_table = "results/supplementary_tables/summary_of_element_gene_categories.tsv",
    k562_guide_targets_supp_table = "results/supplementary_tables/k562_guide_targets_supp_table.tsv",
    wtc11_guide_targets_supp_table = "results/supplementary_tables/wtc11_guide_targets_supp_table.tsv"
  output:
    touch("results/supplementary_tables/create_all_supplementary_tables.done")
