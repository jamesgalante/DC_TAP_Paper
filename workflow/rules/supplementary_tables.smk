
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
    
    
rule create_all_supplementary_tables:
  input:
    summary_of_element_gene_categories_supplementary_table = "results/supplementary_tables/summary_of_element_gene_categories.tsv"
  output:
    touch("results/supplementary_tables/create_all_supplementary_tables.done")
