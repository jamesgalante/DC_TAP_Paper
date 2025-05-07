# DC-TAP Workflow: Analyzing Unbiased TAP-Seq Screens

This repository contains a robust and reproducible workflow for analyzing datasets generated from an unbiased DC TAP-Seq screen on K562 and WTC11. The final output file can be found in `results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_and_merged_elements.tsv`. An explanation of each column in that table can be found in `resources/formatting_dc_tap_results/columns_of_Final_DC_TAP_Seq_Results_table.txt`

## Analysis Summary

Setup the CellRanger outputs for SCEPTRE differential expression analysis
- These steps are summarized in `sceptre_setup.smk` and involve processing the CellRanger outputs, modifying the guide_targets.tsv files and setting up SCEPTRE objects for differential expression

Run differential expression analysis and power simulations using SCEPTRE
- These steps are detailed in `sceptre_power_analysis.smk` where most rules are trying to efficiently run the power simulations based on SCEPTRE's package. The results of the power simulations and differential expression analysis are summarized in `rule format_sceptre_output` and then passed for further formatting

Overlap the tested elements with genomic features and liftOver element coordinates
- In `create_encode_output.smk`, the elements are overlapped with genomic features (promoter, gene body) in order to interpret the differential expression results. The results of K562 and WTC11 are combined at the end of this ruleset.

Format all previous steps, and add different flags based on overlap with genomic features and chromatin features
- In `formatting_dc_tap_results.smk`, the results of the previous analyses are summarized, confidence intervals are added based on a dev branch of SCEPTRE, notes about the original design are added, elements are categorized as "distal" or "promoter" and further categorized based on relationship with the tested gene. Thus a category for each element-gene pair is defined. Finally, chromatin categories are assigned to each element based on Maya Sheth's pipeline.

## Important Files

Screen Result Files (Created in `sceptre_setup.smk`)
  - Raw Gene Counts Matrix: `results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/dge.rds`
  - Raw Guide Counts Matrix: `results/process_validation_datasets/K562_DC_TAP_Seq/raw_counts/perturb_status.rds`
  - Metadata File: `results/process_validation_datasets/K562_DC_TAP_Seq/metadata.rds`
  - Guide Design File: `results/process_validation_datasets/K562_DC_TAP_Seq/guide_targets.tsv`
  - All Pairs Tested for Differential Expression: `results/process_validation_datasets/K562_DC_TAP_Seq/gene_gRNA_group_pairs.rds`
  - Input for SCEPTRE Differential Expression: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/sceptre_diffex_input.rds`

Differential Expression Output Files (Created in `sceptre_power_analysis.smk`)
  - SCEPTRE Differential Expression Results: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/results_run_discovery_analysis.rds`
  - SCEPTRE Calibration Check Results (Negative Control Guide Testing): `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/results_run_calibration_check.rds`
  - SCEPTRE Object post-differential expression analysis: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/final_sceptre_object.rds`
  - SCEPTRE Differential Expression Summary Statistics: `results/process_validation_datasets/K562_DC_TAP_Seq/differential_expression/analysis_summary.txt`

Singleton Differential Expression Output Files (Created in `sceptre_power_analysis.smk`)
  - All Results - See (2) for descriptions of each file: `results/process_validation_datasets/K562_DC_TAP_Seq/singleton_differential_expression`

Power Analysis Output Files (created in `sceptre_power_analysis.smk`)
  - Batched Outputs (intermediate files): `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/effect_size_*/`
  - Combined Batched Outputs (intermediate files): `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/combined_power_analysis_results_es_*.tsv`
  - Analyzed Power Sim Outputs (intermediate files): `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/power_analysis_results_es_*.tsv`
  - Final Power Sim Output: `results/process_validation_datasets/K562_DC_TAP_Seq/power_analysis/output_0.13gStd_Sceptre_perCRE.tsv`

Intermediate Labelling of Elements (Created in `create_encode_output.smk`):
  - All Results: `results/create_encode_output/ENCODE/*`
  - The order of execution for this ruleset is `create_encode_dataset` (hg19), `liftover_enhancers` & `liftover_crispr_dataset`, `filter_crispr_dataset` (hg38), `create_ensemble_encode`, `create_ensemble_epbenchmarking`
  - Essentially everything done in this ruleset is represented in the final file

Final Screen Output Files (Created in `formatting_dc_tap_results.smk`):
  - SCEPTRE Differential Expression Results with Confidence Intervals (See point 2 for file descriptions): `results/formatted_dc_tap_results/differential_expression_w_confidence_intervals_K562_DC_TAP_Seq/*`
  - Final Output w/ EG categories and specific pairs modified: `results/formatted_dc_tap_results/results_with_element_gene_pair_categories_modified.tsv`
  - Summary Statistics of EG categories: `results/formatted_dc_tap_results/summary_K562.tsv`
  - Final Output w/ Added columns for 500bp extension and merging: `results/formatted_dc_tap_results/resized_and_merged_input_for_chromatin_categorization_pipeline.tsv`
  - Final Output w/ Chromatin Categories calculated for resized/merged regions: `results/formatted_dc_tap_results/Final_DC_TAP_Seq_Results_w_Chromatin_Categories_on_resized_elements.tsv`

Supplementary Tables:
  - Summary of element-gene categories: `results/supplementary_tables/summary_of_element_gene_categories.tsv`