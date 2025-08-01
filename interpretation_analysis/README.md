This directory contains code, intermediate files, and results used for the interpretation of DC-TAP-seq results and comparison to previous screens.


- `dc_tap_analysis.functions.R`: function definitions
- `dc_tap_analysis.run.R`: code to produce results for Fig. 5, 6, S4a
- `categorized_pairs`: used for Fig. 5 and 6; produced using results from ENCODE-rE2G, chrom-annotate, and DC-TAP-seq analysis pipelines
  - `dc_tap_gasperini.random_distal_element_pairs.categorized.tsv.gz`: DC-TAP random screen DE-G pairs and Gasperini DE-G pairs with chromatin category, promoter class, direct effect rate, and power annotations
  - `other_crispr_screens.categorized.tsv.gz`: same as above but for CRISPR perturbation data from 8 other studies
  - `all_gw_pairs.categorized.summary.tsv`
  - `thresholds.tsv`: quantiles for various chromatin features for distal elements in 5 cell types used for categorization
 


