# Script: downsampling_reads_for_UMIs.py

### SETUP =====================================================================

# Saving state for debugging
import sys
import pickle

def save_state(filename):
	with open(f"RDA_objects/downsampling_reads_for_UMIs.pkl", "wb") as f:
		pickle.dump(globals(), f, protocol=pickle.HIGHEST_PROTOCOL)
	print("Saved state")

save_state("downsampling_reads_for_UMIs")
# sys.exit()  # Uncomment to stop manually after saving state

# Open log file for messages, warnings, and errors
log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log


### LOADING FILES =============================================================

print("Loading packages")
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
import scanpy as sc
import os

print("Loading input files")
# Get dataset paths from Snakemake input object
datasets = [
	{
		"name": "DC-TAP-seq (MOI6)",
		"mol_info": snakemake.input.moi6_mol_info,
		"filtered_h5": snakemake.input.moi6_filtered_h5
	},
	{
		"name": "CRISPRi Direct Capture",
		"mol_info": snakemake.input.crisprdi_mol_info,
		"filtered_h5": snakemake.input.crisprdi_filtered_h5
	},
	{
		"name": "DC-TAP-seq (MOI3)",
		"mol_info": snakemake.input.moi3_mol_info,
		"filtered_h5": snakemake.input.moi3_filtered_h5
	}
]

print(f"Will process {len(datasets)} datasets:")
for dataset in datasets:
	print(f"  - {dataset['name']}")
	

### ANALYSIS =========================================================

def get_complexity_df(read_df, dataset_name, fracs=None):
	"""
	Process reads data to perform downsampling analysis.
	This follows Maddie's approach but without the TPM related parts.
	"""
	# Get the raw reads
	read_df['reads'] = read_df['count'].astype(int)
	read_df.rename(columns={'barcode_idx': 'cell'}, inplace=True)
	
	# Create a dataframe with each read represented by a row
	print(f"Expanding reads for {dataset_name}...")
	df_repeat = read_df.reindex(read_df[['umi', 'reads']].index.repeat(read_df.reads))
	
	# Define downsampling fractions (0% to 100%)
	if fracs is None:
		fracs = np.arange(0.1, 1.1, 0.1)
	
	# Create empty results dataframe
	results = pd.DataFrame(index=range(len(fracs)),
						  columns=['proportion', 'mean_reads_per_cell', 'mean_umis_per_cell_per_gene', 'dataset'])
	
	# Downsample at each fraction
	for i, frac in enumerate(fracs):
		print(f"Processing fraction {frac} for {dataset_name}")
		
		# Take a sampling of our reads at that fraction
		sample = df_repeat.sample(frac=frac, random_state=42)
		
		# Total reads leftover after sampling
		total_reads = sample.shape[0]
		
		# Drop the duplicate reads now to count unique molecular identifiers per cell per gene
		umis_cell_gene = sample.drop_duplicates().groupby(['cell', 'gene_name']).count()['umi'].reset_index()
		
		# Remake a cell counts matrix (cell x gene with umi counts per gene in that cell)
		cell_counts_matrix = umis_cell_gene.pivot(index='cell', columns='gene_name').fillna(0)
		
		# Get metrics from that matrix
		mean_umis_per_cell_per_gene = cell_counts_matrix.mean().mean()
		
		# Total reads per cell considers all in the sampling
		mean_reads_per_cell = total_reads / len(sample.cell.unique())
		
		# Store results
		results.loc[i, 'proportion'] = frac
		results.loc[i, 'mean_reads_per_cell'] = mean_reads_per_cell
		results.loc[i, 'mean_umis_per_cell_per_gene'] = mean_umis_per_cell_per_gene
		results.loc[i, 'dataset'] = dataset_name
	
	return results

def process_dataset(dataset):
	"""Process a single dataset using Maddie's approach"""
	name = dataset["name"]
	mol_info_path = dataset["mol_info"]
	filtered_h5_path = dataset["filtered_h5"]
	
	print(f"Processing {name}...")
	
	# Check if files exist
	if not os.path.exists(mol_info_path) or not os.path.exists(filtered_h5_path):
		print(f"Files not found for dataset: {name}")
		return None
	
	try:
		# Read in final output of counts
		counts = sc.read_10x_h5(filtered_h5_path)
		ensg_id_map = counts.var.reset_index()
		
		# Read in molecule_info.h5
		mol_dfs = {}
		names = 'barcode_idx umi count feature_idx'.split()
		for name_field in names:
			with h5py.File(mol_info_path, 'r') as f:
				mol_dfs[name_field] = pd.DataFrame(np.array(f[name_field]))
		
		molecule_info = pd.concat(mol_dfs, axis=1)
		molecule_info.columns = molecule_info.columns.droplevel(1)
		
		# Read features and barcodes
		with h5py.File(mol_info_path, "r") as f:
			ds_obj = f['features']
			gene_id_df = pd.DataFrame(ds_obj['id'][()].astype(str))
			
		with h5py.File(mol_info_path, "r") as f:
			ds_obj = f['barcodes']
			barcode_id_df = pd.DataFrame(ds_obj[()].astype(str))
		
		# Get the actual cell names and gene names
		molecule_info['barcode_idx'] = molecule_info['barcode_idx'].map(lambda x: barcode_id_df.iloc[x, 0] if x < len(barcode_id_df) else None)
		molecule_info['feature_idx'] = molecule_info['feature_idx'].map(lambda x: gene_id_df.iloc[x, 0] if x < len(gene_id_df) else None)
		
		# Map to gene names from the counts matrix
		if isinstance(counts.var.index, pd.MultiIndex):
			gene_map = dict(zip(counts.var.index.get_level_values(0), counts.var.index.get_level_values(1)))
		else:
			gene_map = dict(zip(ensg_id_map['gene_ids'], ensg_id_map['index'] if 'index' in ensg_id_map.columns else ensg_id_map.index))
			
		molecule_info['gene_name'] = molecule_info['feature_idx'].map(gene_map)
		molecule_info.dropna(subset=['gene_name'], inplace=True)
		
		# Filter to cells in the counts matrix
		cell_barcodes = counts.obs.index.str.replace('-1$', '', regex=True)
		reads_df = molecule_info[molecule_info['barcode_idx'].isin(cell_barcodes)]
		
		# Run the downsampling analysis
		if not reads_df.empty:
			return get_complexity_df(reads_df, name)
		else:
			print(f"No matching cells found for {name}")
			return None
			
	except Exception as e:
		print(f"Error processing {name}: {str(e)}")
		return None

# Process all datasets
all_results = []

for dataset in datasets:
	results = process_dataset(dataset)
	if results is not None:
		all_results.append(results)
		

### SAVE OUTPUT ===============================================================

print("Saving output files")

# Combine all results
if all_results:
	combined_results = pd.concat(all_results)
	
	# Save the raw results
	combined_results.to_csv(snakemake.output.combined_results, index=False)
	
	# Create plot
	plt.figure(figsize=(10, 6))
	sns.lineplot(data=combined_results, x='mean_reads_per_cell', y='mean_umis_per_cell_per_gene', 
				 hue='dataset', marker='o', markersize=8)
	
	plt.title('DC-TAP-seq Effect Size Analysis', fontsize=16, fontweight='bold')
	plt.xlabel('Mean Total Reads per Cell', fontsize=14)
	plt.ylabel('Mean UMIs per Cell per Gene', fontsize=14)
	plt.grid(alpha=0.3)
	plt.tight_layout()
	
	# Set colors similar to the R plot
	colors = ["#4477AA", "#EE6677", "#228833"]
	for i, line in enumerate(plt.gca().get_lines()):
		if i < len(colors):
			line.set_color(colors[i])
	
	plt.savefig(snakemake.output.plot_pdf)
	plt.savefig(snakemake.output.plot_png)
	print("Analysis complete. Plots saved")
else:
	print("No data to plot. Check file paths and processing.")
	

### CLEAN UP ==================================================================

print("Closing log file")
log.close()
