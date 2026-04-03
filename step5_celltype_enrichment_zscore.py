import pandas as pd
import numpy as np
from scipy.stats import zscore

print("Loading synaptic glycoproteins...")

# -----------------------------------
# FILE PATHS
# -----------------------------------

synaptic_file = r"C:\Users\brob6\ALG13 Bioinformatics\step3_synaptic_glycoproteins.csv"

genes_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_genes-rows.csv"
samples_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_samples-columns.csv"
matrix_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_exon-matrix.csv"

output_file = r"C:\Users\brob6\ALG13 Bioinformatics\step5_celltype_enrichment_zscore.csv"

# -----------------------------------
# LOAD TARGET GENE LIST
# -----------------------------------

syn = pd.read_csv(synaptic_file)

target_genes = syn["gene_name"].str.upper().unique()

print("Target genes:", len(target_genes))

# -----------------------------------
# LOAD ALLEN GENE METADATA
# -----------------------------------

print("Loading Allen gene metadata...")

genes = pd.read_csv(genes_file)

genes["gene"] = genes["gene"].str.upper()

matching_rows = genes[genes["gene"].isin(target_genes)]

print("Matching genes in Allen dataset:", len(matching_rows))

gene_indices = matching_rows.index.tolist()

# -----------------------------------
# LOAD SAMPLE METADATA
# -----------------------------------

print("Loading sample metadata...")

samples = pd.read_csv(samples_file)

cell_ids = samples["sample_name"]
cell_types = samples["class"]

print("Unique cell types:")
print(cell_types.unique())

# -----------------------------------
# LOAD EXPRESSION MATRIX
# -----------------------------------

print("Loading expression rows for matched genes...")

matrix = pd.read_csv(
    matrix_file,
    header=None,
    skiprows=lambda x: x not in gene_indices,
)

# first column is index artifact
matrix = matrix.iloc[:,1:]

matrix.columns = cell_ids

print("Loaded matrix shape:", matrix.shape)

# -----------------------------------
# Z-SCORE NORMALIZATION PER GENE
# -----------------------------------

print("Performing gene-wise Z-score normalization...")

matrix_z = matrix.sub(matrix.mean(axis=1), axis=0)
matrix_z = matrix_z.div(matrix.std(axis=1).replace(0, np.nan), axis=0)

# -----------------------------------
# CELL-TYPE ENRICHMENT
# -----------------------------------

print("Calculating normalized cell-type enrichment...")

results = {}

for cell_type in cell_types.unique():

    mask = (cell_types == cell_type).to_numpy()

    subset = matrix_z.loc[:, mask]

    results[cell_type] = subset.mean().mean()

results_df = pd.DataFrame.from_dict(
    results,
    orient="index",
    columns=["zscore_expression"]
)

results_df = results_df.sort_values(
    "zscore_expression",
    ascending=False
)

print("\nZ-score normalized enrichment:")
print(results_df)

# -----------------------------------
# SAVE RESULTS
# -----------------------------------

results_df.to_csv(output_file)

print("\nSaved results to:")
print(output_file)