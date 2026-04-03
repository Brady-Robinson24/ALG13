import pandas as pd
import numpy as np

# --------------------------------------------------
# FILE PATHS
# --------------------------------------------------

synaptic_file = r"C:\Users\brob6\ALG13 Bioinformatics\step3_synaptic_glycoproteins.csv"

matrix_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_exon-matrix.csv"
genes_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_genes-rows.csv"
samples_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_samples-columns.csv"


# --------------------------------------------------
# LOAD SYNAPTIC GLYCOPROTEINS
# --------------------------------------------------

print("Loading synaptic glycoproteins...")

synaptic = pd.read_csv(synaptic_file)

target_genes = (
    synaptic["gene_name"]
    .dropna()
    .astype(str)
    .str.upper()
    .str.strip()
    .unique()
)

target_genes = set(target_genes)

print("Target genes:", len(target_genes))


# --------------------------------------------------
# LOAD ALLEN GENE METADATA
# --------------------------------------------------

print("Loading Allen gene metadata...")

genes = pd.read_csv(genes_file)

genes["gene"] = (
    genes["gene"]
    .astype(str)
    .str.upper()
    .str.strip()
)

gene_rows = genes[genes["gene"].isin(target_genes)]

print("Matching genes in Allen dataset:", len(gene_rows))

if len(gene_rows) == 0:
    raise ValueError("No gene matches found between synaptic list and Allen dataset.")


# --------------------------------------------------
# LOAD SAMPLE METADATA
# --------------------------------------------------

print("Loading sample metadata...")

samples = pd.read_csv(samples_file)

cell_ids = samples["sample_name"]
cell_types = samples["class"]

print("Unique cell types:")
print(cell_types.unique())


# --------------------------------------------------
# LOAD EXPRESSION MATRIX ROWS
# --------------------------------------------------

print("Loading expression rows for matched genes...")

matrix = pd.read_csv(
    matrix_file,
    header=None,
    index_col=0,
    skiprows=lambda x: x not in gene_rows.index
)

matrix.index = gene_rows["gene"].values
matrix.columns = cell_ids

print("Loaded matrix shape:", matrix.shape)


# --------------------------------------------------
# CALCULATE CELL-TYPE EXPRESSION
# --------------------------------------------------

print("Calculating cell-type enrichment...")

results = {}

for ct in cell_types.unique():

    cells = samples[samples["class"] == ct]["sample_name"]

    subset = matrix[cells]

    results[ct] = subset.mean().mean()


results_df = pd.DataFrame.from_dict(
    results,
    orient="index",
    columns=["mean_expression"]
)

results_df = results_df.sort_values("mean_expression", ascending=False)


# --------------------------------------------------
# SAVE RESULTS
# --------------------------------------------------

output_file = r"C:\Users\brob6\ALG13 Bioinformatics\step5_celltype_enrichment.csv"

results_df.to_csv(output_file)

print("\nCell-type enrichment:")
print(results_df)

print("\nSaved results to:")
print(output_file)
