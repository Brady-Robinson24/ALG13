import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

print("Loading synaptic glycoproteins...")

# -----------------------------
# FILE PATHS
# -----------------------------

synaptic_file = r"C:\Users\brob6\ALG13 Bioinformatics\step3_synaptic_glycoproteins.csv"

genes_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_genes-rows.csv"
samples_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_samples-columns.csv"
matrix_file = r"C:\Users\brob6\ALG13 Bioinformatics\Cell type\human_MTG_2018-06-14_exon-matrix.csv"

output_file = r"C:\Users\brob6\ALG13 Bioinformatics\step6_gene_level_gaba_enrichment_filtered.csv"

plot_file = r"C:\Users\brob6\ALG13 Bioinformatics\step6_top_gaba_genes_plot.png"


# -----------------------------
# LOAD TARGET GENES
# -----------------------------

syn = pd.read_csv(synaptic_file)

target_genes = syn["gene_name"].str.upper().unique()

print("Target genes:", len(target_genes))


# -----------------------------
# LOAD ALLEN GENE METADATA
# -----------------------------

print("Loading Allen gene metadata...")

genes = pd.read_csv(genes_file)

genes["gene"] = genes["gene"].str.upper()

matching_rows = genes[genes["gene"].isin(target_genes)]

print("Matching genes in Allen dataset:", len(matching_rows))

gene_indices = matching_rows.index.tolist()
gene_symbols = matching_rows["gene"].values


# -----------------------------
# LOAD SAMPLE METADATA
# -----------------------------

print("Loading sample metadata...")

samples = pd.read_csv(samples_file)

cell_ids = samples["sample_name"]
cell_types = samples["class"]

gaba_mask = (cell_types == "GABAergic").values
glut_mask = (cell_types == "Glutamatergic").values

print("GABA cells:", gaba_mask.sum())
print("Glut cells:", glut_mask.sum())


# -----------------------------
# LOAD EXPRESSION MATRIX
# -----------------------------

print("Loading expression rows for matched genes...")

matrix = pd.read_csv(
    matrix_file,
    header=None,
    skiprows=lambda x: x not in gene_indices
)

# remove index column artifact
matrix = matrix.iloc[:, 1:]

matrix.columns = cell_ids
matrix.index = gene_symbols

print("Loaded matrix shape:", matrix.shape)


# -----------------------------
# CALCULATE GENE ENRICHMENT
# -----------------------------

results = []

for gene in matrix.index:

    gene_expr = matrix.loc[gene]

    gaba_mean = gene_expr[gaba_mask].mean()
    glut_mean = gene_expr[glut_mask].mean()

    log2_fc = np.log2((gaba_mean + 1e-6) / (glut_mean + 1e-6))

    results.append({
        "gene": gene,
        "mean_GABA": gaba_mean,
        "mean_Glut": glut_mean,
        "log2FC_GABA_vs_Glut": log2_fc
    })

results_df = pd.DataFrame(results)


# -----------------------------
# FILTER LOW EXPRESSION GENES
# -----------------------------

print("Applying expression filter...")

results_df["mean_expression"] = (
    results_df["mean_GABA"] + results_df["mean_Glut"]
) / 2

filtered = results_df[results_df["mean_expression"] > 1]

print("Genes remaining after filter:", len(filtered))


# -----------------------------
# SORT BY ENRICHMENT
# -----------------------------

filtered = filtered.sort_values(
    "log2FC_GABA_vs_Glut",
    ascending=False
)


print("\nTop enriched genes:")
print(filtered.head(15))


# -----------------------------
# SAVE RESULTS
# -----------------------------

filtered.to_csv(output_file, index=False)

print("\nSaved filtered results to:")
print(output_file)


# -----------------------------
# PLOT TOP GENES
# -----------------------------

top = filtered.head(15)

plt.figure(figsize=(8,6))

sns.barplot(
    data=top,
    y="gene",
    x="log2FC_GABA_vs_Glut",
    palette="viridis"
)

plt.title("Top GABA-Enriched Synaptic Glycoproteins")
plt.xlabel("log2(GABA / Glut expression)")
plt.ylabel("Gene")

plt.tight_layout()

plt.savefig(plot_file, dpi=300)

plt.show()

print("\nSaved plot to:")
print(plot_file)