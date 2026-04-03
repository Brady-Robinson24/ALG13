import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

# --------------------------------------------------
# FILE PATHS
# --------------------------------------------------

synaptic_file = r"C:\Users\brob6\ALG13 Bioinformatics\step3_synaptic_glycoproteins.csv"
brain_file = r"C:\Users\brob6\ALG13 Bioinformatics\rna_brain_region_hpa.tsv"

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

print("Loading synaptic glycoproteins...")
synaptic = pd.read_csv(synaptic_file)

print("Loading HPA brain expression dataset...")
brain = pd.read_csv(brain_file, sep="\t")

print("\nColumns in brain dataset:")
print(list(brain.columns))

# --------------------------------------------------
# GET GENE LIST
# --------------------------------------------------

genes = synaptic["gene_name"].dropna().unique()

print("\nSynaptic candidate genes:", len(genes))

# --------------------------------------------------
# FILTER BRAIN DATA FOR THESE GENES
# --------------------------------------------------

brain = brain[brain["Gene name"].isin(genes)].copy()

# --------------------------------------------------
# CREATE LOWERCASE REGION COLUMN
# --------------------------------------------------

brain["region_lower"] = brain["Brain region"].str.lower()

# --------------------------------------------------
# DEFINE EPILEPSY CIRCUIT REGIONS
# --------------------------------------------------

epi_regions = [
    "hippocampus",
    "cerebral cortex",
    "amygdala",
    "thalamus",
    "basal ganglia"
]

# --------------------------------------------------
# ANALYSIS
# --------------------------------------------------

results = []

for gene in genes:

    g = brain[brain["Gene name"] == gene]

    epi = g[g["region_lower"].isin(epi_regions)]["nTPM"].values
    other = g[~g["region_lower"].isin(epi_regions)]["nTPM"].values

    if len(epi) == 0 or len(other) == 0:
        continue

    mean_epi = np.mean(epi)
    mean_other = np.mean(other)

    log2fc = np.log2((mean_epi + 1e-6) / (mean_other + 1e-6))

    stat, p = ttest_ind(epi, other, equal_var=False)

    results.append({
        "gene": gene,
        "mean_epi": mean_epi,
        "mean_other": mean_other,
        "log2FC_epi_vs_other": log2fc,
        "p_value": p
    })

# --------------------------------------------------
# CREATE RESULTS TABLE
# --------------------------------------------------

results_df = pd.DataFrame(results)

# --------------------------------------------------
# SORT BY ENRICHMENT
# --------------------------------------------------

results_df = results_df.sort_values(
    by="log2FC_epi_vs_other",
    ascending=False
)

# --------------------------------------------------
# SAVE RESULTS
# --------------------------------------------------

output_file = r"C:\Users\brob6\ALG13 Bioinformatics\step4_epilepsy_region_enrichment.csv"

results_df.to_csv(output_file, index=False)

print("\nTop candidates enriched in epilepsy regions:")
print(results_df.head(10))

print("\nSaved results to:")
print(output_file)