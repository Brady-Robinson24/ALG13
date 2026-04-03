import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# --------------------------------------------------
# FILES
# --------------------------------------------------

BRAIN_FILE = "transcript_rna_brain.tsv"
GENE_MAP_FILE = "hpa_gene_info.tsv"

# --------------------------------------------------
# GENES OF INTEREST
# --------------------------------------------------

genes = [
    "GABRA1","GABRA2","GABRA3","GABRA5",
    "GABRB1","GABRB2","GABRB3","GABRG2","GABRD",
    "SLC6A1","GABBR1","GABBR2","SLC32A1",
    "NLGN2","GPHN",
    "GRIA1","GRIA2","GRIA3","GRIA4",
    "GRIN1","GRIN2A","GRIN2B",
    "GRIK2","GRIK5",
    "SLC1A2","SLC1A3",
    "NLGN1","NLGN3",
    "NRXN1","NRXN2","NRXN3",
    "CNTNAP2","CNTN2","CNTN4",
    "L1CAM","NCAM1",
    "LRRTM2","LRRTM3",
    "ADAM22","ADAM23",
    "LGI1",
    "SCN1A","SCN2A","SCN8A",
    "CACNA1A","CACNA2D1","CACNA2D2",
    "KCNQ2","KCNQ3",
    "HCN1","HCN2",
    "NPTX1","NPTX2",
    "RELN",
    "EFNB2","NTN1",
    "ALG13"
]

# --------------------------------------------------
# REGIONS
# --------------------------------------------------
# "Epilepsy-relevant" anatomical group
epi_regions = [
    "hippocampus",
    "amygdala",
    "cerebral cortex",
    "thalamus",
    "basal ganglia"
]

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

print("Loading transcript-level brain dataset...")
brain_df = pd.read_csv(BRAIN_FILE, sep="\t")

print("Loading gene mapping table...")
gene_map = pd.read_csv(GENE_MAP_FILE, sep="\t")

# HPA mapping file columns from your printout:
# "Ensembl" and "Gene"
gene_map = gene_map[["Ensembl", "Gene"]].drop_duplicates()
gene_map.columns = ["ensgid", "Gene"]

# --------------------------------------------------
# MERGE TRANSCRIPT DATA WITH GENE SYMBOLS
# --------------------------------------------------

print("Merging transcript data with gene symbols...")
brain_df = brain_df.merge(gene_map, on="ensgid", how="left")

# Keep only genes of interest
brain_df = brain_df[brain_df["Gene"].isin(genes)].copy()

# --------------------------------------------------
# FIND TPM COLUMNS
# --------------------------------------------------

tpm_cols = [c for c in brain_df.columns if c.startswith("TPM.")]
if len(tpm_cols) == 0:
    raise ValueError("No TPM columns found. Check the file structure.")

# --------------------------------------------------
# COLLAPSE TRANSCRIPTS -> GENES
# --------------------------------------------------
# IMPORTANT:
# For transcript-level TPM to gene-level expression, summing transcript TPMs
# within a gene is a reasonable approximation and better than taking the mean.

print("Collapsing transcript isoforms to gene-level TPM by summation...")
gene_expr = brain_df.groupby("Gene")[tpm_cols].sum()

# --------------------------------------------------
# SPLIT COLUMNS INTO EPILEPSY-RELEVANT VS OTHER
# --------------------------------------------------

epi_cols = []
other_cols = []

for c in tpm_cols:
    cl = c.lower()
    if any(region in cl for region in epi_regions):
        epi_cols.append(c)
    else:
        other_cols.append(c)

print(f"Epilepsy-region samples: {len(epi_cols)}")
print(f"Other brain samples: {len(other_cols)}")

if len(epi_cols) == 0:
    raise ValueError("No epilepsy-region TPM columns found. Check region strings.")
if len(other_cols) == 0:
    raise ValueError("No non-epilepsy TPM columns found.")

# --------------------------------------------------
# RUN TWO-TAILED WELCH T-TEST FOR EACH GENE
# --------------------------------------------------

results = []

for gene in genes:
    if gene not in gene_expr.index:
        results.append({
            "Gene": gene,
            "mean_epi": np.nan,
            "mean_other": np.nan,
            "log2FC_epi_vs_other": np.nan,
            "p_value": np.nan,
            "n_epi": len(epi_cols),
            "n_other": len(other_cols)
        })
        continue

    x = gene_expr.loc[gene, epi_cols].astype(float).values
    y = gene_expr.loc[gene, other_cols].astype(float).values

    mean_epi = np.mean(x)
    mean_other = np.mean(y)

    # small pseudocount to avoid log2(0)
    log2fc = np.log2((mean_epi + 1e-6) / (mean_other + 1e-6))

    # two-tailed Welch t-test
    stat, p = ttest_ind(x, y, equal_var=False, nan_policy="omit")

    results.append({
        "Gene": gene,
        "mean_epi": mean_epi,
        "mean_other": mean_other,
        "log2FC_epi_vs_other": log2fc,
        "p_value": p,
        "n_epi": len(x),
        "n_other": len(y)
    })

results_df = pd.DataFrame(results)

# --------------------------------------------------
# MULTIPLE TESTING CORRECTION
# --------------------------------------------------

valid_mask = results_df["p_value"].notna()
adj = multipletests(results_df.loc[valid_mask, "p_value"], method="fdr_bh")

results_df.loc[valid_mask, "p_adj_fdr"] = adj[1]
results_df.loc[valid_mask, "significant_nominal_p05"] = results_df.loc[valid_mask, "p_value"] < 0.05
results_df.loc[valid_mask, "significant_fdr_p05"] = results_df.loc[valid_mask, "p_adj_fdr"] < 0.05

# Top contenders:
# must be significantly different after FDR and higher in epilepsy-relevant regions
results_df["top_contender"] = (
    (results_df["log2FC_epi_vs_other"] > 0) &
    (results_df["p_adj_fdr"] < 0.05)
)

# Sort for convenience
results_df = results_df.sort_values(
    by=["top_contender", "log2FC_epi_vs_other", "p_adj_fdr"],
    ascending=[False, False, True]
)

# --------------------------------------------------
# SAVE OUTPUT
# --------------------------------------------------

results_df.to_csv("alg13_epilepsy_region_vs_other_brain_ttest_results.csv", index=False)

# Also save only the top contenders
top_df = results_df[results_df["top_contender"]].copy()
top_df.to_csv("alg13_top_contenders_fdr05.csv", index=False)

print("\nAnalysis complete.\n")
print(results_df[[
    "Gene", "mean_epi", "mean_other", "log2FC_epi_vs_other",
    "p_value", "p_adj_fdr", "top_contender"
]].head(20))
