import pandas as pd

glyco_file = r"C:\Users\brob6\ALG13 Bioinformatics\N-glycosylation_full_hPROTEIN_list.csv"
hpa_file = r"C:\Users\brob6\ALG13 Bioinformatics\rna_brain_region_hpa.tsv"

print("Loading glycoprotein dataset...")
glyco = pd.read_csv(glyco_file)

print("Columns in glycoprotein dataset:")
print(glyco.columns)

print("Loading HPA brain dataset...")
hpa = pd.read_csv(hpa_file, sep="\t")

glyco_genes = glyco["gene_name"].dropna().unique()
print("Total glycoproteins:", len(glyco_genes))

# --------------------------------------------------
# MATCH GLYCOPROTEINS TO BRAIN EXPRESSION DATA
# --------------------------------------------------

brain_glyco = hpa[hpa["Gene name"].isin(glyco_genes)].copy()

print("Matching glycoproteins found in HPA dataset:", brain_glyco["Gene name"].nunique())

# --------------------------------------------------
# COMPUTE MEAN EXPRESSION ACROSS BRAIN REGIONS
# --------------------------------------------------

brain_expr = (
    brain_glyco.groupby("Gene name")["nTPM"]
    .mean()
    .reset_index()
)

# --------------------------------------------------
# APPLY EXPRESSION THRESHOLD
# --------------------------------------------------

threshold = 1

brain_filtered = brain_expr[
    brain_expr["nTPM"] >= threshold
]

print("Brain-expressed glycoproteins (nTPM ≥ 1):", len(brain_filtered))

# --------------------------------------------------
# SAVE OUTPUT
# --------------------------------------------------

output_file = r"C:\Users\brob6\ALG13 Bioinformatics\step1_brain_expressed_glycoproteins.csv"

brain_filtered.to_csv(output_file, index=False)

print("Saved filtered dataset to:")
print(output_file)