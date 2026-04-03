import pandas as pd

# ----------------------------------------
# FILE PATHS
# ----------------------------------------

brain_file = r"C:\Users\brob6\ALG13 Bioinformatics\step1_brain_expressed_glycoproteins.csv"
glyco_file = r"C:\Users\brob6\ALG13 Bioinformatics\N-glycosylation_full_hPROTEIN_list.csv"
uniprot_file = r"C:\Users\brob6\ALG13 Bioinformatics\uniprot_human_annotations.tsv"

# ----------------------------------------
# LOAD DATA
# ----------------------------------------

brain = pd.read_csv(brain_file)
glyco = pd.read_csv(glyco_file)
uniprot = pd.read_csv(uniprot_file, sep="\t")

# ----------------------------------------
# MERGE GLYCO + BRAIN
# ----------------------------------------

merged = glyco.merge(brain, left_on="gene_name", right_on="Gene name")

print("Brain glycoproteins:", merged["gene_name"].nunique())

# ----------------------------------------
# MERGE WITH UNIPROT
# ----------------------------------------

merged = merged.merge(uniprot, left_on="gene_name", right_on="Gene Names", how="left")

# ----------------------------------------
# SECRETORY PATHWAY FILTER
# ----------------------------------------

secretory = merged[
    merged["Signal peptide"].notna() |
    merged["Transmembrane"].notna() |
    merged["Subcellular location [CC]"].str.contains("Secreted|Cell membrane|Extracellular", na=False)
]

print("Secretory-pathway proteins:", secretory["gene_name"].nunique())

# ----------------------------------------
# SAVE
# ----------------------------------------

output = r"C:\Users\brob6\ALG13 Bioinformatics\step2_secretory_glycoproteins.csv"
secretory.to_csv(output, index=False)

print("Saved to:", output)