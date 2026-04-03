import pandas as pd

# --------------------------------------------------
# FILE PATHS
# --------------------------------------------------

secretory_file = r"C:\Users\brob6\ALG13 Bioinformatics\step2_secretory_glycoproteins.csv"
syngo_file = r"C:\Users\brob6\ALG13 Bioinformatics\Syngo\annotations.xlsx"

# --------------------------------------------------
# LOAD DATA
# --------------------------------------------------

print("Loading secretory glycoproteins...")
secretory = pd.read_csv(secretory_file)

print("Loading SynGO annotations...")
syngo = pd.read_excel(syngo_file)

# --------------------------------------------------
# SHOW COLUMNS (useful for debugging)
# --------------------------------------------------

print("\nColumns in SynGO file:")
print(list(syngo.columns))

# --------------------------------------------------
# EXTRACT SYNAPTIC GENE SYMBOLS
# --------------------------------------------------

if "hgnc_symbol" not in syngo.columns:
    raise ValueError("Expected column 'hgnc_symbol' not found in SynGO dataset.")

syngo_genes = syngo["hgnc_symbol"].dropna().unique()

print("\nTotal SynGO genes:", len(syngo_genes))

# --------------------------------------------------
# INTERSECT WITH SECRETORY GLYCOPROTEINS
# --------------------------------------------------

print("Secretory glycoproteins:", secretory["gene_name"].nunique())

synaptic = secretory[secretory["gene_name"].isin(syngo_genes)].copy()

print("Synaptic glycoproteins:", synaptic["gene_name"].nunique())

# --------------------------------------------------
# SAVE OUTPUT
# --------------------------------------------------

output_file = r"C:\Users\brob6\ALG13 Bioinformatics\step3_synaptic_glycoproteins.csv"

synaptic.to_csv(output_file, index=False)

print("\nSaved synaptic glycoproteins to:")
print(output_file)