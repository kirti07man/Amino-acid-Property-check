# Amino-acid-Property-check
Python tool for classifying IDR amino-acid segments by physicochemical propensity (hydrophobic, polar, charged, aromatic) from Excel sequence data.

# ============================================
# 1) Upload your Excel file
# ============================================
from google.colab import files
import pandas as pd

uploaded = files.upload()  # Upload your .xlsx file

excel_filename = list(uploaded.keys())[0]
print("Loaded file:", excel_filename)

df = pd.read_excel(excel_filename)
print("Columns in your file:", df.columns.tolist())
display(df.head())

# ============================================
# 2) SET THESE COLUMN NAMES CORRECTLY
# ============================================
# Change these to EXACTLY match your column names in Excel
AA_COL = 'Res_name'          # e.g. 'AA', 'AminoAcid', 'Residue', etc.
REGION_COL = 'Region'     # e.g. 'IDR', 'Region', 'IDR_region', etc.

# ============================================
# 3) Mapping: single-letter AA -> category
# ============================================

aa_category = {
    # Non-polar (hydrophobic)
    'A': "Non-polar (hydrophobic)",
    'V': "Non-polar (hydrophobic)",
    'L': "Non-polar (hydrophobic)",
    'I': "Non-polar (hydrophobic)",
    'M': "Non-polar (hydrophobic)",
    'G': "Non-polar (hydrophobic)",
    'P': "Non-polar (hydrophobic)",

    # Aromatic
    'F': "Aromatic",
    'W': "Aromatic",
    'Y': "Aromatic",

    # Polar, hydrophilic (uncharged)
    'S': "Polar, hydrophilic (uncharged)",
    'T': "Polar, hydrophilic (uncharged)",
    'N': "Polar, hydrophilic (uncharged)",
    'Q': "Polar, hydrophilic (uncharged)",
    'C': "Polar, hydrophilic (uncharged)",

    # Polar, hydrophilic (acidic)
    'D': "Polar, hydrophilic (acidic)",
    'E': "Polar, hydrophilic (acidic)",

    # Polar, hydrophilic (basic)
    'K': "Polar, hydrophilic (basic)",
    'R': "Polar, hydrophilic (basic)",
    'H': "Polar, hydrophilic (basic)",
}

# Make sure amino acid codes are uppercase single letters
df['AA_clean'] = df[AA_COL].astype(str).str.strip().str.upper()

# Map to category
df['AA_Category'] = df['AA_clean'].map(aa_category).fillna("Unknown")

print("Annotated residue-level table:")
display(df.head(20))

# ============================================
# 4) Region-wise category fractions + dominant type
# ============================================

# Drop rows with missing region
df_valid = df.dropna(subset=[REGION_COL])

# Count residues per (Region, Category)
region_cat_counts = (
    df_valid
    .groupby([REGION_COL, 'AA_Category'])
    .size()
    .reset_index(name='count')
)

print("Counts per region-category:")
display(region_cat_counts)

# Convert to fractions
region_cat_fraction = (
    region_cat_counts
    .groupby(REGION_COL)
    .apply(lambda x: x.assign(fraction=x['count'] / x['count'].sum()))
    .reset_index(drop=True)
)

print("Fractions per region-category:")
display(region_cat_fraction)

# Pivot to wide table: each IDR = one row
region_cat_table = region_cat_fraction.pivot(
    index=REGION_COL,
    columns='AA_Category',
    values='fraction'
).fillna(0.0)

# Dominant category per IDR region
region_cat_table['Dominant_Category'] = region_cat_table.idxmax(axis=1)

print("Region-level summary (IDR1–IDR5):")
display(region_cat_table)

# ============================================
# 5) Save outputs and download
# ============================================
residue_outfile = "amino_acid_with_category.xlsx"
region_outfile = "idr_region_category_summary.xlsx"

df.to_excel(residue_outfile, index=False)
region_cat_table.to_excel(region_outfile)

print("Saved files:", residue_outfile, "and", region_outfile)

from google.colab import files
files.download(residue_outfile)
files.download(region_outfile)
