from pathlib import Path
import os
import pandas as pd
import numpy as np
import glob

# -------------------------------------------------------------------------
# METADATA / FEATURE-TABLE CONSISTENCY
#
# Sample identifiers and BioProject assignments were extensively curated and
# validated upstream when constructing Supplementary_table1.xlsx and the final
# feature tables. The code therefore assumes:
#   1. Sample ID is unique in SRA_info_index;
#   2. each Sample ID listed in SRA_info_index is present in the merged feature
#      table;
#   3. BioProject Number mappings are internally consistent.
#
# These are treated as properties of the curated input dataset rather than
# re-validated at this preprocessing step.

file_path = 'Data//Additional file 1.xlsx'
df = pd.read_excel(file_path, sheet_name='Table S1', skiprows=1)
df = df[df["Set"] == "Training"] # Unique studies = 27
df["Citation"] = df["Author"].astype(str) + " et al., " + df["Publication Year"].astype(str)

sra_data = pd.read_excel(file_path, sheet_name='Table S2', skiprows=1)
sra_data = sra_data[sra_data["Set"]=="Training"] # samples 1654
sra_data = sra_data[['Sample ID', 'BioProject Number', 'Category', 'Disease']]
# -------------------------------------------------------------------------


# Filtering cut-offs
rich_cutoff     = 5 # richness cutoff per sample
pf              = 0.05 # prevalence filtering

def read_genus_table(path):
    with open(path, 'r', encoding='utf-8') as f: 
        lines = f.readlines()

    header_idx = next((i for i,l in enumerate(lines) if l.lstrip().startswith('#OTU ID')), None)
    header = lines[header_idx].lstrip('#').strip().split('\t')
    table_data = pd.read_csv(path, sep='\t', header=None, names=header, skiprows=header_idx+1, dtype=str)

    for c in table_data.columns:
        if c != 'OTU ID': table_data[c] = pd.to_numeric(table_data[c], errors='coerce')
    return table_data

paths = sorted(glob.glob(os.path.join('Feature_tables', '**', '*_final_table_depth.tsv'), recursive=True))

allowed_projects = set(df[df["Set"] == "Training"]["BioProject Number"].astype(str))

# Keep only the files in paths whose BioProject ID is in the Excel rows labeled Training.
filtered_files = [f for f in paths if os.path.basename(f).split("_")[2] in allowed_projects]
print(f"Found {len(filtered_files)} genus tables in index")

# CREATE ALL GENUS TABLES MERGED
merged = None
for p in filtered_files:
    table_data = read_genus_table(p)
    sample_cols = [c for c in table_data.columns if c != 'OTU ID']
    if merged is not None:
        dupes = set(sample_cols) & set([c for c in merged.columns if c != 'OTU ID'])
        if dupes: raise ValueError(f"Duplicate sample IDs found: {dupes} in {p}")
        merged = merged.merge(table_data, on='OTU ID', how='outer')
    else: merged = table_data.copy()

for c in merged.columns:
    if c != 'OTU ID': merged[c] = pd.to_numeric(merged[c], errors='coerce').fillna(0.0)
sample_cols = [c for c in merged.columns if c != 'OTU ID']
merged = merged[['OTU ID'] + sorted(sample_cols)]

# CLEAN TAXONOMY 
BAD_LABELS = {"__", "uncultured", "Unknown_Family"}

def has_bad_token(x):
    x = x.strip()
    return (
        x == ""
        or x == "__"
        or any(b.lower() in x.lower() for b in BAD_LABELS if b != "__")
    )

def simplify_taxonomy(tax):
    parts = [p.strip() for p in str(tax).split(";")]
    rank_vals = {}
    for p in parts:
        if "__" in p:
            rank, val = p.split("__", 1)
            rank_vals[rank] = val.strip()

    # analyze genus first. There is some mislabeling. e.g some taxa at the genus level have 'family' taxonomy. 
    g_val = rank_vals.get("g", "")
    genus_bad = (has_bad_token(g_val) or g_val.endswith("ceae") or g_val.endswith("ales"))
    genus_known = not genus_bad

    # pass 2: clean
    cleaned = []
    for p in parts:
        if "__" in p:
            rank, val = p.split("__", 1)
            v = val.strip()

            is_bad = (has_bad_token(v) 
                or (rank == "g" and (v.endswith("ceae") or v.endswith("ales")))
                or (rank == "f" and v.endswith("ales") and not genus_known)
            )
            cleaned.append("__" if is_bad else p)
        else:
            cleaned.append("__" if has_bad_token(p) else p)

    return ";".join(cleaned)

merged["OTU ID"] = merged["OTU ID"].apply(simplify_taxonomy)
merged = merged.groupby("OTU ID", as_index=False).sum(numeric_only=True)

# Manual filter of contaminants; identified in environmental control samples across several studies.
OTU_IDS_TO_EXCLUDE = {
    "d__Bacteria;p__Acidobacteriota;c__Acidobacteriae;o__Subgroup_2;f__Subgroup_2;g__Subgroup_2",
    "d__Bacteria;p__Firmicutes;c__Bacilli;o__Mycoplasmatales;f__Mycoplasmataceae;g__Ureaplasma",
    "d__Bacteria;p__Fusobacteriota;c__Fusobacteriia;o__Fusobacteriales;f__Fusobacteriaceae;g__Cetobacterium",
    "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Weeksellaceae;g__Cloacibacterium"
}
if OTU_IDS_TO_EXCLUDE:
    n_before = merged.shape[0]
    merged = merged[~merged["OTU ID"].isin(OTU_IDS_TO_EXCLUDE)].copy()
    

# Samples in SRA / metadata

merged = merged[["OTU ID", *sra_data["Sample ID"]]].copy()

ra_sub = merged.set_index('OTU ID').copy() 
ra_sub_initial = ra_sub.copy() # raw abundance # cound data

# -------------------------------------------------------------------------
# INITIAL SPARSITY FILTERING
#
# Prevalence filtering (5%) is intentionally performed before final cohort/study
# eligibility filtering. At this stage the objective is to remove extremely
# sparse taxa from the harmonized genus table. The final analytical cohort
# is defined subsequently using sample- and study-level inclusion criteria.
# -------------------------------------------------------------------------

prevalence = (ra_sub > 0).sum(axis=1) / ra_sub.shape[1]
keep_taxa = prevalence[prevalence >= pf].index
ra_sub = ra_sub.loc[keep_taxa]

# Richness filter 5 taxa
richness = (ra_sub > 0).sum(axis=0)
keep_samples = richness[richness >= rich_cutoff].index
ra_sub = ra_sub[keep_samples]
ra_raw = ra_sub.copy()
ra_raw.to_csv('Data/Genus_level_raw_count_featureTable.tsv', sep = '\t')

# Convert to relative abundance using the 184 taxa
col_sums = ra_sub.sum(axis=0)
ra_sub = ra_sub.div(col_sums, axis=1).fillna(0.0) # 184 taxa and 1654 samples

# Keep only metadata rows for samples present in ra_sub
meta = sra_data[sra_data['Sample ID'].isin(keep_samples)].copy()

# -------------------------------------------------------------------------
# FINAL COHORT DEFINITION
#
# After sparsity- and sample-level filtering, define the analytical cohort by
# retaining BioProject x Category groups containing >=10 samples.
# This filtering is intentionally performed after the initial prevalence and
# richness filters because cohort eligibility is based on the samples remaining
# after microbiome quality/sparsity filtering.
# -------------------------------------------------------------------------

meta['proj_cat_count'] = meta.groupby(['BioProject Number', 'Category'])['Sample ID'].transform('count')
meta_filtered = meta[meta['proj_cat_count'] >= 10].copy()


# SANITY CHECK; ENFORCE min 10 samples per BioProject overall ----------
proj_counts = meta_filtered['BioProject Number'].value_counts()
keep_bioprojects = proj_counts[proj_counts >= 10].index
meta_filtered = meta_filtered[meta_filtered['BioProject Number'].isin(keep_bioprojects)]

# Samples to keep after both filters
keep_samples = meta_filtered['Sample ID'].unique() #1654 samples; above code was used just as a safety check

# Filter ra_sub to the final set of samples
ra_sub = ra_sub[keep_samples].copy()
filtered_t = ra_sub.T.reset_index()
filtered_t.rename(columns={'index': 'Sample ID'}, inplace=True)

# === Merge metadata ===
merged_meta = filtered_t.merge(sra_data, on='Sample ID', how='left')
merged_meta = merged_meta.merge(df, on='BioProject Number', how='left')


meta_cols = [c for c in merged_meta.columns if not c.startswith("d__")]
taxa_cols = [c for c in merged_meta.columns if c.startswith("d__")]

n_t, n_s = merged_meta.shape
print(f"Number of samples: {n_t}, and number of taxa: {len(taxa_cols)}")
print(f"Number of unique studies: {merged_meta['BioProject Number'].nunique()}")

training_set = merged_meta.set_index(['BioProject Number', 'Sample ID']) # multi indexed table
X = training_set[[c for c in training_set.columns if c.startswith("d__")]]
y = training_set["Category"].map({'case': False, 'control': True})

# Count Sample ID per Category in merged_meta
samples_per_category = (
    merged_meta.groupby("Category")["Sample ID"]
    .nunique()              # use .count() if you want raw row counts instead
    .reset_index(name="n_Sample ID")
    .sort_values("n_Sample ID", ascending=False)
)

print("\nSample ID per Category:")
print(samples_per_category.to_string(index=False))