# Load metadata
file_path_validation = 'Data/Supplementary_table1.xlsx'
df_validation = pd.read_excel(file_path_validation, sheet_name='Validation')

sra_data_validation = pd.read_excel(file_path_validation, sheet_name='SRA_validation')
sra_data_validation = sra_data_validation[['Sample_id', 'BioProject.Number', 'Category', 'Disease']]

def read_genus_table(path):
    with open(path, 'r', encoding='utf-8') as f: 
        lines_validation = f.readlines()

    header_idx_validation = next((i for i,l in enumerate(lines_validation) if l.lstrip().startswith('#OTU ID')), None)
    header_validation = lines_validation[header_idx_validation].lstrip('#').strip().split('\t')
    table_data_validation = pd.read_csv(path, sep='\t', header=None, names=header_validation, skiprows=header_idx_validation+1, dtype=str)

    for c in table_data_validation.columns:
        if c != 'OTU ID': table_data_validation[c] = pd.to_numeric(table_data_validation[c], errors='coerce')
    return table_data_validation

paths_validation = sorted(glob.glob(os.path.join('Feature_tables', '**', '*_final_table_depth.tsv'), recursive=True))
allowed_projects_validation = set(df_validation["BioProject.Number"].dropna().astype(str))
filtered_files_validation = [f for f in paths_validation if os.path.basename(f).split("_")[2] in allowed_projects_validation]

print(f"Found {len(filtered_files_validation)} genus tables in validation")

# CREATE ALL GENUS TABLES MERGED
merged_validation = None
for p in filtered_files_validation:
    table_data_validation = read_genus_table(p)
    sample_cols_validation = [c for c in table_data_validation.columns if c != 'OTU ID']
    if merged_validation is not None:
        dupes_validation = set(sample_cols_validation) & set([c for c in merged_validation.columns if c != 'OTU ID'])
        if dupes_validation: 
            raise ValueError(f"Duplicate sample IDs found: {dupes_validation} in {p}")
        merged_validation = merged_validation.merge(table_data_validation, on='OTU ID', how='outer')
    else: 
        merged_validation = table_data_validation.copy()

for c in merged_validation.columns:
    if c != 'OTU ID': 
        merged_validation[c] = pd.to_numeric(merged_validation[c], errors='coerce').fillna(0.0)
sample_cols_validation = [c for c in merged_validation.columns if c != 'OTU ID']
merged_validation = merged_validation[['OTU ID'] + sorted(sample_cols_validation)]

merged_validation["OTU ID"] = merged_validation["OTU ID"].apply(simplify_taxonomy)
merged_validation = merged_validation.groupby("OTU ID", as_index=False).sum(numeric_only=True)

ra_validation = merged_validation.set_index('OTU ID').copy() 

filtered_t_val = ra_validation.T.reset_index()
filtered_t_val.rename(columns={'index': 'Sample_id'}, inplace=True)

# === Merge metadata ===
validation_meta = filtered_t_val.merge(sra_data_validation, on='Sample_id', how='left')
validation_meta = validation_meta.merge(df_validation, on='BioProject.Number', how='left')
validation_meta = validation_meta[validation_meta["Sample_id"].isin(sra_data_validation["Sample_id"])]

meta_cols = [c for c in validation_meta.columns if not c.startswith("d__")]
taxa_cols = [c for c in validation_meta.columns if c.startswith("d__")]

n_t, n_s = validation_meta.shape
print(f"Number of samples: {n_t}, and number of taxa: {len(taxa_cols)}")
print(f"Number of unique studies: {validation_meta['BioProject.Number'].nunique()}")

validation_set = validation_meta.set_index(['BioProject.Number', 'Sample_id']) # multi indexed table
X_val = validation_set[[c for c in validation_set.columns if c.startswith("d__")]]
# Relative abundance!
X_val = X_val.div(X_val.sum(axis=1), axis=0) # manually checked
y_val = validation_set["Category"].map({'case': False, 'control': True})
