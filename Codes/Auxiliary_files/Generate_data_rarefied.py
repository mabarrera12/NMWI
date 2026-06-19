# These code uses the rarefied data for the diversity index

# Load metadata
file_path_rare = 'Data/Supplementary_table1.xlsx'
df_rare = pd.read_excel(file_path_rare, sheet_name='Index')
df_rare = df_rare[df_rare['Author'] != 'Miao']

def read_genus_table(path):
    with open(path, 'r', encoding='utf-8') as f: 
        lines_rare = f.readlines()

    header_idx_rare = next((i for i,l in enumerate(lines_rare) if l.lstrip().startswith('#OTU ID')), None)
    header_rare = lines_rare[header_idx_rare].lstrip('#').strip().split('\t')
    table_data_rare = pd.read_csv(path, sep='\t', header=None, names=header_rare, skiprows=header_idx_rare+1, dtype=str)

    for c in table_data_rare.columns:
        if c != 'OTU ID': 
            table_data_rare[c] = pd.to_numeric(table_data_rare[c], errors='coerce')
    return table_data_rare

# Esto es lo que hay que cambiar.
paths_rare = sorted(glob.glob(os.path.join('Feature_tables', '**', '*_table_rarefied.tsv'), recursive=True))
allowed_projects_rare = set(df_rare["BioProject.Number"].dropna().astype(str))
filtered_files_rare = [f for f in paths_rare if os.path.basename(f).split("_")[0] in allowed_projects_rare]

print(f"Found {len(filtered_files_rare)} genus tables in index")

# CREATE ALL GENUS TABLES MERGED
merged_rare = None
for p in filtered_files_rare:
    table_data_rare = read_genus_table(p)
    sample_cols_rare = [c for c in table_data_rare.columns if c != 'OTU ID']
    if merged_rare is not None:
        dupes_rare = set(sample_cols_rare) & set([c for c in merged_rare.columns if c != 'OTU ID'])
        if dupes_rare: 
            raise ValueError(f"Duplicate sample IDs found: {dupes_rare} in {p}")
        merged_rare = merged_rare.merge(table_data_rare, on='OTU ID', how='outer')
    else: 
        merged_rare = table_data_rare.copy()

for c in merged_rare.columns:
    if c != 'OTU ID': 
        merged_rare[c] = pd.to_numeric(merged_rare[c], errors='coerce').fillna(0.0)
sample_cols_rare = [c for c in merged_rare.columns if c != 'OTU ID']
merged_rare = merged_rare[['OTU ID'] + sorted(sample_cols_rare)]

merged_rare["OTU ID"] = merged_rare["OTU ID"].apply(simplify_taxonomy)
merged_rare = merged_rare.groupby("OTU ID", as_index=False).sum(numeric_only=True)

ra_sub_rare = merged_rare.set_index('OTU ID').copy() 
ra_sub_initial_rare = ra_sub_rare.copy() # raw abundance

rarefied_merged = ra_sub_initial_rare[keep_samples]
rarefied_t = rarefied_merged.T.reset_index()
rarefied_t.rename(columns={'index': 'Sample_id'}, inplace=True)

rarefied_t = rarefied_t.set_index('Sample_id')
rarefied_t = rarefied_t[taxa_cols].T
m_rarefied = rarefied_t.copy()
m_rarefied.shape
