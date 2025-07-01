import pandas as pd

# === File paths ===
input_file = "/mnt/storage13/ahri/human_genotyping/analysis2/output/filtered_snps.txt"
output_file = "/mnt/storage13/ahri/human_genotyping/analysis2/output/per_sample_snps.csv"

# === Mapping chromosome to gene ===
gene_map = {
    '1': 'ACKR1', '4': 'Dantu', '11': 'HBB', 'X': 'G6PD',
    'NC_000001.11': 'ACKR1', 'NC_000004.12': 'Dantu',
    'NC_000011.10': 'HBB', 'NC_000023.11': 'G6PD'
}

# === Load and clean the input file ===
df = pd.read_csv(input_file, sep='\t')
df = df[df['GT'].notna()]
df = df[~df['SAMPLE'].str.startswith('NC')] 
df = df[~df['RS'].str.startswith('.')] 

# === Create gene-RS label for columns ===
df['Gene'] = df['CHROM'].map(gene_map).fillna(df['CHROM'])
df['Marker'] = df['Gene'] + "_" + df['RS'].astype(str)


# === Pivot into wide format ===
pivot_df = df.pivot(index='SAMPLE', columns='Marker', values='GT')
pivot_df.index.name = 'Sample_ID'
pivot_df.reset_index(inplace=True)

# === Replace missing values with 'NA' ===
pivot_df = pivot_df.fillna("NA")

# === Split columns into two header rows ===
columns = pivot_df.columns.tolist()
header1, header2 = [], []
for col in columns:
    if col == 'Sample_ID':
        header1.append('Sample_ID')
        header2.append('')
    else:
        gene, rsid = col.split('_', 1)
        header1.append(gene)
        header2.append("rs" + rsid if rsid.isdigit() else rsid)

# === Write to CSV with two-line header ===
with open(output_file, 'w') as f:
    f.write(','.join(header1) + '\n')
    f.write(','.join(header2) + '\n')
    pivot_df.to_csv(f, index=False, header=False)
