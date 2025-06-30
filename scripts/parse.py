from collections import defaultdict
import os
import pandas as pd

# === Step 1: Define CHROM → Gene mapping ===
gene_map = {
    '1': 'ACKR1',
    '4': 'Dantu',
    '11': 'HBB',
    'X': 'G6PD',
    'NC_000001.11': 'ACKR1',
    'NC_000004.12': 'Dantu',
    'NC_000011.10': 'HBB',
    'NC_000023.11': 'G6PD'
}

# === Step 2: Parse filtered_snps.txt ===
snp_data = defaultdict(lambda: {
    'chrom': 'NA',
    'pos': 'NA',
    'ref': '',
    'alt': '',
    'total': 0,
    'ref/ref': 0,
    'ref/alt': 0,
    'alt/alt': 0,
    'gene': 'NA'
})

filtered_snps_file = '/mnt/storage13/ahri/human_genotyping/analysis2/output/filtered_snps.txt'
with open(filtered_snps_file, 'r') as infile:
    for line in infile:
        if line.strip() == '':
            continue
        parts = line.strip().split('\t')
        if len(parts) < 12:
            continue

        rsid = parts[11].strip()
        if rsid in {'RS', '.'} or not rsid.isdigit():
            continue

        chrom = parts[1]
        pos = parts[2]
        ref = parts[3]
        alt = parts[4]
        genotype = parts[6]

        gene = gene_map.get(chrom, 'NA')

        snp_data[rsid]['chrom'] = chrom
        snp_data[rsid]['pos'] = pos
        snp_data[rsid]['ref'] = ref
        snp_data[rsid]['alt'] = alt
        snp_data[rsid]['gene'] = gene
        snp_data[rsid]['total'] += 1

        if genotype == '0/0':
            snp_data[rsid]['ref/ref'] += 1
        elif genotype in ('0/1', '1/0'):
            snp_data[rsid]['ref/alt'] += 1
        elif genotype == '1/1':
            snp_data[rsid]['alt/alt'] += 1

# === Step 3: Add missing RS IDs with coverage info ===
coverage_dir = "/mnt/storage13/ahri/human_genotyping/analysis2/output"
missing_rs_file = os.path.join(coverage_dir, "missing_rs_ids_clinvar_positions.txt")
df_missing = pd.read_csv(missing_rs_file, sep="\t", dtype=str)
unique_rsids = df_missing['RS_ID'].unique()

for rsid in unique_rsids:
    if rsid not in snp_data:
        rs_file_path = os.path.join(coverage_dir, f"rs{rsid}.txt")
        total_samples = 'NA'
        if os.path.isfile(rs_file_path):
            try:
                df_cov = pd.read_csv(rs_file_path, sep='\t')
                # Exclude SAMPLEs that start with "NC" or are "unassigned"
                df_cov = df_cov[~df_cov['SAMPLE'].str.startswith('NC')]
                df_cov = df_cov[df_cov['SAMPLE'] != 'unassigned']
                total_samples = (df_cov['COV'] > 50).sum()
            except Exception as e:
                print(f"Warning: Could not read file {rs_file_path}: {e}")
                total_samples = 'NA'

        row = df_missing[df_missing['RS_ID'] == rsid].iloc[0]
        chrom = row['CHROM']
        gene = gene_map.get(chrom, 'NA')
        snp_data[rsid] = {
            'chrom': chrom,
            'pos': row['POS'],
            'ref': row['REF'],
            'alt': row['ALT'],
            'total': total_samples,
            'ref/ref': 'NA',
            'ref/alt': 'NA',
            'alt/alt': 'NA',
            'gene': gene
        }


# === Step 4: Write final summary ===
summary_out = os.path.join(coverage_dir, 'snp_summary.txt')
with open(summary_out, 'w') as outfile:
    outfile.write("RS_ID\tCHROM\tPOS\tREF\tALT\tGene\tTotal_Samples\tRef/Ref\tRef/Alt\tAlt/Alt\n")
    for rsid in sorted(snp_data.keys(), key=int):
        d = snp_data[rsid]
        outfile.write(f"{rsid}\t{d['chrom']}\t{d['pos']}\t{d['ref']}\t{d['alt']}\t{d['gene']}\t{d['total']}\t{d['ref/ref']}\t{d['ref/alt']}\t{d['alt/alt']}\n")
