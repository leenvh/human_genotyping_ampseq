import os
import argparse
import pandas as pd
from collections import defaultdict

# === Argument parser ===
parser = argparse.ArgumentParser(description="Summarize SNP genotype counts by district and G6PD genotype classification.")
parser.add_argument('--filtered-snps', required=True, help='Path to filtered_snps.txt')
parser.add_argument('--missing-rs', required=True, help='Path to missing_rs_ids_clinvar_positions.txt')
parser.add_argument('--coverage-dir', required=True, help='Directory with rs*.txt coverage files')
parser.add_argument('--metadata', required=True, help='Path to metadata.csv with participant sex info')
parser.add_argument('--output-dir', required=True, help='Directory to write output summary files')
args = parser.parse_args()

# === Setup ===
filtered_snps_file = args.filtered_snps
missing_rs_file = args.missing_rs
coverage_dir = args.coverage_dir
metadata_file = args.metadata
output_dir = args.output_dir
os.makedirs(output_dir, exist_ok=True)

# === Mapping keys ===
gene_map = {
    '1': 'ACKR1', '4': 'Dantu', '11': 'HBB', 'X': 'G6PD',
    'NC_000001.11': 'ACKR1', 'NC_000004.12': 'Dantu',
    'NC_000011.10': 'HBB', 'NC_000023.11': 'G6PD'
}
district_map = {
    'AC': 'Abobo Catholic',
    'GM': 'Maksegnit',
    'KD': 'Dembi Dollo',
    'MM': 'Mizan'
}
def get_district(sample_id):
    return district_map.get(sample_id[:2], 'Unknown')

# === Step 1: Parse filtered_snps.txt ===
snp_data = defaultdict(lambda: {
    'chrom': 'NA', 'pos': 'NA', 'ref': '', 'alt': '',
    'gene': 'NA', 'total': 0, 'ref/ref': 0, 'ref/alt': 0, 'alt/alt': 0
})
with open(filtered_snps_file, 'r') as infile:
    for line in infile:
        if not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) < 12:
            continue
        sample = parts[0]
        if sample.startswith("NC") or sample == "unassigned":
            continue
        rsid = parts[11].strip()
        if rsid in {'RS', '.'} or not rsid.isdigit():
            continue
        chrom, pos, ref, alt, genotype = parts[1], parts[2], parts[3], parts[4], parts[6]
        gene = gene_map.get(chrom, 'NA')
        district = get_district(sample)
        key = (rsid, district)
        snp_data[key]['chrom'] = chrom
        snp_data[key]['pos'] = pos
        snp_data[key]['ref'] = ref
        snp_data[key]['alt'] = alt
        snp_data[key]['gene'] = gene
        snp_data[key]['total'] += 1
        if genotype == '0/0':
            snp_data[key]['ref/ref'] += 1
        elif genotype in ('0/1', '1/0'):
            snp_data[key]['ref/alt'] += 1
        elif genotype == '1/1':
            snp_data[key]['alt/alt'] += 1

# === Step 2: Add counts from missing RSIDs ===
df_missing = pd.read_csv(missing_rs_file, sep="\t", dtype=str)
for rsid in df_missing['RS_ID'].unique():
    rs_path = os.path.join(coverage_dir, f"rs{rsid}.txt")
    if not os.path.isfile(rs_path):
        continue
    try:
        df_cov = pd.read_csv(rs_path, sep='\t')
        df_cov = df_cov[~df_cov['SAMPLE'].str.startswith('NC')]
        df_cov = df_cov[df_cov['SAMPLE'] != 'unassigned']
        df_cov = df_cov[df_cov['COV'] > 50]
        df_cov['district'] = df_cov['SAMPLE'].apply(get_district)
    except Exception as e:
        print(f"Warning reading {rs_path}: {e}")
        continue
    for district, sub_df in df_cov.groupby('district'):
        total_samples = sub_df.shape[0]
        row = df_missing[df_missing['RS_ID'] == rsid].iloc[0]
        chrom = row['CHROM']
        gene = gene_map.get(chrom, 'NA')
        key = (rsid, district)
        if key not in snp_data:
            snp_data[key] = {
                'chrom': chrom,
                'pos': row['POS'],
                'ref': row['REF'],
                'alt': row['ALT'],
                'gene': gene,
                'total': total_samples,
                'ref/ref': total_samples, 'ref/alt': '0', 'alt/alt': '0'
            }

# === Step 3: Load metadata ===
df_meta = pd.read_csv(metadata_file)
sample_sex_map = dict(zip(df_meta['participants_code'], df_meta['gender']))

# === Step 4: Write outputs ===
out_g6pd = os.path.join(output_dir, 'snp_summary_g6pd.txt')
out_other = os.path.join(output_dir, 'snp_summary_ackr1_dantu_hbb.txt')
with open(out_other, 'w') as out1, open(out_g6pd, 'w') as out2:
    out1.write("RS_ID\tCHROM\tPOS\tREF\tALT\tGene\tDistrict\tTotal_Samples\tRef/Ref\tRef/Alt\tAlt/Alt\n")
    out2.write("RS_ID\tCHROM\tPOS\tREF\tALT\tGene\tDistrict\tTotal_Samples\tMale_Hemizygote\tFemale_Het\tFemale_Hom\tWT\n")

    for (rsid, district), d in sorted(snp_data.items(), key=lambda x: (int(x[0][0]), x[0][1])):
        if d['gene'] != "G6PD":
            out1.write(f"{rsid}\t{d['chrom']}\t{d['pos']}\t{d['ref']}\t{d['alt']}\t{d['gene']}\t{district}\t"
                       f"{d['total']}\t{d['ref/ref']}\t{d['ref/alt']}\t{d['alt/alt']}\n")
        else:
            male_hemi = female_het = female_hom = wt = 0
            with open(filtered_snps_file, 'r') as infile:
                for line in infile:
                    if not line.strip():
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 12:
                        continue
                    sample = parts[0]
                    if sample.startswith("NC") or sample == "unassigned":
                        continue
                    sample_district = get_district(sample)
                    sex = sample_sex_map.get(sample, 'NA')
                    this_rsid = parts[11].strip()
                    if (this_rsid != rsid) or (sample_district != district):
                        continue
                    genotype = parts[6]
                    if sex == "Male":
                        if genotype == '1/1':
                            male_hemi += 1
                        elif genotype == '0/0':
                            wt += 1
                    elif sex == "Female":
                        if genotype == '0/0':
                            wt += 1
                        elif genotype in ('0/1', '1/0'):
                            female_het += 1
                        elif genotype == '1/1':
                            female_hom += 1
            total = male_hemi + female_het + female_hom + wt
            if total == 0 and d['total'] > 0:
                total = d['total']
                wt = total
                male_hemi = female_het = female_hom = 0
            out2.write(f"{rsid}\t{d['chrom']}\t{d['pos']}\t{d['ref']}\t{d['alt']}\t{d['gene']}\t{district}\t"
                       f"{total}\t{male_hemi}\t{female_het}\t{female_hom}\t{wt}\n")
