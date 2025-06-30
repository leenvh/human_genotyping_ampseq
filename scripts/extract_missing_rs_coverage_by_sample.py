import os
import pandas as pd
from glob import glob
from collections import defaultdict

# Path setup
positions_file = "/mnt/storage13/ahri/human_genotyping/analysis2/output/missing_rs_ids_clinvar_positions.txt"
coverage_dir = "/mnt/storage13/ahri/human_genotyping/analysis2"
output_dir = "/mnt/storage13/ahri/human_genotyping/analysis2/output"
os.makedirs(output_dir, exist_ok=True)

# Mapping CHROM (if needed) — adjust to your reference naming convention
chrom_map = {
    'X': 'NC_000023.11',
    '1': 'NC_000001.11',
    '11': 'NC_000011.10',
    # Add more if needed
}

# Load RS_IDs and positions
df_rs = pd.read_csv(positions_file, sep='\t')
df_rs['CHROM_mapped'] = df_rs['CHROM'].astype(str).map(chrom_map).fillna(df_rs['CHROM'])

# Organize positions by RS_ID
rs_positions = defaultdict(list)
for _, row in df_rs.iterrows():
    rs_positions[row['RS_ID']].append((str(row['CHROM_mapped']), int(row['POS'])))

# Get all coverage files
coverage_files = glob(os.path.join(coverage_dir, "*.coverage.txt"))

# For each RS_ID, collect matching rows
for rsid, pos_list in rs_positions.items():
    rsid_outfile = os.path.join(output_dir, f"rs{rsid}.txt")
    with open(rsid_outfile, 'w') as out:
        out.write("REF\tPOS\tCOV\tA\tC\tG\tT\tDEL\tREFSKIP\tSAMPLE\n")
        for cov_file in coverage_files:
            try:
                df_cov = pd.read_csv(cov_file, sep='\t')
            except Exception as e:
                print(f"Could not read {cov_file}: {e}")
                continue
            for chrom, pos in pos_list:
                match = df_cov[(df_cov['REF'] == chrom) & (df_cov['POS'] == pos)]
                if not match.empty:
                    match.to_csv(out, sep='\t', header=False, index=False)
