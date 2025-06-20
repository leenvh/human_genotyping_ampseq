import os
import argparse
import pandas as pd
from glob import glob

# Argument parser
parser = argparse.ArgumentParser(description="Summarize amplicon coverage across samples")
parser.add_argument("--input", required=True, help="Folder containing *_coverage_mean.txt files")
parser.add_argument("--output", required=True, help="Folder to save output files")
args = parser.parse_args()

# Create output folder if it doesn't exist
os.makedirs(args.output, exist_ok=True)

# Get coverage files
coverage_files = glob(os.path.join(args.input, "*_coverage_mean.txt"))

# Filter out files where the sample name starts with 'NC'
coverage_files = [
    f for f in coverage_files
    if not os.path.basename(f).replace("_coverage_mean.txt", "").startswith(("NC", "unassigned"))
]

print("Coverage files after removing those starting with 'NC':", len(coverage_files))
print(coverage_files)

all_data = []

for filepath in coverage_files:
    sample_name = os.path.basename(filepath).replace("_coverage_mean.txt", "")
    df = pd.read_csv(filepath, sep="\t", header=None, names=["CHROM", "START", "END", "AMPLICON", "DEPTH"])
    df["SAMPLE"] = sample_name
    all_data.append(df)

# Combine all rows
combined_df = pd.concat(all_data)

# Pivot into a matrix: rows = amplicons, columns = samples
matrix_df = combined_df.pivot(index="AMPLICON", columns="SAMPLE", values="DEPTH")

# Save matrix
matrix_outfile = os.path.join(args.output, "amplicon_coverage_matrix.tsv")
matrix_df.to_csv(matrix_outfile, sep="\t")

# Binary success matrix: coverage > threshold (e.g., 10x)
success_matrix = (matrix_df > 50).astype(int)

# % success per amplicon and per sample
amplicon_success_rate = success_matrix.mean(axis=1) * 100
sample_success_rate = success_matrix.mean(axis=0) * 100

# Print to terminal
print("Amplicon success %: % of samples with coverage ≥ 50x for each amplicon")
print(amplicon_success_rate)

print("\nSample success %: % of amplicons within each sample which have coverage ≥ 50x")
print(sample_success_rate)

# Save summary to file
summary_file = os.path.join(args.output, "amplicon_coverage_matrix_summary.txt")
with open(summary_file, "w") as f:
    f.write("Amplicon success (% of samples with coverage ≥ 50x for each amplicon):\n")
    f.write(amplicon_success_rate.sort_index().to_string())
    f.write("\n\nSample success (% of amplicons within each sample which have coverage ≥ 50x):\n")
    f.write(sample_success_rate.sort_index().to_string())

