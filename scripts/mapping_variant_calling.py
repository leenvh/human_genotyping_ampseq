#! /usr/bin/env python

#! /usr/bin/env python
"""
Amplicon Nanopore Variant Calling Pipeline
------------------------------------------

This script aligns Nanopore amplicon reads, calls variants using FreeBayes,
applies filtering and annotation, and summarizes SNPs and indels for each sample.
Edits from the original script include removing samclip as this got rid of amplicon reads, 
changing to minimap2 for mapping as this is faster, and allowing freebayes to do gvcf chunking.

EXAMPLE USAGE:
    python sophie_nanopore_amplicon_script_minimap2.py \
    --index-file sample_file.csv \
    --ref GCF_000001405.40_GRCh38.p14_genomic.fna \
    --gff GCF_000001405.40_GRCh38.p14_genomic.gff \
    --bed GRCh38_amplicon_targets_updated_sorted.bed \
    --clinvar clinvar_GRCh38.vcf.gz \
    --output-dir path/to/outputdir \
    --min-base-qual 30 \
    --threads 40

REQUIRED ARGUMENTS:
    --index-file        CSV file sampleIDs in 'sample' column
    --ref               Reference FASTA file (indexed with faidx, make sure no .dict file already exists in the folder)
    --gff               GFF annotation file for consequence prediction
    --bed               BED file with target amplicon regions

OPTIONAL ARGUMENTS:
    --threads           Number of threads to use (default: 10)
    --min-base-qual     Minimum base quality for variant calling (default: 20)

OUTPUT:
    - Sorted BAMs and coverage stats for each sample
    - Jointly called, filtered, and annotated SNPs and indels (VCF + TSV)
    - Per-base coverage summary across amplicon targets
"""


# ___IMPORTS___
import sys
import argparse
import subprocess as sp
import csv
import fastq2matrix as fm  # Custom module for sequence processing utilities
from fastq2matrix import run_cmd  # Helper function for running shell commands
from collections import defaultdict
import gzip
from datetime import datetime
import os
import textwrap

# ___ LOGGING FUNCTION ___
def log(msg):
    """Prints a timestamped message to stdout (also flushes for real-time logging)."""
    print(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)

# ___ MAIN FUNCTION ___
def main(args):
    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)
    # Read sample IDs from CSV index file
    samples = []
    reader = csv.DictReader(open(args.index_file))
    if "sample" not in reader.fieldnames:
        # Handles BOM issues in some CSVs
        reader = csv.DictReader(open(args.index_file, encoding='utf-8-sig'))
    for row in reader:
        if row["sample"] == "":
            continue
        samples.append(row["sample"])
    log(f"Loaded {len(samples)} samples from index file: {args.index_file}")

    # Create sequence dictionary and FASTA index
    # Manually remove .dict file if already exists before calling GATK to avoid error
    dict_path = os.path.splitext(args.ref)[0] + ".dict"
    if os.path.exists(dict_path):
        log(f"Existing .dict file found at {dict_path}. Deleting it to avoid GATK failure.")
        os.remove(dict_path)
    fm.create_seq_dict(args.ref)
    fm.faidx(args.ref)
    log(f"Reference indexing complete for: {args.ref}")
    # Check if all *_coverage_mean.txt files already exist
    coverage_files_present = all(os.path.exists(f"{sample}_coverage_mean.txt") for sample in samples)

    if coverage_files_present:
        log("All *_coverage_mean.txt files are present. Skipping alignment and coverage calculation steps.")
    else:
        # Alignment & Coverage per Sample
        for sample in samples:
            args.sample = sample
            log(f"Starting alignment for sample: {sample}")

            # Minimap2 alignment + sort BAM
            run_cmd("minimap2 -x map-ont --MD -t %(threads)s -R '@RG\\tID:%(sample)s\\tSM:%(sample)s\\tPL:nanopore' -a %(ref)s %(sample)s.fastq.gz | samtools sort -@ %(threads)s -o %(sample)s.bam -" % vars(args))
            log(f"Finished alignment for sample: {sample}")

            # Index BAM + generate QC stats
            run_cmd("samtools index %(sample)s.bam" % vars(args))
            run_cmd("samtools flagstat %(sample)s.bam > %(sample)s.flagstat.txt" % vars(args))

            # Coverage stats
            run_cmd("mosdepth -x -b %(bed)s %(sample)s --thresholds 1,10,20,30  %(sample)s.bam" % vars(args))
            run_cmd("bedtools coverage -a %(bed)s -b %(sample)s.bam -mean > %(sample)s_coverage_mean.txt" % vars(args))

        
    # Write BAM file list as bam_list.txt for the joint genotyping
    with open("bam_list.txt", "w") as O:
        for s in samples:
            O.write("%s.bam\n" % (s))

    # Ensure bam_list.txt exists and is not empty
    assert os.path.exists("bam_list.txt") and os.path.getsize("bam_list.txt") > 0, "bam_list.txt is empty!"

    # VARIANT CALLING (FreeBayes in GVCF mode)(chunking of gvcf is enabled here for Aedes genome as this is so large. For Anopheles we used to have the additional argument --gvcf-dont-use-chunk true)
    log(f"Starting variant calling with FreeBayes")
    run_cmd("freebayes -f %(ref)s -t %(bed)s -L bam_list.txt --haplotype-length -1 --min-coverage 30 --min-base-quality %(min_base_qual)s --gvcf > combined.genotyped.vcf" % vars(args))
    log(f"Finished variant calling with FreeBayes")

    # Normalize and compress VCF
    run_cmd("bcftools view --threads %(threads)s -T %(bed)s combined.genotyped.vcf | "
            "bcftools norm -f %(ref)s | bcftools sort -Oz -o combined.genotyped.vcf.gz" % vars(args))

    # VARIANT FILTERING
    log(f"Starting filtering for FMT/DP>30")
    run_cmd("bcftools filter -i 'FMT/DP>30' -S . combined.genotyped.vcf.gz | "
            "bcftools view --threads %(threads)s -i 'QUAL>30' | bcftools sort | bcftools norm -m - -Oz -o combined.genotyped_filtered_FMTDP_30.vcf.gz" % vars(args))
    log(f"Finished filtering for FMT/DP>30")

    # SNP ANNOTATION + EXPORT
    run_cmd("bcftools view --threads %(threads)s -v snps combined.genotyped_filtered_FMTDP_30.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o %(output_dir)s/snps.vcf.gz" % vars(args))

    # Index SNP VCF
    run_cmd("tabix -f %(output_dir)s/snps.vcf.gz" % vars(args))

    # Rename chromosomes
    run_cmd(textwrap.dedent(r"""
    zcat %(output_dir)s/snps.vcf.gz | \
    awk '{gsub(/NC_000001.11/, "1");  \
          gsub(/NC_000004.12/, "4");  \
          gsub(/NC_000011.10/, "11"); \
          gsub(/NC_000023.11/, "X");  \
          print;}' \
    > %(output_dir)s/snps_num.vcf
    """ % vars(args)))

    run_cmd("bgzip -c %(output_dir)s/snps_num.vcf > %(output_dir)s/snps_num.vcf.gz" % vars(args))
    run_cmd("bcftools index -f %(output_dir)s/snps_num.vcf.gz" % vars(args))

    run_cmd("SnpSift annotate %(clinvar)s %(output_dir)s/snps_num.vcf.gz > %(output_dir)s/snps_annotated.vcf" % vars(args))
    run_cmd("bcftools query %(output_dir)s/snps_annotated.vcf -f '[%%SAMPLE\\t%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%QUAL\\t%%GT\\t%%TGT\\t%%DP\\t%%AD\\t%%INFO/BCSQ\\t%%RS\\t%%CLNDN\\n]' > %(output_dir)s/combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt" % vars(args))
    run_cmd("sed -i '1iSAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tTGT\tDP\tAD\tBCSQ\tRS\tCLNDN' %(output_dir)s/combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt" % vars(args))

    # Fix RS ID for specific SNP
    run_cmd("awk -F'\t' 'BEGIN{OFS=FS} $2==\"4\" && $3==\"143781321\" && $4==\"A\" && $5==\"G\" {$12=\"186873296\"} {print}' %(output_dir)s/combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt > tmp && mv tmp %(output_dir)s/combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt" % vars(args))

    # INDEL ANNOTATION + EXPORT
    run_cmd("bcftools view --threads %(threads)s -v indels combined.genotyped_filtered_FMTDP_30.vcf.gz | bcftools csq -p a -f %(ref)s -g %(gff)s -Oz -o %(output_dir)s/indels.vcf.gz" % vars(args))

    run_cmd("tabix %(output_dir)s/indels.vcf.gz" % vars(args))

    run_cmd(textwrap.dedent(r"""
    zcat %(output_dir)s/indels.vcf.gz | \
    awk '{gsub(/NC_000001.11/, "1");  \
          gsub(/NC_000004.12/, "4");  \
          gsub(/NC_000011.10/, "11"); \
          gsub(/NC_000023.11/, "X");  \
          print;}' \
    > %(output_dir)s/indels_num.vcf
    """ % vars(args)))

    run_cmd("bgzip -c %(output_dir)s/indels_num.vcf > %(output_dir)s/indels_num.vcf.gz" % vars(args))

    run_cmd("bcftools index -f %(output_dir)s/indels_num.vcf.gz" % vars(args))

    run_cmd("SnpSift annotate %(clinvar)s %(output_dir)s/indels_num.vcf.gz > %(output_dir)s/indels_annotated.vcf" % vars(args))
    run_cmd("bcftools query %(output_dir)s/indels_annotated.vcf -f '[%%SAMPLE\\t%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%QUAL\\t%%GT\\t%%TGT\\t%%DP\\t%%AD\\t%%INFO/BCSQ\\t%%RS\\t%%CLNDN\\n]' > %(output_dir)s/combined.genotyped_filtered_FMTDP_30_formatted.indels.trans.txt" % vars(args))
    run_cmd("sed -i '1iSAMPLE\tCHROM\tPOS\tREF\tALT\tQUAL\tGT\tTGT\tDP\tAD\tBCSQ\tRS\tCLNDN' %(output_dir)s/combined.genotyped_filtered_FMTDP_30_formatted.indels.trans.txt" % vars(args))
   
    # DEPTH EXTRACTION ACROSS AMPLICON POSITIONS
    bedlines = []
    amplicon_positions = []
    for l in open(args.bed):
        row = l.strip().split()
        bedlines.append(row)
        for p in range(int(row[1]), int(row[2])):
            amplicon_positions.append((row[0], p))

    def overlap_bedlines(a, bedlines):
        """Find overlaps between per-base coverage regions and BED intervals."""
        overlaps = []
        for b in bedlines:
            if b[0] == a[0]:
                overlap = max(0, min(int(a[2]), int(b[2])) - max(int(a[1]), int(b[1])))
                if overlap > 0:
                    overlaps.append([b[0], max(int(a[1]), int(b[1])), min(int(a[2]), int(b[2]))])
        return overlaps

    # Collect depth info from mosdepth .bed.gz outputs
    dp = defaultdict(dict)
    for s in samples:
        for l in gzip.open(f"{s}.per-base.bed.gz"):
            row = l.decode().strip().split()
            overlaps = overlap_bedlines(row, bedlines)
            if overlaps:
                for overlap in overlaps:
                    for pos in range(int(overlap[1]), int(overlap[2])):
                        dp[s][(row[0], pos)] = int(row[3])

# ___ ARGUMENT PARSER ___
parser = argparse.ArgumentParser(description='Amplicon sequencing analysis script', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--index-file', type=str, help='CSV file with "sample" column for sample IDs', required=True)
parser.add_argument('--ref', type=str, help='Reference FASTA file', required=True)
parser.add_argument('--gff', type=str, help='GFF3 annotation file for consequence prediction', required=True)
parser.add_argument('--bed', type=str, help='BED file with amplicon or target regions', required=True)
parser.add_argument('--threads', default=10, type=int, help='Threads for parallel processing')
parser.add_argument('--min-base-qual', default=30, type=int, help='Minimum base quality for FreeBayes')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('--clinvar', type=str, help='ClinVar VCF file for annotation', required=True)
parser.add_argument('--output-dir', type=str, default="output", help='Folder to save output files')
parser.set_defaults(func=main)

# ___ RUN MAIN ___
args = parser.parse_args()
args.func(args)

