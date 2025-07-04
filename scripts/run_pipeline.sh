# Basecalling using dorado and demultiplexing by nanopore barcodes. Needs to be run on G1 server which has dorado installed

conda activate dorado

dorado_pipe \
dorado --sample_sheet /mnt/storageG1/data/experiments/exp181_humanamplicons_AHRI/sample_sheet_exp181.csv \
--barcoding_kit SQK-NBD114-96 \
--ID exp181_humanamplicons_AHRI \
--data_directory /mnt/storageG1/data/experiments/exp181_humanamplicons_AHRI/exp181/  \
--output_directory /mnt/storageG1/data/experiments/exp181_humanamplicons_AHRI/exp181/rerun_nanopore_pipeline_out_exp181/ \
--basecalling_model dna_r10.4.1_e8.2_400bps_sup@v5.0.0 \
--ref_seq /mnt/storageG1/data/experiments/exp181_humanamplicons_AHRI/GCF_000001405.40_GRCh38.p14_genomic.fna \
--kraken false \
--threads 20 \
--trimmer none

# Transfer output folder to s13

# Demultiplexing by plate barcodes (barcodes included in the primers, e.g. BC1-BC11 etc)
conda activate fastq2matrix

python /path/to/demux_nanopore_plates_edgesize100.py \
-p plate_layout.csv \
-b internal_barcodes.csv \
-f /mnt/storage13/ahri/human_genotyping/rerun_nanopore_pipeline_out_exp181/fastq/fastq \
-o samples.csv \
--demux-script /path/to/demux_nanopore_amplicon.py \
-m 0 \
-t 8

# bgzip fastq files and make sample_file.csv
for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done
ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv

# Mapping and variant calling
python /path/to/mapping_variant_calling.py \
--index-file /path/to/sample_file.csv \
--ref /path/to/GCF_000001405.40_GRCh38.p14_genomic.fna \
--gff /path/to/GCF_000001405.40_GRCh38.p14_genomic.gff \
--bed /path/to/GRCh38_amplicon_targets_updated_sorted.bed \
--clinvar /path/to/clinvar_GRCh38.vcf.gz \
--output-dir path/to/outputdir \
--min-base-qual 30 \
--threads 40

# Make coverage plots
Rscript /path/to/coverage_plots.R /path/to/outputdir/combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt /path/to/outputdir/

# Create amplicon coverage matrix
python /path/to/create_ampliconxsample_coverage_matrix.py --input /path/to/coverage_files --output /path/to/outputdir

# Compute SNP allele frequencies
python /path/to/compute_SNP_allele_frequencies.py \
  --snp-file path/to/combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt \
  --coverage-matrix path/to/amplicon_coverage_matrix.tsv \
  --bed-file path/to/GRCh38_amplicon_targets_updated_sorted_num.bed \
  --output path/to/outputdir/snp_frequencies_filtered_by_coverage_and_altDP20.tsv


# Filter combined.genotyped_filtered_FMTDP_30_formatted.snps.trans to exclude samples with ./. Or . In GT, and remove unassigned:
awk -F'\t' 'NR==1 || ($1 != "unassigned" && $7 != "./." && $7 != ".")' path/to/combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt > path/to/outputdir/filtered_snps.txt


# Extract coverage files for rsIDs that are missing from data
python /path/to/extract_missing_rs_coverage_by_sample.py \
  --positions path/to/missing_rs_ids_clinvar_positions.txt \
  --coverage-dir path/to/*.coverage.txt files \
  --output-dir path/to/outputdir


# Make snp summary files including missing rsIDS and per district
python /path/to/parse.py \
  --filtered-snps /path/to/filtered_snps.txt \
  --missing-rs /path/to/missing_rs_ids_clinvar_positions.txt \
  --coverage-dir /path/to/rs*.txt coverage files \
  --metadata /path/to/metadata.csv \
  --output-dir /path/to/outputdir

# Make pie chart figures
Rscript /path/to/piecharts.R \
  path/to/snp_summary_ackr1_dantu_hbb.txt \
  path/to/snp_summary_g6pd.txt \
  path/to/snp_genotype_pies_ackr1_dantu_hbb.pdf \
  path/to/snp_genotype_pies_g6pd.pdf
