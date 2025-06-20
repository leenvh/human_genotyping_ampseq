# Basecalling using dorado and demultiplexing by nanopore barcodes. Needs to be run on G1 server which has dorado installed

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

# Demultiplexing by plate barcodes (barocdes included in the primers, e.g. BC1-BC11 etc)
python demux_nanopore_plates_edgesize100.py \
-p plate_layout.csv \
-b internal_barcodes.csv \
-f /mnt/storage13/ahri/human_genotyping/rerun_nanopore_pipeline_out_exp181/fastq/fastq \
-o samples.csv \
-m 0 \
-t 8

# bgzip fastq files and make sample_file.csv
for f in *.fastq ; do bgzip -c "$f" > "${f%.*}.fastq.gz" ; done
ls *.fastq | sed 's/.fastq//' > samples.txt
echo -e "sample" | cat - samples.txt > samples_header.txt
sed 's/ \+/,/g' samples_header.txt > sample_file.csv

# Mapping and variant calling
python mapping_variant_calling.py \
--index-file sample_file.csv \
--ref GCF_000001405.40_GRCh38.p14_genomic.fna \
--gff GCF_000001405.40_GRCh38.p14_genomic.gff \
--bed GRCh38_amplicon_targets_updated_sorted.bed \
--min-base-qual 30 \
--threads 40

# Add rs code for Dantu variant manually, as it is not included in the clinvar file
awk -F'\t' 'BEGIN{OFS=FS} 
$2=="4" && $3=="143781321" && $4=="A" && $5=="G" {$12="186873296"} 
{print}' combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt > tmp && mv tmp combined.genotyped_filtered_FMTDP_30_formatted.snps.trans.txt

# Create amplicon coverage matrix
python create_ampliconxsample_coverage_matrix.py

# Compute SNP allele frequencies
python compute_SNP_allele_frequencies.py

