
Workflow for analysis of nanopore amplicon sequencing data to genotype variants in the G6PD, HBB, and ACKR1 (Duffy) genes, as well as the Dantu variant.<br>

Assay methods adapted from https://pubmed.ncbi.nlm.nih.gov/37495620/ <br>
Scripts adapted from https://github.com/LSHTMPathogenSeqLab/amplicon-seq/tree/main and https://github.com/sophiemoss/smoss_ampseq/tree/main

## Workflow
Follow [run_pipeline.sh](scripts/run_pipeline.sh)

## Input files
- Bed file: [GRCh38_amplicon_targets_updated_sorted.bed](input_files/GRCh38_amplicon_targets_updated_sorted.bed)
- Bed file with numeric chromosomes: [GRCh38_amplicon_targets_updated_sorted_num.bed](input_files/GRCh38_amplicon_targets_updated_sorted_num.bed)
- Internal barcodes example file: [internal_barcodes.csv](input_files/internal_barcodes.csv)
- Plate layout example file: [plate_layout.csv](input_files/plate_layout.csv)
- Reference genome (size too large for upload): download from https://www.ncbi.nlm.nih.gov/genome/guide/human or scp from /mnt/storage13/ahri/human_genotyping/github/input_files/GCF_000001405.40_GRCh38.p14_genomic.fna
- gff file (size too large for upload): download from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/ or scp from /mnt/storage13/ahri/human_genotyping/github/input_files/GCF_000001405.40_GRCh38.p14_genomic.gff
- Clinvar file (size too large for upload): download via below commands or scp from /mnt/storage13/ahri/human_genotyping/github/input_files/clinvar_GRCh38.vcf.gz
<pre>curl -s ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz > clinvar_GRCh38.vcf.gz
tabix -f clinvar_GRCh38.vcf.gz<pre> 
