#!/bin/bash

echo "script start: download and initial sequencing read quality control"
date

# Get run accessions
sqlite3 -batch -noheader -csv /shared/projects/2314_medbioinfo/pascal/central_database/sample_collab.db \
"SELECT run_accession FROM sample_annot spl LEFT JOIN sample2bioinformatician s2b ON spl.patient_code = s2b .patient_code WHERE username='akarlsson';" > akarlsson_run_accessions.txt

# Create dir to save fastq files
mkdir ../data/sra_fastq

module load sra-tools
cat akarlsson_run_accessions.txt | srun --cpus-per-task=1 --time=00:30:00 xargs fastq-dump --readids --gzip \
--outdir ../data/sra_fastq/ --disable-multithreading --split-3

module load seqkit

# In analysis directory:
mkdir ./fastqc

module load fastqc
srun --cpus-per-task=2 --time=00:30:00 xargs -I{} -a akarlsson_run_accessions.txt fastqc --outdir ./fastqc/ --threads 2 --noextract ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz

# MErging paired end reads
module load flash2
mkdir ../data/merged_pairs

# Merge data files
srun --cpus-per-task=2 --time=00:30:00 xargs -a akarlsson_run_accessions.txt -n 1 -I{} flash2 --threads=2 -z \
-d ../data/merged_pairs -o {}.flash ../data/sra_fastq/{}_1.fastq.gz ../data/sra_fastq/{}_2.fastq.gz 2>&1 | tee -a akarlsson_flash2.log 


############################## Align to PhiX

# Directory for reference
mkdir ../data/reference_seqs

# Install edirect toolkit first!
# sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
efetch -db nuccore -id NC_001422 -format fasta > ../data/reference_seqs/PhiX_NC_001422.fna

module load bowtie2

# Create index database
mkdir ../data/bowtie2_DBs
srun bowtie2-build -f ../data/reference_seqs/PhiX_NC_001422.fna ../data/bowtie2_DBs/PhiX_bowtie2_DB

mkdir ./bowtie

# Align to genome
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/PhiX_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz \
 -S bowtie/akarlsson_merged2PhiX.sam --threads 8 --no-unal 2>&1 | tee bowtie/akarlsson_bowtie_merged2PhiX.log

############################## Align to SARS-CoV-2
efetch -db nuccore -id NC_045512 -format fasta > ../data/reference_seqs/SC2_NC_045512.fna

# Create index database
mkdir ../data/SC2_DBs
srun bowtie2-build -f ../data/reference_seqs/SC2_NC_045512.fna ../data/bowtie2_DBs/SC2_bowtie2_DB

# Align to genome
srun --cpus-per-task=8 bowtie2 -x ../data/bowtie2_DBs/SC2_bowtie2_DB -U ../data/merged_pairs/ERR*.extendedFrags.fastq.gz \
 -S bowtie/akarlsson_merged2SC2.sam --threads 8 --no-unal 2>&1 | tee bowtie/akarlsson_bowtie_merged2SC2.log

# Result for SARS-CoV-2: 4291 (0.07%) aligned exactly 1 time

############################## multiQC
srun multiqc --force --title "akarlsson sample sub-set" ../data/merged_pairs/ ./fastqc/ ./akarlsson_flash2.log ./bowtie/

date
echo "script end."


