#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Download fastq files from SRA
#################################################################
module load sratoolkit/2.8.2
fastq-dump SRR034450
fastq-dump SRR034451
fastq-dump SRR034452
fastq-dump SRR034453
