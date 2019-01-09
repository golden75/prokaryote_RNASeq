#!/bin/bash
#SBATCH --job-name=fastqc_trimmed
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --partition=cbcworkshop
#SBATCH --qos=cbcworkshop
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Quality Check of the Reads
#################################################################
module load fastqc/0.11.5

fastqc --outdir . ../sickle_quality_control/SRR034450_trimmed.fastq
fastqc --outdir . ../sickle_quality_control/SRR034451_trimmed.fastq
fastqc --outdir . ../sickle_quality_control/SRR034452_trimmed.fastq
fastqc --outdir . ../sickle_quality_control/SRR034453_trimmed.fastq
