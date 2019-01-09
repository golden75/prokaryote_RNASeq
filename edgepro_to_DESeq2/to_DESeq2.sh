#!/bin/bash
#SBATCH --job-name=to_DESeq2
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500M
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
module load python/2.7.8

python trim_epro2deseq.py deseqFile
