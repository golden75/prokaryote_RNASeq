#!/bin/bash
#SBATCH --job-name=rpkm_to_deseq
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=500M
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Copy RPKM files 
#################################################################

cp ../edgepro_gene_expression/*.rpkm_0 . 

#################################################################
# Create the deseq file
#################################################################
/isg/shared/apps/EDGE_pro/1.3.1/additionalScripts/edgeToDeseq.perl SRR034450.out.rpkm_0 SRR034451.out.rpkm_0 SRR034452.out.rpkm_0 SRR034453.out.rpkm_0
