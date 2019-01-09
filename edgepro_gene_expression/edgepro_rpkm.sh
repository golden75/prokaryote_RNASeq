#!/bin/bash
#SBATCH --job-name=edgepro_rpkm
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=20G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=first.last@uconn.edu
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

echo `hostname`

#################################################################
# Gene Expression with Edge-Pro
#################################################################
module load EDGE_pro/1.3.1

edge.pl -g ../reference_genome/NC_003210.fna \
        -p ../reference_genome/NC_003210.ptt \
        -r ../reference_genome/NC_003210.rnt \
        -u ../sickle_quality_control/SRR034450_trimmed.fastq \
        -o SRR034450.out \
        -s /isg/shared/apps/EDGE_pro/1.3.1 \
        -t 8

edge.pl -g ../reference_genome/NC_003210.fna \
        -p ../reference_genome/NC_003210.ptt \
        -r ../reference_genome/NC_003210.rnt \
        -u ../sickle_quality_control/SRR034451_trimmed.fastq \
        -o SRR034451.out \
        -s /isg/shared/apps/EDGE_pro/1.3.1 \
        -t 8

edge.pl -g ../reference_genome/NC_003210.fna \
        -p ../reference_genome/NC_003210.ptt \
        -r ../reference_genome/NC_003210.rnt \
        -u ../sickle_quality_control/SRR034452_trimmed.fastq \
        -o SRR034452.out \
        -s /isg/shared/apps/EDGE_pro/1.3.1 \
        -t 8

edge.pl -g ../reference_genome/NC_003210.fna \
        -p ../reference_genome/NC_003210.ptt \
        -r ../reference_genome/NC_003210.rnt \
        -u ../sickle_quality_control/SRR034453_trimmed.fastq \
        -o SRR034453.out \
        -s /isg/shared/apps/EDGE_pro/1.3.1 \
        -t 8
