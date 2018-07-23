# prokaryote_RNASeq

This repository is a usable, publicly available tutorial for analyzing differential expression data and creating topological gene networks. All steps have been provided for the UConn CBC Xanadu cluster here with appropriate headers for the Slurm scheduler that can be modified simply to run.  Commands should never be executed on the submit nodes of any HPC machine.  If working on the Xanadu cluster, you should use sbatch scriptname after modifying the script for each stage.  Basic editing of all scripts can be performed on the server with tools such as nano, vim, or emacs.  If you are new to Linux, please use <a href="https://bio Informatics.uconn.edu/unix-basics/">this</a> handy guide for the operating system commands.  In this guide, you will be working with common bio Informatic file formats, such as <a href="https://en.wikipedia.org/wiki/FASTA_format">FASTA</a>, <a href="https://en.wikipedia.org/wiki/FASTQ_format">FASTQ</a>, <a href="https://en.wikipedia.org/wiki/SAM_(file_format)">SAM/BAM</a>, and <a href="https://en.wikipedia.org/wiki/General_feature_format">GFF3/GTF</a>. You can learn even more about each file format <a href="https://bio Informatics.uconn.edu/resources-and-events/tutorials/file-formats-tutorial/">here</a>. If you do not have a Xanadu account and are an affiliate of UConn/UCHC, please apply for one <a href="https://bio Informatics.uconn.edu/contact-us/">here</a>.


<h2 id="Header_1">Introduction and programs</h2>
This tutorial will serve as an introduction to analysis of prokaryote RNASeq data with an associated reference genome. The tools used for this analysis strictly consist of freeware and open-source software - any bioinformatician can perform the following analysis without any licenses!   

<strong>Software download links:</strong>

<a href="https://github.com/ncbi/sra-tools">SRA Toolkit</a>

<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/">FastQC</a>

<a href="https://github.com/najoshi/sickle">Sickle</a>

<a href="http://ccb.jhu.edu/software/EDGE-pro/">EDGE-pro</a>

<a href="http://www.r-project.org/">R</a>

<a href="https://www.rstudio.com/">RStudio</a>

<a href="http://bioconductor.org/packages/release/bioc/html/DESeq2.html">DESeq2</a>



<h2 id="Header_2"> Download raw reads in fastq format: SRA Toolkit</h2>
The first step is to retrieve the biological sequence data from the <a href="https://www.ncbi.nlm.nih.gov/sra">Sequence Read Archive</a>. The data that we will be using is from this <a href="https://www.ncbi.nlm.nih.gov/bioproject/PRJNA116667">experiment</a>, which analyzes two strains of Listeria monocytogenes (10403S and DsigB) with two replicates per strain, resulting in a total of four raw read files. We will use the <span style="color: #339966;">fastq-dump</span> utility from the SRA toolkit to download the raw files by accession number into fastq format. This format is necessary because the software used to perform quality control, Sickle, requires fastq files as input.

<pre style="color: silver; background: black;">module load sratoolkit/2.8.1
fastq-dump SRR034450
fastq-dump SRR034451
fastq-dump SRR034452
fastq-dump SRR034453</pre>

The converted reads are located in the <code>/UCHC/LABS/CBC/Tutorials/Listeria/raw_data</code>. The script used is located at: <code>/UCHC/LABS/CBC/Tutorials/Listeria/raw_data/sra_download.sh</code> Please note that this directory is write protected, meaning that the script will need to be run from a different directory. This is also the case for all scripts referenced throughout this tutorial - copy the scripts to a different directory before running them!

Since the fastqdump operation can be resource intensive, we will be submitting this script as a job to the cluster. Before doing so, replace the email in the line: <code><span style="color: #0000ff;">#SBATCH --mail-user= firstname.lastname@uconn.edu</span></code> with your own email. Then submit the script with the <code><span style="color: #339966;">sbatch</span></code> command:
<pre style="color: silver; background: black;"> sbatch sra_download.sh </pre>

For more information on using the Xanadu cluster and SLURM please refer to our tutorial on <a href="https://bioinformatics.uconn.edu/resources-and-events/tutorials/xanadu/">Understanding the Xanadu HPC Resource</a>



<h2 id="Header_3"> Checking the quality of the reads using FASTQC</h2>
FastQC can be used to give an impression of the quality of the data before any further analysis such as quality control. We will run FastQC over the command line on just one of the .fastq files for demonstration purposes.

<pre style="color: silver; background: black;">
fastqc [-o output dir] seqfile1

    -o --outdir     Create all output files in the specified output directory
</pre>

<pre style="color: silver; background: black;">
module load fastqc/0.11.5
fastqc --outdir . ../raw_data/SRR034450.fastq
fastqc --outdir . ../raw_data/SRR034451.fastq
fastqc --outdir . ../raw_data/SRR034452.fastq
fastqc --outdir . ../raw_data/SRR034453.fastq
</pre>


Once the program finishes running, the output will produce a series of html formatted statistics file with the name<code><span style="color: #339966;">SRR034450_fastqc.html</span></code>. 

<pre style="color: silver; background: black;">
fastqc_raw_data/
├── SRR034450_fastqc.html
├── SRR034451_fastqc.html
├── SRR034452_fastqc.html
└── SRR034453_fastqc.html</pre>

Copy this file to your desktop and open it with a web browser to view the contents, which will contain summary graphs and data such as the 'Per base sequence quality graph' below.

![](images/SRR034450_PerBaseQuality.png)


<h2 id="Header_4"> Sickle: Quality Control on raw reads</h2>

The next step is to perform quality control on the reads using sickle. Since our reads are all unpaired reads, we indicate this with the se option in the sickle command. 
<pre style="color: silver; background: black;">
Usage: sickle se [options] -f [fastq sequence file] -t [quality type] -o [trimmed fastq file]

Command:
se	single-end sequence trimming

Options:
-f  flag designates the input file 
-o  output file  
-q  the minimum quality (sickle defaults to 20) and 
-l  the minimum read length
-t  flag designates the type of read</pre> 

Unfortunately, despite the reads being Illumina reads, the average quality did not meet sickle's minimum for Illumina reads, hence the sanger option.

<pre style="color: silver; background: black;">
module load sickle/1.33

sickle se -f ../raw_data/SRR034450.fastq -t sanger -o SRR034450_trimmed.fastq
sickle se -f ../raw_data/SRR034451.fastq -t sanger -o SRR034451_trimmed.fastq
sickle se -f ../raw_data/SRR034452.fastq -t sanger -o SRR034452_trimmed.fastq
sickle se -f ../raw_data/SRR034453.fastq -t sanger -o SRR034453_trimmed.fastq</pre>

The script is called `sickle_qc.sh` which can be found in `/UCHC/LABS/CBC/Tutorials/Listeria/sickle_quality_control` folder. Following the sickle run it will produce the trimmed read files:
<pre style="color: silver; background: black;">
sickle_quality_control/
├── SRR034450_trimmed.fastq
├── SRR034451_trimmed.fastq
├── SRR034452_trimmed.fastq
└── SRR034453_trimmed.fastq</pre>


<h2 id="Header_5"> Checking the quality of the trimmed reads using FASTQC</h2>
Now we will use the FASTQC tools to check the quality of reads after trimming.
<pre style="color: silver; background: black;">
module load fastqc/0.11.5

fastqc --outdir . ../sickle_quality_control/SRR034450_trimmed.fastq
fastqc --outdir . ../sickle_quality_control/SRR034451_trimmed.fastq
fastqc --outdir . ../sickle_quality_control/SRR034452_trimmed.fastq
fastqc --outdir . ../sickle_quality_control/SRR034453_trimmed.fastq</pre>

The above script is called `fastqc_trimmed.sh` and it is located at `/UCHC/LABS/CBC/Tutorials/Listeria/fastqc_trimmed_reads` This will produce the html files with statistics, which can be downloaded and viewed.
<pre style="color: silver; background: black;">
fastqc_trimmed_reads
├── SRR034450_trimmed_fastqc.html
├── SRR034451_trimmed_fastqc.html
├── SRR034452_trimmed_fastqc.html
└── SRR034453_trimmed_fastqc.html</pre>

![](images/SRR034450_trimmed_perbase_quality.png)


<h2 id="Header_6"> EDGE-pro: Gene expression</h2>

Before we get started with EDGE-pro, we need to retrieve the Listeria reference genome and its protein and rna tables. By searching the NCBI genome database, we learn that the <a href="http://www.ncbi.nlm.nih.gov/genome/159">EGD-e strain is the reference genome</a>. We will use NCBI's ftp website: <a href="ftp://ftp.ncbi.nih.gov/">ftp://ftp.ncbi.nih.gov/</a> to download the files. Since Listeria is a bacterial genome, navigate to the genome directory then bacteria directory. Note that there are multiple genomes for Listeria -- navigate to the EGD-e assembly: <a href="ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Listeria_monocytogenes_EGD_e_uid61583/">ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Listeria_monocytogenes_EGD_e_uid61583/</a>

The reference genome is stored in the .fna file, protein table in the .ptt file, and rna table in the .rnt file. We will use the wget command line utility to download these files.

<pre style="color: silver; background: black;">
wget http://genome2d.molgenrug.nl/Bacteria/Listeria_monocytogenes_EGD_e_uid61583/NC_003210.fna
wget http://genome2d.molgenrug.nl/Bacteria/Listeria_monocytogenes_EGD_e_uid61583/NC_003210.ptt
wget http://genome2d.molgenrug.nl/Bacteria/Listeria_monocytogenes_EGD_e_uid61583/NC_003210.rnt</pre>

Once downloaded, these files are located in `/UCHC/LABS/CBC/Tutorials/Listeria/reference_genome` folder in Xanadu server.
<pre style="color: silver; background: black;">
reference_genome
├── NC_003210.fna
├── NC_003210.ptt
└── NC_003210.rnt</pre>

Now we can run EDGE-pro to generate gene expression levels. The -g option takes a reference genome, -p flag a protein table, -r flag an rnt table, -u flag a fastq file (the trimmed file from sickle). The -o flag will be the prefix before each of the EDGE-pro output files, and the -s flag is the directory where the EDGE-pro executables are located.

<pre style="color: silver; background: black;">
useage:
[OMP_NUM_THREADS=n] PATH/edge.perl [-g genome] [-p ptt] [-r rnt] [-u reads] [options]

MANDATORY FILES:
-g genome: fasta file containing bacterial genome
-p ptt: ptt file with coordinates of coding genes, in Genbank format
-r rnt: rnt file with coordinates of rRNAs and tRNAs, in Genbank format
-u reads: fastq file of reads

OPTIONAL FILES/PARAMETERS:
-s source_dir: It is a string specifying the absolute of relative path to the directory that contains all scripts. Default: working directory.
-o prefix: It is a string specifying the prefix of all output files</pre>

The Edge-Pro program can be called using following method, and it will create rpkm files for each of the trimmed input fasta files.
<pre style="color: silver; background: black;">
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
        -t 8</pre>

The above script is called <code>edgepro_rpkm.sh</code> and it can be located at <code>/UCHC/LABS/CBC/Tutorials/Listeria/edgepro_gene_expression</code>

Edge-Pro program will execute and will be producing bunch of files for each run, and we are interested in the rpkm files it will produced. 
<pre style="color: silver; background: black;">
edgepro_gene_expression/
├── SRR034450.out.alignments
├── SRR034450.out.chrNum
├── SRR034450.out.counts_0
├── SRR034450.out.genes_0
├── SRR034450.out.genomeSize
├── SRR034450.out.index.1.bt2
├── SRR034450.out.index.2.bt2
├── SRR034450.out.index.3.bt2
├── SRR034450.out.index.4.bt2
├── SRR034450.out.index.rev.1.bt2
├── SRR034450.out.index.rev.2.bt2
├── SRR034450.out.multipleCounts_0
├── SRR034450.out.numberMultiReads
├── SRR034450.out.numberReads
├── SRR034450.out.numberUniqueReads
├── SRR034450.out.readLength
├── SRR034450.out.rpkm_0
├── SRR034450.out.rpkm.numberReads
├── SRR034450.out.rRNA.numberReads
└── SRR034450.out.uniqueCounts_0</pre>


<h2 id="Header_7"> EDGE-pro to DESeq</h2>
