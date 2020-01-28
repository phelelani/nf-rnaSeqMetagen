# nf-rnaSeqMetagen
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/728)

*nf-rnaSeqMetagen* is a [Nextflow](http://nextflow.io/)


To use the `nf-rnaSeqMetagen` pipeline, the following dependencies are required:
   1. Installed softwares:
      - [`Nextflow`](https://www.nextflow.io/)
      - [`Singularity`](http://singularity.lbl.gov/)
   2. `Singularity` [containers](https://www.singularity-hub.org/collections/770) with the required applications/programs for executing the workflow:
      - `nf-rnaSeqMetagen-kraken2.sif`
      - `nf-rnaSeqMetagen-multiqc.sif`
      - `nf-rnaSeqMetagen-star.sif`
      - `nf-rnaSeqMetagen-trinity.sif`
      - `nf-rnaSeqMetagen-upset.sif`
   3. Reference genome, annotation and indexes
      - Reference genome (`.fa`/`.fasta`) and genome annotation (`.gtf`) files.
      - Reference genome indexes (`bowtie2` & `STAR` - see *1.3.* below on how to generate the indexes).
 
---

<p align="center">
  <img width="600" src="nf-rnaSeqMetagen.png">
</p>

## 1. Obtaining the `nf-rnaSeqMetagen` pipeline and preparing data
First, you need to clone the `nf-rnaSeqMetagen` repository onto you machine. You can either use `git` or `nextflow` (see the two methods below). I recommend using `nextflow` and creating you own `config` file (will explain later) for executing the workflow in the directory of your choosing. The rest of this documentation assumes that you have used `nextflow` to clone this workflow - If your're an expert and have used `git` to clone the workflow - you know what to do :)
```bash
## Using nextflow
nextflow pull https://github.com/phelelani/nf-rnaSeqMetagen
```
Content of the repository (will be in "$HOME/.nextflow/assets/phelelani/nf-rnaSeqCount"):
```bash
nf-rnaSeqMetagen
|--containers                       ## Folder for Singularity images and recipes (in case you want to build yourself). All downloaded images go here!
|  |--Singularity.kraken2           ## Singularity recipe file for
|  |--Singularity.multiQC           ## Singularity recipe file for 
|  |--Singularity.star              ## Singularity recipe file for 
|  |--Singularity.trinity           ## Singularity recipe file for 
|  |--Singularity.upset             ## Singularity recipe file for
|--templates                        ## Folder for extra scripts for the pipeline.
|  |--create_matrix.R               ## Script for 
|  |--get_taxons.sh                 ## Script for 
|--LICENSE                          ## Duh!
|--main.config                      ## User configuration file! All inputs, outputs and options GO HERE!! ONLY file that SHOULD be modified by user!
|--main.nf                          ## Main nf-rnaSeqMetagen nextflow scripts.
|--nextflow.config                  ## Pipeline configuration file! DO NOT EDIT!!!
|--nf-rnaSeqMetagen.png             ## Pipeline flow diagram
|--README.md                        ## Duh!
```
To get the `help menu` for the workflow, execute the following from anywherre on your system aftercloning the repository:
```
nextflow run nf-rnaSeqMetagen --help
```
The command above will give you the following usage information and options for running the `nf-rnaSeqMetagen` workflow:
```
====================================================================================================
#####################################  nf-rnaSeqMetagen v0.2   #####################################
====================================================================================================

USAGE:
nextflow run nf-rnaSeqMetagen -profile "slurm" --data "/path/to/data" --genome "/path/to/genome.fa" --genes "/path/to/genes.gtf"

HELP:
nextflow run nf-rnaSeqMetagen --help

MANDATORY ARGUEMENTS:
-profile     STRING    Executor to be used. Available options:
				"standard"          : Local execution (no job scheduler).
                "slurm"             : SLURM scheduler.
--data       FOLDER    Path to where the input data (FASTQ files) is located. Supported FASTQ files:
				[ fastq | fastq.gz | fastq.bz2 | fq | fq.gz | fq.bz2 ]
--genome     FILE      The whole genome FASTA sequence. Supported FASTA files:
    			[ fasta | fa | fna ]
--genes      FILE      The genome annotation GFT file. Supported GTF file:
				[ gtf ]
--db         FOLDER    Path to where the Kraken2 database will be saved (or where it is located if already created).
                       Default: $PWD/kraken2db
--mode       STRING    To specify which step of the workflow you are running (see https://github.com/phelelani/nf-rnaSeqMetagen).
                       Available options:
				"prep.Containers"   : For downloading Singularity containers used in this workflow.
                "prep.STARIndex"    : For indexing your reference genome using STAR.
                "prep.BowtieIndex"  : For indexing your reference genome using Bowtie2.
                "prep.KrakenDB"     : For building the Kraken2 database.
                "run.FilterClassify": For performing metagenomics analysis, i.e., filtering and classification.

OPTIONAL ARGUEMENTS:
--help                 To show this menu.
--out        FOLDER    Path to where the output should be directed.
                       Default: $PWD/results_nf-rnaSeqMetagen
--pairedEnd            If working with paired-end FASTQ files (default).
--singleEnd            If working with single-end FASTQ files.
--max_memory STRING    Maximum memory you have access to.
                       Default: "200.GB"
--max_cpus   STRING    Maximum CPUs you have access to.
                       Default: "24"
--max_time   STRING    Maximum time you have access to.
                       Default: "24.h"
====================================================================================================
```

---

### 1.1. Download test datasets (optional)
We will now download the reference genome (along with its annotation file) from Ensembl. We will also download the FASTQ files from the H3ABioNet site, which we will analyse using the `nf-rnaSeqMetagen` workflow. *__NB__: Skip this section if you have your own data to analyse using this workflow! This section is only for getting data to practice using the `nf-rnaSeqMetagen` workflow!* 

- [x] Download and decompress the mouse reference genome along with its annotation:
```
## Make a directory for the reference genome:
mkdir reference

## Download the reference genome (FASTA) and annotation file (GTF) files and put them into the newlly created directory:
wget -c -O reference/genome.fa.gz ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz
wget -c -O reference/genes.gtf.gz ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
gunzip reference/genome.fa.gz
gunzip reference/genes.gtf.gz
```

- [x] Download RNA-seq test dataset from H3ABioNet:
```
## Make a directory for the data:
mkdir data

## Download the data:
for sample in sample{37..42}_R{1,2}.fastq.gz; do wget -c -O data/$sample http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/$sample; done
```
### 1.2. Download the `Singularity` containers (required to execute the pipeline):
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.Containers
```

### 1.3. Generating genome indexes.
To generate the `STAR` and `Bowtie2` genome indexes, run the following commands:
```bash
## Generate STAR indexes
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.STARIndex --genome "reference/genome.fa" --genes "reference/genes.gtf"

## Generate Bowtie2 indexes:
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.BowtieIndex --genome "reference/genome.fa" --genes "reference/genes.gtf"
```

### 1.4. Creating the Kraken2 database:
To create the Kraken2 database, run the following command:
```bash
## Create Kraken2 database
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.KrakenDB
```

We are now ready to execute the workflow!

---

## 2. Executing the main `nf-rnaSeqMetagen` pipeline
As seen on the `help menu` above, there are a couple of options that you can use with this workflow. It can become a bit tedious and confusing having to specify these commands everytime you have to execute the each section for the analysis. To make your life easier, we will create a configuration script that we will use in this tutorial (we will pass this using the `-c` option of `nextflow`). You can name it whatever you want, but for now, lets call it `myparams.config`. We will add the mandatory arguements for now, but as you become more farmiliar with the workflow - you can experiment with other options. You can use your favourite text editor to create the `myparams.config` file. Copy and paste the the parameters below:
```
params {
    data = "data/"
    genome = "reference/genome.fa"
    genes = "reference/genes.fa"
}
```
Obviously - the above `myparams.config` assumes that you have been following this tutorial. If you have your data lying around somewhere in your system, you need to put the full path to where your the `data`, `genome` and `genes` files are. Since the `--mode` will keep changing, we will add this on the command as we do the analysis. Now that we have the mandatory arguements in our `myparams.config`, lets do some analysis

### 2.1. Read Filtering and Classification:
To perform filtering of host reads and classification of exogeneous reads, use this command:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode run.FilterClassify -c myparams.config
```

CONGRATULATIONS for getting this far!! :) You can now explore the results and use the read counts to perform differential expression analysis!

---

## 3. Explore `nf-rnaSeqMetagen` results

<!-- ``` -->
<!-- - [1] Read QC (optional)         =>    `<output_directory>/1_RQC` -->
<!-- - [2] Read Trimming (optional)   =>    `<output_directory>/2_Read_Trimming` -->
<!-- - [3] Read Alignment             =>    `<output_directory>/3_Read_Alignment` -->
<!-- - [4] Read Counting              =>    `<output_directory>/4_Read_Counts` -->
<!-- - [5] MultiQC                    =>    `<output_directory>/5_MultiQC -->
<!-- - [6] Workflow tracing           =>    `<output_directory>/workflow-tracing -->
<!-- ``` -->
<!-- In addition to the 5 directories created for each step in the results directory, a directory `workflow-tracing` is created to monitor the resources used in each step. This directory will contain 4 files for each step (--mode) of the workflow: -->
<!-- - `nf-rnaSeqMetagen_<mode>_report.html` -->
<!-- - `nf-rnaSeqMetagen_<mode>_timeline.html` -->
<!-- - `nf-rnaSeqMetagen_<mode>_trace.txt` -->

<!-- These files contain detailed information on the resources (CPU, MEMORY and TIME) usage of each of the process in the different pipeline steps. The `<output_directory>` directory structure is summarized below: -->

<!-- ```bash -->
<!-- <output_directory> -->
<!-- |--1_Read_QC -->
<!-- |  |--<sample_1>_R1.fastqc.html .. <sample_N>_R1.fastqc.html -->
<!-- |  |--<sample_1>_R2.fastqc.html .. <sample_N>_R1.fastqc.html -->
<!-- |--2_Read_Trimming -->
<!-- |  |--<sample_1>.1P.fastq.gz .. <sample_N>.1P.fastq.gz -->
<!-- |  |--<sample_1>.2P.fastq.gz .. <sample_N>.2P.fastq.gz -->
<!-- |--3_Read_Alignment -->
<!-- |  |--<sample_1>_Aligned.out.bam .. <sample_N>_Aligned.out.bam -->
<!-- |  |--<sample_1>_Log.final.out .. <sample_N>_Log.final.out -->
<!-- |  |--<sample_1>_Log.out .. <sample_N>_Log.out -->
<!-- |  |--<sample_1>_Log.progress.out .. <sample_N>_Log.progress.out -->
<!-- |  |--<samplle_1>_SJ.out.tab .. <sample>_SJ.out.tab -->
<!-- |--4_Read_Counts -->
<!-- |  |--featureCounts -->
<!-- |  |  |--gene_counts_final.txt -->
<!-- |  |  |--gene_counts.txt -->
<!-- |  |  |--gene_counts.txt.jcounts -->
<!-- |  |  |--gene_counts.txt.summary -->
<!-- |  |--htseqCounts -->
<!-- |  |  |--gene_counts_final.txt -->
<!-- |  |  |--<sample>.txt .. <sample>.txt -->
<!-- |--5_MultiQC -->
<!-- |  |--multiqc_data -->
<!-- |  |--multiqc_report.html -->
<!-- |--workflow-tracing -->
<!-- |  |--nf-rnaSeqMetagen_run.MultiQC_{report.html,timeline.html,trace.txt} -->
<!-- |  |--nf-rnaSeqMetagen_run.ReadAlignment_{report.html,timeline.html,trace.txt} -->
<!-- |  |--nf-rnaSeqMetagen_run.ReadCounting_{report.html,timeline.html,trace.txt} -->
<!-- |  |--nf-rnaSeqMetagen_run.ReadTrimming_{report.html,timeline.html,trace.txt} -->
<!-- |  |--nf-rnaSeqMetagen_run.ReadQC_{report.html,timeline.html,trace.txt} -->
<!-- ``` -->
<!-- **NB:** I am working on further improving the pipleine and the associated documentation, feel free to share comments and suggestions! -->

<!-- --- -->
