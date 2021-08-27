[![GitHub license](https://img.shields.io/github/license/phelelani/nf-rnaSeqCount)](https://github.com/phelelani/nf-rnaSeqCount/blob/master/LICENSE) [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)

[biotools:nf-rnaseqmetagen](https://bio.tools/nf-rnaseqmetagen)

`nf-rnaSeqMetagen` is a [Nextflow](http://nextflow.io/)
        
To use the `nf-rnaSeqMetagen` pipeline, the following are required:
1. Software dependencies:
   - [`Nextflow`](https://www.nextflow.io/)
   - [`Singularity`](http://singularity.lbl.gov/)
2. RNA-seq data (paired-end for now - support for single-ended reads to follow)
3. Reference genome (FASTA sequences) and its annotation file (GFT)

---

<p align="center">
  <img width="832px" src="assets/images/nf-rnaSeqMetagen.png">
</p>

---

## 1. Obtaining the `nf-rnaSeqMetagen` Pipeline and Preparing Data
First, you need to clone the `nf-rnaSeqMetagen` repository onto you machine. You can either use `git` or `nextflow`. I recommend you use `nextflow`. The rest of this documentation assumes that you have used `nextflow` to clone this workflow.
```bash
nextflow pull https://github.com/phelelani/nf-rnaSeqMetagen
```
<script id="asciicast-308777" src="https://asciinema.org/a/308777.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="6" data-speed="3" data-loop="0"></script>

Content of the repository (located in `$HOME/.nextflow/assets/phelelani/nf-rnaSeqCount`):
<script id="asciicast-308808" src="https://asciinema.org/a/308808.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="12" data-speed="3" data-loop="0"></script>

To get the `help menu` for the workflow, execute the following command from anywherre on your system:
```
nextflow run nf-rnaSeqMetagen --help
```
<script id="asciicast-308804" src="https://asciinema.org/a/308804.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>

---

### 1.1. Download test datasets (optional)
We will now download the reference genome (along with its annotation file) from Ensembl. We will also download the FASTQ files from the H3ABioNet site, which we will analyse using the `nf-rnaSeqMetagen` workflow. *__NB__: Skip this section if you have your own data to analyse using this workflow! This section is only for getting data to practice using the `nf-rnaSeqMetagen` workflow!* 

Make directories:
```
mkdir example
cd example
mkdir reference
mkdir data
```
<script id="asciicast-308945" src="https://asciinema.org/a/308945.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="8" data-speed="3" data-loop="0"></script>

Download and decompress the mouse reference genome along with its annotation:
```
wget -c -O reference/genome.fa.gz ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz
```
<script id="asciicast-308949" src="https://asciinema.org/a/308949.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="18" data-speed="3" data-loop="0"></script>


```
wget -c -O reference/genes.gtf.gz ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
```
<script id="asciicast-308953" src="https://asciinema.org/a/308953.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="18" data-speed="3" data-loop="0"></script>

```
gunzip reference/genome.fa.gz
gunzip reference/genes.gtf.gz
```
<script id="asciicast-308955" src="https://asciinema.org/a/308955.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="8" data-speed="3" data-loop="0"></script>

Download RNA-seq test dataset from H3ABioNet: <a href="examples/data/get_data.sh" target="_blank">script</a>.
```
cd data
wget https://phelelani.github.io/nf-rnaSeqMetagen/examples/data/get_data.sh
```
<script id="asciicast-309156" src="https://asciinema.org/a/309156.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="14" data-speed="3" data-loop="0"></script>

```
sh get_data.sh
ls -l 
cd ..
```
<script id="asciicast-309177" src="https://asciinema.org/a/309177.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="12" data-speed="7" data-loop="0"></script>


### 1.2. Download the `singularity` containers (required to execute the pipeline):
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.Containers
```
<script id="asciicast-308816" src="https://asciinema.org/a/308816.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>

### 1.3. Generating genome indexes.
To generate the `STAR` genome indexes, run the following commands:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.GenomeIndexes --genome "$PWD/reference/genome.fa" --genes "$PWD/reference/genes.gtf"
```
<script id="asciicast-309150" src="https://asciinema.org/a/309150.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>


### 1.4. Creating the Kraken2 database:
To create the Kraken2 database, run the following command:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.KrakenDB --db $PWD/K2DB
```
<script id="asciicast-309125" src="https://asciinema.org/a/309125.js" async data-autoplay="false" data-size="small" data-cols="160" data-rows="43" data-speed="3" data-loop="0"></script>

We are now ready to execute the workflow!

---

## 2. Executing the Main `nf-rnaSeqMetagen` Pipeline
As seen on the `help menu` above, there are a couple of options that you can use with this workflow. It can become a bit tedious and confusing having to specify these commands everytime you have to execute the each section for the analysis. To make your life easier, we will create a configuration script that we will use in this tutorial (we will pass this using the `-c` option of `nextflow`). You can name it whatever you want, but for now, lets call it `myparams.config`. We will add the mandatory arguements for now, but as you become more farmiliar with the workflow - you can experiment with other options. You can use your favourite text editor to create the `myparams.config` file. Copy and paste the the parameters below:
```
params {
    data   = $PWD/data
    out    = $PWD/myresults
    genome = $PWD/reference/genome.fa
    genes  = $PWD/reference/gene.gtf
    db     = $PWD/K2DB
}

```

To perform filtering of host reads and classification of exogeneous reads, use this command:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode run.FilterClassify -c myparams.config
```

---

## 3. Exploring `nf-rnaSeqMetagen` Results

```
- [1] Sample analysis directories  =>    `<output_directory>/<sample_1> .. <sample_N>`
- [2] MultiQC                      =>    `<output_directory>/MultiQC`
- [3] Upset tool                   =>    `<output_directory>/upset`
- [4] Workflow tracing             =>    `<output_directory>/workflow-tracing
```
### 3.1. MultiQC
View the full MultiQC report <a href="examples/output/MultiQC/multiqc_report.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/MultiQC/multiqc_report.html"></iframe>

### 3.2. Sample analysis directories
#### 3.2.1. Krona report: raw reads (SRR5074528)
View full Krona chart for raw reads <a href="examples/output/SRR5074528/SRR5074528_reads.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="no" src="examples/output/SRR5074528/SRR5074528_reads.html"></iframe>

#### 3.2.2 Krona report: assembled reads (SRR5074528)
View full Krona chart for assembled reads <a href="examples/output/SRR5074528/SRR5074528_fasta.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="no" src="examples/output/SRR5074528/SRR5074528_fasta.html"></iframe>

### 3.3. UpSet visualisation tool
View full UpSet plot <a href="examples/output/upset/index.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/upset/index.html"></iframe>

### 3.4. Workflow tracing
#### 3.4.1. Report
View full Nextflow report <a href="examples/output/workflow-tracing/nf-rnaSeqMetagen_report.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/workflow-tracing/nf-rnaSeqMetagen_report.html"></iframe>

#### 3.4.2. Timeline
View full timeline report <a href="examples/output/workflow-tracing/nf-rnaSeqMetagen_timeline.html" target="_blank">here</a>.
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/workflow-tracing/nf-rnaSeqMetagen_timeline.html"></iframe>

---
