[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/728)
    
`nf-rnaSeqMetagen` is a [Nextflow](http://nextflow.io/)
        
To use the `nf-rnaSeqMetagen` pipeline, the following dependencies are required:
1. Installed softwares:
   - [`Nextflow`](https://www.nextflow.io/)
   - [`Singularity`](http://singularity.lbl.gov/)
2. Reference genome, annotation and indexes
   - Reference genome (`.fa`/`.fasta`) and genome annotation (`.gtf`) files.
- Reference genome indexes (`STAR` - see *1.3.* below on how to generate the indexes).
 
---

<p align="center">
  <img width="600px" src="assets/images/nf-rnaSeqMetagen.png">
</p>

---

## 1. Obtaining the `nf-rnaSeqMetagen` pipeline and preparing data
First, you need to clone the `nf-rnaSeqMetagen` repository onto you machine. You can either use `git` or `nextflow`. I recommend you use `nextflow`. The rest of this documentation assumes that you have used `nextflow` to clone this workflow.
```bash
nextflow pull https://github.com/phelelani/nf-rnaSeqMetagen
```
<script id="asciicast-308777" src="https://asciinema.org/a/308777.js" async data-autoplay="false" data-size="small" data-cols="150" data-rows="6" data-speed="1.5" data-loop="0"></script>

Content of the repository (located in `$HOME/.nextflow/assets/phelelani/nf-rnaSeqCount`):
<script id="asciicast-308808" src="https://asciinema.org/a/308808.js" async data-autoplay="false" data-size="small" data-cols="150" data-rows="12" data-speed="1.5" data-loop="0"></script>

To get the `help menu` for the workflow, execute the following command from anywherre on your system:
```
nextflow run nf-rnaSeqMetagen --help
```
<script id="asciicast-308804" src="https://asciinema.org/a/308804.js" async data-autoplay="false" data-size="small" data-cols="150" data-rows="43" data-speed="1.5" data-loop="0"></script>

---

### 1.1. Download test datasets (optional)
We will now download the reference genome (along with its annotation file) from Ensembl. We will also download the FASTQ files from the H3ABioNet site, which we will analyse using the `nf-rnaSeqMetagen` workflow. *__NB__: Skip this section if you have your own data to analyse using this workflow! This section is only for getting data to practice using the `nf-rnaSeqMetagen` workflow!* 

Make directories:
```
mkdir example
cd example
mkdir reference
```
<script id="asciicast-308945" src="https://asciinema.org/a/308945.js" async data-autoplay="false" data-size="small" data-cols="150" data-rows="12" data-speed="1.5" data-loop="0"></script>


Download and decompress the mouse reference genome along with its annotation:
```
wget -c -O reference/genome.fa.gz ftp://ftp.ensembl.org/pub/release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz
wget -c -O reference/genes.gtf.gz ftp://ftp.ensembl.org/pub/release-68/gtf/mus_musculus/Mus_musculus.GRCm38.68.gtf.gz
gunzip reference/genome.fa.gz
gunzip reference/genes.gtf.gz
```

Download RNA-seq test dataset from H3ABioNet:
```
mkdir data
for sample in sample{37..42}_R{1,2}.fastq.gz; do wget -c -O data/$sample http://h3data.cbio.uct.ac.za/assessments/RNASeq/practice/dataset/$sample; done
```

### 1.2. Download the `Singularity` containers (required to execute the pipeline):
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.Containers
```
<script id="asciicast-308816" src="https://asciinema.org/a/308816.js" async data-autoplay="false" data-size="small" data-cols="150" data-rows="43" data-speed="1.5" data-loop="0"></script>

### 1.3. Generating genome indexes.
To generate the `STAR` genome indexes, run the following commands:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.STARIndex --genome "$PWD/reference/genome.fa" --genes "$PWD/reference/genes.gtf"
```

### 1.4. Creating the Kraken2 database:
To create the Kraken2 database, run the following command:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode prep.KrakenDB --db $PWD/K2DB
```

We are now ready to execute the workflow!

---

## 2. Executing the main `nf-rnaSeqMetagen` pipeline
As seen on the `help menu` above, there are a couple of options that you can use with this workflow. It can become a bit tedious and confusing having to specify these commands everytime you have to execute the each section for the analysis. To make your life easier, we will create a configuration script that we will use in this tutorial (we will pass this using the `-c` option of `nextflow`). You can name it whatever you want, but for now, lets call it `myparams.config`. We will add the mandatory arguements for now, but as you become more farmiliar with the workflow - you can experiment with other options. You can use your favourite text editor to create the `myparams.config` file. Copy and paste the the parameters below:
```
params {
    data    = "$PWD/data"
    db      = "$PWD/K2DB"
    genome  = "$PWD/reference/genome.fa"
    genes   = "$PWD/reference/genes.fa"
}
```
Obviously - the above `myparams.config` assumes that you have been following this tutorial. If you have your data lying around somewhere in your system, you need to put the full path to where your the `data`, `genome` and `genes` files are. Since the `--mode` will keep changing, we will add this on the command as we do the analysis. Now that we have the mandatory arguements in our `myparams.config`, lets do some analysis

### 2.1. Read Filtering and Classification:
To perform filtering of host reads and classification of exogeneous reads, use this command:
```bash
nextflow run nf-rnaSeqMetagen -profile slurm --mode run.FilterClassify -c myparams.config
```

---

## 3. Explore `nf-rnaSeqMetagen` results

```
- [1] Sample analysis directories  =>    `<output_directory>/<sample_1> .. <sample_N>`
- [2] MultiQC                      =>    `<output_directory>/MultiQC`
- [3] Upset tool                   =>    `<output_directory>/upset`
- [4] Workflow tracing             =>    `<output_directory>/workflow-tracing
```
### MultiQC
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/MultiQC/multiqc_report.html"></iframe>

### Sample Analysis Directories
#### Krona Report: Raw Reads (SRR5074528)
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="no" src="examples/output/SRR5074528/SRR5074528_reads.html"></iframe>

#### Krona Report: Assembled Reads (SRR5074528)
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="no" src="examples/output/SRR5074528/SRR5074528_fasta.html"></iframe>

### UpSet Visualisation Tool
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/upset/index.html"></iframe>

### Workflow Tracing
#### Report
<div>
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/workflow-tracing/nf-rnaSeqMetagen_report.html"></iframe>
</div>

#### Timeline
<iframe style="overflow: hidden; margin: 0px; border: 1px solid grey; display: inline-block; width: 832px; float: none; visibility: visible; height: 723px;" seamless scrolling="yes" src="examples/output/workflow-tracing/nf-rnaSeqMetagen_timeline.html"></iframe>

---
