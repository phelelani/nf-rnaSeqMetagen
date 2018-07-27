# nf-rnaSeqMetagen
*nf-rnaSeqMetagen* is a [Nextflow](http://nextflow.io/) 

<p align="center">
  <img width="1000" src="nf-rnaSeqMetagen.png">
</p>

# 1. Pipeline Dependencies
To use the rnaSeqCount pipeline, the following dependencies are required:
### 1.1. Softwares
- [x] [Nextflow](https://www.nextflow.io/)
- [x] [Singularity](http://singularity.lbl.gov/)

### 1.2. Singularity Containers
- [x] https://www.singularity-hub.org/collections/728

## 1.3. Reference Genome and Indexes
- [x] Reference Genome (.fa) and Genome Annotation (.gtf) files
- [x] Reference Genome Indexes (```bowtie2``` & ```STAR``` - see below on how to generate)
- [x] Kraken database (Installation instructions: http://ccb.jhu.edu/software/kraken/MANUAL.html#installation)

# 2. Optaining the ```nf-rnaSeqCount``` pipeline
The ```nf-rnaSeqMetagen``` pipeline can be obtain using any of the following methods:

### 2.1. Using the ```git``` command:
- [x] ```git clone https://github.com/phelelani/nf-rnaSeqMetagen.git```

### 2.2. Using the ```nextflow``` command:
- [x] ```nextflow pull phelelani/nf-rnaSeqMetagen```
- [x] ```nextflow pull https://github.com/phelelani/nf-rnaSeqMetagen.git```
- [x] ```nextflow clone phelelani/nf-rnaSeqMetagen <target-dir>```

# 3. Generating genome indexes.
To generate the ```STAR``` and ```bowtie2``` indexes for the reference genome, the ```Singularity``` containers first need to be downloaded from [Singularity Hub](ttps://www.singularity-hub.org). The ```prepareData.nf``` script can be used to download and prepare data (generate indexes) to be be used with the ```nf-rnaSeqCount``` pipeline. The ```prepareData.nf``` can be run in three different modes:
- [x] ```getContainers```: for downloading the required ```Singularity``` containers.
- [x] ```generateStarIndex```: for generating ```STAR``` indexes.
- [x] ```generateBowtieIndex```: for generating ```bowtie2``` indexes.
- [x] ```generateKrakenDB``` : for downloading and generating ```kraken``` database.

To generate the genome indexes, run the following commands in the pipeline directory:

### 3.1 Download ```Singularity``` containers:
```
nextflow run prepareData.nf --mode getContainers -profile pbsPrepare
```

### 3.2. Generate ```STAR``` index
```
nextflow run prepareData.nf --mode generateStarIndex -profile pbsPrepare
```

### 3.3. Generate ```bowtie2``` index
```
nextflow run prepareData.nf --mode generateBowtieIndex -profile pbsPrepare
```

NB: The ```-profile``` option can either be one of depending on the scheduler.

# 4. Pipeline Execution
The ```nf-rnaSeqCount``` pipeline can be run in one of two ways:

### 4.1 By editing the ```parameters.config``` file and specifying the parameters (recommended)
Edit ```main.config```:
```
/*
 *  USE THIS FILE TO SPECIFY YOUR PARAMETERS. ALLOWED PARAMETERS ARE AS FOLLOWS:
 *  ============================================================================
 *  data     : Path to where the input data is located (where fastq files are located).
 *  out      : Path to where the output should be directed.
 *  db       :
 *  taxonomy :
 *  genome   : The whole genome sequence.
 *  index    : Path to where the STAR index files are locaded.
 *  bind     : Paths to be passed onto the singularity image (Semi-colon separated).
 *  help     : Print out help menu. Passed as "--help" to the "main.nf" script for detailed information
 */
params {
    data     = "/path/to/data"
    out      = "/path/to/output"
    filetype = "fastq.gz"
    db       = "/path/to/kraken_db/"
    taxonomy = "/path/to/taxonomy"
    genome   = "/path/to/genome.fa"
    index    = "/path/to/STARIndex"
    bind     =  "/path/to/bind_1;/path/to/bind_2"
    help     = null
}

```

To run the pipeline:
```
nextflow run main.nf
```

### 4.1. Directly from the command line by supplying the required parameters
```
nextflow run main.nf --data '/path/to/data' \
    --out '/path/to/output' \
    --genome '/path/to/genome.fa' \
    --index '/path/to/STARIndex' \
    --genes '/path/to/genes.gtf' \
    --bind '/path/to/bind;/another/path/to/bind'
```


# References
