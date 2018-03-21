#!/usr/bin/env nextflow

// PIPELINE PARAMETERS - Edit if brave... Else, specify options on command line
params.data      = "/spaces/phelelani/ssc_data/data_trimmed/inflated"                                                   // Path to where the input data is located (where fastq files are located).
params.out       = "/spaces/phelelani/ssc_data/nf-rnaSeqMetagen"                                                          // Path to where the output should be directed.
params.db        = "/global/blast/kraken_db/kraken_std"
params.taxonomy  = "/global/blast/taxonomy"
params.genome    = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa"    // The whole genome sequence
params.index     = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex"                     // Path to where the STAR index files are locaded
params.bind      = '/global/;/spaces/'                                                                                  // Paths to be passed onto the singularity image

// DO NOT EDIT FROM HERE
data_path        = file(params.data, type: 'dir')                                                                       // Path to where the input data is located (where fastq files are located). 
out_path         = file(params.out, type: 'dir')                                                                        // Path to where the output should be directed.
db               = file(params.db, type: 'dir')
taxonomy         = file(params.taxonomy, type: 'dir')
genome           = file(params.genome, type: 'file')                                                                    // The whole genome sequence
index            = file(params.index, type: 'dir')                                                                      // Path to where the STAR index files are locaded 
bind             = params.bind.split(';')          
//======================================================================================================
//
//
//
//======================================================================================================
// HELP MENU
if (params.help) {
    log.info ''
    log.info "===================================="
    log.info "         nf-rnaSeqMetagen v0.1        "
    log.info "===================================="
    log.info ''
    log.info 'USAGE: '
    log.info 'nextflow run main.nf --data /path/to/data --out /path/to/output --db /path/to/kraken-db --taxonomy /path/to/taxonomy --genome /path/to/genome.fa --index /path/to/STARIndex --bind /path/to/bind1;/path/to/bind2'
    log.info ''
    log.info 'HELP: '
    log.info 'nextflow run main.nf --help'
    log.info ''
    log.info 'MANDATORY ARGUEMENTS:'
    log.info '    --data     FOLDER    Path to where the input data is located (fastq | fq)'
    log.info '    --out      FOLDER    Path to where the output should be directed (will be created if it does not exist).'
    log.info '    --db       FOLDER    
    log.info '    --taxonomy FOLDER   
    log.info '    --genome   FILE      The whole genome sequence (fasta | fa | fna)'
    log.info '    --index    FOLDER    Path to where the STAR index files are locaded'
    log.info '    --bind     FOLDER(S) Paths to be passed onto the singularity image'
    log.info ''
    log.info "====================================\n"
    exit 1
}

// RUN INFO
log.info "===================================="
log.info "           nf-rnaSeqCount           "
log.info "===================================="
log.info "Input data          : ${data_path}"
log.info "Output path         : ${out_path}"
log.info "Kraken database     : ${db}"
log.info "Taxonomy database   : ${taxonomy}"
log.info "Genome              : ${genome}"
log.info "Genome Index (STAR) : ${index}"
log.info "Paths to bind       : ${bind}"
log.info "====================================\n"
//======================================================================================================
//
//
//
//======================================================================================================
// PIPELINE START
//Create output directory
out_path.mkdir()

// Get input reads
read_pair = Channel.fromFilePairs("${data_path}/*R[1,2].fastq", type: 'file') 
                   .ifEmpty { error "ERROR - Data input: \nOooops... Cannot find any '.fastq' or '.fq' files in ${data_path}. Please specify a folder with '.fastq' or '.fq' files."}

// 1. 
process runSTAR_process {
    cpus 6
    memory '40 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, file(reads) from read_pair
    
    output:
    set sample, file("${sample}_*") into star_results
    set sample, file("${sample}_Unmapped*") into unmapped_kraken, unmapped_trinity
    
    """	
    STAR --runMode alignReads \
       --genomeDir ${index} \
       --readFilesIn ${reads.get(0)} ${reads.get(1)} \
       --runThreadN 5 \
       --outSAMtype BAM SortedByCoordinate \
       --outReadsUnmapped Fastx \
       --outFileNamePrefix ${sample}_
       
    sed -i 's|\\s.[0-9]\$|\\/1|g' ${sample}_Unmapped.out.mate1 
    sed -i 's|\\s.[0-9]\$|\\/2|g' ${sample}_Unmapped.out.mate2

    """ 
}

// 2. 
process runKrakenClassifyReads_process {
    cpus 6
    memory '150 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, file(reads) from unmapped_kraken
    
    output:
    set sample, file("kraken_report_reads.krak") into kraken_classified_reads 
    
    """	
    kraken --db ${db} \
        --fastq-input \
        --paired ${reads.get(0)} ${reads.get(1)} \
        --threads 5 \
        --output kraken_report_reads.krak
    """ 
}

// 3. 
process runTrinityAssemble_process {
     cpus 6
     memory '150 GB'
     time '50h'
    scratch '$HOME/tmp'
     tag { sample }
     publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
     input:
     set sample, file(reads) from unmapped_trinity
    
     output:
     set sample, "trinity_${sample}/Trinity.fasta" into trinity_assembled_reads

     """
     Trinity --seqType fq \
        --max_memory 150G \
        --left ${reads.get(0)} --right ${reads.get(1)} \
        --SS_lib_type RF \
        --CPU 5 \
        --output  trinity_${sample}
     """
 }

// 4.
process runKrakenClassifyFasta_process{
    cpus 6
    memory '150 GB'
    time '100h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(fasta) from trinity_assembled_reads

    output:
    set sample, file("kraken_report_fasta.krak") into kraken_classified_fasta 

    """	
    kraken --db ${db} \
        --fasta-input ${fasta} \
        --threads 5 \
        --output kraken_report_fasta.krak
    """ 
}

all_classified = kraken_classified_reads.merge( kraken_classified_fasta ) { listA, listB -> [ listA[0], [listA[1], listB[1]] ] }

// 5. 
process runKronareport{
    cpus 8
    memory '5 GB'
    time '10h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, file(kraken) from all_classified
    
    output:
    set sample, file("*") into html

    """
    function createChart {
    cut -f 2,3 \$1 > \$(sed 's/.krak/.kron/' <<< "\$1")
    ktImportTaxonomy \$(sed 's/.krak/.kron/' <<< "\$1") \
        -tax ${taxonomy} \
        -o \$(sed 's/.krak/.html/' <<< "\$1")
    }
    createChart ${kraken.get(0)}
    createChart ${kraken.get(1)}
    """
}

// 6a. Collect files for STAR QC
star_results.collectFile () { item -> [ 'qc_star.txt', "${item.get(1).find { it =~ 'Log.final.out' } }" + ' ' ] }
.set { qc_star }

// 6. Get QC for STAR, HTSeqCounts and featureCounts
process runMultiQC_process {
    cpus 1
    memory '5 GB'
    time '10h'
    scratch '$HOME/tmp'
    tag { sample }
    publishDir "$out_path/report_QC", mode: 'copy', overwrite: false

    input:
    file(star) from qc_star

    output:
    file('*') into multiQC

    """
    multiqc `< ${star}` --force
    """
}
//======================================================================================================
//
//
//
//======================================================================================================
// WORKFLOW SUMMARY
workflow.onComplete {
    println "===================================="
    println "Pipeline execution summary:"
    println "===================================="
    println "Execution command   : ${workflow.commandLine}"
    println "Execution name      : ${workflow.runName}"
    println "Workflow start      : ${workflow.start}"
    println "Workflow end        : ${workflow.complete}"
    println "Workflow duration   : ${workflow.duration}"
    println "Workflow completed? : ${workflow.success}"
    println "Work directory      : ${workflow.workDir}"
    println "Project directory   : ${workflow.projectDir}"
    println "Execution directory : ${workflow.launchDir}"
    println "Configuration files : ${workflow.configFiles}"
    println "Workflow containers : ${workflow.container}"
    println "exit status : ${workflow.exitStatus}"
    println "Error report: ${workflow.errorReport ?: '-'}"
    println "===================================="
}

workflow.onError {
    println "Oohhh DANG IT!!... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
//======================================================================================================
