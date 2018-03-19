#!/usr/bin/env nextflow

params.data = "/global/blast/ssc_data"
params.taxonomy = "/global/blast/taxonomy"
params.out = "/spaces/phelelani/ssc_data/nf-rnaSeqMetagen"
params.db = "/global/blast/kraken_db/kraken_std"
params.ref = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
params.index = "/global/blast/reference_genomes/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex" 

out_path = file(params.out)
taxonomy = params.taxonomy
data_path = params.data
ref = params.ref
index = params.index
db = params.db

out_path.mkdir()

// Create a new channel from input data/reads
read_pair = Channel.fromFilePairs("${data_path}/*_R{1,2}.fastq", type: 'file')

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
    /bin/hostname
    multiqc `< ${star}` --force
    """
}
