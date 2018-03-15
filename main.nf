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
read_pair = Channel.fromFilePairs("${data_path}/*_R{1,2}.fq", type: 'file')

// 1. 
process runSTAR_process {
    cpus 6
    memory '40 GB'
    time '100h'
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
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, file(reads) from unmapped_kraken
    
    output:
    set sample, file("kraken_output.txt") into kraken_classified_reads 
    
    """	
    kraken --db ${db} \
        --fastq-input \
        --paired ${reads.get(0)} ${reads.get(1)} \
        --threads 5 \
        --output kraken_output.txt
    """ 
}

// 3. 
process runTrinityAssemble_process {
     cpus 6
     memory '150 GB'
     time '50h'
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
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

    input:
    set sample, file(fasta) from trinity_assembled_reads

    output:
    set sample, file("kraken_outputFasta.txt") into kraken_classified_fasta 

    """	
    kraken --db ${db} \
        --fasta-input ${fasta} \
        --threads 5 \
        --output kraken_outputFasta.txt
    """ 
}

all_classified = kraken_classified_reads.merge( kraken_classified_fasta ) { listA, listB -> [ listA[0], [listA[1], listB[1]] ] }

// 5. 
process runKronareport{
    cpus 8
    memory '200 GB'
    time '50h'
    tag { sample }
    publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
    input:
    set sample, file(kraken) from all_classified
    
    output:
    set sample, file("*") into html

    """
    function doKrona {
        cut -f 2,3 \$1 > \$(sed 's/.kraken/.krona/' <<< "\$1")
        ktImportTaxonomy \$(sed 's/.kraken/.krona/' <<< "\$1") \
        -tax ${taxonomy} \
        -o \$(sed 's/.kraken/.html/' <<< "\$1")
    }
    doKrona ${kraken.get(0)}
    doKrona ${kraken.get(1)}
    """
}


html.subscribe {println it}

