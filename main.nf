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
    // set sample, file("${sample}_*.bam") into bam_file
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

// // 2. 
// process runKrakenClassifyReads_process {
//     cpus 6
//     memory '150 GB'
//     time '100h'
//     tag { sample }
//     publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
//     input:
//     set sample, file(reads) from unmapped_kraken
    
//     output:
//     set sample, file("kraken_output.txt") into kraken_classified_reads 
    
//     """	
//     kraken --db ${db} \
//         --fastq-input \
//         --paired ${reads.get(0)} ${reads.get(1)} \
//         --threads 5 \
//         --output kraken_output.txt
//     """ 
// }

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

// // 4.
// process runKrakenClassifyFasta_process{
//     cpus 6
//     memory '150 GB'
//     time '100h'
//     tag { sample }
//     publishDir "$out_path/${sample}", mode: 'copy', overwrite: false

//     input:
//     set sample, file(fasta) from trinity_assembled_reads

//     output:
//     set sample, file("kraken_outputFasta.txt") into kraken_classified_fasta 

//     """	
//     kraken --db ${db} \
//         --fasta-input ${fasta} \
//         --threads 5 \
//         --output kraken_outputFasta.txt
//     """ 
// }


// star_results.subscribe { println it }
// kraken_classified_reads.subscribe { println it }
// kraken_classified_fasta.subscribe { println it }

// all_classified = classified_reads.merge( classified_fasta ) { listA, listB -> [ listA[0], [listA[1], listB[1]] ] }
// all_classified.println { it }

// // 5. 
// process runKronareport{
//     cpus 8
//     memory '200 GB'
//     time '50h'
//     tag { sample }
//     publishDir "$out_path/${sample}", mode: 'copy', overwrite: false
    
//     input:
//     set sample, file(kraken) from classified_reads
    
//     output:
//     set sample, file("results.html") into html

//     """
    
//     cut -f 2,3 ${kraken} > columns.txt
//     ktImportTaxonomy columns.txt                \
//        --tax ${taxonomy}                        \
//        -o results.html

//        //cut -f 2,3 caskiSubset_01/kraken_output.txt caskiSubset_02/kraken_output.txt caskiSubset_03/kraken_output.txt > data.txt
//        //cut -f 2,3 caskiSubset_01/kraken_outputFasta.txt caskiSubset_02/kraken_outputFasta.txt caskiSubset_03/kraken_outputFasta.txt > dataFasta.txt

//        //odds  = Channel.from([1, 3, 5, 7, 9]);
//        //evens = Channel.from([2, 4, 6]);

// //	 odds
// //	  .merge( evens ) { o, e -> [o, e] }
// //   	  .subscribe { println it }
//     """
// }


// html.subscribe {println it}

