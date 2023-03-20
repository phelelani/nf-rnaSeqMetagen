#!/usr/bin/env nextflow
nextflow.enable.dsl=2

genome     = file(params.genome, type: 'file', checkIfExists: true)
genes      = file(params.genes, type: 'file', checkIfExists: true)
index_dir  = file(params.genome, type: 'file', checkIfExists: true).getParent()
out_dir    = file(params.out, type: 'dir')
out_dir.mkdir()

// 1.  ALIGN READS TO REFERENCE GENOME
process run_STAR {
    label 'maxi'
    tag { sample }
    publishDir "${out_dir}/${sample}", mode: 'copy', overwrite: true, pattern: "${sample}*.{out,tab}"
    
    input:
    tuple val(sample), path(reads)
    
    output:
    tuple val(sample), path("${sample}*.{out,tab}"), emit: star_results
    tuple val(sample), path("${sample}_Unmapped*"), emit: unmapped_reads

    """
    STAR --runMode alignReads \
        --genomeDir ${index_dir} \
        --readFilesCommand gunzip -c \
        --readFilesIn ${reads.findAll().join(' ')} \
        --runThreadN ${task.cpus} \
        --outSAMtype BAM Unsorted \
        --outReadsUnmapped Fastx \
        --outFileNamePrefix ${sample}_
    """
}

process run_FixSeqNames {
    label 'mini'
    tag { sample }
    publishDir "${out_dir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(unmapped)
            
    output:
    tuple val(sample), path("${sample}_unmapped*"), emit: unmapped_reads

    """
    sed 's| \\(.*\\)\$|\\/1|g' ${unmapped.find { it =~ 'mate1' } } > ${sample}_unmapped_R1.fastq
    sed 's| \\(.*\\)\$|\\/2|g' ${unmapped.find { it =~ 'mate2' } } > ${sample}_unmapped_R2.fastq
    """
}

// 2. Run KRAKEN to classify the raw reads that aren't mapped to the reference genome.
process run_KrakenClassifyReads {
    label 'maxi'
    tag { sample }
    publishDir "${out_dir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("${sample}_reads.krak"), emit: kraken_reads_report
    tuple val(sample), path("${sample}_classified_*.fastq"), emit: kraken_classified_reads
    tuple val(sample), path("${sample}_unclassified_*.fastq"), emit: kraken_unclassified_reads

    """	
    kraken2 --db ${db} \
        --paired ${reads.findAll().join(' ')} \
        --threads ${task.cpus} \
        --classified-out ${sample}_classified#.fastq \
        --unclassified-out ${sample}_unclassified#.fastq \
        --output ${sample}_reads.krak
    """ 
}

// 3. Assemble the reads, emit: longer contigs/sequences for classification.
process run_TrinityAssemble {
    label 'maxi'
    tag { sample }
    publishDir "${out_dir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("trinity_${sample}/Trinity.fasta"), emit: trinity_assembled_reads

    """
    Trinity --seqType fq \
        --max_memory 150G \
        --left ${reads.find { it =~ 'R1' } } \
        --right ${reads.find { it =~ 'R2' } } \
        --SS_lib_type RF \
        --CPU ${task.cpus} \
        --output trinity_${sample}
    """
}

// 4. Run KRAKEN to classify the assembled FASTA sequences.
process run_KrakenClassifyFasta {
    label 'maxi'
    tag { sample }
    publishDir "${out_dir}/${sample}", mode: 'copy', overwrite: true, pattern: "*{_fasta.krak,_classified.fasta}"

    input:
    tuple val(sample), path(fasta)

    output:
    tuple val(sample), path("${sample}_fasta.krak"), emit: kraken_fasta_report
    tuple val(sample), path("${sample}_classified.fasta"), emit: kraken_classified_fasta
    tuple val(sample), path("${sample}_unclassified.fasta"), emit: kraken_unclassified_fasta

    """	
    kraken2 --db ${db} \
        ${fasta} \
        --threads ${task.cpus} \
        --classified-out ${sample}_classified.fasta \
        --unclassified-out ${sample}_unclassified.fasta \
        --output ${sample}_fasta.krak
    """ 
}

// // For each val(sample), create a list with [ SAMPLE_NAME, READ, FASTA ] by merging the classified outputs (reads and fasta) from KRAKEN
// kraken_reads_report.join(kraken_fasta_report)
//     .map { it -> [ it[0], [ it[1], it[2] ] ] }
//     .set { all_kraken_reports }

//5. Create that pretty KRONA report for all samples (reads and fasta)
process run_KronaReport {
    label 'mini'
    tag { sample }
    publishDir "${out_dir}/${sample}", mode: 'copy', overwrite: true

    input:
    tuple val(sample), path(kraken)

    output:
    tuple val(sample), path("*.html"), emit: html
    tuple val(sample), path("*reads.kron"), emit: reads_krona
    tuple val(sample), path("*fasta.kron"), emit: fasta_krona

    """
    function createChart {
        awk '\$1=="C" { print \$2"\t"\$3 }' \$1 > \$(sed 's/.krak/.kron/' <<< "\$1")
        ktImportTaxonomy \$(sed 's/.krak/.kron/' <<< "\$1") \
            -tax \$2 \
            -o \$(sed 's/.krak/.html/' <<< "\$1")
        }
    createChart ${kraken.get(0)} ${taxonomy}
    createChart ${kraken.get(1)} ${taxonomy}
    """
}

// fasta_krona_seqs.join(kraken_classified_fasta)
//     .map { it -> [ it[0], [ it[1], it[2] ] ] }
//     .set { krona_fasta_pair }

process run_CollectTaxSeqs {
    label 'mini'
    tag { sample }
    publishDir "${out_dir}/${sample}/taxon_sequences", mode: 'copy', overwrite: true
    
    input: 
    tuple val(sample), path(test)
            
    output:
    tuple val(sample), path("taxid_*.fasta"), emit: taxon_sequences
            
    """
    for id in \$(awk '{ print \$2 }' ${test.get(0)} | sort -gu)
    do
        grep -A1 --no-group-separator "kraken:taxid|\$id\$" ${test.get(1)} | sed 's/path=\\[\\(.*\\)\\] //' > "taxid_"\$id".fasta"
    done
    """
}
        
// // 6a. Collect files for STAR QC
// star_results.collectFile () { item -> [ 'qc_star.txt', "${item.get(1).find { it =~ 'Log.final.out' } }" + ' ' ] }
//     .set { qc_star }

// 6. Get QC for STAR, HTSeqCounts and featureCounts
process run_MultiQC {
    label 'mini'
    tag { "Get QC Information" }
    publishDir "${out_dir}/MultiQC", mode: 'copy', overwrite: true
    
    input:
    path(star)
    
    output:
    path('*'), emit: multiQC
    
    """
    multiqc `< ${star}` --force
    """
}

// // 7a. Collect all the krona file locations and put them in a text file
// fasta_krona.collectFile () { item -> [ 'fasta_krona_files.txt', "${item.get(1)}" + '\n' ] }
//     .set { fasta_krona_list } 

// fasta_krona_report.collectFile () { item -> [ 'fasta_krona_report.txt', "${item.get(1).find { it =~ 'fasta.kron' } }" + '\n' ] }

process run_CopyUpsetDir {
    label 'mini'
    tag { "Copy UpSet Tool" }
    publishDir "${out_dir}/upset", mode: 'copy', overwrite: true
    
    output:
    path("*"), emit: upset_dir
    
    """
    /bin/hostname
    rsync -avhP /opt/upset/* .
    """
}
        
// 7b. Prepare data for creating the matrix for UpSet: json file, taxonomy file, and sample taxids 
process run_PrepareMatrixData {
    label 'mini'
    tag { "Prepare Matrix Data" }

    input:
    path(list)
    
    output:
    path("*.{dmp,taxon,json}"), emit: matrix_files

    """
    get_taxons.sh ${taxonomy} ${list}
    """
}

// 7 Create the UpSet matrix
process run_CreateMatrix {
    label 'mini'
    tag { "Create UpSet Matrix" }
    publishDir "${out_dir}/upset/data/nf-rnaSeqMetagen", mode: 'copy', overwrite: true
    
    input:
    path(list)

    output:
    tuple val("upset_files"), path(list), path("*"), emit: the_matrix

    """
    create_matrix.R
    """
}
