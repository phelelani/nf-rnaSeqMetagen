#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_GenerateSTARIndex {
    label 'maxi'
    tag { "Generate Star Index" }
    publishDir "${index_dir}", mode: 'copy', overwrite: true

    output:
    set val("starIndex"), path("*") into star_index
    
    """
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir . \
        --genomeFastaFiles ${genome} \
        --sjdbGTFfile ${genes} \
        --sjdbOverhang 99
    """
}

process run_GenerateBowtieIndex {
    label 'maxi'
    tag { "Generate Bowtie2 Index" }
    publishDir "$index_dir", mode: 'copy', overwrite: true
    
    output:
    set val("bowtieIndex"), path"*") into bowtie_index
    
    """
    bowtie2-build --threads ${task.cpus} ${genome} genome
    """
}        
