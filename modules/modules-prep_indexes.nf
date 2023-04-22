#!/usr/bin/env nextflow
nextflow.enable.dsl=2

genome     = file(params.genome, type: 'file', checkIfExists: true)
genes      = file(params.genes, type: 'file', checkIfExists: true)
index_dir  = file(params.genome, type: 'file', checkIfExists: true).getParent()

process run_GenerateSTARIndex {
    label 'maxi'
    tag { "Generate Star Index" }
    publishDir "${index_dir}", mode: 'copy', overwrite: true

    output:
    tuple val("starIndex"), path("*"), emit: star_index
    
    """
    STAR --runThreadN ${task.cpus} \
        --runMode genomeGenerate \
        --genomeDir . \
        --genomeFastaFiles ${genome} \
        --sjdbGTFfile ${genes} \
        --sjdbOverhang 99
    """
}
