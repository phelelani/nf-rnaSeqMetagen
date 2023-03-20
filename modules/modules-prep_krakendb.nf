#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_GenerateKrakenDB {
    label 'maxi'
    tag { "Generate Kraken DB" }
    publishDir "${db}", mode: 'copy', overwrite: false
    
    output:
    path("*.k2d"), emit:kraken_db
    path("taxonomy/taxdump.tar.gz"), emit taxonomy_dump
    
    """
    kraken2-build --standard --threads ${task.cpus} --db .
    """
}

process run_UpdateTaxonomy {
    label 'mini'
    tag { "Update NCBI Taxonomy" }
    publishDir "${taxonomy}", mode: 'copy', overwrite: true
    
    input:
    path(dmp)
    
    output:
    path("*"), emit: taxonomy_update
    
    """
    /opt/KronaTools-2.8/updateTaxonomy.sh --only-build --preserve .
    """
}

