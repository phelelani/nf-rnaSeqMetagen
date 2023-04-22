#!/usr/bin/env nextflow
nextflow.enable.dsl=2

db       = file(params.db, type: 'dir')

// https://benlangmead.github.io/aws-indexes/k2
process run_DownloadK2DBIndexes {
    label 'mini'
    tag { "Download_K2DB" }
    memory '1 GB'
    publishDir "${db}", mode: 'copy', overwrite: true
    
    output:
    path("*"), emit: k2db
    
    """
    lftp -e 'pget -n10 https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz; exit'
    tar -vxf k2_standard_20230314.tar.gz
    rm -rf k2_standard_20230314.tar.gz
    """
}

process run_DownloadTaxonomy {
    label 'mini'
    tag { "Download_Taxonomy" }
    memory '1 GB'
    publishDir "${db}", mode: 'copy', overwrite: true
    
    output:
    path("taxonomy"), emit: taxonomy
    path("taxonomy/taxdump.tar.gz"), emit: taxonomy_dump
    
    """
    kraken2-build --download-taxonomy --db .
    """
}

process run_UpdateTaxonomy {
    label 'mini'
    tag { "Update_Taxonomy" }
    memory '1 GB'
    publishDir "${db}/taxonomy", mode: 'copy', overwrite: true

    input:
    path(taxonomy_dump)
    
    output:
    path("*"), emit: k2db
    
    """
    /opt/KronaTools-2.8/updateTaxonomy.sh --only-build --preserve .
    """
}
