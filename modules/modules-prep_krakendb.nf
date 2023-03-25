#!/usr/bin/env nextflow
nextflow.enable.dsl=2

db       = file(params.db, type: 'dir')
taxonomy = file("${db}/taxonomy", type: 'dir')

process run_DownloadK2DBLibs {
    label 'mini'
    tag { "Download ${k2db_libs}" }
    memory '1 GB'
    publishDir "${db}", mode: 'copy', overwrite: true
    
    input:
    each k2db_libs

    output:
    path("${k2db_libs}_download.log"), emit: download_log
    path("*"), emit: genomic_libraries
    
    """
    kraken2-build --download-library ${k2db_libs} --db . --use-ftp |& tee ${k2db_libs}_download.log
    """
}

process run_BuildK2DB {
    label 'maxi'
    tag { "Build K2DB" }
    publishDir "${db}", mode: 'copy', overwrite: true
    
    input:
    path(logs)
    path(genomic_libraries)

    output:
    path("${k2db_libs}_log"), emit: build_log
    path("taxonomy/taxdump.tar.gz"), emit: taxonomy_dump
    
    """
    kraken2-build --build --db . |& tee k2db_build.log
    """
}

process run_UpdateTaxonomy {
    label 'mini'
    tag { "Update NCBI Taxonomy" }
    publishDir "${taxonomy}", mode: 'copy', overwrite: true
    
    input:
    path(build_log)
    
    output:
    path("*"), emit: taxonomy_update
    
    """
    /opt/KronaTools-2.8/updateTaxonomy.sh --only-build --preserve .
    """
}
