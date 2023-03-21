#!/usr/bin/env nextflow
nextflow.enable.dsl=2

db       = file(params.db, type: 'dir')
taxonomy = file("${db}/taxonomy", type: 'dir')

process run_DownloadK2DBLibs {
    label 'mini'
    tag { "Download ${k2db_libs}" }

    input:
    each k2db_libs

    output:
    path("${k2db_libs}_download.log"), emmit: download_log
    
    """
    kraken2-build --download-library ${k2db_libs} --db ${db} --use-ftp
    cp .command.log > ${k2db_libs}_download.log
    """
}

process run_BuildK2DB {
    label 'maxi'
    tag { "Build K2DB" }
    publishDir "${db}/logs", mode: 'copy', overwrite: true
    
    input:
    path(logs)

    output:
    path("${k2db_libs}_log"), emmit: build_log
    
    """
    kraken2-build --build --db ${db}
    cp .command.log > k2db_build.log
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

