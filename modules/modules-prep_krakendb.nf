#!/usr/bin/env nextflow
nextflow.enable.dsl=2

db       = file(params.db, type: 'dir')

// https://benlangmead.github.io/aws-indexes/k2
process run_DownloadK2DBIndexes {
    publishDir "${db}", mode: 'copy', overwrite: true
    label 'mini'
    tag { "Download_K2DB" }
    memory '1 GB'

    output:
    path("*"), emit: k2db
    
    """
    lftp -e 'pget -n10 https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230314.tar.gz; exit'
    tar -vxf k2_standard_20230314.tar.gz
    kraken2-build --download-taxonomy --db .
    /opt/KronaTools-2.8/updateTaxonomy.sh --only-build --preserve ./taxonomy/
    """
}
