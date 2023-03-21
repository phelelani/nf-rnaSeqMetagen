#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// k2db_libs = ['archaea', 'bacteria', 'plasmid', 'viral', 'human', 'fungi', 'plant', 'protozoa', 'nr', 'nt', 'UniVec', 'UniVec_Core']
db = file(params.db, type: 'dir')

process run_DownloadK2DBLibs {
    label 'mini'
    tag { "Download ${k2db_libs}" }
    errorStrategy 'ignore'

    input:
    each k2db_libs
    
    """
    kraken2-build --download-library ${k2db_libs} --db ${db} --use-ftp
    """
}

// process run_UpdateTaxonomy {
//     label 'mini'
//     tag { "Update NCBI Taxonomy" }
//     publishDir "${taxonomy}", mode: 'copy', overwrite: true
    
//     input:
//     path(dmp)
    
//     output:
//     path("*"), emit: taxonomy_update
    
//     """
//     /opt/KronaTools-2.8/updateTaxonomy.sh --only-build --preserve .
//     """
// }

