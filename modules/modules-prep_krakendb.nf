#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process run_DownloadK2DBLibs {
    label 'mini'
    tag { "Download ${k2db_libs}" }

    input:
    each k2db_libs
    
    """
    kraken2-build --download-library ${k2bd_libs} --db ${db}
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

