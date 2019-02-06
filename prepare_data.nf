#!/usr/bin/env nextflow

/* PARAMETERS NEEDED:
 --genome      : Can be specified in prepare_data.config
 --genes       : Can be specified in prepare_data.config
 --kraken_dir  : Can be specified in prepare_data.config
 --bind        : Can be specified in prepare_data.config
 --mode        : < getContainers | generateStarIndex | generateBowtieIndex | generateKrakenDB >
 -profile      : < prepare | pbsPrepare >
*/

def checkGenome() {
    if(params.genome == null) {
        exit 1, "Please provide a FASTA sequence of the reference genome."
    } else{
        genome = file(params.genome, type: 'file')  // The whole genome sequence.
    }
    return genome
}

def checkGenes() {    
    if(params.genes == null) {
        exit 1, "Please provide an annotation GTF file."
    } else{
        genes = file(params.genes, type: 'file')  // The genome annotation file.
    }
    return genes
}

def checkKrakenDir() {
    if(params.kraken_dir == null) {
        exit 1, "Please provide a path to save the Kraken database"
    } else{
        out_path = file(params.kraken_dir, type: 'dir') 
    }
    return out_path
}


switch (params.mode) {
    case ['getContainers']:
        link_base = "shub://phelelani/nf-rnaSeqMetagen:"
        shub_images = Channel.from( ["${link_base}star", "${link_base}kraken2", "${link_base}trinity", "${link_base}multiqc"] )
        
        process downloadContainers_process {
            cpus 1
            memory '2 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { "Downloading: $link" }
            publishDir "$baseDir/containers", mode: 'copy', overwrite: true, pattern: "*.simg"
            
            input:
            each link from shub_images
            
            output:
            file("*.simg") into containers
        
            """
            singularity pull ${link}
            """
        }
        break
        
    case ['generateStarIndex']:
        
        checkGenome()
        checkGenes()
        out_path = genome.getParent()
        
        process generateSTARIndex {
            cpus 13
            memory '100 GB'
            time '24h'
            scratch '$HOME/tmp'
            tag { "Generate Star Index" }
            publishDir "$out_path", mode: 'copy', overwrite: true
            
            output:
                file("*") into star_index
            
            """
            STAR --runThreadN 12 \
                --runMode genomeGenerate \
                --genomeDir . \
                --genomeFastaFiles ${genome} \
                --sjdbGTFfile ${genes} \
                --sjdbOverhang 99
            """
        }
        break

    case ['generateBowtieIndex']:
        
        checkGenome()
        out_path = genome.getParent()
        
        process generateBowtie2Index {
            container "$baseDir/containers/phelelani-nf-rnaSeqMetagen-master-trinity.simg"
            cpus 13
            memory '100 GB'
            time '24h'
            scratch '$HOME/tmp'
            tag { "Generate Bowtie2 Index" }
            publishDir "$out_path", mode: 'copy', overwrite: false
        
            output:
            file("*") into bowtie_index
            
            """
            bowtie2-build --threads 12 ${genome} genome
            """
        }   
        break

    case ['generateKrakenDB']:
        
        checkKrakenDir()
        
        process generateKrakenDB {
            cpus 7
            memory '200 GB'
            time '48h'
            scratch '$HOME/tmp'
            tag { "Generate Kraken DB" }
            publishDir "$out_path", mode: 'copy', overwrite: true
            
            output:
            file("*") into kraken_db
        
            """
            kraken2-build --standard --threads 6 --db kraken2_std
            kraken2-build --cleanup --db kraken2_std
            """
        }   

        process generateKrakenDB {
            cpus 1
            memory '5 GB'
            time '48h'
            scratch '$HOME/tmp'
            tag { "Generate Taxonomy" }
            publishDir "$out_path", mode: 'copy', overwrite: true

            output:
            file("*") into taxonomy
            
            """
            /opt/KronaTools-2.7/updateTaxonomy.sh taxonomy
            """
        }   
        break
}
