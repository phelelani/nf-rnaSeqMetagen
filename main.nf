#!/usr/bin/env nextflow
println "clear".execute().text
//  DO NOT EDIT FROM HERE!! - Unless you brave like King Shaka of course!
/*  ======================================================================================================
 *  HELP MENU
 *  ======================================================================================================
 */
line="=".multiply(100)
ver="nf-rnaSeqMetagen v0.2"
if (params.help) {
    println "\n${line}"
    println "#".multiply(48 - ("${ver}".size() / 2 )) + "  ${ver}   " + "#".multiply(48 - ("${ver}".size() / 2 ))
    println "${line}\n"
    // log.info ''
    // log.info "===================================="
    // log.info "         nf-rnaSeqMetagen v0.1        "
    // log.info "===================================="
    // log.info ''
    // log.info 'USAGE: '
    // log.info 'nextflow run main.nf --data /path/to/data --out /path/to/output --db /path/to/kraken-db --taxonomy /path/to/taxonomy --genome /path/to/genome.fa --index /path/to/STARIndex --bind /path/to/bind1;/path/to/bind2'
    // log.info ''
    // log.info 'HELP: '
    // log.info 'nextflow run main.nf --help'
    // log.info ''
    // log.info 'MANDATORY ARGUEMENTS:'
    // log.info '    --data     FOLDER    Path to where the input data is located (fastq | fq)'
    // log.info '    --out      FOLDER    Path to where the output should be directed (will be created if it does not exist).'
    // log.info '    --db       FOLDER    Path to where the Kraken database is installed'
    // log.info '    --taxonomy FOLDER    Path to where the taxonomy database is installed'
    // log.info '    --genome   FILE      The whole genome sequence (fasta | fa | fna)'
    // log.info '    --index    FOLDER    Path to where the STAR index files are locaded'
    // log.info '    --bind     FOLDER(S) Paths to be passed onto the singularity image'
    // log.info ''
    println "${line}\n"
    exit 1
}
pairedEnd  = null
singleEnd  = null
max_memory = 200.GB
max_cpus   = 24
max_time   = 24.h
//
//
/*  ======================================================================================================
 *  CHECK ALL USER INPUTS
 *  ======================================================================================================
 */

// USER PARAMETER INPUT: DATA DIRECTORY
// ---- THESE DO NOT REQUIRE DATA!!
if (params.mode in [ "prep.Containers", "prep.STARIndex", "prep.BowtieIndex", "prep.KrakenDB" ]) {
    if (params.data == null) {
        data_dir = "YOU HAVEN'T SPECIFIED THE DITA DIRECTORY YET! PLEASE SPECIFY BEFORE RUNNING THE WORKFLOW"
    } else{
        data_dir = file(params.data, type: 'dir')
    }
} else if (params.data == null && prarams.mode == "run.FilterClassify") {
    exit 1, "$data_error"
} else{
    data_dir = file(params.data, type: 'dir')
}

// USER PARAMETER INPUT: OUTPUT DIRECTORY
if(params.out == null) {
    out_dir = file("${PWD}/results_nf-rnaSeqMetagen", type: 'dir')
} else{
    out_dir = file(params.out, type: 'dir')
}

// USER PARAMETER INPUT: KRAKEN2 DB DIRECTORY
if(params.db== null) {
    db = file("${PWD}/kraken2_db", type: 'dir')
} else{
    db = file(params.db, type: 'dir')
}

// USER PARAMETER INPUT: GENOME FASTA FILE
if(params.genome == null) {
    exit 1, "$genome_error"
} else{
    genome = file(params.genome, type: 'file')
    index_dir = genome.getParent()
}

// USER PARAMETER INPUT: GENOME ANNOTATION FILE (GFT/GFF)
if(params.genes == null) {
    exit 1, "$genes_error"
} else{
    genes = file(params.genes, type: 'file') 
}

// USER INPUT MODE: WHICH ANALYSIS TO RUN!
if(params.mode == null) {
    exit 1, "$mode_error"
} else {
    mode = params.mode
}

// USER STRANDED MODE: ARE WE DOING PAIRED- OR SINGLE-END?
if(params.singleEnd == null && params.pairedEnd == null) {
    stranded = "paired-end"
} else if(params.singleEnd) {
    stranded = "single-end"
} else if(params.pairedEnd){
    stranded = "paired-end"
} else {}

// USER PARAMETER INPUT: PATHS TO BE BINDED TO THE IMAGE
bind_dir      = [ params.data, out_dir, db, new File("${params.genome}").getParent(), new File("${params.genes}").getParent() ]
    .unique()
    .collect { it -> "-B ${it}"}
    .join("\n" + ' '.multiply(26))
    .toString()

// OUTPUT DIRECTORIES
out_dir.mkdir()
filter_dir    = file("${out_dir}/filtering", type: 'dir')
multiqc_dir   = file("${out_dir}/report_MultiQC", type: 'dir')
ext           = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"

// case ['fastq','fq']:
//     ext = params.filetype
//     read_file_cmd = ''
//     break
// case ['fastq.gz','fq.gz']:
//     ext = params.filetype
//     read_file_cmd = '--readFilesCommand gunzip -c'
//     break
// case ['fastq.bz2','fq.bz2']:
//     ext = params.filetype
//     read_file_cmd = '--readFilesCommand bunzip2 -c'
//     break
// case null:
//     ext = 'fastq.gz'
//     read_file_cmd = '--readFilesCommand gunzip -c'
//     break
// }



//
/*  ======================================================================================================
 *  RUN INFO
 *  ======================================================================================================
 */
// log.info "===================================="
// log.info "           nf-rnaSeqMetagen         "
// log.info "===================================="
// log.info "Input data          : ${data_path}"
// log.info "Output path         : ${out_path}"
// log.info "Input file extension  : $ext"
// log.info "STAR readFile command : $read_file_cmd"
// log.info "Kraken database     : ${db}"
// log.info "Taxonomy database   : ${taxonomy}"
// log.info "Genome              : ${genome}"
// log.info "Genome Index (STAR) : ${index}"
// log.info "Paths to bind       : ${bind}"
// log.info "====================================\n"
//
//
/*  ======================================================================================================
 *  PIPELINE START
 *  ======================================================================================================
 */

// DATA ABSENT ERRORS
main_data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your input data! There are no FASTQ file in the directory:
\t${data_dir}
Please ensure that you have given me the correct directory for you FASTQ input reads using the \'--data\' option!
${line}
"""

// GET INPUT DATA DEPENDING ON USER "MODE"
if(mode in ["prep.Containers", "prep.STARIndex", "prep.BowtieIndex", "prep.KrakenDB"]) {
    // OPTIONS FOR PREPARING DATA
    switch (mode) {
        case ["prep.Containers"]:
            
            break
        case ["prep.STARIndex","prep.BowtieIndex"]:
            
            break
        case ["prep.KrakenDB"]:
            
            break
    }
} else if(mode == "run.FilterClassify") {
    // OPTIONS FOR PERFORMING THE ANALYSES
    if(stranded == "paired-end") {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{${ext}}", type: 'file')
            .ifEmpty { exit 1, "$main_data_error" }
    } else if(stranded == "single-end") {
        read_pairs = Channel.fromFilePairs("${data_dir}/*.{${ext}}", type: 'file', size:1)
            .ifEmpty { exit 1, "$main_data_error" }
    }
} else {
    exit 1, "$mode_error"
}

switch (mode) {
        // ========== THIS SECTION IS FOR PREPPING DATA (SINGULARITY IMAGES, STAR INDEXES AND BOWTIE INDEXES)
    case ['prep.Containers']: 
        base = "shub://phelelani/nf-rnaSeqMetagen:"
        images = Channel.from( ["${base}star", "${base}kraken2", "${base}upset", "${base}multiqc", "${base}trinity"] )
        
        process run_DownloadContainers {
            label 'mini'
            scratch '$HOME/tmp'
            tag { "Downloading: ${link}" }
            maxForks 2
            publishDir "$PWD/containers", mode: 'copy', overwrite: true
            
            input:
            each link from images
            
            output:
            file("*.sif") into containers
            
            """
            singularity pull nf-rnaSeqMetagen-${link.substring(32,)}.sif ${link}
            """
        }
        break
        // ==========
        
        //
    case ['prep.STARIndex']:
        process run_GenerateSTARIndex {
            label 'maxi'
            tag { "Generate Star Index" }
            publishDir "$index_dir", mode: 'copy', overwrite: true
            
            output:
            set val("starIndex"), file("*") into star_index
            
            """
            STAR --runThreadN ${task.cpus} \
                --runMode genomeGenerate \
                --genomeDir . \
                --genomeFastaFiles ${genome} \
                --sjdbGTFfile ${genes} \
                --sjdbOverhang 99
            """
        }

        star_index.subscribe { 
            println "\nSTAR index files generated:"
            it[1].each { 
                item -> println "\t${item}" 
            }
            println " "
        }
        break
        // ==========
        
        //
    case ['prep.BowtieIndex']:
        process run_GenerateBowtie2Index {
            label 'maxi'
            tag { "Generate Bowtie2 Index" }
            publishDir "$index_dir", mode: 'copy', overwrite: true
        
            output:
            set val("bowtieIndex"), file("*") into bowtie_index
            
            """
            bowtie2-build --threads ${task.cpus} ${genome} genome
            """
        }   
    
        bowtie_index.subscribe { 
            println "\nBowtie2 index files generated:"
            it[1].each { 
                item -> println "\t${item}" 
            }
            println " "
        }
        break

        
    case ['prep.KrakenDB']:
        process run_GenerateKrakenDB {
            cpus 7
            memory '200 GB'
            time '48h'
            scratch '$HOME/tmp'
            tag { "Generate Kraken DB" }
            publishDir "$db", mode: 'copy', overwrite: true
           
            output:
            file("*") into kraken_db
            
            """
            kraken2-build --standard --threads 6 --db kraken2DB
            """
        }
        break
        // ========== PREPPING STEPS/OPTIONS END HERE!

        // ========== THIS SECTION IS FOR THE MAIN WORKFLOW!
    case ['run.FilterClassify']:
        // 1.  ALIGN READS TO REFERENCE GENOME
        process run_STAR {
            cpus 13
            memory '50 GB'
            time '24h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$filter_dir/${sample}", mode: 'copy', overwrite: true, pattern: "${sample}*.{out,tab}"
    
            input:
            set sample, file(reads) from read_pair
    
            output:
            set sample, file("${sample}*.{out,tab}") into star_results
            set sample, file("${sample}_Unmapped*") into unmapped_kraken, unmapped_trinity
    
            """	
            /bin/hostname
            STAR --runMode alignReads \
                --genomeDir ${index} ${read_file_cmd} \
                --readFilesIn ${reads.get(0)} ${reads.get(1)} \
                --runThreadN 12 \
                --outSAMtype BAM Unsorted \
                --outReadsUnmapped Fastx \
                --outFileNamePrefix ${sample}_
       
            sed -i 's|\\s.[0-9]\$|\\/1|g' ${sample}_Unmapped.out.mate1 
            sed -i 's|\\s.[0-9]\$|\\/2|g' ${sample}_Unmapped.out.mate2

            """ 
        }

        // 2. Run KRAKEN to classify the raw reads that aren't mapped to the reference genome.
        process run_KrakenClassifyReads {
            cpus 13
            memory '50 GB'
            time '24h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$filter_dir/${sample}", mode: 'copy', overwrite: true, pattern: "${sample}_*.fastq"

            input:
            set sample, file(reads) from unmapped_kraken

            output:
            set sample, file("${sample}_reads.krak") into kraken_reads_report
            set sample, file("${sample}_classified_*.fastq") into kraken_classified_reads
            set sample, file("${sample}_unclassified_*.fastq") into kraken_unclassified_reads

            """	
            /bin/hostname
            kraken2 --db ${db} \
                --paired ${reads.get(0)} ${reads.get(1)} \
                --threads 12 \
                --classified-out ${sample}_classified#.fastq \
                --unclassified-out ${sample}_unclassified#.fastq \
                --output ${sample}_reads.krak
            """ 
        }

        // 3. Assemble the reads into longer contigs/sequences for classification.
        process run_TrinityAssemble {
            cpus 13
            memory '150 GB'
            time '48h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$filter_dir/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(reads) from unmapped_trinity

            output:
            set sample, "trinity_${sample}/Trinity.fasta" into trinity_assembled_reads

            """
            /bin/hostname
            Trinity --seqType fq \
               --max_memory 150G \
               --left ${reads.get(0)} --right ${reads.get(1)} \
               --SS_lib_type RF \
               --CPU 12 \
               --output  trinity_${sample}
            """
         }

        // 4. Run KRAKEN to classify the assembled FASTA sequences.
        process run_KrakenClassifyFasta {
            cpus 13
            memory '50 GB'
            time '24h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$filter_dir/${sample}", mode: 'copy', overwrite: true, pattern: "*.fasta"

            input:
            set sample, file(fasta) from trinity_assembled_reads

            output:
            set sample, file("${sample}_fasta.krak") into kraken_fasta_report
            set sample, file("${sample}_classified.fasta") into kraken_classified_fasta
            set sample, file("${sample}_unclassified.fasta") into kraken_unclassified_fasta

            """	
            kraken2 --db ${db} \
                ${fasta} \
                --threads 12 \
                --classified-out ${sample}_classified.fasta \
                --unclassified-out ${sample}_unclassified.fasta \
                --output ${sample}_fasta.krak
            """ 
        }

        // For each sample, create a list with [ SAMPLE_NAME, READ, FASTA ] by merging the classified outputs (reads and fasta) from KRAKEN 
        kraken_reads_report.join(kraken_fasta_report)
        .map { it -> [ it[0], [ it[1], it[2] ] ] }
        .set { all_kraken_reports }

        //5. Create that pretty KRONA report for all samples (reads and fasta)
        process run_KronaReport {
            cpus 1
            memory '1 GB'
            time '6h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$filter_dir/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(kraken) from all_kraken_reports

            output:
            set sample, file("*") into html
            set sample, file("*.kron") into krona_report

            """
            function createChart {
                cut -f 2,3 \$1 > \$(sed 's/.krak/.kron/' <<< "\$1")
                ktImportTaxonomy \$(sed 's/.krak/.kron/' <<< "\$1") \
                    -tax \$2 \
                    -o \$(sed 's/.krak/.html/' <<< "\$1")
                }
            createChart ${kraken.get(0)} ${taxonomy}
            createChart ${kraken.get(1)} ${taxonomy}
            """
        }

        //-tax ${taxonomy} \
        // 6a. Collect files for STAR QC
        star_results.collectFile () { item -> [ 'qc_star.txt', "${item.get(1).find { it =~ 'Log.final.out' } }" + ' ' ] }
        .set { qc_star }

        // 6. Get QC for STAR, HTSeqCounts and featureCounts
        process run_MultiQC {
            cpus 1
            memory '2 GB'
            time '1h'
            scratch '$HOME/tmp'
            tag { "Get QC Information" }
            publishDir "$multiqc_dir", mode: 'copy', overwrite: true

            input:
                file(star) from qc_star

            output:
                file('*') into multiQC

            """
            multiqc `< ${star}` --force
            """
        }

        // 7a. Collect all the krona file locations and put them in a text file
        krona_report.collectFile () { item -> [ 'fasta_krona_report.txt', "${item.get(1).find { it =~ 'fasta.kron' } }" + '\n' ] }
        .set { krona_report_list } 

        // 7b. Prepare data for creating the matrix for UpSet: json file, taxonomy file, and sample taxids 
        process run_PrepareMatrixData {
            cpus 1
            memory '1 GB'
            time '1h'
            scratch '$HOME/tmp'
            tag { "Prepare Matrix Data" }
            publishDir "$filter_dir", mode: 'copy', overwrite: true
            echo 'true'

            input:
            file(list) from krona_report_list

            output:
            file("*") into whole_list
            file("upset/data/nf-rnaSeqMetagen/*") into final_list

            shell:
            database = "${db}"
            file_list = "${list}"
            json_file = "nf-rnaSeqMetagen.json"
            names_file = "names_table.dmp"
            template 'get_taxons.sh'
        }

        // 7 Create the UpSet matrix
        process run_CreateMatrix {
            cpus 1
            memory '1 GB'
            time '1h'
            scratch '$HOME/tmp'
            tag { "Create UpSet Matrix" }
            publishDir "$multiqc_dir/upset/data/nf-rnaSeqMetagen", mode: 'copy', overwrite: true
            echo 'true'

            input:
            file(list) from final_list

            output:
            file("*") into the_matrix

            shell:
            template 'create_matrix.R'
        }
        break
        // ==========
}

/*  ======================================================================================================
 *  WORKFLOW SUMMARY 
 *  ======================================================================================================
 */
summary="nf-rnaSeqMetagen v0.2 - Execution Summary:"
workflow.onComplete {
    println "\n${line}"
    println "#".multiply(48 - ("${summary}".size() / 2 )) + "  ${summary}  " + "#".multiply(48 - ("${summary}".size() / 2 ))    
    println "Execution command   : ${workflow.commandLine}"
    println "Execution name      : ${workflow.runName}"
    println "Workflow start      : ${workflow.start}"
    println "Workflow end        : ${workflow.complete}"
    println "Workflow duration   : ${workflow.duration}"
    println "Workflow completed? : ${workflow.success}"
    println "Work directory      : ${workflow.workDir}"
    println "Project directory   : ${workflow.projectDir}"
    println "Execution directory : ${workflow.launchDir}"
    println "Configuration files : ${workflow.configFiles}"
    println "Workflow containers : ${workflow.container}"
    println "exit status         : ${workflow.exitStatus}"
    println "Error report        : ${workflow.errorReport ?: '-'}"
    println "${line}\n"
    println "\n"
}

workflow.onError {
    println "Oohhh DANG IT!!... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
//======================================================================================================
