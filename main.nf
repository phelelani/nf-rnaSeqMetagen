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
} else if (params.data == null && params.mode == "run.FilterClassify") {
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
    db = file("${PWD}/kraken2db", type: 'dir')
    taxonomy = file("$db/taxonomy", type: 'dir')
} else{
    db = file(params.db, type: 'dir')
    taxonomy = file("$db/taxonomy", type: 'dir')
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
// filter_dir    = file("${out_dir}/filtering", type: 'dir')
// multiqc_dir   = file("${out_dir}/report_MultiQC", type: 'dir')
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


//  ======================================================================================================
//  RUN INFO
//  ======================================================================================================
options="nf-rnaSeqMetagen v0.2 - Input/Output and Parameters:"
println "\n" + "=".multiply(100)
println "#".multiply(48 - ("${options}".size() / 2 )) + "  ${options}  " + "#".multiply(48 - ("${options}".size() / 2 ))
println "=".multiply(100)
println "Input data              : $data_dir"
println "Input data type         : $stranded"
println "Output directory        : $out_dir"
println "Kraken2 DB directory    : $db"
// println ' '.multiply(26) + "- ${filter_dir.baseName}"
// println ' '.multiply(26) + "- ${multiqc_dir.baseName}"
println "Genome                  : $genome"
println "Genome annotation       : $genes"
println "Paths to bind           : $bind_dir"
println "=".multiply(100)
println " "

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

mode_error = """
${line}
Oooh no!! Looks like there's an serious issue in your command! 
I do not recognise the \'--mode ${params.mode}\' option you have given me, or you have not given me any \'--mode\' option at all!
The allowed options for \'--mode\' are:
\tprep.Containers\t\t: For downloading Singularity containers used in this workflow.
\tprep.STARIndex\t\t: For indexing your reference genome using STAR.
\tprep.BowtieIndex\t: For indexing your reference genome using Bowtie2.
\trun.ReadQC\t\t: For performing general QC on your reads using FastQC. 
\trun.ReadTrimming\t: For trimming low quality bases and removing adapters from your reads using Trimmmomatic.
\trun.ReadAlignment\t: For aligning your reads to your reference genome using STAR.
\trun.ReadCounting\t: For counting features in your reads using HTSeq-count and featureCounts.
\trun.MultiQC\t\t: For getting a summary of QC through the analysis using MultiQC.
\nPlease use one of the above options with \'--mode\' to run the nf-rnaSeqCount workflow!
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
            tag { "Downloading: ${link}" }
            publishDir "$PWD/containers", mode: 'copy', overwrite: true
            
            input:
            each link from images
            
            output:
            file("*.sif") into containers
            
            """
            singularity pull nf-rnaSeqMetagen-${link.substring(34,)}.sif ${link}
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
        process run_GenerateBowtieIndex {
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
            label 'maxi'
            tag { "Generate Kraken DB" }
            publishDir "$db", mode: 'copy', overwrite: true
           
            output:
            file("*") into kraken_db
            
            """
            kraken2-build --standard --threads ${task.cpus} --db .
            """
        }
        break
        // ========== PREPPING STEPS/OPTIONS END HERE!

        // ========== THIS SECTION IS FOR THE MAIN WORKFLOW!
    case ['run.FilterClassify']:
        // 1.  ALIGN READS TO REFERENCE GENOME
        process run_STAR {
            label 'maxi'
            tag { sample }
            // publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true, pattern: "${sample}*.{out,tab}"
    
            input:
            set sample, file(reads) from read_pairs
    
            output:
            set sample, file("${sample}*.{out,tab}") into star_results
            set sample, file("${sample}_Unmapped*") into unmapped_reads
    
            """	
            /bin/hostname
            STAR --runMode alignReads \
                --genomeDir ${index_dir} \
                --readFilesCommand gunzip -c \
                --readFilesIn ${reads.findAll().join(' ')} \
                --runThreadN ${task.cpus} \
                --outSAMtype BAM Unsorted \
                --outReadsUnmapped Fastx \
                --outFileNamePrefix ${sample}_
            """ 
        }

        process run_FixSeqNames {
            label 'mini'
            tag { sample }
            // publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(unmapped) from unmapped_reads
            
            output:
            set sample, file("${sample}_unmapped*") into unmapped_kraken, unmapped_trinity

            """
            sed 's| \\(.*\\)\$|\\/1|g' ${unmapped.find { it =~ 'mate1' } } > ${sample}_unmapped_R1.fastq
            sed 's| \\(.*\\)\$|\\/2|g' ${unmapped.find { it =~ 'mate2' } } > ${sample}_unmapped_R2.fastq
            """
        }
        
        // 2. Run KRAKEN to classify the raw reads that aren't mapped to the reference genome.
        process run_KrakenClassifyReads {
            label 'maxi'
            tag { sample }
            // publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(reads) from unmapped_kraken

            output:
            set sample, file("${sample}_reads.krak") into kraken_reads_report
            set sample, file("${sample}_classified_*.fastq") into kraken_classified_reads
            set sample, file("${sample}_unclassified_*.fastq") into kraken_unclassified_reads

            """	
            /bin/hostname
            kraken2 --db ${db} \
                --paired ${reads.findAll().join(' ')} \
                --threads ${task.cpus} \
                --classified-out ${sample}_classified#.fastq \
                --unclassified-out ${sample}_unclassified#.fastq \
                --output ${sample}_reads.krak
            """ 
        }

        // 3. Assemble the reads into longer contigs/sequences for classification.
        process run_TrinityAssemble {
            label 'maxi'
            memory '150 GB'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(reads) from unmapped_trinity

            output:
            set sample, "trinity_${sample}/Trinity.fasta" into trinity_assembled_reads

            """
            /bin/hostname
            Trinity --seqType fq \
               --max_memory 150G \
               --left ${reads.find { it =~ 'R1' } } \
               --right ${reads.find { it =~ 'R2' } } \
               --SS_lib_type RF \
               --CPU ${task.cpus} \
               --output trinity_${sample}
            """
         }

        // 4. Run KRAKEN to classify the assembled FASTA sequences.
        process run_KrakenClassifyFasta {
            label 'maxi'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true, pattern: "*{_fasta.krak,_classified.fasta}"

            input:
            set sample, file(fasta) from trinity_assembled_reads

            output:
            set sample, file("${sample}_fasta.krak") into kraken_fasta_report
            set sample, file("${sample}_classified.fasta") into kraken_classified_fasta
            set sample, file("${sample}_unclassified.fasta") into kraken_unclassified_fasta

            """	
            kraken2 --db ${db} \
                ${fasta} \
                --threads ${task.cpus} \
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
            label 'mini'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true

            input:
            set sample, file(kraken) from all_kraken_reports

            output:
            set sample, file("*.html") into html
            set sample, file("*reads.kron") into reads_krona
            set sample, file("*fasta.kron") into fasta_krona, fasta_krona_seqs

            """
            function createChart {
                awk '\$1=="C" { print \$2"\t"\$3 }' \$1 > \$(sed 's/.krak/.kron/' <<< "\$1")
                ktImportTaxonomy \$(sed 's/.krak/.kron/' <<< "\$1") \
                    -tax \$2 \
                    -o \$(sed 's/.krak/.html/' <<< "\$1")
                }
            createChart ${kraken.get(0)} ${taxonomy}
            createChart ${kraken.get(1)} ${taxonomy}
            """
        }

        fasta_krona_seqs.join(kraken_classified_fasta)
        .map { it -> [ it[0], [ it[1], it[2] ] ] }
        .set { krona_fasta_pair }
        
        process run_CollectTaxSeqs {
            label 'mini'
            tag { sample }
            publishDir "$out_dir/${sample}/taxon_sequences", mode: 'copy', overwrite: true

            input: 
            set sample, file(test) from krona_fasta_pair
            
            output:
            set sample, file("taxid_*.fasta") into taxon_sequences
            
            """
            for id in \$(awk '{ print \$2 }' ${test.get(0)} | sort -gu)
            do
                grep -A1 --no-group-separator "kraken:taxid|\$id\$" ${test.get(1)} | sed 's/path=\\[\\(.*\\)\\] //' > "taxid_"\$id".fasta"
            done
            """
        }
        
        // 6a. Collect files for STAR QC
        star_results.collectFile () { item -> [ 'qc_star.txt', "${item.get(1).find { it =~ 'Log.final.out' } }" + ' ' ] }
        .set { qc_star }

        // 6. Get QC for STAR, HTSeqCounts and featureCounts
        process run_MultiQC {
            label 'mini'
            tag { "Get QC Information" }
            publishDir "$out_dir/report_multiqc", mode: 'copy', overwrite: true

            input:
            file(star) from qc_star

            output:
            file('*') into multiQC

            """
            multiqc `< ${star}` --force
            """
        }

        // 7a. Collect all the krona file locations and put them in a text file
        fasta_krona.collectFile () { item -> [ 'fasta_krona_files.txt', "${item.get(1)}" + '\n' ] }
        .set { fasta_krona_list } 

        // fasta_krona_report.collectFile () { item -> [ 'fasta_krona_report.txt', "${item.get(1).find { it =~ 'fasta.kron' } }" + '\n' ] }

        process run_CopyUpsetDir {
            label 'mini'
            tag { "Copy UpSet Tool" }
            publishDir "$out_dir/upset", mode: 'copy', overwrite: true
            
            output:
            file("*") into upset_dir

            """
            /bin/hostname
            rsync -avhP /opt/upset/* .
            """
        }
        
        // 7b. Prepare data for creating the matrix for UpSet: json file, taxonomy file, and sample taxids 
        process run_PrepareMatrixData {
            label 'mini'
            tag { "Prepare Matrix Data" }

            input:
            file(list) from fasta_krona_list

            output:
            file("*.{dmp,taxon,json}") into matrix_files

            shell:
            tax_names = "${taxonomy}"
            file_list = "${list}"
            json_file = "nf-rnaSeqMetagen.json"
            names_file = "names_table.dmp"
            template 'get_taxons.sh'
        }

        // 7 Create the UpSet matrix
        process run_CreateMatrix {
            label 'mini'
            tag { "Create UpSet Matrix" }
            publishDir "$out_dir/upset/data/nf-rnaSeqMetagen", mode: 'copy', overwrite: true

            input:
            file(list) from matrix_files

            output:
            set val("upset_files"), file(list), file("*") into the_matrix

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
    println "${line}\n"
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
