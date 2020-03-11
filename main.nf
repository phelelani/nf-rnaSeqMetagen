#!/usr/bin/env nextflow
println "clear".execute().text

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
    println "USAGE:"
    println "nextflow run nf-rnaSeqMetagen -profile slurm --data /path/to/data --genome path/to/genome.fa --genes /path/to/genes.gtf --db /path/to/db\n" 
    println "HELP:"
    println "nextflow run nf-rnaSeqMetagen --help\n"
    println "MANDATORY ARGUEMENTS:"
    println "-profile     STRING    Executor to be used. Available options:"
    println "\t\t\t\t\"standard\"          : Local execution (no job scheduler)."
    println "\t\t\t\t\"slurm\"             : SLURM scheduler."
    println "--mode       STRING    To specify which step of the workflow you are running (see https://github.com/phelelani/nf-rnaSeqMetagen)."
    println "                       Availeble options:"
    println "\t\t\t\t\"prep.Containers\"   : For downloading Singularity containers used in this workflow."
    println "\t\t\t\t\"prep.GenomeIndexes\"    : For indexing your reference genome using STAR."
    println "\t\t\t\t\"prep.KrakenDB\"     : For building the Kraken2 database."
    println "\t\t\t\t\"run.FilterClassify\": For performing metagenomics analysis, i.e., filtering and classification.\n"
    println "--data       FOLDER    Path to where the input data (FASTQ files) is located. Supported FASTQ files:"
    println "\t\t\t\t[ fastq | fastq.gz | fastq.bz2 | fq | fq.gz | fq.bz2 ]"
    println "--genome     FILE      The whole genome FASTA sequence. Supported FASTA files:"
    println "\t\t\t\t[ fasta | fa | fna ]"
    println "--genes      FILE      The genome annotation GFT file. Supported GTF file:"
    println "\t\t\t\t[ gtf ]"
    println "--db         FOLDER    Path to where the Kraken2 database will be saved (or where it is located if already created)."
    println "                       Default: \$PWD/kraken2db"
    println "OPTIONAL ARGUEMENTS:"
    println "--help                 To show this menu."
    println "--out        FOLDER    Path to where the output should be directed."
    println "                       Default: \$PWD/results_nf-rnaSeqMetagen"
    println "--max_memory STRING    Maximum memory you have access to."
    println "                       Default: \"200.GB\""
    println "--max_cpus   STRING    Maximum CPUs you have access to."
    println "                       Default: \"24\""
    println "--max_time   STRING    Maximum time you have access to."
    println "                       Default: \"24.h\""
    println "${line}\n"
    exit 1
}

/*  ======================================================================================================
 *  CHECK ALL USER INPUTS
 *  ======================================================================================================
 */

// MAIN USER INPUT ERRORS
mode_error = """
${line}
Oooh no!! Looks like there's an serious issue in your command! 
I do not recognise the \'--mode ${params.mode}\' option you have given me, or you have not given me any \'--mode\' option at all!
The allowed options for \'--mode\' are:
\tprep.Containers\t\t: For downloading Singularity containers used in this workflow.
\tprep.GenomeIndexes\t\t: For indexing your reference genome using STAR.
\tprep.KrakenDB\t\t: For building the Kraken2 database.
\trun.FilterClassify\t: For performing metagenomics analysis, i.e., filtering and classification.
\nPlease use one of the above options with \'--mode\' to run the nf-rnaSeqMetagen workflow!
${line}
"""

data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--data ${params.data}\' option you have given me, or you have not given me any \'--data\' option at all!
Please provide a valid directory with you input FASTQ reads with the \'--data\' option to run the nf-rnaSeqMetagen workflow! 
${line}
"""

genome_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--genome ${params.genome}\' option you have given me, or you have not given me any \'--genome\' option at all!
Please provide a valid FASTA file (.fasta or .fa) for your reference genome with the \'--genome\' option to run the nf-rnaSeqMetagen workflow! 
${line}
"""

genes_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--genes ${params.genes}\' option you have given me, or you have not given me any \'--genes\' option at all!
Please provide a valid GTF annotation file (.gtf) for your reference genome with the \'--genes\' option to run the nf-rnaSeqMetagen workflow! 
${line}
"""

db_error = """
${line}
Oooh no!! Looks like there's a serious issue in your command! 
I do not recognise the \'--db ${params.db}\' option you have given me, or you have not given me any \'--db\' option at all!
Please provide a valid directory with your Kraken database with the \'--db\' option to run the nf-rnaSeqMetagen workflow! 
${line}
"""

// EMPTY LIST FOR COLLECTING ALL THE PATHS TO BIND TO SINGULARITY IMAGE
bind_dirs = []

// USER PARAMETER INPUT: DATA DIRECTORY
switch (params.data) {
    case [null]:
        data_dir = "NOT SPECIFIED!"
        break
    default:
        data_dir = file(params.data, type: 'dir')
        bind_dirs.add(data_dir)
}

// USER PARAMETER INPUT: GENOME FASTA FILE
switch (params.genome) {
    case [null]:
        genome = "NOT SPECIFIED!"
        break
    default:
        genome = file(params.genome, type: 'file', checkIfExists: true)
        index_dir = genome.getParent()
        bind_dirs.add(genome.getParent())
}

// USER PARAMETER INPUT: GENOME ANNOTATION FILE (GFT/GFF)
switch (params.genes) {
    case [null]:
        genes = "NOT SPECIFIED!"
        break
    default:
        genes = file(params.genes, type: 'file', checkIfExists: true)
        bind_dirs.add(genes.getParent())
}

// USER PARAMETER INPUT: OUTPUT DIRECTORY
switch (params.out) {
    case [null]:
        out_dir = file("${PWD}/results_nf-rnaSeqMetagen", type: 'dir')
        break
    default:
        out_dir = file(params.out, type: 'dir')
        bind_dirs.add(out_dir)
}

// USER PARAMETER INPUT: KRAKEN2 DB DIRECTORY
switch (params.db) {
    case [null]:
        db = "NOT SPECIFIED!"
        break
    default:
        db = file(params.db, type: 'dir')
        db.mkdir()
        taxonomy = file("$db/taxonomy", type: 'dir')
        bind_dirs.add(db)
}

def breakIfNull(parameter,error) {
    if (parameter == null) {
        exit 1, error
    } else {}
}

// DATA ABSENT ERRORS
main_data_error = """
${line}
Oooh no!! Looks like there's a serious issue in your input data! There are no FASTQ file in the directory:
\t${data_dir}
Please ensure that you have given me the correct directory for you FASTQ input reads using the \'--data\' option!
${line}
"""

// USER INPUT MODE: WHICH ANALYSIS TO RUN!
switch (params.mode) {
    case [null]:
        exit 1, "$mode_error"
        
    // DATA NOT REQUIRED FOR THESE MODES!
    case [ "prep.Containers", "prep.GenomeIndexes", "prep.KrakenDB" ]:
        mode = params.mode
        
        switch (mode) {
            case ["prep.Containers"]:
                break
                
            case ["prep.GenomeIndexes"]:
                breakIfNull(params.genome,"$genome_error")
                // breakIfNull(params.genes,"$genes_error")
                break
                
            case ["prep.KrakenDB"]:
                breakIfNull(params.db,"$db_error")
                break
        }
        break
        

    // MAIN WORKFLOW PARAMETER CHECKS
    case ["run.FilterClassify"]:
        mode = params.mode
        
        // BREAK THE WORKFLOW IF THE FOLLOWING PARAMETERS AREN'T SPECIFIED!
        breakIfNull(params.data,"$data_error")
        breakIfNull(params.genome,"$genome_error")
        breakIfNull(params.genes,"$genes_error")
        breakIfNull(params.db,"$db_error")

        // GET THE INPUT DATA!
        ext = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"
        
        // GET DATA
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read,_}[1,2]*.{${ext}}", type: 'file')
            .ifEmpty { exit 1, "$main_data_error" }
        
        // OUTPUT DIRECTORIES
        out_dir.mkdir()
        break
    default:
        exit 1, "$mode_error"
        break
}

// USER PARAMETER INPUT: PATHS TO BE BINDED TO THE IMAGE
bind_dirs = bind_dirs
    .unique()
    .collect { it -> "-B ${it}"}
    .join("\n" + ' '.multiply(26))
    .toString()

/*  ======================================================================================================
 *  RUN INFO
 *  ======================================================================================================
 */
options="nf-rnaSeqMetagen v0.2 - Input/Output and Parameters:"
println "\n" + "=".multiply(100)
println "#".multiply(48 - ("${options}".size() / 2 )) + "  ${options}  " + "#".multiply(48 - ("${options}".size() / 2 ))
println "=".multiply(100)
println "Input data              : $data_dir"
println "Output directory        : $out_dir"
println "Genome                  : $genome"
println "Genome annotation       : $genes"
println "Kraken2 DB directory    : $db"
println "Paths to bind           : $bind_dirs"
println "=".multiply(100)
println " "

/*  ======================================================================================================
 *  PIPELINE START
 *  ======================================================================================================
 */

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
        
    case ['prep.GenomeIndexes']:
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


    //     star_index.subscribe { 
    //         println "\nSTAR index files generated:"
    //         it[1].each { 
    //             item -> println "\t${item}" 
    //         }
    //         println " "
    //     }
    //     break
    //     // ==========
        
    // case ['prep.BowtieIndex']:
        
        // process run_GenerateBowtieIndex {
        //     label 'maxi'
        //     tag { "Generate Bowtie2 Index" }
        //     publishDir "$index_dir", mode: 'copy', overwrite: true
        
        //     output:
        //     set val("bowtieIndex"), file("*") into bowtie_index
            
        //     """
        //     bowtie2-build --threads ${task.cpus} ${genome} genome
        //     """
        // }   
    
        // bowtie_index.subscribe { 
        //     println "\nBowtie2 index files generated:"
        //     it[1].each { 
        //         item -> println "\t${item}" 
        //     }
        //     println " "
        // }
        break
        // ==========
        
    case ['prep.KrakenDB']:
        process run_GenerateKrakenDB {
            label 'maxi'
            tag { "Generate Kraken DB" }
            publishDir "$db", mode: 'copy', overwrite: false
           
            output:
            file("*.k2d") into kraken_db
            file("taxonomy/taxdump.tar.gz") into taxonomy_dump
            
            """
            kraken2-build --standard --threads ${task.cpus} --db .
            """
        }

        process run_UpdateTaxonomy {
            label 'mini'
            tag { "Update NCBI Taxonomy" }
            publishDir "$taxonomy", mode: 'copy', overwrite: true
            
            input:
            file(dmp) from taxonomy_dump
            
            output:
            file("*") into taxonomy_update
            
            """
            /opt/KronaTools-2.7/updateTaxonomy.sh --only-build --preserve .
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


            // name="${reads.get(0)}"
            // if [[ "\$name" =~ ".fastq.gz" || "\$name" =~ ".fq.gz" ]]
            // then
            //     read_file_cmd = '--readFilesCommand gunzip -c'
            // elif [[ "\$name" =~ ".fastq.bz2" || "\$name" =~ ".fq.bz2" ]]
            // then 
            //     read_file_cmd = '--readFilesCommand bunzip2 -c'
            // else
            //     read_file_cmd = ""
            // fi


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
            cpus 20
            maxForks 5
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
            publishDir "$out_dir/MultiQC", mode: 'copy', overwrite: true

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
