#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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

include { run_GenerateSTARIndex } from './modules/modules-prep_indexes.nf'
include { run_DownloadK2DBIndexes; run_DownloadTaxonomy; run_UpdateTaxonomy } from './modules/modules-prep_krakendb.nf'
include { run_STAR; run_FixSeqNames; run_KrakenClassifyReads;
         run_TrinityAssemble; run_KrakenClassifyFasta; run_KronaReport;
         run_CollectTaxSeqs; run_MultiQC; run_CopyUpsetDir;
         run_PrepareMatrixData; run_CreateMatrix } from './modules/modules-filter_classify.nf'

// 
workflow PREP_INDEXES {
    main:
    run_GenerateSTARIndex()
}

// 
workflow PREP_KRAKENDB {
    main:
    run_DownloadK2DBIndexes()
    run_DownloadTaxonomy()
    run_UpdateTaxonomy(run_DownloadTaxonomy.out.taxonomy)
}

// 
workflow FILTER_CLASSIFY {
    take:
    read_pairs        

    main:
    run_STAR(read_pairs)
    run_FixSeqNames(run_STAR.out.unmapped_reads)
    run_KrakenClassifyReads(run_FixSeqNames.out.unmapped_reads)
    run_TrinityAssemble(run_FixSeqNames.out.unmapped_reads)
    run_KrakenClassifyFasta(run_TrinityAssemble.out.trinity_assembled_reads)
    run_KrakenClassifyReads.out.kraken_reads_report
        .join(run_KrakenClassifyFasta.out.kraken_fasta_report)
        .map { it -> [ it[0], [ it[1], it[2] ] ] }
        .set { all_kraken_reports }
    run_KronaReport(all_kraken_reports)
    run_KronaReport.out.fasta_krona
        .map { it -> [ it[0], [ it[1], it[2] ] ] }
        .set { krona_fasta_pair }
    krona_fasta_pair.view()
    // run_CollectTaxSeqs(krona_fasta_pair)
    // run_STAR.out.star_results
    //     .collectFile() { item -> [ 'qc_star.txt', "${item.get(1).find { it =~ 'Log.final.out' } }" + ' ' ] }
    //     .set { qc_star }
    // run_MultiQC(qc_star)
    // run_KronaReport.out.fasta_krona
    //     .collectFile() { item -> [ 'fasta_krona_files.txt', "${item.get(1)}" + '\n' ] }
    //     .set { fasta_krona_list }
    // run_CopyUpsetDir()
    // run_PrepareMatrixData(fasta_krona_list)
    // run_CreateMatrix(run_PrepareMatrixData.out.matrix_files)
}

workflow {
    switch (mode) {
        case ['prep.GenomeIndexes']:
            PREP_INDEXES()
            break
        case ['prep.KrakenDB']:
            PREP_KRAKENDB()
            break
        case ['run.FilterClassify']:
            FILTER_CLASSIFY(read_pairs)
            break
        default:
            exit 1, "NO WORKFLOW GIVEN!"
            break
    }
}

//         summary="nf-rnaSeqMetagen v0.2 - Execution Summary:"
//         workflow.onComplete {
//             println "\n${line}"
//             println "#".multiply(48 - ("${summary}".size() / 2 )) + "  ${summary}  " + "#".multiply(48 - ("${summary}".size() / 2 ))    
//             println "${line}\n"
//             println "Execution command   : ${workflow.commandLine}"
//             println "Execution name      : ${workflow.runName}"
//             println "Workflow start      : ${workflow.start}"
//             println "Workflow end        : ${workflow.complete}"
//             println "Workflow duration   : ${workflow.duration}"
//             println "Workflow completed? : ${workflow.success}"
//             println "Work directory      : ${workflow.workDir}"
//             println "Project directory   : ${workflow.projectDir}"
//             println "Execution directory : ${workflow.launchDir}"
//             println "Configuration files : ${workflow.configFiles}"
//             println "Workflow containers : ${workflow.container}"
//             println "exit status         : ${workflow.exitStatus}"
//             println "Error report        : ${workflow.errorReport ?: '-'}"
//             println "${line}\n"
//             println "\n"
//         }

// }
// workflow.onError {
//     println "Oohhh DANG IT!!... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
// }
