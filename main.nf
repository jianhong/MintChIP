#!/usr/bin/env nextflow
/*
========================================================================================
                         qiubio-nf/mintchip
========================================================================================
 qiubio-nf/mintchip Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/jianhong/MintChIP
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def json_schema = "$projectDir/nextflow_schema.json"

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run jianhong/MintChIP -r dev --barcode barcode.txt --genome GRCh37 -profile conda"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help("$projectDir/nextflow_schema.json", command)
    exit 0
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
run_name = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    run_name = workflow.runName
}

if(params.barcode){ch_barcode = file(params.barcode, checkIfExists: true)}else{exit 1, "barcode is required!"}

// Open channel for left and right files and merge it into triples, the
// first entry is the LCS of the file names that can be used as a read
// pair identifier.
ch_input = Channel.fromPath("${params.input}").splitCsv(header: true)
ch_fastq = ch_input.map{ 
    row -> 
        fastq1 = row.remove("fastq_1")
        fastq2 = row.remove("fastq_2")
        meta = row
        meta.id = row.group + "_R" + row.replicate
        [meta, [fastq1, fastq2]]
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
/*
 * Print parameter summary
 */
def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --    IMPORT LOCAL MODULES/SUBWORKFLOWS     -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////
def getParam(modules, module) {
    return modules[module]?:[:]
}
include { JO_CHECKSUMS   } from './modules/checksum/checksum' addParams(options: getParam(modules, 'checksum'))
include { JE_DEMULTIPLEX } from './modules/demultiplex/je' addParams(options: getParam(modules, 'je'))
include { JO_TRIM        } from './modules/trimmomatic/trimmomatic' addParams(options: getParam(modules, 'trimmomatic'))
include { JO_CHIPSEQ     } from './modules/chipseq/chipseq' addParams(options: getParam(modules, 'jo_chipseq'))

workflow {
    /*
     * checksum
     */
    JO_CHECKSUMS(ch_fastq) 
    /*
     * Demultiplex
     */
    JE_DEMULTIPLEX(ch_fastq, ch_barcode)
    
    /*
     * Split the channel
     */
    JE_DEMULTIPLEX.out.fastq
        .transpose()
        .map{
            meta, fastq1, fastq2 ->
                info = meta.clone()
                info.remove("md5_1")
                info.remove("md5_2")
                info.sample_name = fastq1.simpleName
                info.id = info.id + "_" + fastq1.simpleName
                [info, fastq1, fastq2]
        }
        .set{paired_reads}
    
    /*
     * Trimm
     */
    JO_TRIM(paired_reads)
    
    /*
     * Run ChIPseq pipeline
     */
    JO_TRIM.out.paired
        .map{
            meta, fastq1, fastq2 ->
                [
                 group:meta.group+"_"+meta.sample_name, 
                 replicate:meta.replicate, 
                 fastq_1:fastq1,
                 fastq_2:fastq2,
                 antibody:meta.antibody,
                 control:meta.control
                 ]
        }
        .collect()
        .set{chipseq_input}
    JO_CHIPSEQ(chipseq_input)
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    def multiqc_report = []
    Completion.email(workflow, params, summary_params, run_name, projectDir, multiqc_report, log)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////