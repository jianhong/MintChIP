/*
 * This file holds several functions used to perform standard checks for the nf-core pipeline template.
 */

class Checks {

    static void aws_batch(workflow, params) {
        if (workflow.profile.contains('awsbatch')) {
            assert !params.awsqueue || !params.awsregion : "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
            // Check outdir paths to be S3 buckets if running on AWSBatch
            // related: https://github.com/nextflow-io/nextflow/issues/813
            assert !params.outdir.startsWith('s3:') : "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
            // Prevent trace files to be stored on S3 since S3 does not support rolling files.
            assert params.tracedir.startsWith('s3:') :  "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
        }
    }

    static void hostname(workflow, params, log) {
        Map colors = Headers.log_colours(params.monochrome_logs)
        if (params.hostnames) {
            def hostname = "hostname".execute().text.trim()
            params.hostnames.each { prof, hnames ->
                hnames.each { hname ->
                    if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                        log.info "=${colors.yellow}====================================================${colors.reset}=\n" +
                                  "${colors.yellow}WARN: You are running with `-profile $workflow.profile`\n" +
                                  "      but your machine hostname is ${colors.white}'$hostname'${colors.reset}.\n" +
                                  "      ${colors.yellow_bold}Please use `-profile $prof${colors.reset}`\n" +
                                  "=${colors.yellow}====================================================${colors.reset}="
                    }
                }
            }
        }
    }

    // Show a big warning message if we're not running MACS
    static void macs2_warn(params, log) {
        if (!params.macs_gsize) {
            def warnstring = params.genome ? "supported for '${params.genome}'" : 'supplied'
            log.warn "=================================================================\n" +
                     "  WARNING! MACS genome size parameter not $warnstring.\n" +
                     "  Peak calling, annotation and differential analysis will be skipped.\n" +
                     "  Please specify value for '--macs_gsize' to run these steps.\n" +
                     "======================================================================="
        }
    }
    
    // Exit pipeline if incorrect --genome key provided
    static void genome_exists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "============================================================================="
            System.exit(0)
        }
    }
    
     // Get attribute from genome config file e.g. fasta
    static String get_genome_attribute(params, attribute) {
        def val = ''
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }else{
            if(params.containsKey(attribute)){
                val = params[attribute]
            }
        }
        return val
    }
}
