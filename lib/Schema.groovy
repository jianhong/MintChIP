/*
 * This file holds several functions used to perform JSON parameter validation, help and summary rendering for the nf-core pipeline template.
 */

import groovy.json.JsonSlurper

class Schema {
    /*
     * This method tries to read a JSON params file
     */
    private static LinkedHashMap params_get(String path) {
        def params_map = new LinkedHashMap()
        try {
            params_map = params_try(path)
        } catch (Exception e) {
            println "Could not read parameters settings from JSON. $e"
            params_map = new LinkedHashMap()
        }
        return params_map
    }

    /*
    Method to actually read in JSON file using Groovy.
    Group (as Key), values are all parameters
        - Parameter1 as Key, Description as Value
        - Parameter2 as Key, Description as Value
        ....
    Group
        -
    */
    private static LinkedHashMap params_try(String path) throws Exception {

        def json = new File(path).text
        def Map json_params = (Map) new JsonSlurper().parseText(json).get('properties')

        /* Tree looks like this in nf-core schema
        *  properties <- this is what the first get('properties') gets us
             group 1
               properties
               description
             group 2
               properties
               description
             group 3
               properties
               description
        */
        def params_map = new LinkedHashMap()
        json_params.each { key, val ->
            def Map group = json_params."$key".properties // Gets the property object of the group
            def sub_params = new LinkedHashMap()
            group.each { innerkey, value ->
                sub_params.put("$innerkey", [ "$value.type", "$value.description" ])
            }
            params_map.put("$key", sub_params)
        }
        return params_map
    }

    private static Integer params_max_chars(params_map) {
        Integer max_chars = 0
        for (group in params_map.keySet()) {
            def params = params_map.get(group)  // This gets the parameters of that particular group
            for (par in params.keySet()) {
                if (par.size() > max_chars) {
                    max_chars = par.size()
                }
            }
        }
        return max_chars
    }

    private static String params_beautify(params_map) {
        String output = ""
        def max_chars = params_max_chars(params_map) + 1
        for (group in params_map.keySet()) {
            output += group + "\n"
            def params = params_map.get(group)  // This gets the parameters of that particular group
            for (par in params.keySet()) {
                def type = "[" + params.get(par)[0] + "]"
                def description = params.get(par)[1]
                output+= "    \u001B[1m" +  par.padRight(max_chars) + "\u001B[1m" + type.padRight(10) + description + "\n"
            }
            output += "\n"
        }
        return output
    }

    private static String params_help(path, command) {
          String output = "Typical pipeline command:\n\n"
          output += "    ${command}\n\n"
          output += params_beautify(params_get(path))
    }

    private static LinkedHashMap params_summary(workflow, params, run_name) {
        def Map summary = [:]
        if (workflow.revision) summary['Pipeline Release'] = workflow.revision
        summary['Run Name']         = run_name ?: workflow.runName
        summary['Design File']            = params.input
        summary['Genome']                 = params.genome ?: 'Not supplied'
        summary['Fasta File']             = params.fasta
        summary['GTF File']               = params.gtf
        if (params.gene_bed)              summary['Gene BED File'] = params.gene_bed
        if (params.bwa_index)             summary['BWA Index'] = params.bwa_index
        if (params.blacklist)             summary['Blacklist BED'] = params.blacklist
        if (params.bwa_min_score)         summary['BWA Min Score'] = params.bwa_min_score
        summary['MACS2 Genome Size']      = params.macs_gsize ?: 'Not supplied'
        summary['Min Consensus Reps']     = params.min_reps_consensus
        if (params.macs_fdr)              summary['MACS2 FDR'] = params.macs_fdr
        if (params.macs_pvalue)           summary['MACS2 P-value'] = params.macs_pvalue
        if (params.skip_trimming) {
            summary['Trimming Step']      = 'Skipped'
        } else {
            summary['trimmomatic']        = "$params.modules.jo_trimmomatic.args"
        }
        if (params.seq_center)            summary['Sequencing Center'] = params.seq_center
        summary['Fragment Size']          = "$params.fragment_size bp"
        summary['Fingerprint Bins']       = params.fingerprint_bins
        if (params.keep_dups)             summary['Keep Duplicates'] = 'Yes'
        if (params.keep_multi_map)        summary['Keep Multi-mapped'] = 'Yes'
        summary['Save Genome Index']      = params.save_reference ? 'Yes' : 'No'
        if (params.save_trimmed)          summary['Save Trimmed'] = 'Yes'
        if (params.save_align_intermeds)  summary['Save Intermeds'] =  'Yes'
        if (params.save_macs_pileup)      summary['Save MACS2 Pileup'] = 'Yes'
        if (params.skip_peak_qc)          summary['Skip MACS2 Peak QC'] = 'Yes'
        if (params.skip_peak_annotation)  summary['Skip Peak Annotation'] = 'Yes'
        if (params.skip_consensus_peaks)  summary['Skip Consensus Peaks'] = 'Yes'
        if (params.deseq2_vst)            summary['Use DESeq2 vst Transform'] = 'Yes'
        if (params.skip_diff_analysis)    summary['Skip Differential Analysis'] = 'Yes'
        if (params.skip_fastqc)           summary['Skip FastQC'] = 'Yes'
        if (params.skip_picard_metrics)   summary['Skip Picard Metrics'] = 'Yes'
        if (params.skip_preseq)           summary['Skip Preseq'] = 'Yes'
        if (params.skip_plot_profile)     summary['Skip plotProfile'] = 'Yes'
        if (params.skip_plot_fingerprint) summary['Skip plotFingerprint'] = 'Yes'
        if (params.skip_spp)              summary['Skip spp'] = 'Yes'
        if (params.skip_igv)              summary['Skip IGV'] = 'Yes'
        if (params.skip_multiqc)          summary['Skip MultiQC'] = 'Yes'
        summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
        if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
        summary['Output dir']       = params.outdir
        summary['Launch dir']       = workflow.launchDir
        summary['Working dir']      = workflow.workDir
        summary['Script dir']       = workflow.projectDir
        summary['User']             = workflow.userName
        if (workflow.profile.contains('awsbatch')) {
            summary['AWS Region']   = params.awsregion
            summary['AWS Queue']    = params.awsqueue
            summary['AWS CLI']      = params.awscli
        }
        summary['Config Profile'] = workflow.profile
        if (params.config_profile_description) summary['Config Profile Descr']   = params.config_profile_description
        if (params.config_profile_contact)     summary['Config Profile Contact'] = params.config_profile_contact
        if (params.config_profile_url)         summary['Config Profile URL']     = params.config_profile_url
        summary['Config Files'] = workflow.configFiles.join(', ')
        if (params.email || params.email_on_fail) {
            summary['E-mail Address']    = params.email
            summary['E-mail on failure'] = params.email_on_fail
            summary['MultiQC maxsize']   = params.max_multiqc_email_size
        }
        return summary
    }

    static String params_mqc_summary(summary) {
        String yaml_file_text  = """
        id: 'nf-core-MintChIP-summary'
        description: " - this information is collected when the pipeline is started."
        section_name: 'nf-core/MintChIP Workflow Summary'
        section_href: 'https://github.com/jianhong/MintChIP'
        plot_type: 'html'
        data: |
            <dl class=\"dl-horizontal\">
            ${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
            </dl>
        """.stripIndent()

        return yaml_file_text
    }
    
    /*
     * Groovy Map summarising parameters/workflow options used by the pipeline
     */
    private static LinkedHashMap params_summary_map(workflow, params, json_schema) {
        // Get a selection of core Nextflow workflow options
        def Map workflow_summary = [:]        
        if (workflow.revision) {
            workflow_summary['revision'] = workflow.revision
        }
        workflow_summary['runName']      = workflow.runName
        if (workflow.containerEngine) {
            workflow_summary['containerEngine'] = "$workflow.containerEngine"
        }
        if (workflow.container) {
            workflow_summary['container']       = "$workflow.container"
        }
        workflow_summary['launchDir']    = workflow.launchDir
        workflow_summary['workDir']      = workflow.workDir
        workflow_summary['projectDir']   = workflow.projectDir
        workflow_summary['userName']     = workflow.userName
        workflow_summary['profile']      = workflow.profile
        workflow_summary['input']        = params.input
        workflow_summary['barcode']      = params.barcode
        workflow_summary['configFiles']  = workflow.configFiles.join(', ')
        
        // Get pipeline parameters defined in JSON Schema
        def Map params_summary = [:]
        def blacklist  = ['hostnames']
        def params_map = params_get(json_schema)
        for (group in params_map.keySet()) {
            def sub_params = new LinkedHashMap()
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            for (param in group_params.keySet()) {
                if (params.containsKey(param) && !blacklist.contains(param)) {
                    def params_value = params.get(param)
                    def schema_value = group_params.get(param).default
                    def param_type   = group_params.get(param).type
                    if (schema_value == null) {
                        if (param_type == 'boolean') {
                            schema_value = false
                        }
                        if (param_type == 'string') {
                            schema_value = ''
                        }
                        if (param_type == 'integer') {
                            schema_value = 0
                        }
                    } else {
                        if (param_type == 'string') {
                            if (schema_value.contains('$projectDir') || schema_value.contains('${projectDir}')) {
                                def sub_string = schema_value.replace('\$projectDir','')
                                sub_string     = sub_string.replace('\${projectDir}','')
                                if (params_value.contains(sub_string)) {
                                    schema_value = params_value
                                }
                            }
                            if (schema_value.contains('$params.outdir') || schema_value.contains('${params.outdir}')) {
                                def sub_string = schema_value.replace('\$params.outdir','')
                                sub_string     = sub_string.replace('\${params.outdir}','')
                                if ("${params.outdir}${sub_string}" == params_value) {
                                    schema_value = params_value
                                }
                            }
                        }
                    }

                    if (params_value != schema_value) {
                        sub_params.put("$param", params_value)
                    }
                }
            }
            params_summary.put(group, sub_params)
            
        }
        return [ 'Core Nextflow options' : workflow_summary ] << params_summary
    }
    
    /*
     * Beautify parameters for summary and return as string
     */
    private static String params_summary_log(workflow, params, json_schema) {
        String output  = Headers.nf_core(workflow, params.monochrome_logs) + "\n"
        def params_map = params_summary_map(workflow, params, json_schema)
        def max_chars  = params_max_chars(params_map)
        for (group in params_map.keySet()) {
            def group_params = params_map.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                output += group + "\n"
                for (param in group_params.keySet()) {
                    output += "    \u001B[1m" +  param.padRight(max_chars) + ": \u001B[1m" + group_params.get(param) + "\n"
                }
                output += "\n"
            }
        }
        output += Headers.dashed_line(params.monochrome_logs)
        return output
    }
}
