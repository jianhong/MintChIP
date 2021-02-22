// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]

/*
 * Demultiplex by Je
 */
process JO_CHIPSEQ {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    
    input:
    val meta

    output:
    tuple val(meta), path("results/*")
    
    script:
    def ioptions      = initOptions(params.options)
    def args          = ioptions.args
    def out = "${workflow.workDir}/design.csv"
    def FILE_HEADER = meta[0].keySet().join(',')+'\n'
    new File(out).withWriter {
        target ->
            target << FILE_HEADER
            meta.each {
                target << it.values().join(',')+'\n'
            }
    }
    def cmd = workflow.commandLine.replaceFirst(/run .*MintChIP/, "run jianhong/chipseq") +
              " --input \"$out\""
    """
    #nextflow run 
    mkdir results
    cp ${workflow.workDir}/design.csv results/
    echo $cmd > results/cmd.txt
    $cmd
    """
}