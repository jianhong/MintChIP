// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]

/*
 * Demultiplex by Je
 */
process JO_TRIM {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
    
    conda (params.conda ? "${params.conda_softwares.trimmomatic}" : null)
    
    input:
    tuple val(meta), path(readsL), path(readsR)

    output:
    tuple val(meta), path("${meta.id}.R1.paired.fq.gz"), path("${meta.id}.R2.paired.fq.gz"), emit: paired
    tuple val(meta), path("${meta.id}.R1.unpaired.fq.gz"), path("${meta.id}.R2.unpaired.fq.gz"), emit: unpaired
    path "*.version.txt", emit: version
    
    script:
    def ioptions      = initOptions(params.options)
    def args          = ioptions.args2
    if(! ioptions.args instanceof String){
        if(ioptions.args[meta.id]){
            args += " " + ioptions.args[meta.id]
        }
    }
    """
    trimmomatic PE -threads ${task.cpus} \\
                   ${readsL} ${readsR} \\
                   ${meta.id}.R1.paired.fq.gz ${meta.id}.R1.unpaired.fq.gz \\
                   ${meta.id}.R2.paired.fq.gz ${meta.id}.R2.unpaired.fq.gz \\
                   ${args}
    trimmomatic -version 2>&1 > trimmomatic.version.txt
    """
}