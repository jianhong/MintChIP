// Import generic module functions
include { initOptions; saveFiles } from '../functions'

/*
 * Demultiplex by Je
 */
process JO_DEMULTIPLEX {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.tokenize(':')[-1].tokenize('_')[1].toLowerCase(), publish_id:meta.id) }
    
    container (params.universalContainer? "${process.container}":"${process.container}")
    conda (params.conda ? "${params.conda_softwares.je}" : null)
    
    input:
    tuple val(meta), path(reads)
    path barcode
    val options

    output:
    tuple val(meta), path("${meta.id}/jemultiplexer_out_stats.txt"), emit: jemultiplexerStats
    tuple val(meta), path("${meta.id}/*_1.txt.gz"), path("${meta.id}/*_2.txt.gz"), emit: fastq
    path "*.version.txt", emit: version
    
    script:
    def ioptions      = initOptions(options)
    """
    je demultiplex F1=${reads[0]} F2=${reads[1]} BF=${barcode} \\
                   O="${meta.id}" BPOS=READ_2 ADD=true $ioptions.args
    je -v 2>&1 | grep -v JAVA  > je.version.txt
    """
}