// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

process PICARD_MERGESAMFILES {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container (params.universalContainer? "${process.container}":"quay.io/biocontainers/picard:2.23.2--0")
    //container "https://depot.galaxyproject.org/singularity/picard:2.23.2--0"

    conda (params.conda ? "${params.conda_softwares.picard}" : null)

    input:
    tuple val(meta), path(bams)
    val options

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    def software  = getSoftwareName(task.process)
    def ioptions  = initOptions(options)
    def prefix    = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    def bam_files = bams.sort()
    def avail_mem = 3
    if (!task.memory) {
        log.info '[Picard MergeSamFiles] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    if (bam_files.size() > 1) {
        """
        picard \\
            -Xmx${avail_mem}g \\
            MergeSamFiles \\
            $ioptions.args \\
            ${'INPUT='+bam_files.join(' INPUT=')} \\
            OUTPUT=${prefix}.bam
        echo \$(picard MergeSamFiles --version 2>&1) | awk -F' ' '{print \$NF}' > ${software}.version.txt
        """
    } else {
        """
        ln -s ${bam_files[0]} ${prefix}.bam
        echo \$(picard MergeSamFiles --version 2>&1) | awk -F' ' '{print \$NF}' > ${software}.version.txt
        """
    }
}
