// Import generic module functions
include { initOptions; saveFiles } from '../functions'

process JO_TRIMMOMATIC {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:task.process.tokenize(':')[-1].tokenize('_')[1].toLowerCase(), publish_id:meta.id) }

    container "quay.io/biocontainers/trimmomatic:0.3.9--0"
    //container "https://depot.galaxyproject.org/singularity/trimmomatic:0.3.9--0"

    conda (params.conda ? "${params.conda_softwares.trimmomatic}" : null)

    input:
    tuple val(meta), path(reads)
    val options

    output:
    tuple val(meta), path("*.trimmed.fq.gz"), emit: reads
    tuple val(meta), path("*.unpaired.fq.gz"), emit: unpaired_reads
    tuple val(meta), path("*report.txt"), emit: log
    path "*.version.txt", emit: version

    script:
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) cores = (task.cpus as int) - 3
        if (cores < 1) cores = 1
        if (cores > 4) cores = 4
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def ioptions = initOptions(options)
    def prefix   = ioptions.suffix ? "${meta.id}${ioptions.suffix}" : "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        trimmomatic \\
            SE -threads $cores \\
            ${prefix}.fastq.gz \\
            ${prefix}.trimmed.fq.gz \\
            -trimlog trimmomatic.report.txt \\
            $ioptions.args
        trimmomatic -version  > trimmomatic.version.txt
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trimmomatic \\
            PE -threads $cores \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz \\
            ${prefix}_R1.trimmed.fq.gz \\
            ${prefix}_R1.unpaired.fq.gz \\
            ${prefix}_R2.trimmed.fq.gz \\
            ${prefix}_R2.unpaired.fq.gz \\
            -trimlog trimmomatic.report.txt \\
            $ioptions.args
        trimmomatic -version  > trimmomatic.version.txt
        """
    }
}
