/*
 * unzip files
 */
process GUNZIP_FILES {
    tag "${gz.baseName}"
    label 'process_low'

    input:
    path gz

    output:
    path "${gz.baseName}" 

    script: // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    """
    gunzip -k --verbose --stdout --force ${gz} > ${gz.baseName}
    """
}
