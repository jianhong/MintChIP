// Import generic module functions
include { saveFiles } from './functions'

/*
 * Prepare genome intervals for filtering by removing regions in blacklist file
 */
process MAKE_GENOME_FILTER {
    tag "$sizes"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:"genome", publish_id:'') }

    container (params.universalContainer? "${params.container}":"quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0")
    //container "https://depot.galaxyproject.org/singularity/bedtools:2.29.2--hc088bd4_0"

    conda (params.conda ? "${params.conda_softwares.bedtools}" : null)

    input:
    path sizes
    path blacklist
    val options

    output:
    path '*.bed', emit: bed
    path "*.version.txt", emit: version

    script:
    def software = 'bedtools'
    def file_out = "${sizes.simpleName}.include_regions.bed"
    if (params.blacklist) {
        """
        sortBed -i $blacklist -g $sizes | complementBed -i stdin -g $sizes > $file_out
        bedtools --version | sed -e "s/bedtools v//g" > ${software}.version.txt
        """
    } else {
        """
        awk '{print \$1, '0' , \$2}' OFS='\t' $sizes > $file_out
        bedtools --version | sed -e "s/bedtools v//g" > ${software}.version.txt
        """
    }
}
