/*
 * Read QC and trimming
 */

include { JO_DEMULTIPLEX    } from '../process/demultiplex/je'
include { FASTQC            } from '../../nf-core/software/fastqc/main'
include { FASTQC as TRIMQC  } from '../../nf-core/software/fastqc/main'
include { JO_TRIMMOMATIC    } from '../process/trimmomatic/trimmomatic'

workflow JO_FASTQC_DEMULTIPLEX_TRIMMOMATIC {
    take:
    ch_reads           // channel: [ val(meta), [ reads ] ]
    ch_barcode         // channel: [ path(barcode) ]
    skip_fastqc        // boolean: true/false
    skip_trimming      // boolean: true/false
    fastqc_options     //     map: options for FastQC module
    je_options         //     map: options for Je suite demultiplex module
    trimgalore_options //     map: options for TrimGalore! module

    main:
    fastqc_html = Channel.empty()
    fastqc_zip = Channel.empty()
    fastqc_version = Channel.empty()
    if (!skip_fastqc) {
        FASTQC(ch_reads, fastqc_options).html.set { fastqc_html }
        fastqc_zip = FASTQC.out.zip
        fastqc_version = FASTQC.out.version
    }

    je_multiplexer_stats = Channel.empty()
    je_version           = Channel.empty()
    JO_DEMULTIPLEX(
        ch_reads,
        ch_barcode,
        je_options
    )
    je_multiplexer_stats = JO_DEMULTIPLEX.out.jemultiplexerStats
    je_version           = JO_DEMULTIPLEX.out.version
    
    //JO_DEMULTIPLEX.out.fastq.view()
    
    ch_trim_reads_p1 = JO_DEMULTIPLEX.out.fastq
                                  .map{it[1]}
                                  .flatten()
                                  .map{
                                      it ->
                                             entry = it.simpleName.toString()
                                             [it.getParent().getName().toString(), entry.replaceAll("_[ACGTN]+_.*\$", ""), entry.replaceAll("^.*_([ACGTN]+)_.*\$", "\$1"), it]
                                  }
    //ch_trim_reads_p1.view()
    ch_trim_reads_p2 = JO_DEMULTIPLEX.out.fastq
                                  .map{it[2]}
                                  .flatten()
                                  .map{
                                      it ->
                                             def entry = it.simpleName.toString()
                                             [it.getParent().getName().toString(), entry.replaceAll("_[ACGTN]+_.*\$", ""), entry.replaceAll("^.*_([ACGTN]+)_.*\$", "\$1"), it]
                                  }
    //ch_trim_reads_p2.view()
    ch_trim_fq = ch_trim_reads_p1.join(ch_trim_reads_p2, by:[0, 1, 2])
                                 .filter { !(it[3].simpleName.toString() ==~ /_unassigned_/) && !(it[4].simpleName.toString() ==~ /_unassigned_/) }
    //ch_trim_fq.view()
    
    ch_trim_reads = JO_DEMULTIPLEX.out.fastq.map{[it[0].id, it[0]]}
                            .cross(ch_trim_fq)
                            .map{
                                origin_group, split_fqs ->
                                def (group, meta) = origin_group
                                def (group1, name, barcode, fq_1, fq_2) = split_fqs
                                [meta, barcode, name, [fq_1, fq_2]]
                            }
    //ch_trim_reads.view()
    
    unpaired_reads = Channel.empty()
    trim_html = Channel.empty()
    trim_zip = Channel.empty()
    trim_log = Channel.empty()
    trimgalore_version = Channel.empty()
    if (!skip_trimming) {
        JO_TRIMMOMATIC(ch_trim_reads, trimgalore_options).reads.set { ch_trim_reads }
        unpaired_reads = JO_TRIMMOMATIC.out.unpaired_reads
        trim_log = JO_TRIMMOMATIC.out.log
        if (!skip_fastqc){
            TRIMQC(ch_trim_reads, fastqc_options).html.set { trim_html }
            trim_zip = TRIMQC.out.zip
        }
    }

    emit:
    fastqc_html           // channel: [ val(meta), [ html ] ]
    fastqc_zip            // channel: [ val(meta), [ zip ] ]
    fastqc_version        //    path: *.version.txt
    
    je_multiplexer_stats  // channel: [val(meta), [stats.txt] ]
    je_version            //    path: *.version.txt

    reads = ch_trim_reads // channel: [ val(meta), [ reads ] ]
    unpaired_reads        // channel: [ val(meta), [ unpaired_reads ] ]
    trim_html             // channel: [ val(meta), [ html ] ]
    trim_zip              // channel: [ val(meta), [ zip ] ]
    trim_log              // channel: [ val(meta), [ txt ] ]
    trimgalore_version    //    path: *.version.txt
}
