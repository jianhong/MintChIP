#!/usr/bin/env nextflow

//  main.nf
//
//
//  Created by jianhong ou on 11/7/18.
//  modified from holtgrewe/ngs_pipelines


import ParamsHelper


// ---------------------------------------------------------------------------
// Read preprocessing: demultiplex, and
// ---------------------------------------------------------------------------

if (params.verbose)
echo true

ParamsHelper.checkNonEmptyParam(params.dataDir, "dataDir");
ParamsHelper.checkNonEmptyParam(params.barcode, "barcode");

log.info """\
R E A D   P R E P R O C E S S I N G   P I P E L I N E
=====================================================
dataDir : ${params.dataDir}
barcode : ${params.barcode}
outdir  : ${params.outdir}
"""
.stripIndent()


// Open channel for left and right files and merge it into triples, the
// first entry is the LCS of the file names that can be used as a read
// pair identifier.
Channel
  .fromFilePairs("${params.dataDir}/*_{R,}{1,2}.f{ast,}q.gz")
  .set {readPairs}

// barcode file.
barcode = file(params.barcode)

// Duplicate the read pairs into one queue for runFastQCOriginal
// and runTrimming.
readPairs.into{readPairsFastQCOriginal; readPairsDemultiplex}

// --------------------------------------------------------------------------
// Step 1a) Run FastQC
//
// - yields report
// --------------------------------------------------------------------------

process runFastQCOriginal {
cpus params.fastqc.cpus

input:
set runID, file(reads) from readPairsFastQCOriginal

storeDir "${params.outdir}/reports/fastqc-original/${runID}"

output:
set file('*.zip'), file('*.html') into fastqcOutputOriginal

script:
"""
${params.fastqc.path} -t ${params.fastqc.cpus} -o . ${reads}
"""
}

// --------------------------------------------------------------------------
// Step 1b) Run demultiplex
//
// - yields demultiplexed read, used as downstream input
// --------------------------------------------------------------------------

process runDemultiplex {
cpus params.je.cpus

input:
set runID, file(reads) from readPairsDemultiplex
file barcode

storeDir "${params.outdir}/demultiplexed"

output:
file("${runID}/${runID}__jemultiplexer_out_stats.txt") into jemultiplexerStats
file("${runID}/*_1.txt.gz") into readPairsDemultiedTrimmingL
file("${runID}/*_2.txt.gz") into readPairsDemultiedTrimmingR

script:
"""
# call je
${params.je.path} demultiplex \\
F1=${reads[0]} \\
F2=${reads[1]} \\
BF=$barcode \\
O=${runID} \\
BPOS=READ_2 \\
ADD=true \\
Q=$params.je.MIN_BASE_QUALITY
for i in ${runID}/* 
do
i=`basename \$i`
mv ${runID}/\$i ${runID}/${runID}__\$i
done
"""
}

readPairsDemultiedTrimmingL.flatten()
  .merge(readPairsDemultiedTrimmingR.flatten())
  .into{readPairsDemultiedTrimPF; readPairsDemultiedTrimFlat}

// --------------------------------------------------------------------------
// Step 2) Run adapter trimming
//
// - yields trimmed read, used as downstream input
// --------------------------------------------------------------------------

process runTrimming {
cpus params.trimmomatic.cpus

input:
set file(readL), file(readR) from readPairsDemultiedTrimFlat

storeDir "${params.outdir}/trimmed"

output:
set file("${readL.simpleName}.R1.paired.fq.gz"), file("${readL.simpleName}.R2.paired.fq.gz") into readPairsTrimmed
set file("${readL.simpleName}.R1.unpaired.fq.gz"), file("${readL.simpleName}.R2.unpaired.fq.gz") into readUnpairesTrimmed

script:
"""
# call Trimmomatic
echo "PE: ${readL} ${readR}"
${params.trimmomatic.path} PE -phred33 -threads ${params.trimmomatic.cpus} \\
${readL} ${readR} \\
${readL.simpleName}.R1.paired.fq.gz ${readL.simpleName}.R1.unpaired.fq.gz ${readL.simpleName}.R2.paired.fq.gz ${readL.simpleName}.R2.unpaired.fq.gz \\
ILLUMINACLIP:${params.trimmomatic.adapters}:2:30:10 \\
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:18
"""
}


readPairsTrimmed.into{readPairsFastQCTrimmed; readPairsRunMapping}

// --------------------------------------------------------------------------
// Step 3a) Run FastQC on trimmed
//
// - yields report
// --------------------------------------------------------------------------

process runFastQCTrimmed {
cpus params.fastqc.cpus

input:
set file(readL), file(readR) from readPairsFastQCTrimmed

storeDir "${params.outdir}/reports/fastqc-trimmed"

output:
set file('*.zip'), file('*.html') into fastqcOutputTrimmed

script:
"""
${params.fastqc.path} -t ${params.fastqc.cpus} -o . ${readL} ${readR}
"""
}

// --------------------------------------------------------------------------
// Step 2b) Align reads using BWA-MEM
//
// - align reads
// - sort
// - mark duplicates
// - yields alignment for downstream processing
// --------------------------------------------------------------------------

// The alignments are written to the temporary files alignment.bam. These
// BAM files are already sorted.

// bwa mem -M -t 24 /work/yk170/Reference/mm10/Sequence/BWAIndex/genome.fa R1.paired.fq.gz R2.paired.fq.gz > bwa.sam
process runBWA {
  cpus params.bwa.cpus
  
  input:
    set file(readL), file(readR) from readPairsRunMapping
  
  storeDir "${params.outdir}/bwa"
  
  output:
    set file("${readL.simpleName}.bam"), file("${readL.simpleName}.bam.bai") into mappedFiles
  
  script:
    """
    ${params.bwa.path} mem -M -t ${params.bwa.cpus} ${params.bwa.index} ${readL} ${readR} > ${readL.simpleName}.sam
    ${params.samtools.path} sort -@ ${params.bwa.cpus} -n ${readL.simpleName}.sam -o ${readL.simpleName}.bam
    ${params.samtools.path} fixmate -@ ${params.bwa.cpus} -m ${readL.simpleName}.bam ${readL.simpleName}.fixout.bam
    ${params.samtools.path} sort ${readL.simpleName}.fixout.bam -@ ${params.bwa.cpus} -o ${readL.simpleName}.sort.bam
    ${params.samtools.path} markdup ${readL.simpleName}.sort.bam ${readL.simpleName}.rem.bam -@ ${params.bwa.cpus} -r
    mv ${readL.simpleName}.rem.bam ${readL.simpleName}.bam
    ${params.samtools.path} index ${readL.simpleName}.bam
    """
}

// --------------------------------------------------------------------------
// Step 2c) run Homer
//
// - makeTagDirectory
// - findPeaks
// - annotatePeaks
// - pos2bed
// --------------------------------------------------------------------------

process runHOMER {
  cpus params.homer.cpus
  
  input:
    set file(bam), file(bamIndex) from mappedFiles
  
  storeDir "${params.outdir}/homer"
  
  output:
    set file("${bam.simpleName}*") into homerFiles
  
  script:
    """
    ${params.homer.makeTagDirectory} ${bam.simpleName}_Tagdir ${bam} -sspe
    ${params.homer.makeUCSCfile} ${bam.simpleName}_Tagdir -name ${bam.simpleName}_Chr1-10 \\
    -skipChr chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \\
    -o ${bam.simpleName}_Chr1-10.bedgraph -color 0,0,204 -norm 1e7
    ${params.homer.makeUCSCfile} ${bam.simpleName}_Tagdir -name ${bam.simpleName}_Chr11 \\
    -skipChr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \\
    -o ${bam.simpleName}_Chr11-.bedgraph -norm 1e7
    ${params.homer.findPeaks} ${bam.simpleName}_Tagdir -region -size 1000 -minDist 2000 -C 0 -L 50 \\
    -o ${bam.simpleName}_Calledpeaks.txt
    ${params.homer.annotatePeaks} ${bam.simpleName}_Calledpeaks.txt mm10 > ${bam.simpleName}_Annotatedlist.txt
    ${params.homer.pos2bed} ${bam.simpleName}_Calledpeaks.txt -o ${bam.simpleName}.bed -track ${bam.simpleName}
    """
}
