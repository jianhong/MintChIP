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

String pwd = System.getProperty("user.dir")

// input files
def inputs = []
if (params.inputs){
  new File(params.inputs).eachLine{
    line -> 
      inputs << line.split("\\s")
  }
}
def withInput = []
if (params.inputs){
  new File(params.inputs).eachLine{
    line -> 
      withInput << line.split("\\s")[0]
  }
}

//additional parameters for trimmomatic
Map<String,String> trimmomatic = new HashMap<String,String>()
if (params.trim){
  new File(params.trim).eachLine{
    line -> 
      def p=line.split("\\t")
      trimmomatic.put(p[0], p[1])
  }
}
//additional parameters for findPeak
Map<String,String> findPeaks = new HashMap<String,String>()
if (params.findPeaks){
  new File(params.findPeaks).eachLine{
    line -> 
      def p=line.split("\\t")
      findPeaks.put(p[0], p[1])
  }
}

log.info """\
R E A D   P R E P R O C E S S I N G   P I P E L I N E
=====================================================
dataDir : ${params.dataDir}
barcode : ${params.barcode}
outdir  : ${params.outdir}
input   : ${inputs}
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
file("${runID}/${runID}%jemultiplexer_out_stats.txt") into jemultiplexerStats
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
mv ${runID}/\$i ${runID}/${runID}%\$i
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

storeDir "${params.outdir}/trimmed/${GROUP}"

output:
set GROUP, read0, file("${read0}.R1.paired.fq.gz"), file("${read0}.R2.paired.fq.gz") into readPairsTrimmed
set GROUP, read0, file("${read0}.R1.unpaired.fq.gz"), file("${read0}.R2.unpaired.fq.gz") into readUnpairesTrimmed

script:
GROUP = readL.simpleName
GROUP = GROUP.tokenize('%')[0]
read0 = readL.simpleName
read0 = read0.replace("${GROUP}%", "")
read0 = read0.replaceAll("_[ACGTN]+_1\$", "")
additional = trimmomatic.get(read0)
if(!additional){
  additional = ""
}else{
println "${read0} additional parameter for trimmomatic: ${additional}"
}
"""
# call Trimmomatic
echo "PE: ${readL} ${readR}"
${params.trimmomatic.path} PE -threads ${params.trimmomatic.cpus} \\
${readL} ${readR} \\
${read0}.R1.paired.fq.gz ${read0}.R1.unpaired.fq.gz ${read0}.R2.paired.fq.gz ${read0}.R2.unpaired.fq.gz \\
${params.trimmomatic.options} \\
${additional}
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
set runID, sampleName, file(readL), file(readR) from readPairsFastQCTrimmed

storeDir "${params.outdir}/reports/fastqc-trimmed/${runID}"

output:
set runID, sampleName, file('*.zip'), file('*.html') into fastqcOutputTrimmed

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

process runBWA {
  cpus params.bwa.cpus
  
  input:
    set runID, sampleName, file(readL), file(readR) from readPairsRunMapping
  
  storeDir "${params.outdir}/bwa/${runID}"
  
  output:
    set sampleName, runID, file("${runID}%${readL.simpleName}.bam"), file("${runID}%${readL.simpleName}.bam.bai") into mappedFiles
  
  script:
    """
    ${params.bwa.path} mem -M -t ${params.bwa.cpus} ${params.bwa.index} ${readL} ${readR} > ${readL.simpleName}.sam
    ${params.samtools.path} sort -@ ${params.bwa.cpus} -n ${readL.simpleName}.sam -o ${readL.simpleName}.bam
    ${params.samtools.path} fixmate -@ ${params.bwa.cpus} -m ${readL.simpleName}.bam ${readL.simpleName}.fixout.bam
    ${params.samtools.path} sort ${readL.simpleName}.fixout.bam -@ ${params.bwa.cpus} -o ${readL.simpleName}.sort.bam
    ${params.samtools.path} markdup ${readL.simpleName}.sort.bam ${readL.simpleName}.rem.bam -@ ${params.bwa.cpus} -r
    mv ${readL.simpleName}.rem.bam ${runID}%${readL.simpleName}.bam
    ${params.samtools.path} index ${runID}%${readL.simpleName}.bam
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

// pair exp with inputs

Channel.from(inputs).into{inputCh; inputChcp}
mappedFiles.into{mappedFiles1; mappedFiles2; mappedFiles3; mappedFiles4}
//inputChcp.println()
//mappedFiles2.println()

inputChFile=inputChcp.map{it-> "${it[2]}%${it[3]}"}.flatten().cross(mappedFiles2.map{it -> ["${it[1]}%${it[0]}", it]})
  .map{it -> it[1]}
//inputChFile.println{it->"inputChFile: $it"}

mappedPair=mappedFiles1.map{it -> ["${it[1]}%${it[0]}", it]}.cross(inputCh.map{it-> ["${it[0]}%${it[1]}", "${it[2]}%${it[3]}"]})
  .map{it -> [it[1][1], it[0][1]]}

//mappedPair.println{it->"mappedPair: $it"}

// first item input, second item exp
forHomerWithInput = inputChFile.cross(mappedPair)
  .map{it -> [it[0][1], it[1][1]]}.map{it->it.flatten()}

//forHomerWithInput.println()

// without input
forHomerWithoutInput = mappedFiles3.filter{it->!(it[0] in withInput)}
//forHomerWithoutInput.println()

process runHOMERwithoutInput {
  cpus params.homer.cpus
  
  input:
    set expSampleName, expGroup, file(expbam), file(expbamIndex) from forHomerWithoutInput
  
  storeDir "${params.outdir}/homer/${expGroup}/${expSampleName}"
  
  output:
    file("${expSampleName}*") into homerFiles0
  
  script:
    additional = findPeaks.get(expSampleName)
    if(!additional){
      additional = ""
    }else{
      println "${expSampleName} additional parameter for findPeaks: ${additional}"
    }
    
    """
    ${params.homer.makeTagDirectory} ${expSampleName}_Tagdir ${expbam} -sspe
    ${params.homer.makeUCSCfile} ${expSampleName}_Tagdir -name ${expSampleName}_Chr1-10 \\
    -skipChr chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \\
    -o ${expSampleName}_Chr1-10.bedgraph -color 0,0,204 -norm 1e7
    ${params.homer.makeUCSCfile} ${expSampleName}_Tagdir -name ${expSampleName}_Chr11 \\
    -skipChr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \\
    -o ${expSampleName}_Chr11-.bedgraph -norm 1e7
    
    ${params.homer.findPeaks} ${expSampleName}_Tagdir ${params.homer.findPeaksOptions} \\
    -o ${expSampleName}_Calledpeaks.txt \\
    ${additional}

    ${params.homer.annotatePeaks} ${expSampleName}_Calledpeaks.txt mm10 > ${expSampleName}_Annotatedlist.txt
    ${params.homer.pos2bed} ${expSampleName}_Calledpeaks.txt -o ${expSampleName}.bed -track ${expSampleName}
    """
}


process runHOMERwithInput {
  cpus params.homer.cpus
  
  input:
    set inputSampleName, inputGroup, file(inputbam), file(inputbamIndex), expSampleName, expGroup, file(expbam), file(expbamIndex) from forHomerWithInput
  
  storeDir "${params.outdir}/homer/${expGroup}/${expSampleName}"
  
  output:
    file("${expSampleName}*") into homerFiles1
    file("${inputSampleName}*") into homerFiles2
    file("${inputGroup}%${inputSampleName}_*") into homerFiles3
    file("${expGroup}%${expSampleName}_*") into homerFiles4
  
  script:
    additional = findPeaks.get(expSampleName)
    if(!additional){
      additional = ""
    }else{
      println "${expSampleName} additional parameter for findPeaks: ${additional}"
    }
    """
    ${params.homer.makeTagDirectory} ${inputGroup}%${inputSampleName}_Tagdir ${inputbam} -sspe
    ${params.homer.makeUCSCfile} ${inputGroup}%${inputSampleName}_Tagdir -name ${inputSampleName}_Chr1-10 \\
    -skipChr chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \\
    -o ${inputSampleName}_Chr1-10.bedgraph -color 0,0,204 -norm 1e7
    ${params.homer.makeUCSCfile} ${inputGroup}%${inputSampleName}_Tagdir -name ${inputSampleName}_Chr11 \\
    -skipChr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \\
    -o ${inputSampleName}_Chr11-.bedgraph -norm 1e7
    
    ${params.homer.makeTagDirectory} ${expGroup}%${expSampleName}_Tagdir ${expbam} -sspe
    ${params.homer.makeUCSCfile} ${expGroup}%${expSampleName}_Tagdir -name ${expSampleName}_Chr1-10 \\
    -skipChr chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY \\
    -o ${expSampleName}_Chr1-10.bedgraph -color 0,0,204 -norm 1e7
    ${params.homer.makeUCSCfile} ${expGroup}%${expSampleName}_Tagdir -name ${expSampleName}_Chr11 \\
    -skipChr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \\
    -o ${expSampleName}_Chr11-.bedgraph -norm 1e7
    
    ${params.homer.findPeaks} ${expGroup}%${expSampleName}_Tagdir ${params.homer.findPeaksOptions} \\
    -i ${inputGroup}%${inputSampleName}_Tagdir \\
    -o ${expSampleName}_Calledpeaks.txt \\
    ${additional}

    ${params.homer.annotatePeaks} ${expSampleName}_Calledpeaks.txt mm10 > ${expSampleName}_Annotatedlist.txt
    ${params.homer.pos2bed} ${expSampleName}_Calledpeaks.txt -o ${expSampleName}.bed -track ${expSampleName}
    """

}
