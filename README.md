# ![qiubio-nf-core/Mint-ChIP](assets/mintchiplogo.png)

## Introduction

**qiubio/MintChIP** is a bioinformatics analysis pipeline used for Mint-ChIP-seq data based on [nfcore/chipseq](https://nf-co.re/chipseq).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

This pipeline will generate the UCSC genome browser track hub and metagene analysis
resuls in addition to original output of **nfcore/chipseq** pipeline.

The most important change for this pipeline from **nfcore/chipseq** pipeline is 
that I tried to improve the reproducibility of the pipeline depend on conda but 
not docker for following 2 reasons:

1. I can not use docker in our cluster. 

2. The memory required for the pipeline is too heavy for personal computer if using docker.

However, conda always throw errors when create environment even I use modules. 
I add conda_softwares section in module.conf setting to make the pipeline more flexible 
to figure out this issue.
I also changed the R/Biocondactor package installation methods from conda installation 
to BiocManager installation, which will be much slower than conda installation. 
The reason for that is because lots of package in conda is malfunction. By using
BiocManager to avoid the dependece issues.



## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Split the reads by [Je](https://gbcs.embl.de/portal/tiki-index.php?page=Je)
2. Adapter trimming ([`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic))
3. Alignment ([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
5. Merge alignments from multiple libraries of the same sample ([`picard`](https://broadinstitute.github.io/picard/))
    1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    2. Filtering to remove:
        * reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
        * reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that arent marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
        * reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); *paired-end only*)
        * reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
    3. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
    4. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    5. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
    6. Calculate genome-wide IP enrichment relative to control ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
    7. Calculate strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC ([`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools))
    8. Call broad/narrow peaks ([`MACS2`](https://github.com/taoliu/MACS))
    9. Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
    10. Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
    11. Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
    12. Differential binding analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/release/bioc/html/DESeq2.html))
6. Visualisation the tracks.
    1. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
    2. Create UCSC genome browser track hub for bigWig tracks [trackhub](https://daler.github.io/trackhub/quickstart.html).
7. Present QC for raw read, alignment, peak-calling and differential binding results ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))
     


## Installation by conda

```bash
conda update conda
nextflow pull jianhong/MintChiP -r dev
srun --mem 60G -c 2 nextflow run jianhong/MintChiP -profile test -r dev --conda
```

## Update

```bash
conda activate MintChiP
nextflow pull jianhong/MintChiP
```

## Remove

```bash
conda activate MintChiP
nextflow drop jianhong/MintChiP
conda deactivate
conda remove --name MintChiP --all
conda info --envs
```


## design table

| group | replicate | fastq_1 | fastq_2 | antibody | control | track_color | track_group |
|-------|-----------|---------|---------|----------|---------|-------------|-------------|
| WT | 1 | fastq/WT1.fastq.gz| | ANT1 | Input | #E69F00 | SAMPLE |
| WT | 2 | fastq/WT2.fastq.gz| | ANT1 | Input | #E69F00 | SAMPLE |
| KD | 1 | fastq/KD1.fastq.gz| | ANT1 | Input | #0000FF | SAMPLE |
| KD | 2 | fastq/KD2.fastq.gz| | ANT1 | Input | #0000FF | SAMPLE |
| Input | 1 | fastq/KD1.fastq.gz| |  |  | #000000 | SAMPLE |
| Input | 2 | fastq/KD2.fastq.gz| |  |  | #000000 | SAMPLE |

## metagene analysis

```
nextflow run jianhong/chipseq -profile test -resume --genomicElements beds/*.bed
```

## Get help

Please create an issue to submit your questions.


## Errors


Unexpected error [InvocationTargetException]

 -- Check script 'MintChIP/./modules/local/subworkflow/input_check.nf' at line: 17 or see '.nextflow.log' file for more details

solusion: check the input design table. make sure the fastq file exist.
