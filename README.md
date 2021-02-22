# ![qiubio-nf-core/Mint-ChIP](assets/mintchiplogo.png)

## Introduction

**qiubio/MintChIP** is a bioinformatics analysis pipeline used for Mint-ChIP-seq data based on [jianhong/chipseq](https://github.com/jianhong/chipseq).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.



## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Split the reads by [Je](https://gbcs.embl.de/portal/tiki-index.php?page=Je)
3. Run [QiuBio-nf/chipseq pipeline](https://github.com/jianhong/chipseq)

## Installation

Step1. Install nextflow

```bash
wget -qO- https://get.nextflow.io | bash
```

Step2. Pull jianhong/MintChIP. Option '-r dev' will clone dev brance of the pipeline.

```bash
nextflow pull jianhong/MintChiP -r dev
```

Step3. Test the pipeline

```bash
nextflow run jianhong/MintChIP -profile test -r dev --conda
```

## Update

```bash
nextflow pull jianhong/MintChiP
```

## Remove

```bash
nextflow drop jianhong/MintChiP
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
nextflow run jianhong/MintChiP -profile test -resume --genomicElements beds/*.bed
```

## Get help

Please create an issue to submit your questions.


## Errors


Unexpected error [InvocationTargetException]

 -- Check script 'MintChIP/./modules/local/subworkflow/input_check.nf' at line: 17 or see '.nextflow.log' file for more details

solusion: check the input design table. make sure the fastq file exist.
