# MintChIP
pipeline for MintChIP

## install

### requirements

[nextflow](https://www.nextflow.io/)

[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

[Je](https://gbcs.embl.de/portal/tiki-index.php?page=Je)

[trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

[bwa](http://bio-bwa.sourceforge.net/)

[samtools](http://samtools.sourceforge.net/)

[Homer](http://homer.ucsd.edu/homer/)

[R](https://www.r-project.org/), ChIPQC[http://bioconductor.org/packages/release/bioc/html/ChIPQC.html]

### install

```
git clone https://github.com/jianhong/MintChIP.git
```

## configuration

edit the nextflow.config file depend on your local setting.

## prepare barcode file

barcode file should be a text file with two columns: sample_Name and barcode. sample file is bar.txt

## prepare fastq files

put the original fastq files in the folder of ${dataDir}/fastq/original/

## prepare inputs

put the inputs for the samples with inputs. sample file is inputs.txt. The column should be 
treatment fastq file name, treatment barcode file name, input fastq file name, input barcode file name. 

## prepare additional parameters

put additional parameters for findPeaks or trimmomatic in a text file. sample file is findPeaks.txt. 
The delimiter must be tab.

## run pipeline

```
nextflow run main.nf --dataDir fastq/original --barcode bar.txt
```

### with inputs

```
nextflow run main.nf --dataDir fastq/original --barcode bar.txt --inputs inputs.txt
```

### with additional parameters for findPeaks

```
nextflow run main.nf --dataDir fastq/original --barcode bar.txt --findPeaks findPeaks.txt
```

### with additional parameters for trimmomatic

```
nextflow run main.nf --dataDir fastq/original --barcode bar.txt --trim trimmomatic.txt
```

## output

The output files will be in the folder of results be default. You can change it in the nextflow.config file.
