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

## run pipeline

```
nextflow run main.nf --dataDir fastq/original --barcode bar.txt
```

## output

The output files will be in the folder of results be default. You can change it in the nextflow.config file.
