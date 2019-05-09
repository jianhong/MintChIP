#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
# args[1] = samples
# args[2] = annotation
# args[3] = reportFolder

# Load ChIPQC package
library(ChIPQC)

samples = read.csv(args[1])
annotation = args[2]
reportFolder = args[3]

sampleExp = ChIPQC(samples, annotation=annotation)
ChIPQCreport(sampleExp, reportFolder=reportFolder)
