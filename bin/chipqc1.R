#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = T)
# args[1] = bam
# args[2] = narrow peaks
# args[3] = id

# Load ChIPQC package
library(ChIPQC)
library(rtracklayer)

bam = args[1]
peak_file = args[2]
id = args[3]
peaks <- import(peak_file, format = "bed")
if(length(peaks)>0){
  tryCatch({
    sample = ChIPQCsample(bam, peaks = peaks)
    ChIPQCreport(sample, reportName=id, reportFolder=id)
  }, error=function(e){
    dir.create(id, recursive = TRUE)
    writeLines(as.character(e), con = file.path(id, "error.txt"))
  })
}else{
  dir.create(id, recursive = TRUE)
}

