library(ExomeDepth)
library(methods)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
input = args[1]  #input folder
output = args[2]  #output file
name <- (substr(output, 0 ,str_length(output)-4 ))
target.file = args[3]  #Probe target bed file  
reference.file = args[4]  #Reference genome (fasta)   
exons = paste(args[5],sep = "")  #Exon target tsv file

refbam.files <- paste0(input, dir(input, "bam$"))
exons.hg19 = read.table(exons,sep = "\t", header = TRUE)
my.refcounts <- getBamCounts(bed.frame = exons.hg19, bam.files = refbam.files, include.chr = FALSE, referenceFasta = reference.file)
save(my.refcounts,file = output)
