.libPaths("/hpc/diaggen/software/production/R_libs/ExomeDepthv2.0.0/3.5.1/")

library(ExomeDepth)
library(methods)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
input=args[1]  #input folder
output=args[2]  #output folder  
name<-(substr(output, 0 ,str_length(output)-4 ))
target.file=args[3]  #Probe target bed file  
reference.file=args[4]  #Reference genome (fasta)   
exons=paste(args[5],sep="")  #Exon target tsv file

pathtorefbams <- input
refbam.files <- paste0(pathtorefbams, dir(pathtorefbams, "bam$"))
exons.hg19= read.table(exons,sep="\t", header=TRUE)
refcounts <- getBamCounts(
    bed.frame = exons.hg19,
    bam.files = refbam.files,
    include.chr = FALSE,
    referenceFasta = reference.file
)
save(refcounts,file = output)

