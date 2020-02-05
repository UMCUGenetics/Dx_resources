.libPaths("/hpc/cog_bioinf/diagnostiek/production/R_libs/3.5.1")

library(ExomeDepth)
library(methods)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
name<-(substr(args[2], 0 ,str_length(args[2])-4 ))
input=args[1]
output=args[2]
target.file=args[3]
reference.file=args[4]
exons=paste(args[5],sep="")

pathToRefBams <- input
refbam.files <- paste0(pathToRefBams, dir(pathToRefBams, "bam$"))
exons.hg19= read.table(exons,sep="\t", header=TRUE)
my.refcounts <- getBamCounts(bed.frame = exons.hg19,
                         bam.files = refbam.files,
                          include.chr = FALSE,
                         referenceFasta = reference.file)
save(my.refcounts,file = output)

