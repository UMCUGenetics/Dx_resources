library(ExomeDepth)
library(methods)
library(stringr)

options(scipen = 50)

args = commandArgs(trailingOnly=TRUE)
auto_ref_file = args[1]  # argument ExomeDepth reference file 
target.file = args[2]  # argument target BED file
reference.file = args[3]  # argument Reference genome
input_exon.hg19 = paste(args[4], sep="")  # argument target Exon file
probability = as.numeric(args[5])  # argument transition probability
bam.files = args[6]  # input bam file (including full path)
model = args[7]  # used calling model (eg HC/UMCU)
refset = args[8]  # used reference set (eg Jan2020)
expectedCNVLength = as.numeric(args[9])  # expected length CNV (eg 50000)
run_id = args[10]  # runID. 
log_extension = args[11]  # extension for CNV log file 
igv_extension = args[12]  # extension for igv file
sample_id = args[13]  # SampleID

exons.hg19 = read.table(input_exon.hg19, sep="\t", header=TRUE)
data(Conrad.hg19)

################## Call autosomal events  ###################################
print("Calling CNVs")
load(file=auto_ref_file)
my.counts <- getBamCounts(
    bed.frame = exons.hg19,
    bam.files = bam.files,
    include.chr = FALSE,
    referenceFasta = reference.file
)

save(my.counts, file = "my.counts")

my.counts.dafr <- as(my.counts[, colnames(my.counts)], 'data.frame')
my.counts.dafr$chromosome <- gsub(as.character(my.counts.dafr$space),
    pattern = 'chr',
    replacement = ''
) ##remove the annoying chr letters

samplecounts.mat <- as.matrix(my.counts.dafr[, grep(names(my.counts.dafr), pattern='*.bam')])
nsamples <- ncol(samplecounts.mat)

my.refcounts.dafr <- as(my.refcounts[, colnames(my.refcounts)], 'data.frame')

my.refcounts.dafr$chromosome <- gsub(as.character(my.refcounts.dafr$space),
    pattern = 'chr',
    replacement = ''
) ##remove the annoying chr letters

my.ref.samples <- colnames(my.refcounts.dafr)[7:(ncol(my.refcounts.dafr)-1)]

#loop over samples in my.counts
for (i in 1:nsamples) {
  my.current.samplename <- colnames(my.counts.dafr[6+i])
  message(my.current.samplename)
  my.reference.set <- as.matrix(my.refcounts.dafr[, my.ref.samples])
  my.choice <- select.reference.set(test.counts=samplecounts.mat[, i],
      reference.counts=(my.reference.set),
      bin.length=(my.counts.dafr$end - my.counts.dafr$start)/1000,
      n.bins.reduced = 10000
  )

  file_name=paste(model, refset, sample_id, run_id, log_extension, sep="_")
  write_file = file(file_name, "w")
  line1=paste("Number of selected reference samples ",toString(length(my.choice[[1]])), sep="\t")
  write(line1, write_file, append=T)
  line2=paste("Selected reference samples ", toString(my.choice[[1]]), sep="\t")
  write(line2, write_file, append=T)
  close(write_file)

  my.matrix <- as.matrix( my.refcounts.dafr[, my.choice$reference.choice, drop = FALSE])
  my.reference.selected <- apply(X = my.matrix, MAR = 1, FUN = sum)

  #CNV calling
  all.exons <- new('ExomeDepth',
      test = samplecounts.mat[, i],
      reference = my.reference.selected,
      formula = 'cbind(test, reference) ~ 1' 
  )

  all.exons <- CallCNVs(x = all.exons,
      transition.probability = probability,
      chromosome = my.counts.dafr$space,
      start = my.counts.dafr$start,
      end = my.counts.dafr$end,
      name = my.counts.dafr$names,
      expected.CNV.length = expectedCNVLength
  )

  all.exons <- AnnotateExtra(x = all.exons,
      reference.annotation = Conrad.hg19.common.CNVs,
      min.overlap = 0.5,
      column.name = 'Conrad.hg19'
  )

  str(all.exons)
  output.file <- paste(model, refset, sample_id, run_id,'exome_calls.csv', sep = "_")


  save(all.exons, file=paste(model, refset, sample_id, run_id, "all.exons", sep = "_"))
  refsize <- toString(length(my.choice[[1]]))
  correlation <- all.exons@refcorrelation
  print_array <- cbind(all.exons@CNV.calls, correlation, refsize)
  write.csv(file = output.file,
      x = print_array,
      row.names = FALSE
  )
}

print(all.exons@refcorrelation)

## Select reference ratio + boundries
x <- all.exons
anno <- x@annotations
selected <- which(anno$start >= 0)		# select all exons
anno <- anno[selected,]				# select all exons
anno$expected <- x@expected[ selected ]
anno$freq <- x@test[ selected ]/ (x@reference[selected ] + x@test[selected])
anno$middle <- 0.5*(anno$start + anno$end)
anno$ratio <- anno$freq/ anno$expected
anno$ratio[is.na(anno$ratio)] <- 0
anno$log2ratio <- log2(anno$ratio)
anno$log2ratio[is.na(anno$log2ratio)] <- 0
anno$test <- x@test[ selected ]
anno$reference <- x@reference[ selected ]
anno$total.counts <- anno$test + anno$reference
if (length( x@phi ) == 1) anno$phi <- x@phi else anno$phi <- x@phi [ selected ]

message(str(anno))
for (i in 1:nrow(anno)) {
    anno$my.min.norm[ i ] <- qbetabinom (p = 0.01, size = anno$total.counts[ i ], phi = anno$phi[ i ], prob = anno$expected[ i ])
    anno$my.max.norm[ i ] <- qbetabinom (p = 0.99, size = anno$total.counts[ i ], phi = anno$phi[ i ], prob = anno$expected[ i ])
}

anno$my.min.norm.prop <- anno$my.min.norm / anno$total.counts
anno$my.max.norm.prop <- anno$my.max.norm / anno$total.counts
head(anno)
chroms <- anno$chromosome
starts <- anno$start - 1
ends <- anno$end
name <- anno$name
freq <- anno$freq
observed <- anno$ratio
log2ratio<-anno$log2ratio
expected <- anno$expected
min <- anno$my.min.norm.prop/anno$expected
max <- anno$my.max.norm.prop/anno$expected
correl <- all.exons@refcorrelation

ref_df = data.frame(chroms, starts, ends, name, observed, log2ratio, freq, expected, min, max, correl)
colnames(ref_df) <- c("chr", "start", "end", "locusID", "ratio_test", "log2ratio_test", "frequency_test", "ratio_expected", "ref_min.ratio", "ref_max.ratio", "refset_correlation")
write.table(ref_df, paste(model, refset, sample_id, run_id, igv_extension, sep = "_"), sep="\t", row.names=FALSE, quote=FALSE)
