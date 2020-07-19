# The yeastRNASeq experiment data package
# contains FASTQ files from an RNA seq 
# experiment in yeast. When the package 
# is installed, you can access one of the
# FASTQ files by the path given by
# 
# Question 1: What fraction of reads in this 
# file has an A nucleotide in the 5th base
# of the read?

# BiocManager::install('yeastRNASeq')
# BiocManager::install('ShortRead')


library(yeastRNASeq)
library(ShortRead)
library(magrittr)

fastqFilePath <- 
  system.file("reads", 
              "wt_1_f.fastq.gz",
              package = "yeastRNASeq")

reads <- readFastq(fastqFilePath)

log <- 
  sread(narrow(reads, start = 5, end = 5)) == 'A'

length(which(log)) / length(log)

# Answer: [1] 0.363841


# Question 2
# 
# This is a continuation of Question 1.
# 
# Question: What is the average numeric 
# quality value of the 5th base of these reads?

x <- reads %>% 
  narrow(start = 5, end = 5) %>% 
  quality() %>% 
  as('matrix') %>% 
  as.numeric() %>% 
  mean()

# Answer: 28.93346



# Question 3
# 
# The leeBamViews experiment data package contains 
# aligned BAM files from an RNA seq experiment in 
# yeast (the same experiment as in Questions 1 and
# 2, but that is not pertinent to the question).
# You can access one of the BAM files by the path 
# given by
# 
# These reads are short reads (36bp) and have been
# aligned to the genome using a standard aligner, 
# ie. potential junctions have been ignored (this 
# makes some sense as yeast has very few junctions 
# and the reads are very short).
# 
# A read duplicated by position is a read where at 
# least one more read shares the same position.
# 
# We will focus on the interval from 800,000 to 
# 801,000 on yeast chromosome 13.
# 
# Question: In this interval, how many reads are 
# duplicated by position?

BiocManager::install('Rsamtools')
library(Rsamtools)

# BiocManager::install('leeBamViews') 
library(leeBamViews)
bamFilePath <- 
  system.file("bam", "isowt5_13e.bam", package="leeBamViews")

bamFile <- BamFile(bamFilePath)

seqinfo(bamFile)

aln <- scanBam(bamFile)

length(aln)
class(aln)

aln <- aln[[1]]

names(aln)

gr <- GRanges(seqnames = 'Scchr13', 
              ranges = IRanges(start = 800000, end = 801000))


params <- ScanBamParam(which = gr, what = scanBamWhat())

aln <- scanBam(bamFile, param = params)[[1]]                         
   

aln$pos %>% 
  table() %>% 
  .[. > 1] %>% 
  sum()

# Answer: 129



# Question 4
# 
# This is a continuation of Question 3.
# 
# The package contains 8 BAM files in total, 
# representing 8 different samples from 4 groups. 
# A full list of file paths can be had as
# 
# An objective of the original paper was the discovery 
# of novel transcribed regions in yeast. One such region 
# is Scchr13:807762-808068.
# 
# Question: What is the average number of reads across
# the 8 samples falling in this interval?


bpaths <- 
  list.files(
    system.file("bam", package="leeBamViews"), 
    pattern = "bam$", 
    full=TRUE)                         

bpaths                         

# TODO: get the number of reads in the desired interval
# from each file and average them 

gr <- GRanges(seqnames = 'Scchr13', 
              ranges = IRanges(start = 807762, end = 808068))

params <- ScanBamParam(which = gr, what = scanBamWhat())


lapply(bpaths, 
       function(path) {
        bam <- scanBam(path, param = params)[[1]]$pos %>% length
      }) %>% unlist %>% mean()
                       
# Answer: 90.25

# In the lecture on the oligo package an ExpressionSet with 18 samples 
# is constructed, representing normalized data from an Affymetrix gene
# expression microarray. The samples are divided into two groups given 
# by the group\verb|group|group variable.
# 
# Question: What is the average expression across samples in the control 
# group for the “8149273” probeset (this is a character identifier, not 
#                                   a row number).

# BiocManager::install('oligo')

library(GEOquery)
library(oligo)

supp_files <- 
  getGEOSuppFiles('GSE38792')

list.files('GSE38792')
untar("GSE38792/GSE38792_RAW.tar", 
      exdir = 'GSE38792/CEL')


cel_files <- list.files('GSE38792/CEL', full.names = TRUE)

raw_data <- read.celfiles(cel_files)

exprs(raw_data)[1:4, 1:4]

filename <- sampleNames(raw_data)

pData(raw_data)$filename <- filename 

sample_names <- sub(".*_", "", filename)
sample_names <- sub(".CEL.gz$", "", sample_names)
sampleNames(raw_data) <- sample_names

pData(raw_data)$group <-
  ifelse(grepl('^OSA', sample_names), 
                                      'OSA', 'Control')

pData(raw_data)

norm_data <- rma(raw_data)

exprs(norm_data["8149273",]) %>% mean()

# Answer: 7.039974
# Question 5: NOT 11.434


# Question 6
# 
# This is a continuation of Question 5.
# 
# Use the limma package to fit a two group 
# comparison between the control group and 
# the OSA group, and borrow strength across
# the genes using eBayes()\verb|eBayes()|eBayes().
# Include all 18 samples in the model fit.
# 
# Question: What is the absolute value of the log 
# foldchange (logFC\verb|logFC|logFC) of the gene
# with the lowest P.value\verb|P.value|P.value.

# BiocManager::install('limma')

library(limma)

design <- model.matrix( ~ pData(norm_data)$group)
fit <- lmFit(norm_data, design)
fit <- eBayes(fit)

topTable(fit, n = 1)


# Question 6: 0.7126
# Question 6: NOT 0.38

# 7.Question 7
# 
# This is a continuation of Question 6.
# 
# Question: How many genes are differentially 
# expressed between the two groups at an adj.P.value\
# verb|adj.P.value|adj.P.value cutoff of 0.05?

design <- model.matrix( ~ pData(norm_data)$group - 1)

colnames(design) <- c('Control', 'OSA')

fit <- lmFit(norm_data, design)

contast.matrix <- makeContrasts('OSA - Control', levels = design)

contast.matrix

fit_contrast <- contrasts.fit(fit, contast.matrix)

fit_contrast <- eBayes(fit_contrast)

topTable(fit_contrast)

# Question 7: 0
                       
                       
                       
                       
                       
                       
                       
                       