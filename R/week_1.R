library(IRanges)

# IRanges -- Integer Intervals
# IRanges form the base of a more
# useful tool in Gemomic Analysis -- GRanges. 
# GRanges are simply a wrapper around GRanges that 
# have more attributes more closely associated with 
# Genomic Analysis

# To form an IRanges you only must supply 
# 2 of the 3 arguments, and thus the third
# can be inferred (start, end, width). Then 
# as many intervals as were supplied will be 
# populated into an IRanges list (NOT A MATRIX) 
ir1 <- IRanges(start = 1:3, end = 3:5)
print(ir1)

ir2 <- IRanges(start = 1:3, width=3)
print(ir2)

all.equal(ir1, ir2)

start(ir1)
end(ir1)
width(ir1)

# we can also assign attributes in this same way
width(ir1) <- 7
print(ir1)

# IRanges may have names
names(ir1) <- paste0('A', 1:3)
print(ir1)

# Notice how and IRange is a 1-d vector!
dim(ir1)
length(ir1)
ir1[1]
ir1[2]
ir1[3]

# Thus we can concatenate two ranges with c()
c(ir1[1], ir1[2])

# Reduce will take the union of the IRanges Set, sort of.
# What this really does is calculcates the normal irange
# of the set, such that each integer only occurs once in the
# set, and their are as few ranges as possible. This result
# is acheived from the `reduce()` function, and many 
# IRanges / GRanges functions return normal IRanges.
reduce(ir1)

# In a Genomic context, this can be seen as giving a set 
# with overlapping exons  

# Disjoin can be said to do the oposite of reduce, giving 
# the disjoint sets of IRanges. Again each integer will only
# appear once in the returned set of IRanges, however now
# each original IRanges will be returned only with its 
# distict sections. 

# `disjoin()` can be a bit confusing, but drawing it out seems
# to help the understanding

print(ir1)
disjoin(ir1)

# Functions that allow for intra-range mainpulation where 
# each range gets mapped to a new range
# shift(), narrow(), flank(), resize(), restrict()
?shift
?narrow
?flank
?resize
?restrict

# There are different from the intra-range functions that we
# saw before such as reduce() and disjoin()

print(ir1)
resize(ir1, width = 1, fix = 'start')
resize(ir1, width = 1, fix = 'center')

# resizing fixing on center is useful for barcoded DNA strands
# where the barcoding is equidistant form the center and not removed
# in pre-analysis steps

# IRanges as sets
# available functions: 
# union(), intersect(), setdiff(), gaps()

ir1 <- IRanges(start = seq.int(1, 5, 2), width = 1)
ir2 <- IRanges(start = seq.int(4, 6), width = 1)
print(ir1)
print(ir2)

union(ir1, ir2)
intersect(ir1, ir2)
reduce(c(ir1, ir2))

?"intra-range-methods"

# There are also element-wise (pair-wise) versions of these
# functions with the "p" pre-fix such as: 
# punion(), pintersect(), psetdiff(), pgap()

# findOverlaps

ir1 <- IRanges(start = c(1,4,8), end = c(3,7,10))
ir2 <- IRanges(start = c(3,4), width = 3)
print(ir1)
print(ir2)
ov <- findOverlaps(ir1, ir2)
print(ov)

# This returns an incidence matrix between
# the query and subject ranges, being ir1 and 
# ir2 respectively 

queryHits(ov)
subjectHits(ov)

# Another way of saying that there is hit between
# the subject and the query is that they have a 
# non-empty intersection

intersect(ir1[queryHits(ov)[1]], 
          ir2[subjectHits(ov)[2]])

# unique(queryHits) gives you the indices of the query 
# ranges which actually had an overlap. unique is needed as
# a query range may overlap multiple subject ranges

queryHits(ov)
unique(queryHits(ov))

# There are many arguments to find overlaps and some might provide
# backpocket tricks that are needed at times
# Ex. minoverlap can limit overlaps returned to only 
# overlaps that provide a certain number of bases

args(findOverlaps)

# countOverlaps returns the number of overlaps, this is more efficient than
# calling length(findOverlaps(ir1, ir2))

countOverlaps(ir1, ir2)

## Finding nearest IRanges()
# nearest, precede and follow

nearest(ir1, ir2)


## GenomicRanges -- GRanges

# install.packages('BiocManager')
# BiocManager::install('GenomicRanges')
  
library(GenomicRanges)


gr <- GRanges(seqnames = 'chr1', 
              strand = c('+', '-', '+'), 
              ranges = IRanges(start = c(1,3,5), width = 3))

# accessor functions: strand(), seqnames(), ranges(), 
# start(), end(), width()

# Also since we have the strand we have access to operations that 
# are strand dependent such as upstream() and downstream()

# flank() will retieve the GRanges of x upstream nucleotides
flank(gr, 2, start = FALSE)

## GRanges, seqinfo

seqinfo(gr)

seqlengths(gr) <- c('chr1' = 10)

seqlevels(gr)

seqlengths(gr)

# gaps() returns the stretches of the genome
# not covered by the GRanges

gaps(gr)

# lets add another chromosome to the GRanges

seqlevels(gr) <- c('chr1', 'chr2')
seqnames(gr) <- c('chr1', 'chr2', 'chr1')

# seqnames in combination with seqlevels are
# used in sorting GRanges
sort(gr)

seqlevels(gr) <- c('chr2', 'chr1')

sort(gr)

# We can also associate a genome with a GRanges
# this is stored in the seqinfo of the GRanges
genome(gr) <- 'hg19'

seqinfo(gr)

# This is very useful to alert you if you are 
# analyzing data from multiple genomes to not 
# get the readings confused during analysis

gr2 <- gr

genome(gr2) <- 'hg18'

findOverlaps(gr, gr2)

## DataFrame

# The S4Vectors package introduced the DataFrame class, 
# a very similar object to data.frame, though allows 
# columns of any class

ir <- IRanges(start = 1:2, width = 3)

df1 <- DataFrame(iranges = ir)

df1

df1$iranges

# Notice what happens when we try to make 
# a data.frame of IRanges

df2 <- data.frame(iranges = ir)

df2

# GRanges can have associated metadata

gr <- GRanges(seqname = 'chr1',
              strand = c('+', '-', '+'), 
              ranges = IRanges(start = c(1,3,5), 
                               width = 3))

values(gr) <- DataFrame(score = c(0.1, 0.5, 0.3))

gr              

values(gr)           
              
gr$score
              
gr$score2 <- gr$score * 0.2
              
# Note that the strand information is used when calling
# findOverlap with two GRanges, and the * strand will match
# both + and - strands

gr2 <- GRanges(seqnames = c('chr1', 'chr2', 'chr1'), 
               strand = '*', 
               ranges = IRanges(start = c(1,3,5), width = 3))

gr2              
gr              

findOverlaps(gr, gr2)

# A useful utility function to only return 
# the overlapping regions
subsetByOverlaps(gr, gr2)

# Another useful utility for conversion

df <- data.frame(chr = 'chr1', start = 1:3, end = 4:6, score = 7:9)
df

makeGRangesFromDataFrame(df)              

# Biology usecase I 
# Suppose we want to identify transcription factor (TF)
# binding sites that overlaps known SNPs
# 
# Input objects are 
# snps: a GRanges (of width 1)
# TF: a GRanges
#
# findOverlaps(snps, TF) -- watch out for strand

# Biology usecase II
#
# Suppose we have a set of differntially methylation regions (DMRs)
# (think genomic regions) and a set of CpG Islands and we want to find
# all DMRs withink 10kb of a CpG Island
#
# Input objects are 
# dmrs: a GRanges
# islands: a GRanges
#
# big_islands <- resize(islands, width = 20e3 + width(islands), fix = 'center')
# findOverlaps(dmrs, big_islands) -- watch out for strand

## Advanced Genomic Ranges

BiocManager::install('GenomeInfoDb')

library(GenomeInfoDb)
library(GenomicRanges)

# drop/keep seqlevels()

gr <- GRanges(seqnames = c('chr1', 'chr2'), 
              ranges = IRanges(start = 1:2, end = 4:5))

dropSeqlevels(gr, 'chr1', pruning.mode = 'coarse')

keepSeqlevels(gr, 'chr2', pruning.mode = 'coarse')

# We can also get rid of strange looking chromosomes with 

keepStandardChromosomes(gr)

gr <- GRanges(seqnames = c('chr1', 'chrU345'), 
              ranges = IRanges(start = 1:2, end = 4:5))

keepStandardChromosomes(gr, pruning.mode = 'coarse')

# Change chromosome naming style

gr <- GRanges(seqnames = 'chr1',
              ranges = IRanges(start = 1:2, width = 2))

new_style <- mapSeqlevels(seqlevels(gr), 'NCBI')

gr <- renameSeqlevels(gr, new_style)

seqlevels(gr)

## AnnotationHub

library(AnnotationHub)

# First we create an AnnotationHub instance

ah <- AnnotationHub()
ah

# use the [] operator to see information about the 
# listing

ah[1]

# use [[]] to retrieve the information from the list

ah[[1]]

unique(ah$dataprovider)

unique(ah$rdataclass)

ah <- subset(ah, species == 'Homo sapiens')

query(ah, 'H3K4me3')

hist <- display(ah)

##  Usecase -- Basic GRanges and AnnotationHub

library(GenomicRanges)
library(rtracklayer)
library(AnnotationHub)

# First use AnnotationHub to get data on 
# homo sapiens

ah <- AnnotationHub()
ah <- subset(ah, species == 'Homo sapiens')

# Next we search for keyworkds

qhs <- query(ah, 'H3K4me3')
qhs <- query(qhs, 'Gm12878')

qhs

# Note that we have some searches that dont have Gm12878
# in their title, this illustrates how query searches over 
# multiple fields

qhs$title

# This result is a great illustration of the mess of public data.
# It turns our thatE116is a RoadmapEpigenetics code for the Gm12878 cell 
# line. The first 5 hits are from ENCODE, hosted at UCSC and the last 6 
# hits are from Roadmap Epigenomics hosted at the Broad Institute. 
# The Roadmapdata is different representation (and peaks) from the 
# same underlying data. For the ENCODE data,two different centers did 
# the same experiment in the same cell line (Broad, hit 1) and (Uw, hit 2-5),
# whereUwexposed data on two replicates (whatever that means). 
# These two experiments seems tobe analyzed using different algorithms. 
# It is even possible that the Roadmap data is from the sameraw data but 
# just analyzed using different algorithms

# Lets take a look at the narrowPeak data

gr1 <- subset(qhs, title == "wgEncodeUwHistoneGm12878H3k4me3StdPkRep1.narrowPeak.gz")[[1]]
gr1

gr2 <- subset(qhs, title == 'E116-H3K4me3.narrowPeak.gz')[[1]]
gr2

summary(width(gr1))
summary(width(gr2))

# we will stick with gr1 for now.

# Lets get some promoter coordinates

qhs <- query(ah, 'RefSeq')
qhs

# There are so many RefSeq datasets as there is a RefSeq
# Genes and RefSeq other for each genome

library(magrittr)
genome(gr1) %>% unique()
qhs$genome

refseq <- qhs[qhs$genome == 'hg19' & qhs$title == 'RefSeq Genes']

refseq

refseq <- refseq[[1]]

refseq

# lets look at the number of isoforms per gene name

table(table(refseq$name))

# this shows that almost all genes have a single transcript

# these ranges do not include introns

promoters <- promoters(refseq)
table(width(promoters))

# Now we compute which promoters have H3K4me3 peak in them 

ov <- findOverlaps(promoters, gr1)
ov

# lets compute the percent of peaks that are in a promoter region 

length(unique(queryHits(ov))) / length(gr1)

# and the percent of promoters that have a peak in them 

length(unique(subjectHits(ov))) / length(promoters)

# Lets take look at the size of these results vs the 
# size of the human genome 

sum(width(reduce(gr1))) / 1e6

sum(width(reduce(promoters))) / 1e6

# the human genome is 300e6 megabases, since these numbers
# are so small compared to the size of the human genome, this
# should convince us that the overlap is highly unlikely to happen
# by chance

sum(width(intersect(gr1, promoters))) / 1e6
# 0

# we need to ignore the strand
sum(width(intersect(gr1, promoters, ignore.strand = TRUE))) / 1e6
# 3.02

# lets computer which bases are in promotors / peaks

prom <- reduce(promoters, ignore.strand = TRUE)
peaks <- reduce(gr1)
both <- intersect(prom, peaks)
only_prom <- setdiff(prom, both)
only_peaks <- setdiff(peaks, both)
overlap_mat <- matrix(0, ncol = 2, nrow = 2)
colnames(overlap_mat) <- c('in_peaks', 'out_peaks')
rownames(overlap_mat) <- c('in_promoters', 'out_promoters')
overlap_mat[1,1] <- sum(width(both))
overlap_mat[1,2] <- sum(width(only_prom))
overlap_mat[2,1] <- sum(width(only_peaks))
overlap_mat[2,2] <- 3e9 - sum(overlap_mat)
overlap_mat <- round(overlap_mat / 1e6, 2)
overlap_mat

# This calculation is very back-of-the-envelope, though 
# lets compute and odds-ratio for this table

odds_ratio <- overlap_mat[1,1] * overlap_mat[2,2] / 
             (overlap_mat[2,1] * overlap_mat[1,2])

odds_ratio 
# 18.2

# This odds-ratio shows an enrichment of peaks in promoters
# We can get a feel of how the genome size effects this 

overlap_mat[2,2] <- 1.5e9 / 1e6

odds_ratio <- overlap_mat[1,1] * overlap_mat[2,2] / 
             (overlap_mat[2,1] * overlap_mat[1,2])
odds_ratio

# The odds-ratio (determinant) is still smaller, but 
# is still bigger than 1