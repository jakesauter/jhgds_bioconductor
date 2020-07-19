# Question 1: 
# What is the GC content of “chr22” 
# in the “hg19” build of the human genome?
# Tip: The reference genome includes “N” bases; 
# you will need to exclude those.

library(BSgenome)

all_genomes <- available.genomes()

# Find all genomes with hsapiens
grep('Hsapiens', all_genomes, value = TRUE)
genome <- "BSgenome.Hsapiens.UCSC.hg19"
  
# BiocManager::install(genome)

library(BSgenome.Hsapiens.UCSC.hg19)

genome <- BSgenome.Hsapiens.UCSC.hg19
chr22 <- genome$chr22

x <- alphabetFrequency(chr22, baseOnly = TRUE)

sum(x[c('G', 'C')]) / sum(x[c('A', 'C', 'G', 'T')])

# Answer: 0.4798807

# Background: In the previous assessment we studied H3K27me3 
# “narrowPeak” regions from the H1 cell line (recall that the 
# Roadmap ID for this cell line is “E003”). We want to examine 
# whether the GC content of the regions influence the signal;
# in other words wether the reported results appear biased by 
# GC content.
# 
# Question 2: What is mean GC content of H3K27me3 “narrowPeak” 
# regions from Epigenomics Roadmap from the H1 stem cell line 
# on chr 22.
# 
# Clarification: Compute the GC content for each peak region 
# as a percentage and then average those percentages to compute 
# a number between 0 and 1.

library(AnnotationHub)
library(magrittr)
ah <- AnnotationHub()

q <- ah %>% 
  query("H3K27me3") %>%  
  .[.$genome == 'hg19'] %>% 
  .[.$dataprovider == "BroadInstitute"] %>% 
  .[stringr::str_detect(.$title, 'E003')] %>% 
  .[stringr::str_detect(.$title, 'narrowPeak')] 

data <- q[[1]]

# filter out only chr22

chr22 <- 
  data %>% 
  seqnames(.) %>% 
  grepl('^chr22$', .) %>% 
  subset(data, .)


view <- Views(genome, chr22)

gc_probs <- letterFrequency(view, 'GC', as.prob = TRUE)

gc_probs %>% mean

# Answer: 0.528866

# Question 3: 
# 
# The “narrowPeak” regions includes information 
# on a value they call “signalValue”.
# 
# Question: What is the correlation between GC 
# content and “signalValue” of these regions 
# (on chr22)?

signal_value <- DataFrame(view)$signalValue

cor(gc_probs, signal_value)

# Answer: 0.004467924


# The “narrowPeak” regions are presumably reflective of a 
# ChIP signal in these regions. To confirm this, we want 
# to obtain the “fc.signal” data from AnnotationHub package 
# on the same cell line and histone modification. This data 
# represents a vector of fold-change enrichment of ChIP signal 
# over input.
# 
# Question 4: what is the correlation between the “signalValue” 
# of the “narrowPeak” regions and the average “fc.signal” across 
# the same regions?
#   
# Clarification: First compute the average “fc.signal” for 
# across each region, for example using “Views”; this yields 
# a single number of each region. Next correlate these numbers
# with the “signalValue” of the “narrowPeaks”.

q <- ah %>% 
  query("E003-H3K27me3") 

q

# E003-H3K27me3.fc.signal.bigwig

fc.signal.data <- q[['AH32033']]

# Now we need to compute the average fc.signal 
# across each narrow_peak region 

data = import(fc.signal.data, 
              which = GRanges('chr22', 
                              ranges = IRanges(1, 10^8)), 
              as = 'RleList')$chr22

chr22_ranges <- ranges(chr22)

fc.signal <- Views(data, chr22_ranges)

fc.signal.mean <- mean(fc.signal)

cor(fc.signal.mean, signal_value)

# Answer:  0.9149614


# Referring to the objects made and defined 
# in the previous question.
# 
# Question 5: How many bases on chr22 have an 
# fc.signal greater than or equal to 1?

# Assuming each signal value corresponds
# to one base
length(which(data > 1))

# Answer: 10914671



# The H1 stem cell line is an embryonic stem cell line, 
# a so-called pluripotent cell. Many epigenetic marks 
# change upon differentiation. We will examine this. 
# We choose the cell type with Roadmap ID “E055” which 
# is foreskin fibroblast primary cells.
# 
# We will use the “fc.signal” for this cell type for the 
# H3K27me3 mark, on chr22. We now have a signal track for 
# E003 and a signal track for E055. We want to identify 
# regions of the genome which gain H3K27me3 upon 
# differentiation. These are regions which have a 
# higher signal in E055 than in E003. To do this 
# properly, we would need to standardize (normalize) 
# the signal across the two samples; we will ignore 
# this for now.
# 
# Question 6: Identify the regions of the genome where the
# signal in E003 is 0.5 or lower and the signal in E055 
# is 2 or higher.
# 
# Tip: If you end up with having to intersect two different Views,
# note that you will need to convert the Views to IRanges or GRanges 
# first with 
# ir <- as(vi, "IRanges")\verb|ir <- as(vi, "IRanges")|ir <- as(vi, "IRanges").                        
                                       
                                       

# H3K27me3 E055 | E003 chr22

# E003-H3K27me3.fc.signal.bigwig
q <- ah %>% 
  query("E003-H3K27me3") 

e003_signal_file <- q[['AH32033']]

e003_signal_data <- 
  import(e003_signal_file, 
         which = GRanges('chr22', 
                         ranges = IRanges(1, 10^8)), 
         as = 'RleList')$chr22

# E005-H3K27me3.fc.signal.bigwig
q <- ah %>% 
  query("E055-H3K27me3") 

e055_signal_file <- q[['AH32470']]

e055_signal_data <- 
  import(e055_signal_file, 
         which = GRanges('chr22', 
                         ranges = IRanges(1, 10^8)), 
         as = 'RleList')$chr22


length(which(e003_signal_data <= 0.5 & e055_signal_data >= 2))

# Answer: 1869937

# CpG Islands are dense clusters of CpGs. 
# The classic definition of a CpG Island 
# compares the observed to the expected 
# frequencies of CpG dinucleotides as well 
# as the GC content.
# 
# Specifically, the observed CpG frequency
# is just the number of “CG” dinucleotides
# in a region. The expected CpG frequency is
# defined as the frequency of C multiplied by
# the frequency of G divided by the length of 
# the region.
# 
# Question 7: What is the average observed-to-expected 
# ratio of CpG dinucleotides for CpG Islands on chromosome 22?

q <- query(ah, "CpG Islands")

cpg_islands <- q[q$genome == 'hg19']
cpg_islands %<>% .[[1]]

chr22_cpg_islands <- 
  cpg_islands %>% 
  subset(seqnames == 'chr22')


chr22_cpg_view <- 
  Views(genome, chr22_cpg_islands)


expected_cg_freq <- 
  (letterFrequency(chr22_cpg_view, 'G')  * 
  letterFrequency(chr22_cpg_view, 'C'))  / 
  width(chr22_cpg_islands)

observed_cg_freq <- 
  dinucleotideFrequency(chr22_cpg_view)[,'CG'] 

(observed_cg_freq / expected_cg_freq) %>% mean()

# Answer: 0.8340929


# A TATA box is a DNA element of the form “TATAAA”. 
# Around 25% of genes should have a TATA box in their 
# promoter. We will examine this statement.
# 
# Question 8 : How many TATA boxes are there on chr 22 
# of build hg19 of the human genome?
#   
#   Clarification: You need to remember to search both 
# forward and reverse strands.

genome <- BSgenome.Hsapiens.UCSC.hg19
chr22 <- genome$chr22

countPattern('TATAAA', chr22) + 
countPattern('TATAAA', reverseComplement(chr22)) 

# Answer: 27263s


# Question 9: How many promoters of transcripts on 
# chromosome 22 containing a coding sequence, 
# contains a TATA box on the same strand as the transcript?
#   
# Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene 
# package to define transcripts and coding sequence. Here, we
# defined a promoter to be 900bp upstream and 100bp downstream 
# of the transcription start site.

# BiocManager::install('TxDb.Hsapiens.UCSC.hg19.knownGene')
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

chr22 <- GRanges(seqnames = 'chr22', ranges = IRanges(1, 10^8))

chr22_cds <- 
  subsetByOverlaps(cdsBy(txdb, by = 'tx'), chr22)

chr22_transcripts <- 
  subsetByOverlaps(transcripts(txdb), chr22)

chr22_tx_with_cds <- 
  subsetByOverlaps(chr22_transcripts, chr22_cds)

prom <- promoters(chr22_tx_with_cds, upstream = 900, downstream = 100)

x = matchPattern('TATAAA', genome$chr22)
y = matchPattern('TATAAA', reverseComplement(genome$chr22))


x_range = IRanges(start = start(x), end(x))
y_range = IRanges(start = start(y), end(y))
x_gr = GRanges('chr22', ranges = x_range, strand = '+')
# y_gr = GRanges('chr22', ranges = y_range, strand = '-')
# tata_gr <- c(x_gr, y_gr)
tata_gr <- x_gr

subsetByOverlaps(prom, tata_gr) %>% length() 

# help for question 9
# https://www.coursera.org/learn/bioconductor/discussions/weeks/2/threads/b5xWh91EEea2mQ7cbimPbg


# It is possible for two promoters from different transcripts to 
# overlap, in which case the regulatory features inside the overlap 
# might affect both transcripts. This happens frequently in bacteria.
# 
# Question 10: How many bases on chr22 are part of more than one promoter 
# of a coding sequence?
#   
# Clarification: Use the TxDb.Hsapiens.UCSC.hg19.knownGene package to 
# define transcripts and coding sequence. Here, we define a promoter to 
# be 900bp upstream and 100bp downstream of the transcription start site. 
# In this case, ignore strand in the analysis.



























