library(magrittr)

# Question 1: What is the mean expression across all 
# features for sample 5 in the ALL dataset (from the ALL package)?

# BiocManager::install('ALL')
library(ALL)
data(ALL)
ALL

# Genes / Features are rows, 
# samples either across time 
# or subjects are in columns
exprs(ALL)[1:4, 1:4]

exprs(ALL) %>% 
  .[,5] %>% 
  mean()

# Answer: 5.629627


# Question 2
# 
# We will use the biomaRt package to annotate an Affymetrix 
# microarray. We want our results in the hg19 build of the 
# human genome and we therefore need to connect to Ensembl 
# 75 which is the latest release on this genome version. 
# How to connect to older versions of Ensembl is described 
# in the biomaRt package vignette; it can be achived with 
# the command 
# 
# Question: Using this version of Ensembl, annotate each feature 
# of the ALL dataset with the Ensembl gene id. How many probesets
# (features) are annotated with more than one Ensembl gene id?

# BiocManager::install('biomaRt')
library(biomaRt)


mart <- useMart(host='feb2014.archive.ensembl.org', 
                biomart = "ENSEMBL_MART_ENSEMBL")


ensembl <- useDataset("hsapiens_gene_ensembl", mart)
feature_ids <- featureNames(exprs(ALL))

bm <- 
  getBM(attributes = c('ensembl_gene_id', 'affy_hg_u95av2'), 
        filters = 'affy_hg_u95av2',
        values = feature_ids, 
        mart = ensembl)

dim(bm)


# We have the Affymetrix "Probe IDs" and want the ensemble 
# gene name 

bm %>% 
  .[,2] %>% 
  table() %>% 
  table() %>% 
  tail(-1) %>% 
  sum()

# Answer: 1045


# Question 3: How many probesets (Affymetrix IDs) are annotated 
# with one or more genes on the autosomes (chromosomes 1 to 22).


ensembl %>% 
  listAttributes() %>% 
  head()

bm <- 
  getBM(attributes = c('ensembl_gene_id', 'affy_hg_u95av2', 'chromosome_name'), 
        filters = 'affy_hg_u95av2',
        values = feature_ids, 
        mart = ensembl)

bm %<>% as.data.frame()

autosome_bm <- 
  bm %>% 
  dplyr::filter(chromosome_name %in%
           as.character(seq.int(1,22)))


autosome_bm %>% 
  .$affy_hg_u95av2 %>% 
  table() %>% table() %>% 
  sum()


# Answer: 11016



# Use the MsetEx dataset from the minfiData package. 
# Part of this question is to use the help system to 
# figure out how to address the question.
# 
# Question 4: What is the mean value of the Methylation 
# channel across the features for sample “5723646052_R04C01”?

# BiocManager::install('minfiData')
library(minfiData)
data(MsetEx)

meth_data <- 
  getMeth(MsetEx) %>% 
  .[, '5723646052_R04C01']

meth_data %>% 
  mean()

# Answer: 7228.277


# Question 5: Access the processed data from 
# NCBI GEO Accession number GSE788. What is
# the mean expression level of sample GSM9024?

# BiocManager::install('GEOquery')
library(GEOquery)

eList <- getGEO('GSE788')

length(eList)

data <- eList[[1]]

pData(data)

exprs_data <- exprs(data)

exprs_data[, 2] %>% mean()

# Answer: 756.432



# We are using the airway dataset from the airway package.
# 
# Question 6: What is the average of the average length across 
# the samples in the expriment?

# BiocManager::install('airway')
library(airway)
data(airway)

# assay(airway, 'counts')

# airway %>%
#   rowRanges() %>%
#   width() %>%
#   mean()

airway %>% 
  colData() %>% 
  .$avgLength %>% 
  mean()

# Answer: 113.75


# We are using the airway dataset from the airway package.
# The features in this dataset are Ensembl genes.
# 
# Question 7: What is the number of Ensembl genes which have 
# a count of 1 read or more in sample SRR1039512?

data <- 
  assay(airway, 'counts') %>% 
  .[, 'SRR1039512']

length(which(data >= 1))

# Answer: 25699




# Question 8: The airway dataset contains more than 
# 64k features. How many of these features overlaps 
# with transcripts on the autosomes (chromosomes 1-22) 
# as represented by the TxDb.Hsapiens.UCSC.hg19.knownGene package?
#   
# Clarification: A feature has to overlap the actual transcript,
# not the intron of a transcript. So you will need to make sure that 
# the transcript representation does not contain introns.

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene















# The expression measures of the airway dataset are 
# the number of reads mapping to each feature. In 
# the previous question we have established that 
# many of these features do not overlap autosomal
# transcripts from the TxDb.Hsapiens.UCSC.hg19.knownGene.
# But how many reads map to features which overlaps 
# these transcripts?
#   
# Question 9: For sample SRR1039508, how big a percentage 
# (expressed as a number between 0 and 1) of the total 
# reads in the airway dataset for that sample, are part 
# of a feature which overlaps an autosomal 
# TxDb.Hsapiens.UCSC.hg19.knownGene transcript?






# Consider sample SRR1039508 and only consider features which 
# overlap autosomal transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene.
# We should be able to very roughly divide these transcripts
# into expressed and non expressed transcript. Expressed transcripts
# should be marked by H3K4me3 at their promoter. The airway dataset
# have assayed “airway smooth muscle cells”. In the Roadmap 
# Epigenomics data set, the E096 is supposed to be “lung”. 
# Obtain the H3K4me3 narrowPeaks from the E096 sample using 
# the AnnotationHub package.
# 
# Question: What is the median number of counts per feature 
# (for sample SRR1039508) containing a H3K4me narrowPeak in 
# their promoter (only features which overlap autosomal
# transcripts from TxDb.Hsapiens.UCSC.hg19.knownGene 
# are considered)?
#   
# Clarification: We are using the standard 2.2kb default 
# Bioconductor promotor setting.
# 
# Conclusion Compare this to the median number of counts 
# for features without a H3K4me3 peak. Note that this short 
# analysis has not taken transcript lengths into account 
# and it compares different genomic regions to each other;
# this is highly suscepticle to bias such as sequence bias.
# 





