
library(GenomicRanges)
library(AnnotationHub)
library(magrittr)

# 1.Question 1

# Use the AnnotationHub package to obtain data on "CpG Islands" in the human genome.
# 
# Question: How many islands exists on the autosomes?

ah <- AnnotationHub::AnnotationHub()

q <- query(ah, "CpG Islands")

cpg_islands <- q[q$genome == 'hg19']

# download the cpg_islands data
cpg_islands %<>% .[[1]]

# Question: How many islands exists on the autosomes (chr 1 to 22)?

autosome <- 
  cpg_islands[seqnames(cpg_islands) %>% 
                stringr::str_detect('^chr[0-9]+$')] 

length(autosome)

cpg_islands <- autosome

# Question 2: How many CpG Islands exists on chromosome 4.

chr4 <- autosome[seqnames(autosome) == 'chr4']

length(chr4)

# Question 3.
# Obtain the data for the H3K4me3 histone modification for the H1 cell line 
# rom Epigenomics Roadmap, using AnnotationHub. Subset these regions to only
# keep regions mapped to the autosomes (chromosomes 1 to 22).


q <- ah %>% 
  query("H3K4me3") %>%  
  .[.$genome == 'hg19'] %>% 
  .[.$dataprovider == "BroadInstitute"] %>% 
  .[stringr::str_detect(.$title, 'E003')] %>% 
  .[stringr::str_detect(.$title, 'narrowPeak')] 


# narrow peak?

data <- q[[1]]

# filter out non-autosome

autosome <- 
  data %>% 
  seqnames(.) %>% 
  stringr::str_detect('^chr[0-9]+$') %>% 
  data[.]

autosome %>% 
  width() %>% 
  sum()

H3K4me3 <- autosome

# Question 4
# 
# Obtain the data for the H3K27me3 histone modification for the H1 cell line
# from Epigenomics Roadmap, using the AnnotationHub package. Subset these 
# regions to only keep regions mapped to the autosomes. In the return data, 
# each region has an associated "signalValue".
# 
# Question: What is the mean signalValue across all regions on the standard chromosomes?


# The H1 cell line as assayed and quantified by the Roadmap Epigenomics project. 
# The Roadmap Epigenomics project code for the H1 cell line is E003

q <- ah %>% 
  query("H3K27me3") %>% 
  .[.$genome == 'hg19'] %>% 
  .[.$dataprovider == "BroadInstitute"] %>% 
  .[stringr::str_detect(.$title, 'E003')] %>% 
  .[stringr::str_detect(.$title, 'narrowPeak')] 
  

data <- q[[1]] 

# Only keep autosome

autosome <- 
  data %>% 
  seqnames(.) %>% 
  stringr::str_detect('^chr[0-9]+$') %>% 
  data[.]

mean(autosome$signalValue)

H3K27me3 <- autosome

# Question 5
# 
# Bivalent regions are bound by both H3K4me3 and H3K27me3.
# 
# Question: Using the regions we have obtained above, how many bases 
# on the standard chromosomes are bivalently marked?

# Thus, bivalent regions are marked by the intesection of the regions
# we have just found


bivalent <- intersect(H3K4me3, H3K27me3)

# How many bases?
bivalent %>% 
  width() %>% 
  sum()

# How many bases?
bivalent %>% 
  reduce() %>% 
  width() %>% 
  sum()

# Question 6
# 
# We will examine the extent to which bivalent regions overlap CpG Islands.
# 
# Question: how big a fraction (expressed as a number between 0 and 1) 
# of the bivalent regions, overlap one or more CpG Islands?

ov <- findOverlaps(bivalent, cpg_islands)

length(unique(queryHits(ov))) / length(bivalent)

# Question 7
# 
# Question: How big a fraction (expressed as a number between 0 and 1) 
# of the bases which are part of CpG Islands, are also bivalent marked.
# ov <- findOverlaps(cpg_islands, bivalent)

x <- intersect(cpg_islands, bivalent)

sum(width(x)) /
  sum(width(cpg_islands))



# Question 8
# 
# Question: How many bases are bivalently marked within 10kb of CpG Islands?
# 
# Tip: consider using the "resize()"" function

expanded_cpg_islands <- 
  cpg_islands %>% 
  resize(width = width(.+10e3), fix = 'center')

x <- intersect(expanded_cpg_islands, bivalent)

sum(width(x))

#==================================================

expanded_cpg_islands <- 
  cpg_islands %>% 
  resize(width = width(.) + 20e3, fix = 'center')

x <- intersect(expanded_cpg_islands, bivalent)

sum(width(x))




# width(range + 2) --> range + 2 (adds 2 to both ends) --> width

# why is this differnce thatn resize with width(range) + 20e3

resize(range, width = width(range+2), fix = 'center')

resize(range, width = width(range) + 4, fix = 'center')


# it is the same so maybe I just selected the wrong answer the first time?


# Question 9
# 
# Question: How big a fraction (expressed as a number between 0 and 1) 
# of the human genome is contained in a CpG Island?
#   
# Tip 1: the object returned by AnnotationHub contains "seqlengths".
# 
# Tip 2: you may encounter an integer overflow. As described in the session 
# on R Basic Types, you can address this by converting integers to numeric 
# before summing them, "as.numeric()".

lens <- seqlengths(cpg_islands) 

lens <- lens[stringr::str_detect(names(lens), '^chr[0-9]+$')]

lens %<>% as.numeric()

total_len <- sum(lens)

sum(width(cpg_islands)) / total_len

# Question 10
# 
# Question: Compute an odds-ratio for the overlap of 
# bivalent marks with CpG islands.

mat <- matrix(0, nrow = 2, ncol = 2)
rownames(mat) <- c('in.bivalent', 'out.bivalent')
colnames(mat) <- c('in.cpg_island', 'out.cpg_island')
mat[1,1] <- sum(width(intersect(bivalent, cpg_islands)))
mat[1,2] <- sum(width(setdiff(bivalent, cpg_islands, ignore.strand = TRUE)))
mat[2,1] <- sum(width(setdiff(cpg_islands, bivalent, ignore.strand = TRUE)))
mat[2,2] <- total_len - sum(width(union(cpg_islands, bivalent)))

odds_ratio <- mat[1,1] * mat[2,2] / 
             (mat[2,1] * mat[1,2])

odds_ratio 






