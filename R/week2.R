library(Biostrings)

dna1 <- DNAString('ACGT-N')
dna1

DNAStringSet('ACG')

dna2 <- DNAStringSet(c('ACGT', 'GTCA', 'GCTA'))
dna2

dna1[2:4]
dna2[2:3]

# To not return a DNAStringSet
# we must use double bracket [[
dna2[2]
dna2[[2]]

names(dna2) <- paste0('seq', 1:3)
dna2

# The fulls set of string classes are: 
#   DNAString: DNA sequences
#   RNAString: RNA sequences
#   AAString Amino Acid sequences (proteins)
#   BString "Big" sequences, any alphabet
#   "XString" -- any of the above


# Basic Functionality: 
width(dna2)
sort(dna2)
rev(dna2)# reverses order of set
rev(dna1)# reverses order of string

# Biological Functionality: 
translate(dna2)

# Computing GC / Alphabet content
alphabetFrequency(dna1)
alphabetFrequency(dna2)

letterFrequency(dna2, 'GC')

consensusMatrix(dna2, as.prob = TRUE)

library(BSgenome)

all_genomes <- available.genomes()
grep("Hsapiens", all_genomes, value = TRUE)
grep("Scerevisiae", all_genomes, value = TRUE)

library(BSgenome.Scerevisiae.UCSC.sacCer2)
Scerevisiae

seqlengths(Scerevisiae)
seqnames(Scerevisiae)

# load a chromosome
Scerevisiae$chrI

# compute the GC content
letterFrequency(Scerevisiae$chrI, "GC", as.prob = TRUE)

# We could use lapply to iterate over the whole genome, 
# being all of the chromosomes, though bsapply from the
# BSgenomes library exists to load, compute and unload
# each chromosome to keep memory overhead down

# Notice that we pass the chromosomal data and function that we 
# wish to apply in the BSParams type, but any options that we would like
# each of the function calls to have we can pass as parameters to bsapply
param <- new('BSParams', X = Scerevisiae, FUN = letterFrequency)

# passing the letters param to bsapply 
# as opposed to in param provies a list
head(bsapply(param, letters = 'GC'))

# The output of bsapply can be simplified as well 
param <- new('BSParams', 
             X = Scerevisiae, 
             FUN = letterFrequency, 
             simplify = TRUE)
head(bsapply(param, letters = 'GC'))

                      
# To conclude, the GC percentage of the genome is 
sum(bsapply(param, letters = 'GC')) / sum(seqlengths(Scerevisiae))
                      
##########################
#  Biostrings - Matching
##########################

library(Biostrings)
library(BSgenome)
library(BSgenome.Scerevisiae.UCSC.sacCer2)

# Biostrings package is very useful for searching the genome

# Global vs local alignment -- matching whole genomes vs genes

# (v)matchPattern -- match a single sequence against one or many
#                    other sequencues

# (v)matchPDict -- match a possibly large set of sequnces against
#                  one sequence or many sequencues

# Both of these functions allows for a small set of mismatches and indels. 
# The term Dict here is used because the function builds a dictionary
# over the sequences

# There are also similar functions using count instead of match
# ex: countPatterns, which return the number of matches instead of the 
# matches

# In many ways these functions are similar to using short read
# aligners like Bowtie, though these functions are designed to be 
# comprehensive (retrun all matches satidfying certain criteria)

dnaseq <- DNAString('ACGTACGT')
matchPattern(dnaseq, Scerevisiae$chrI)

countPattern(dnaseq, Scerevisiae$chrI)

vmatchPattern(dnaseq, Scerevisiae)                      

# We can now use vcountPattern to examine
# matches of a gene sequency across all 
# chromosomes
head(vcountPattern(dnaseq, Scerevisiae))

# Not how the return object of vmatchPattern
# is a GRanges given the exact information of 
# where the string matches. 

# Since the sequency we were searching for was its 
# own reverse compliment we see matches on the forward
# and reverse strands

dnaseq == reverseComplement(dnaseq)


##########################
# Specialized Alignments
##########################