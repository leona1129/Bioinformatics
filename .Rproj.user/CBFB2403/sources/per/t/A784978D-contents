library(dplyr)
install.packages("here")
library(here)
library(seqinr)
library(rentrez) 
library(stringr)
library(Biostrings)

?entrez_search
?entrez_fetch

seq <- entrez_fetch(db="nucleotide",
                    id="NC_001477", 
                    rettype= "fasta")
seq
substr(seq, 1, 100)
seq_split <- stringr::str_split(seq, "")
seq_split 

write(seq, "Dengue_seq.fasta")

dengueseq <- read.fasta("Dengue_seq.fasta")

class(dengueseq)
length(dengueseq)
length(dengueseq$NC_001477.1)



Dengueseq<- readDNAStringSet("Dengue_seq.fasta")

count<- letterFrequency(Dengueseq, letters = c("A","T","G","C"))
count

total<- sum(count)

A <- (count[, "A"] / total) * 100
T <- (count[, "T"] / total) * 100
G <- (count[, "G"] / total) * 100
C <- (count[, "C"] / total) * 100

GC<- G+C

A
T
G
C
GC

library(seqinr)

class(seq)


head(Dengueseq)


