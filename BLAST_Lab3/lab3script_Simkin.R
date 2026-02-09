library(Biostrings)
Refseq<- readDNAStringSet("Lab_3_RefSeq.fasta")
narvavirus<- readDNAStringSet("Narnavirus.fasta")

library(pwalign)
align <- pairwiseAlignment(Refseq, narvavirus)

pid(align)
