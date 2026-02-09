library(Biostrings)
Refseq<- readDNAStringSet("Lab_3_RefSeq.fasta")
narvavirus<- readDNAStringSet("Narnavirus.fasta")

library(pwalign)
align <- pairwiseAlignment(Refseq, narvavirus)

pid(align)













identity(align)
 
class(Refseq)

library(pwalign)

Refseq
narvavirus

aln <- pwalign::aligned(Refseq, narvavirus, type= "global", scoreOnly=FALSE)

aln<- pairwiseAlignment(Refseq, narvavirus)

aligned1 <- aln$a
aligned2 <- aln$b

score(aln)

pattern <- AAStringSet(c("Refseq", "narvavirus"))
subject <- AAString("narvavirus")
pa1 <- pairwiseAlignment(pattern, subject, substitutionMatrix="BLOSUM50",
                         gapOpening=3, gapExtension=1)

pid(aln)
