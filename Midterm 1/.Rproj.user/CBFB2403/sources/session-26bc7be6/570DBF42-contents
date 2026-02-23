library(Biostrings)
library(msa)

Deino.87<- readDNAStringSet("16sRNA_Deino_87seq.fasta")

nchar(Deino.87) [[1]]

msaClustalW(Deino.87)

myalignment <- msa(Deino.87)

consensusseq <- msaConsensusSequence(myalignment)

nchar(consensusseq)
