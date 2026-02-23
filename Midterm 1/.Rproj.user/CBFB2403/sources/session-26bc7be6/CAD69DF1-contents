#set appropriate working directory

#install and load all necessary/useful packages
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
#install.packages("msa")
#install.packages("seqinr")

library(Biostrings)
library(msa)
library(seqinr)

#import and read sequence data
gene_seqs <- "sequences.fasta"
seqs <- readDNAStringSet(gene_seqs)

#check the loaded sequences 
seqs
width(seqs)
names(seqs)

#Align the sequences - multiple sequence alignment 
myalignment <- msa(seqs, method = "Muscle")

#Check the alignment and count the gaps 
print(myalignment, show = "complete") 
aln <- as(myalignment, "DNAStringSet")
sum(letterFrequency(aln, "-"))

#calculate the consensus sequence
consensus <- msaConsensusSequence(myalignment)
nchar(consensus)

#calculate the GC content (%)
counts <- letterFrequency(aln, letters = c("A","T","G","C"))
totalA <- sum(counts[,"A"])
totalT <- sum(counts[,"T"])
totalG <- sum(counts[,"G"])
totalC <- sum(counts[,"C"])
total <- totalA + totalT + totalG + totalC
GC <- (totalG + totalC) / total * 100
GC

#compute distance matrix
aln_seqinr <- msaConvert(myalignment, type = "seqinr::alignment")
dist_mat <- dist.alignment(aln_seqinr, "identity")
dist_mat

#find the average distances (variances) between samples
avg_dist <- apply(dist_matrix, 1, mean)
avg_dist

#Find the most distant pair 
m <- as.matrix(dist_mat)
diag(m) <- NA 

#max distance
idx_max <- which(m == max(m, na.rm=TRUE), arr.ind=TRUE)[1,]
most_distant <- c(rownames(m)[idx_max[1]], colnames(m)[idx_max[2]])
most_distant

#BLASTn result for consensus sequence (DNA)
#Gene: Hemoglobin Beta (HBB)
#Best match accession: LC121775.1
#Query coverage: 100%
#Identity: 100% 
#E-value: 0.0 

#Translate outlier (most different) sequence to protein
outlier <- seqs["Homo_sapiens_6"]
protein <- Biostrings::translate(outlier)
protein

#save protein translation as a FASTA file
names(protein) <- "Homo_sapiens_6_protein"
writeXStringSet(protein, "Homo_sapiens_6_protein.fasta")



