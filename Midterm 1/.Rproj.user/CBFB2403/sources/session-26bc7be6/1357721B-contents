#install any required packages
install.packages("phangorn")


#load the libraries
library(Biostrings)
library(msa)
library(seqinr)
library(phangorn)

#Read the gene sequences 
human <- readDNAStringSet("SMN1_human.fasta")
catfish <- readDNAStringSet("smn1_catfish.fasta")
macaque <- readDNAStringSet("SMN1_macaque.fasta")
mouse <- readDNAStringSet("Smn1_mouse.fasta")
pig <- readDNAStringSet("SMN1_pig.fasta")

#Combine sequences 
seqs <- c(human, catfish, macaque, mouse, pig)
names(seqs) <- c("Human","Catfish","Macaque","Mouse","Pig")

#Run MSA with MUSCLE 
Alignment <- msa(seqs, method = "Muscle")

#Check the alignment gaps
print(Alignment, show = "complete") 
aln <- as(Alignment, "DNAStringSet")
sum(letterFrequency(aln, "-"))

#How long is the alignment 
width(aln)

#Calculate GC%
counts <- letterFrequency(aln, letters = c("A","T","G","C"))
totalA <- sum(counts[,"A"])
totalT <- sum(counts[,"T"])
totalG <- sum(counts[,"G"])
totalC <- sum(counts[,"C"])
total <- totalA + totalT + totalG + totalC
GC <- (totalG + totalC) / total * 100
GC

#Compute distance matrix 
aln_seqinr <- msaConvert(Alignment, type = "seqinr::alignment")
dist_mat <- dist.alignment(aln_seqinr, "identity")
dist_mat


#Find the most distant and closest pair **extra step**
m <- as.matrix(dist_mat)
diag(m) <- NA 

#max distance
idx_max <- which(m == max(m, na.rm=TRUE), arr.ind=TRUE)[1,]
most_distant <- c(rownames(m)[idx_max[1]], colnames(m)[idx_max[2]])
most_distant

#min distance
idx_min <- which(m == min(m, na.rm=TRUE), arr.ind=TRUE)[1,]
closest <- c(rownames(m)[idx_min[1]], colnames(m)[idx_min[2]])
closest

#translate sequence to amino acid 
aln <- as(Alignment, "DNAStringSet")
human_seq <- aln["Human"]
human_nogap <- DNAStringSet(gsub("-", "", as.character(human_seq)))
human_protein <- Biostrings::translate(human_nogap)
human_protein

Alignment_phyDat <- msaConvert(Alignment, type="phangorn::phyDat")

write.phyDat(Alignment_phyDat, file="alignment.fasta", format="fasta")
