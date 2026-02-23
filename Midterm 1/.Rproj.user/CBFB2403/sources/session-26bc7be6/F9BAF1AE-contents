#1. First step is to set the Working Directory to where the FASTA file is located
setwd("C:/Users/jacob/Downloads/BioinformaticsClass")

#2. The next step will be to load in any libraries that may be used for this code.
# As a baseline, I will load in the Biostrings,rentrez, msa, seqinr, ape, 
# and pwalign libraries.

library(Biostrings)
library(pwalign)
library(seqinr)
library(msa)
library(ape)
library(rentrez)

#There are many codes which overlap across the packages, but I feel like having 
# all downloaded allows for more options when performing code for alignments.

#3. (Q1) Import and Align the DNA sequences
# To start this, I will make a vector so that the DNA sequences may be viewed as 
# a string set

Midterm.seq <- readDNAStringSet("Midterm.Take.Home.Sequences.fasta")

# Next, I will use msaMuscle to align the 20 DNA sequences. 
#This alignment will be put into a vector for later use-
#but to check the alignment you can just type the name into the Console

Midterm.Muscle <- msaMuscle(Midterm.seq)

#4 (Q2) and #5 (Q3)
# For this section, I am combining Questions 2 and 3, as I think that if showing 
# that the consensus sequence has limited gaps, then the alignment 
# can be considered a "good" alignment.

# To start, a vector can be made to grab the consensus sequence from the alignment
 
Midterm.Con <- msaConsensusSequence(Midterm.Muscle)

# Next, We can count the letter frequency within the consensus sequence to analyze
# the amount of gaps present, with the lower the amount the better the alignment.
# Before the frequency can be done, however, the consensus sequence must be put into a string set

Midterm.Aln <- as(Midterm.Con, "DNAStringSet")
sum(letterFrequency(Midterm.Aln, "-"))

# The consensus sequence resulted in 0 gaps formed, meaning the overall alignment was
# pretty good and that the gene samples are relatively the same.

#6 (Q4) Calculating GC content:
# For this, I have a function already made from Lab 2, though it is a few extra steps
# I was very proud when I made it and will use it again here.

calculate_gc <- function(sequence){
  sequence <- toupper(sequence)
  gc_count <- sum(charToRaw(sequence) %in% charToRaw("GC"))
  gc_percent <- (gc_count / nchar(sequence)) * 100
  return(gc_percent)
}
gc_content <- calculate_gc(Midterm.Muscle)
cat("GC Content (%):", gc_content, "\n")

# The GC content for the alignment came out to roughly 51.9%


#7 (Q5) Determine any differences between the samples
# To do this, the distance function in seqinr will be used in order to determine how closely
# related each gene is based on the base pair makeup.

# Before we get the distance, however, the alignment must be converted to match the
# alignment for the Seqinr package

Midterm.Muscle.2 <- msaConvert(Midterm.Muscle, type = "seqinr::alignment")

# Next, the distance can be calculated from the new seqinr alignment

Midterm.Distance <- dist.alignment(Midterm.Muscle.2, "identity")

# To illustrate the distance across all samples, the sequences can be made into a table

Midterm.Matrix <- as.matrix(Midterm.Distance)
Midterm.Table <- as.data.frame(Midterm.Matrix)
View(Midterm.Table)

# According to the table created, all of the samples are the same except for Samples 4, 6,
# and 10, which could indicate that those samples are mutated versions of the sampled gene,
# more than likely frameshift mutations since each number of bp is the same across all samples.

#8 (Q6) Determine the gene

# For this section, I used the GenBank database to determine the gene which is being sampled.
# In the database, I used the BLAST tool on Sequence 1, since there are no mutation in that 
# sequence to determine which gene the samples are from. Based on the BLAST search, the gene from 
# which the samples were taken from are the Homo sapiens hbb gene for beta globin, LC121775.1. 
# I knew this was the gene since not only was the query coverage 100%, but also the E-value was at 0.0, 
# meaning that there was no random chance the sequences matched.

#9 (Q7) Translate the DNA sequence to protein

# For this, I will take the most different sample, that being sample 6 with almost 0.1 bp 
# difference, and translate it into the protein sequence to discover any mutations.
# To do this, I will use the alignment tool MEGA and translate the DNA sequence to Protein.
# Once translated, I will export the data to a FASTA file and create a new vector in the 
# R program.

Midterm.Translation <- readAAStringSet("Midterm_Translation.fas")

#10 (Q8) Determine protein match from sample.

# For this, I will use the database GenBank to compare the Sample 6 translated protein
# to what is already known in the database. In this database, I used the protein BLAST feature
# rather than the BLASTn feature I used previously. According to the BLAST, the protein of 
# Sample 6 is potentially a truncated, or shortened, hemoglobin beta chain, with the 
# accession number of AAZ22545.1.

#11 (Q9) Determine diseases from Protein match.

#To determine diseases from the discovered protein, I will use the database PubMed 
# to research any scientific literature related to the protein. According to PubMed, the 
# main disease caused by this mutation is beta thalassemia, which is a genetic disorder where a section of the
# hemoglobin sub unit, in this case the beta portion, is not fully developed and therefore cannot 
# transport enough oxygen within the blood cells. This disease is genetic and likely 
# due to frameshift and nonsense mutations.

# To determine if the subject has the associated disease, I will run a pairwise 
# alignment between the two sequences
# and retrieve the percent identity. To run the pairwise alignment, 
# I downloaded the sample from the database and read the sample as a string set.

Midterm.Protein <- readAAStringSet("Midterm.Protein.fasta")

# Next, I will run a global pairwise alignment to determine the similarity and percent identity.
# This will confirm if the subject is positive for the disease.

Midterm.Protein <- readAAStringSet("Midterm.Protein.fasta")
Midterm.PWaln <- pairwiseAlignment(Midterm.Translation, Midterm.Protein, type = "global")
pid(Midterm.PWaln)

# According to the pairwise alignment, the sample only had a 
# 27.8 percent identity match globally, which could mean that since the mutated 
# protein does not match the protein found in the database, the subject from 
# which sample 6 was taken does in fact have the genetic disease previously identified.

# All Information and script used in this portion have been pushed to GitHub.
