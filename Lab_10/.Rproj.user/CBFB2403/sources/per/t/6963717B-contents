install.packages("UniprotR")
install.packages("protti")
BiocManager::install("GenomicAlignments")
install.packages("r3dmol")

library(BiocManager)
library(Biostrings)
library(msa)
library(seqinr)
library(phangorn)
library(UniprotR)
library(protti)
library(r3dmol)

#load in and read DNA fasta file; translate into amino acid/protein sequence
human <- readDNAStringSet("SMN1_human.fasta")
human_protein <- Biostrings::translate(human)
human_protein

#save translation as .fasta file 
Biostrings::writeXStringSet(human_protein,
                            filepath = "human_protein.fasta", 
                            format = "fasta")

#BLAST protein sequence through UniProt; save top five matches' accession numbers to a txt file
#each accession number should be on a new line 
#hit ENTER at the end of txt file to eliminate errors/warnings

#Load in UniProt accession numbers - read in text file list
accessionstring <- readLines("fiveUniprot.txt")
accessionstring

#Get Gene Ontology (GO) information 
# BP, CC, MF
proteininfo <- GetProteinGOInfo(accessionstring)
PlotGoInfo(proteininfo)

#Visualize Gene Ontology categories - creates bar graph for BP, CC, MF
PlotGOAll(GOObj = proteininfo, Top = 10, 
          directorypath = getwd(), width = 8,
          height = 5)
#Pathology and disease information associated with gene 
#Use accession number string
pathologyinfo <- GetPathology_Biotech(accessionstring)
pathologyinfo
Get.diseases(pathologyinfo)

#Access structural information for each protein from UniProt
uniprot_info <- fetch_uniprot(accessionstring)
uniprot_info

#Accession IDs and PBD Structure IDs
uniprot_info[, c("accession", "xref_pdb")]

#Get PBD structural information
pdb_info <- fetch_pdb(c("1G5V","1MHN"))
pdb_info

#Find Alphafold predicted 3D structure 
af_info <- fetch_alphafold_prediction("Q16637")
af_info
