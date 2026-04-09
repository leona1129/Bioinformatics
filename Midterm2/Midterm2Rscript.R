#load all needed packages 
library(ape)
library(Biostrings)
library(UniprotR)
library(protti)
library(r3dmol)

#read/load in created tree file 
tree <- read.tree("metazoa_raxml.raxml.support")

#root tree using specified outgroup
rooted_tree <- root(tree, outgroup = c("Plakina_jani", "Grantia_compressa"),
                    resolve.root = TRUE)

#save rooted tree as a tree file for further analysis
write.tree(rooted_tree, file = "metazoa_raxml_rooted.treefile")


# export tree as a pdf for visualization
pdf("metazoa_raxml_rooted_tree.pdf", width = 12, height = 10)

#make sure pdf has root and is in a legible size 
plot(rooted_tree, cex = 0.7)

#adds labels (bootstrap values) to the nodes for easier interpretation 
nodelabels()

#closes/saves pdf so it has all specified attributes
dev.off()




# read/load provided DNA sequences
dna <- read.dna("metazoa_alignment.gene.fasta", format = "fasta")

# get sequence names
labs <- rownames(dna)

# find the Homo sapiens sequence 
hs_idx <- grep("Homo_sapiens", labs)
hs_name <- labs[hs_idx[1]]
hs_dna <- dna[hs_idx[1], ]

# convert the sequence to characters and remove gaps 
hs_chars <- as.character(hs_dna)
hs_chars_nogap <- hs_chars[!hs_chars %in% c("-", "n", "N", "?")]

# join the nucleotide sequence together so its one continuous string with no gaps
hs_seq <- paste(hs_chars_nogap, collapse = "")

#trim the sequence so is divisible by 3 to accurately translate
#since there are 3 charcters in eachcodon
trim_len <- nchar(hs_seq) - (nchar(hs_seq) %% 3)
hs_seq_trim <- substr(hs_seq, 1, trim_len)

#convert to DNAString object 
dna_string <- DNAString(hs_seq_trim)

# translate nucleotide sequence into protein
protein <- translate(dna_string)
protein_seq <- as.character(protein)

# save protein sequence as fasta file
writeLines(c(
  paste0(">", hs_name),
  protein_seq
), con = "Homo_sapiens_protein.fasta")

#Load in UniProt accession numbers - read in text file list
accessionstring <- readLines("UniProt.txt")
accessionstring

#Get Gene Ontology (GO) information 
# BP, CC, MF
proteininfo <- GetProteinGOInfo(accessionstring)

#Visualize Gene Ontology categories - creates bar graph for BP, CC, MF
PlotGOAll(GOObj = proteininfo, Top = 10, 
          directorypath = getwd(), width = 8,
          height = 5)
