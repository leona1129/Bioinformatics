#install packages (only once if not already done) and load libraries
#DESeq2 preforms the RNAseq analysis and calculates differential expression
#ggplot2 creates a PCA scatter plot for visualizing the relationships
BiocManager::install("DESeq2")
BiocManager::install("GEOquery")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages(GEOquery)

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(GEOquery)

#load and read RNAseq counts file  
counts <- read.table("GSE121212_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                     header = TRUE,
                     sep = "\t",
                     row.names = 1,
                     check.names = FALSE)

#check the file - number of rows and columns
#sample names and gene IDs
dim(counts)
colnames(counts)[1:20]
head(rownames(counts))

#download GEO metadata for the chosen series 
#Same series (GSE121212) for which RNAseq counts was downloaded
gse <- getGEO("GSE121212", GSEMatrix = TRUE)
meta <- pData(gse[[1]])

#determine column names and view desired data 
#viewing title helps determine where each sample came from (healthy vs lesional vs non-lesional)
colnames(meta)
View(meta[, "title", drop = FALSE])

#filter to strictly select healthy control samples based on the title 
#assign/label the group to healthy_meta
healthy_meta <- meta[
  grepl("^CTRL_.*_healthy$", meta$title, ignore.case = TRUE),]

#filter data to strictly select lesional samples based on the title 
#I did not want to include non-lesional samples so a more direct comparison could be made
#assign/label the group to psoriasis_meta
psoriasis_meta <- meta[
  grepl("^PSO_.*_lesional$", meta$title, ignore.case = TRUE),]

#isolate the GSM ID numbers for both groups for further analysis
healthy_samples <- rownames(healthy_meta)
psoriasis_samples <- rownames(psoriasis_meta)

#combine the isolated/selected sample IDs into a single list  
selected_samples <- c(healthy_samples, psoriasis_samples)

#make sure selected IDs are also in the RNAseq count file 
#create a new list/subset containing only the intersecting/coinciding samples 
#this was done to avoid errors as some IDs were not included in the count 
selected_samples2 <- intersect(selected_samples, colnames(counts))
counts_subset <- counts[, selected_samples2]

#create a table that labels the selected samples as healthy or psoriasis (two columns)
#this allows for cleaner visualization of expression analysis in downstream processing 
sample_info <- data.frame(
  sample = selected_samples2,
  condition = ifelse(selected_samples2 %in% psoriasis_samples,
                     "Psoriasis",
                     "Healthy"))

#set row names as sample IDs
rownames(sample_info) <- sample_info$sample

#convert condition into a factor with levels 
#this is to allow DESeq2 to properly calculate fold changes by defining Healthy as the reference 
sample_info$condition <- factor(sample_info$condition,
                                levels = c("Healthy", "Psoriasis"))

#count the number of samples in each column/category 
table(sample_info$condition)

#Visualize a portion (the head) of the table to ensure it is properly formatted
head(sample_info)

#exclude/filter genes with less than 10 counts (assign to object) to run through DESeq
#this allows for better results/analysis 
dds <- DESeqDataSetFromMatrix(
  countData = counts_subset,
  colData = sample_info,
  design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

#MAIN ANALYSIS 
#run data through DESeq
dds <- DESeq(dds)

#show the names of the model coefficients 
resultsNames(dds)
res <- results(dds, contrast = c("condition", "Psoriasis", "Healthy"))

#results check
head(res)

#summary of unregulated, down regulated, and insignificant genes
summary(res)

#Convert DESeq2 results for easier handling
res_df <- as.data.frame(res)

#save gene IDs as a column
res_df$gene <- rownames(res_df)

#sort genes by significance (most to least significant)
res_df <- res_df[order(res_df$padj), ]

#check data - highest ranked genes 
head(res_df)

#Export the differential expression results table as a CSV file
write.csv(res_df, "psoriasis_vs_healthy_results.csv", row.names = FALSE)

#Keep only statistically significant differentially expressed genes
sig_res <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

#Count how many significant genes were found
nrow(sig_res)

#Count genes upregulated in psoriasis
sum(sig_res$log2FoldChange > 1, na.rm = TRUE)

#Count genes downregulated in psoriasis
sum(sig_res$log2FoldChange < -1, na.rm = TRUE)

#Select the top 10 significant genes 
#includes gene ID, fold change, and adjusted p-value
top10 <- head(sig_res[, c("gene", "log2FoldChange", "padj")], 10)

#Export the top 10 significant genes to a CSV file
write.csv(top10, "top10_psoriasis_genes.csv", row.names = FALSE)

#normalize count data- - varience stabilizing
vsd <- vst(dds, blind = FALSE)

#Run PCA, group samples by condition
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

#PCA scatter plot
p1 <- ggplot(pca_data, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: Psoriasis vs Healthy")

#save PCA plot as an image
ggsave("PCA_psoriasis_vs_healthy.png", p1, width = 6, height = 5)

#create a volcano plot and save as image
png("Volcano_psoriasis_vs_healthy.png", width = 900, height = 700)
plot(res_df$log2FoldChange,
     -log10(res_df$padj),
     pch = 16,
     cex = 0.6,
     xlab = "log2 Fold Change",
     ylab = "-log10 adjusted p-value",
     main = "Volcano Plot: Psoriasis vs Healthy")
dev.off()

#Identify top 30 genes with expression values -averages
topgenes <- na.omit(res_df$gene[1:30])
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)

#create a heat map for the top 30 genes and save as image
png("Heatmap_top30_genes.png", width = 1400, height = 900)
pheatmap(mat,
         annotation_col = sample_info["condition"],
         main = "Top 30 Differentially Expressed Genes")
dev.off()

#read the annotation file (downloaded from GEO) 
annot <- read.delim("Human.GRCh38.p13.annot.tsv.gz",
                    header = TRUE,
                    sep = "\t",
                    quote = "",
                    fill = TRUE,
                    comment.char = "")

#Merge differential expression results with annotation table 
res_annot <- merge(res_df,
                   annot,
                   by.x = "gene",
                   by.y = "GeneID",
                   all.x = TRUE)

#check to make sure gene ID, gene symbol, fold change, and adjusted p-value are shown
head(res_annot[, c("gene", "Symbol", "log2FoldChange", "padj")])

#Identify top 10 most significant annotated genes 
top10_annot <- head(res_annot[order(res_annot$padj),
                              c("gene", "Symbol", "log2FoldChange", "padj")], 10)
#export to a CVS file
write.csv(top10_annot, "top10_psoriasis_gene_names.csv", row.names = FALSE)
