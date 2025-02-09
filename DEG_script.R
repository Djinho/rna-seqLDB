# Set working directory to the location of the dataset
setwd("//wsl$/Ubuntu/home/djinho/rna-seqLDB/datasets/alignment")

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install DESeq2 package
BiocManager::install("DESeq2", force = TRUE)

# Load required libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(viridis)
library(ReactomePA)

# Read gene counts data
counts <- read.delim("gene_counts_cleaned.csv", header = TRUE, row.names = 1, sep = ",")

# Read metadata and align it with the counts data
metadata <- read.csv("RNA-Seq_Meta_Data.csv", header = TRUE)
rownames(metadata) <- metadata$Run
ordered_metadata <- metadata[match(colnames(counts), rownames(metadata)), ]

# Check if metadata rows match counts columns
if (!all(rownames(ordered_metadata) == colnames(counts))) {
  stop("The metadata rows do not match the counts columns! Please check the data.")
}

# Filter genes with low counts
filtered_counts <- counts[rowSums(counts) >= 50, ]

# Save aligned metadata and filtered counts to CSV files
write.csv(ordered_metadata, "Aligned_Meta_Data.csv", row.names = TRUE)
write.csv(filtered_counts, "Filtered_Gene_Counts.csv", row.names = TRUE)

print("Metadata and counts are aligned, and genes with a sum less than 50 have been filtered out.")

# Reload filtered counts and metadata
counts <- read.csv("Filtered_Gene_Counts.csv", header = TRUE, row.names = 1)
metadata <- read.csv("Aligned_Meta_Data.csv", header = TRUE, row.names = 1)

# Set reference level for disease state
metadata$disease_state <- as.factor(metadata$disease_state)
metadata$disease_state <- relevel(metadata$disease_state, ref = "No dementia")

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ disease_state
)

# Filter low-count genes and run DESeq2
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

# Extract results for LBD vs No dementia comparison
res <- results(dds, contrast = c("disease_state", "LBD", "No dementia"))
write.csv(as.data.frame(res), "DESeq2_Results_final.csv")

# Prepare data for volcano plot
res_df <- as.data.frame(res)
res_df$logP <- -log10(res_df$pvalue)
res_df$Significant <- ifelse(res_df$padj < 0.1 & abs(res_df$log2FoldChange) > 0.5, "Yes", "No")

# Generate volcano plot
ggplot(res_df, aes(x=log2FoldChange, y=logP, color=Significant)) +
  geom_point(alpha=0.6, size=1.5) +
  theme_minimal() +
  labs(title="Volcano Plot of DEGs", x="Log2 Fold Change", y="-Log10 P-Value") +
  scale_color_manual(values=c("No"="gray", "Yes"="red")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed", color="blue") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="blue") +
  theme(plot.title = element_text(hjust=0.5))

# Perform rlog transformation for PCA
rld <- rlog(dds, blind=TRUE)
pcaData <- plotPCA(rld, intgroup="disease_state", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Generate PCA plot
ggplot(pcaData, aes(x=PC1, y=PC2, color=disease_state, label=name)) +
  geom_point(size=3, alpha=0.8) +
  geom_text_repel(size=3) +
  theme_minimal() +
  labs(
    title="PCA Plot of RNA-Seq Data",
    x=paste0("PC1: ", percentVar[1], "% variance"),
    y=paste0("PC2: ", percentVar[2], "% variance"),
    color="Group"
  ) +
  theme(
    plot.title = element_text(hjust=0.5),
    legend.position="right"
  ) +
  stat_ellipse(aes(group=disease_state), type="norm", linetype="dashed", size=0.7)

# Install and load AnnotationDbi and org.Hs.eg.db if not already installed
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("AnnotationDbi")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

library(AnnotationDbi)
library(org.Hs.eg.db)

# Map Ensembl IDs to gene symbols
top_genes <- head(order(res$padj), 20)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(dds),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols <- gene_symbols[!is.na(gene_symbols)]
mat <- assay(rld)[top_genes, ]
rownames(mat) <- gene_symbols[rownames(mat)]
mat <- mat[!is.na(rownames(mat)), ]

# Prepare annotation for heatmap
annotation_col <- data.frame(Group = metadata$disease_state)
rownames(annotation_col) <- rownames(metadata)

# Generate heatmap of top 20 DEGs
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  annotation_col = annotation_col,
  color = colorRampPalette(c("blue", "white", "red"))(50),
  main = "Heatmap of Top 20 DEGs",
  fontsize_row = 10,
  fontsize_col = 10,
  annotation_legend = TRUE
)

# Install and load clusterProfiler for GO enrichment
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

library(clusterProfiler)
library(org.Hs.eg.db)

# Perform GO enrichment analysis
sig_ensembl <- rownames(res_df)[res_df$Significant == "Yes" & !is.na(res_df$padj)]
sig_symbols <- gene_symbols[sig_ensembl]
sig_symbols <- na.omit(sig_symbols)

enriched_GO <- enrichGO(
  gene         = sig_symbols,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  universe     = gene_symbols
)

# Visualize GO enrichment results
dotplot(enriched_GO, showCategory = 10)

# Convert gene symbols to Entrez IDs for Reactome analysis
all_entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = rownames(counts),
  keytype = "SYMBOL",
  column = "ENTREZID"
)
all_entrez_ids <- na.omit(all_entrez_ids)

# Perform Reactome pathway enrichment analysis
enriched_Reactome <- enrichPathway(
  gene = all_entrez_ids,
  organism = "human",
  pvalueCutoff = 0.05,
  universe = all_entrez_ids
)

# Check and visualize Reactome enrichment results
if (!is.null(enriched_Reactome)) {
  print("Reactome enrichment successful!")
  print(head(enriched_Reactome))
  dotplot(enriched_Reactome, showCategory = 10)
} else {
  print("No enriched pathways found in Reactome.")
}

# Print full Reactome enrichment results
print(enriched_Reactome)

