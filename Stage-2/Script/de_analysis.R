##############################################
# Differential expression analysis – Arabidopsis
# File Required: counts.txt  +  metadata.csv
##############################################

# 1. Set working directory
setwd("/Users/chairunnisanuramanda/Hackbio-2025/Hackbio-Internship-NGS-2025/Stage-2") # my working directory

# 2. Load libraries
library(DESeq2) # to perform differential expression analysis
library(pheatmap) # to visualise upregulated and downregulated genes in heatmap
library(dplyr) # to perform data manipulation
library(ggplot2) # to create plot for enrichment
library(clusterProfiler) # for enrichment analysis
library(org.At.tair.db)   # for Arabidopsis ID mapping


# 3. Read in data
counts <- read.delim("counts.txt", header = TRUE, comment.char = "#") # result from RNA-Seq Data Analysis
meta   <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = TRUE) # metadata files that we create

# use regular expression to rename the column names to make it clear
colnames(counts) <- gsub("_trimmed_Aligned\\.sortedByCoord\\.out$", "", colnames(counts))

# Preview
head(counts)
head(meta)

# 4. Prepare count matrix
# Keep only columns of counts from samples in metadata
raw_counts <- counts[ , c(1, 7:12) ]

# Add gene IDs as rownames
rownames(raw_counts) <- counts$Geneid

# Remove the Gene id column
raw_counts <- raw_counts[ , -1 ]

# preview raw_counts
head(raw_counts)

# 5. Create DESeq2 dataset
# Design formula compares conditions (Control vs Treated)
meta$condition <- factor(meta$condition, levels = c("control", "treated"))
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData   = meta,
                              design    = ~ condition)

dds$condition <- relevel(dds$condition, ref = "control")

# 6. Run differential expression
dds <- DESeq(dds)
res <- results(dds)    

# Inspect results
head(res)
summary(res)

# 7. Volcano plot
plot(res$log2FoldChange,
     -log10(res$padj),
     pch = 19, cex = 0.4,
     col = "grey",
     xlab = "Log2 Fold Change (UV-C Treated vs Control)",
     ylab = "-log10 Adjusted P-value",
     main = "Volcano Plot")
abline(v = c(-1, 1), lty = 2, col = "blue")
abline(h = -log10(0.05), lty = 2, col = "red")

# Highlight significant genes
# log2(fold change) > 1 and FDR < 0.05 (according to the paper)
up   <- subset(res, padj < 0.05 & log2FoldChange > 1)
down <- subset(res, padj < 0.05 & log2FoldChange < -1)
points(up$log2FoldChange,   -log10(up$padj),   col="salmon",    pch=19, cex=0.5)
points(down$log2FoldChange, -log10(down$padj), col="lightblue", pch=19, cex=0.5)


# 8. Heatmap of DE genes
# combine up & down for visualization
de_genes <- c(rownames(up), rownames(down))

# create heatmap
if (length(de_genes) > 2) {
  pheatmap(raw_counts[de_genes, ],
           scale = "row",
           show_rownames = FALSE,
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           main = "DE Genes (scaled)")
}

# Create a clear heatmap for Top 50 DE Genes
# reorder columns
ordered_cols <- c("ctl_1","ctl_2","ctl_3","uvc_1","uvc_2","uvc_3")

# Sample treatment information
sample_group <- data.frame(
  Treatment = c("Control","Control","Control","UV-C","UV-C","UV-C")
)
ann_colors <- list(
  Treatment = c("Control" = "skyblue", "UV-C" = "tomato")
)

# treatment to sample names
row.names(sample_group) <- colnames(raw_counts) 

# Create heatmap
top <- head(order(res$padj), 50)
pheatmap(raw_counts[top, ordered_cols],
         scale = "row",
         cluster_cols = FALSE,
         fontsize_row = 6, 
         annotation_col = sample_group[ordered_cols, , drop=FALSE],
         annotation_colors = ann_colors,
         main = "Top 50 DE Genes (grouped by treatment)")


# 9. Filter Top 100 Differentially Expressed Genes
# Filter with |log2FC| > 2.5 and FDR < 0.05
res_filtered <- res %>%
  as.data.frame() %>%
  filter(!is.na(padj)) %>%
  filter(abs(log2FoldChange) > 2.5, padj < 0.05)

# Sort by adjusted p-value (or |log2FC|)
top100 <- res_filtered %>%
  arrange(padj) %>%
  head(100)

head(top100)


# 10. Enrichment Analysis
# Get the gene IDs DE result
de_genes <- rownames(top100)

# Map to Entrez IDs for KEGG and GO
entrez <- mapIds(org.At.tair.db,
                 keys = gsub("gene:", "", de_genes), # strip "gene:"
                 column = "ENTREZID",
                 keytype = "TAIR",
                 multiVals = "first")

# drop NAs
entrez_ids <- na.omit(entrez)
length(entrez_ids)

# GO: Biological Process
ego_bp <- enrichGO(gene       = entrez_ids,
                OrgDb         = org.At.tair.db,
                keyType       = "ENTREZID",
                ont           = "BP",         # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE) 

# Plot
dotplot(ego_bp, showCategory = 20) +
  ggtitle("GO Enrichment: Biological Process") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# GO: Cellular Component
ego_cc <- enrichGO(gene       = entrez_ids,
                OrgDb         = org.At.tair.db,
                keyType       = "ENTREZID",
                ont           = "CC",         # Biological Process
                 pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE) 

# Plot
dotplot(ego_cc, showCategory = 20) +
  ggtitle("GO Enrichment: Cellular Component") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 10),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# GO: Molecular Function
ego_mf <- enrichGO(gene       = entrez_ids,
                   OrgDb         = org.At.tair.db,
                   keyType       = "ENTREZID",
                   ont           = "MF",         # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE) 

# Plot
dotplot(ego_mf, showCategory = 20) +
  ggtitle("GO Enrichment: Molecular Function") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 8),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# Write csv for top 20 GO result
# Convert to data frames and take top 20
bp_df <- head(as.data.frame(ego_bp), 20)
cc_df <- head(as.data.frame(ego_cc), 20)
mf_df <- head(as.data.frame(ego_mf), 20)

# Add a column to indicate the ontology
bp_df$Ontology <- "BP"
cc_df$Ontology <- "CC"
mf_df$Ontology <- "MF"

# Combine into one table
go_top20 <- rbind(bp_df, cc_df, mf_df)


# KEGG Pathway Analysis
# Clean TAIR IDs
tair_ids <- sub("^gene:", "", rownames(res))
de_tair  <- tair_ids[res$padj < 0.05 & abs(res$log2FoldChange) > 1]

# Run KEGG enrichment using TAIR IDs
kegg_res <- enrichKEGG(gene         = de_tair,
                       organism     = "ath",
                       keyType      = "kegg", 
                       pAdjustMethod= "BH",
                       qvalueCutoff = 0.05)

# Create plot
dotplot(kegg_res, showCategory = 20) +
  ggtitle("KEGG Pathway Enrichment") +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 5),
    legend.text  = element_text(size = 12),
    legend.title = element_text(size = 14)
  )

# View top 5 pathways
top5_pathways <- kegg_res@result %>%
  arrange(p.adjust) %>%
  head(5)


# 11. Save outputs
write.csv(res,         "all_results.csv")
write.csv(up,          "upregulated.csv")
write.csv(down,        "downregulated.csv")
write.csv(raw_counts,  "raw_counts.csv")
write.csv(top100, "top100_DE_genes.csv", row.names = TRUE)
write.csv(top5_pathways, "top5_pathways.csv", row.names = FALSE)
write.csv(go_top20, "GO_top20_BP_CC_MF.csv", row.names = FALSE)

