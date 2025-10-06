# ---- Install/load DESeq2 ----
if(!"DESeq2" %in% installed.packages()[,"Package"]) {
  if(!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
  BiocManager::install("DESeq2")
}
library(DESeq2)

count_mat <- as.matrix(read.csv("results/count_matrix_filtered.csv", row.names=1))
coldata   <- read.csv("results/coldata.csv", row.names=1)

x_samples <- grep("\\.x$", colnames(count_mat), value = TRUE)
count_mat <- count_mat[, x_samples]

# Clean sample names (remove long suffix)
colnames(count_mat) <- gsub("_gene_quantifications_GRCh38\\.x$", "", colnames(count_mat))

coldata <- coldata[x_samples, , drop = FALSE]
rownames(coldata) <- gsub("_gene_quantifications_GRCh38\\.x$", "", rownames(coldata))

stopifnot(all(rownames(coldata) == colnames(count_mat)))

dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData   = coldata,
                              design    = ~ condition)

dds <- dds[rowSums(counts(dds)) > 1, ]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, alpha=0.05)  
summary(res)

# Significant genes (FDR < 0.05, |log2FC| > 1)
sig <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ]
cat("Significant genes:", nrow(sig), "\n")

# Save results
out_dir <- "results"
write.csv(as.data.frame(res), file=file.path(out_dir, "DESeq2_all_results.csv"))
write.csv(as.data.frame(sig), file=file.path(out_dir, "DESeq2_significant_genes.csv"))
cat("DESeq2 results saved in:", out_dir, "\n")

pdf(file.path(out_dir, "Volcano_plot.pdf"))
plot(res$log2FoldChange, -log10(res$pvalue),
     pch=20, main="Volcano Plot",
     xlab="log2 Fold Change", ylab="-log10(p-value)")
abline(h=-log10(0.05), col="red", lty=2)
dev.off()

pdf(file.path(out_dir, "MA_plot.pdf"))
plotMA(res, main="MA Plot", ylim=c(-5,5))
dev.off()

rld <- rlog(dds, blind=FALSE)
pdf(file.path(out_dir, "PCA_plot.pdf"))
plotPCA(rld, intgroup="condition")
dev.off()

if(!"pheatmap" %in% installed.packages()[,"Package"]) {
  install.packages("pheatmap", repos="https://cloud.r-project.org")
}
library(pheatmap)

top_genes <- head(order(res$padj), 30)
mat <- assay(rld)[top_genes, ]
mat <- mat - rowMeans(mat)

clean_names <- gsub("_gene_quantifications_GRCh38\\.x$", "", colnames(mat))
colnames(mat) <- clean_names

col_annot <- data.frame(condition = coldata$condition)
rownames(col_annot) <- clean_names

pdf(file.path(out_dir, "Heatmap_top30_genes_clean.pdf"))
pheatmap(mat,
         annotation_col = col_annot,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_col = 8,
         fontsize_row = 7)
dev.off()

cat("All plots saved as PDF in:", out_dir, "\n")