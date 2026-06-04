library(GEOquery)

# Download GSE68719
gse <- getGEO("GSE68719", GSEMatrix = TRUE, getGPL = FALSE)
gse


pheno <- pData(gse[[1]])
table(pheno$title)

# Get supplementary files list
getGEOSuppFiles("GSE68719", makeDirectory = TRUE, fetch_files = FALSE)



getGEOSuppFiles("GSE68719", makeDirectory = TRUE, fetch_files = TRUE)

list.files("GSE68719/")

counts <- read.table(
  "GSE68719/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz",
  header = TRUE, row.names = 1, sep = "\t"
)

dim(counts)
counts[1:3, 1:5]

library(DESeq2)


samples <- colnames(counts)
condition <- ifelse(grepl("^C_", samples), "Control", "PD")

coldata <- data.frame(
  row.names = samples,
  condition = factor(condition, levels = c("Control", "PD"))
)

# Check balance
table(coldata$condition)

dim(counts_mat)
dim(coldata)


colnames(counts_mat)
rownames(coldata)

# Remove "symbol" from coldata
coldata <- coldata[rownames(coldata) != "symbol", , drop = FALSE]

# Align counts to coldata
counts_int <- round(counts_mat[, rownames(coldata)])

# Verify
dim(counts_int)
dim(coldata)


dds <- DESeqDataSetFromMatrix(
  countData = counts_int,
  colData = coldata,
  design = ~ condition
)

dds <- dds[rowSums(counts(dds) >= 10) >= 5, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "PD", "Control"))
summary(res)


library(ggplot2)
library(ggrepel)

res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$symbol <- gene_symbols[rownames(res_df)]

# Significant genes (strict cutoff)
res_df$sig <- "NS"
res_df$sig[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "UP"
res_df$sig[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "DOWN"

table(res_df$sig)



# Top genes to label
top_genes <- res_df[res_df$sig != "NS", ]
top_genes <- top_genes[order(top_genes$padj), ][1:20, ]

# Volcano plot
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c("UP" = "#d62728", "DOWN" = "#1f77b4", "NS" = "grey70")) +
  geom_text_repel(data = top_genes, aes(label = symbol), size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot — PD vs Control (GSE68719)",
    subtitle = "Prefrontal Cortex Bulk RNA-seq",
    x = "log2 Fold Change", y = "-log10 adjusted p-value",
    color = "Expression"
  ) +
  theme_bw()

ggsave("GSE68719_volcano.png", width = 8, height = 6, dpi = 300)



# Fix symbol mapping
res_df$symbol <- gene_symbols[match(rownames(res_df), names(gene_symbols))]

# Check
head(res_df[res_df$sig == "UP", c("symbol", "log2FoldChange", "padj")], 5)



head(gene_symbols)
head(names(gene_symbols))

# Rebuild symbol lookup with rownames
symbol_map <- data.frame(
  gene_id = rownames(counts),
  symbol = counts$symbol,
  stringsAsFactors = FALSE
)

# Map to res_df
res_df$symbol <- symbol_map$symbol[match(rownames(res_df), symbol_map$gene_id)]

# Check
head(res_df[res_df$sig == "UP", c("symbol", "log2FoldChange", "padj")], 5)


# Top 20 genes to label
top_genes <- res_df[res_df$sig != "NS" & !is.na(res_df$padj), ]
top_genes <- top_genes[order(top_genes$padj), ][1:20, ]

# Volcano with labels
ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.4, size = 1) +
  scale_color_manual(values = c("UP" = "#d62728", "DOWN" = "#1f77b4", "NS" = "grey70")) +
  geom_text_repel(data = top_genes, aes(label = symbol), size = 3, max.overlaps = 20) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot — PD vs Control (GSE68719)",
    subtitle = "Prefrontal Cortex Bulk RNA-seq",
    x = "log2 Fold Change", y = "-log10 adjusted p-value",
    color = "Expression"
  ) +
  theme_bw()

ggsave("GSE68719_volcano_labeled.png", width = 9, height = 6, dpi = 300)


top30 <- res_df[res_df$sig != "NS", ]
top30 <- top30[order(top30$padj), ][1:30, ]
nrow(top30)

mat <- counts_int[rownames(top30), ]
rownames(mat) <- top30$symbol
dim(mat)
# Step 3
mat_scaled <- t(scale(t(mat)))

# Step 4
ann <- data.frame(Condition = coldata$condition)
rownames(ann) <- colnames(mat)

# Step 5 - save to file
png("GSE68719_heatmap.png", width = 1000, height = 800)
pheatmap(mat_scaled,
         annotation_col = ann,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = colorRampPalette(c("#1f77b4","white","#d62728"))(100),
         main = "Top 30 DE Genes — PD vs Control"
)
dev.off()

getwd()

C:\\Users\\GANAPATHIRAJAN\\OneDrive\\Documents\\GSE68719_heatmap.png


pheatmap(mat_scaled,
         annotation_col = ann,
         show_colnames = FALSE,
         fontsize_row = 8,
         color = colorRampPalette(c("#1f77b4","white","#d62728"))(100),
         main = "Top 30 DE Genes — PD vs Control",
         filename = "GSE68719_heatmap.png",
         width = 10,
         height = 8
)



library(clusterProfiler)
library(org.Hs.eg.db)


up_genes <- res_df$symbol[res_df$sig == "UP" & !is.na(res_df$symbol)]
down_genes <- res_df$symbol[res_df$sig == "DOWN" & !is.na(res_df$symbol)]

# GO enrichment - UP
go_up <- enrichGO(gene = up_genes,
                  OrgDb = org.Hs.eg.db,
                  keyType = "SYMBOL",
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05
)

# Plot
barplot(go_up, showCategory = 15, title = "GO BP — Upregulated in PD")



ggsave("GSE68719_GO_up.png", width = 10, height = 7, dpi = 300)


# Save all results
write.csv(res_df, "GSE68719_DESeq2_results.csv", row.names = TRUE)
save(dds, res, res_df, counts_int, coldata, gene_symbols, symbol_map,
     file = "GSE68719_P01_workspace.RData")
