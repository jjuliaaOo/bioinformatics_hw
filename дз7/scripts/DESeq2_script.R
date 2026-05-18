library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)


expr_raw <- read.csv('/projects/mipt_dbmp_biotechnology/rna_seq_diff_exp/raw_counts_ici_samples.tsv', sep='\t', row.names = 1)

meta = data.frame(RNA = colnames(expr_raw))
data = expr_raw
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ 1)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="~/classes/class_14/normalized_counts_ici_samples.tsv", sep="\t", quote=F, col.names=NA)



sample_metadata <- read.csv("/projects/mipt_dbmp_biotechnology/rna_seq_diff_exp/meta_responses.tsv", row.names = 1, sep='\t')
sample_metadata <- sample_metadata %>%
  filter(X0 %in% c('R', 'NR'))
rownames(sample_metadata) <- gsub("-", ".", rownames(sample_metadata))

normalized_counts <- read.csv('~/classes/class_14/normalized_counts_ici_samples.tsv', sep='\t', row.names = 1)

normalized_counts <- normalized_counts[row.names(sample_metadata)]

normalized_counts <- round(normalized_counts)

dds <- DESeqDataSetFromMatrix(countData = normalized_counts,
                              colData = sample_metadata,
                              design = ~ X0)

dds$X0 <- relevel(dds$X0, ref = "NR")
dds <- DESeq(dds)

res <- results(dds)
summary(res)

plotMA(res, main="DESeq2 MA Plot")

pvalue_threshold <- 10**(-2)
log2fc_threshold <- 3

res$significance <- ifelse(res$padj < pvalue_threshold & abs(res$log2FoldChange) > log2fc_threshold, 
                           "Significant", "Not Significant")

res$log2FoldChange <- as.numeric(res$log2FoldChange)
res$padj <- as.numeric(res$padj)


hgnc_table <- read.csv("/projects/mipt_dbmp_biotechnology/rna_seq_diff_exp/hgnc_complete_set.txt", row.names = 1, sep='\t')
symbol_map <- setNames(hgnc_table$symbol, hgnc_table$ensembl_gene_id)

rownames(res) <- sapply(rownames(res), function(id) {
  if (id %in% names(symbol_map) && !is.na(symbol_map[id])) {
    symbol_map[id]
  } else {
    id
  }
})

significant_genes <- subset(res, padj < pvalue_threshold & abs(log2FoldChange) > log2fc_threshold)

ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("grey", "red")) +  # Non-significant genes in grey, significant in red
  labs(title = "Control Subculture", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal() +
  theme(legend.position = "none") +
  geom_text_repel(data = significant_genes, 
                  aes(label = rownames(significant_genes)), size = 3, max.overlaps = 10) +
  xlim(-10, 10)

write.csv(as.data.frame(res), "~/classes/class_14/differential_expression_results_ICI.csv")
