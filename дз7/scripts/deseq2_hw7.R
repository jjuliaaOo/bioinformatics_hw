setwd("C:/Users/trety/OneDrive/Рабочий стол/hw7")

library(DESeq2)
library(ggplot2)
library(dplyr)
library(ggrepel)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Загружаем подготовленные в ноутбуке таблицы
counts <- read.table(
  "results/counts_for_deseq.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

meta <- read.table(
  "results/meta_for_deseq.tsv",
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

cat("Размер counts:\n")
print(dim(counts))

cat("Размер meta:\n")
print(dim(meta))

cat("\nГруппы:\n")
print(table(meta$response))

# Проверяем совпадение порядка образцов
cat("\nСовпадает ли порядок образцов counts и meta:\n")
print(all(colnames(counts) == rownames(meta)))

# На всякий случай приводим counts к целым числам
counts <- round(as.matrix(counts))
storage.mode(counts) <- "integer"

# Задаём порядок уровней:
# NR — базовая группа, R — сравниваемая группа
meta$response <- factor(meta$response, levels = c("NR", "R"))

# Создаём объект DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = meta,
  design = ~ response
)

# Убираем совсем низкоэкспрессируемые гены
dds <- dds[rowSums(counts(dds)) >= 10, ]

cat("\nРазмер dds после фильтрации низкоэкспрессируемых генов:\n")
print(dim(dds))

# Запускаем DESeq2
dds <- DESeq(dds)

# Нормализованные counts
normalized_counts <- counts(dds, normalized = TRUE)
write.table(
  normalized_counts,
  file = "results/normalized_counts_ici_samples.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

# Результаты: R против NR
res <- results(dds, contrast = c("response", "R", "NR"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

# Упорядочим по padj
res_df <- res_df[order(res_df$padj), ]

# Добавляем колонку значимости
res_df$significant <- ifelse(
  !is.na(res_df$padj) &
    res_df$padj < 0.05 &
    abs(res_df$log2FoldChange) >= 1,
  "yes",
  "no"
)

write.table(
  res_df,
  file = "results/deseq2_results_R_vs_NR.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\nСводка DESeq2:\n")
print(summary(res))

cat("\nКоличество значимых генов при padj < 0.05 и |log2FC| >= 1:\n")
print(table(res_df$significant))

cat("\nПервые строки результата:\n")
print(head(res_df, 10))

# Таблица только значимых генов
sig_df <- res_df[res_df$significant == "yes", ]

write.table(
  sig_df,
  file = "results/deseq2_significant_genes_R_vs_NR.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Volcano plot
plot_df <- res_df
plot_df$minus_log10_padj <- -log10(plot_df$padj)
plot_df$minus_log10_padj[is.infinite(plot_df$minus_log10_padj)] <- NA

plot_df$category <- "not significant"
plot_df$category[
  !is.na(plot_df$padj) &
    plot_df$padj < 0.05 &
    plot_df$log2FoldChange >= 1
] <- "up in R"

plot_df$category[
  !is.na(plot_df$padj) &
    plot_df$padj < 0.05 &
    plot_df$log2FoldChange <= -1
] <- "down in R"

top_genes <- plot_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(10)

volcano <- ggplot(plot_df, aes(x = log2FoldChange, y = minus_log10_padj, color = category)) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_id),
    size = 3,
    max.overlaps = 20
  ) +
  theme_minimal() +
  labs(
    title = "DESeq2 volcano plot: R vs NR",
    x = "log2FoldChange",
    y = "-log10(adjusted p-value)",
    color = "Category"
  )

print(volcano)

ggsave(
  filename = "figures/03_deseq2_volcano_R_vs_NR.png",
  plot = volcano,
  width = 8,
  height = 6,
  dpi = 300
)

cat("\nГотово. Сохранены файлы:\n")
cat("results/normalized_counts_ici_samples.tsv\n")
cat("results/deseq2_results_R_vs_NR.tsv\n")
cat("results/deseq2_significant_genes_R_vs_NR.tsv\n")
cat("figures/03_deseq2_volcano_R_vs_NR.png\n")
