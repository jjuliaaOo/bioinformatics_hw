setwd("C:/Users/trety/OneDrive/Рабочий стол/hw7")

library(limma)
library(ggplot2)
library(dplyr)
library(ggrepel)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

cat("Рабочая папка:\n")
print(getwd())

cat("\nФайлы GSE63885:\n")
print(list.files("data/GSE63885"))

# =========================
# 1. Загрузка данных
# =========================

exp_raw <- read.csv(
  "data/GSE63885/expression_for_limma2.csv",
  header = TRUE,
  check.names = TRUE
)

ann <- read.csv(
  "data/GSE63885/annotation_for_limma.csv",
  header = TRUE,
  row.names = 1,
  check.names = TRUE
)

cat("\nРазмер исходной expression table:\n")
print(dim(exp_raw))

cat("\nРазмер annotation table:\n")
print(dim(ann))

# Первая колонка expression table — названия генов
gene_col <- colnames(exp_raw)[1]

cat("\nКолонка с названиями генов:\n")
print(gene_col)

gene_symbols <- exp_raw[[gene_col]]

exp <- exp_raw[, -1]
rownames(exp) <- make.unique(as.character(gene_symbols))

# Приводим экспрессию к числам
exp <- as.data.frame(
  lapply(exp, function(x) as.numeric(as.character(x))),
  check.names = FALSE
)

rownames(exp) <- make.unique(as.character(gene_symbols))

cat("\nРазмер exp после подготовки:\n")
print(dim(exp))

cat("\nПервые названия генов:\n")
print(head(rownames(exp)))

# =========================
# 2. Поиск колонки pCR / pNC
# =========================

candidate_cols <- c()

for (col in colnames(ann)) {
  values <- trimws(as.character(ann[[col]]))
  if (any(values %in% c("pCR", "pNC"))) {
    candidate_cols <- c(candidate_cols, col)
  }
}

cat("\nКолонки, где есть pCR / pNC:\n")
print(candidate_cols)

if (length(candidate_cols) == 0) {
  stop("Не найдена колонка с pCR / pNC. Нужно проверить annotation_for_limma.csv.")
}

status_col <- candidate_cols[1]

cat("\nВыбранная колонка со статусом:\n")
print(status_col)

cat("\nЗначения выбранной колонки:\n")
print(table(ann[[status_col]], useNA = "ifany"))

# =========================
# 3. Согласование образцов
# =========================

common_samples <- intersect(colnames(exp), rownames(ann))

cat("\nКоличество общих образцов между exp и ann:\n")
print(length(common_samples))

if (length(common_samples) == 0) {
  if (ncol(exp) == nrow(ann)) {
    cat("\nНазвания образцов не совпали, но число образцов одинаковое. Используем порядок из файлов.\n")
    rownames(ann) <- colnames(exp)
    common_samples <- colnames(exp)
  } else {
    stop("Не удалось сопоставить образцы между expression table и annotation.")
  }
}

exp <- exp[, common_samples]
ann <- ann[common_samples, ]

# Оставляем только pCR и pNC
status_values <- trimws(as.character(ann[[status_col]]))
keep_samples <- status_values %in% c("pCR", "pNC")

exp <- exp[, keep_samples]
ann <- ann[keep_samples, ]

status_values <- trimws(as.character(ann[[status_col]]))

cat("\nРазмер exp после фильтрации pCR/pNC:\n")
print(dim(exp))

cat("\nРазмер ann после фильтрации pCR/pNC:\n")
print(dim(ann))

cat("\nГруппы после фильтрации:\n")
print(table(status_values))

if (ncol(exp) == 0 || nrow(ann) == 0) {
  stop("После фильтрации pCR/pNC не осталось образцов. Значит, выбрана неверная колонка.")
}

# =========================
# 4. LIMMA
# =========================

# pNC — базовая группа, pCR — сравниваемая
group <- factor(status_values, levels = c("pNC", "pCR"))

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "pCR_vs_pNC")

cat("\nDesign matrix:\n")
print(head(design))

fit <- lmFit(as.matrix(exp), design)
fit <- eBayes(fit)

limma_res <- topTable(
  fit,
  coef = "pCR_vs_pNC",
  adjust.method = "BH",
  number = Inf
)

limma_res$Gene.Symbol <- rownames(limma_res)

# =========================
# 5. Значимость при logFC 1 / 2 / 3
# Используем обычный P.Value, как требуется в задании
# =========================

limma_res$significant_logFC_1 <- ifelse(
  limma_res$P.Value < 0.05 & abs(limma_res$logFC) >= 1,
  "yes",
  "no"
)

limma_res$significant_logFC_2 <- ifelse(
  limma_res$P.Value < 0.05 & abs(limma_res$logFC) >= 2,
  "yes",
  "no"
)

limma_res$significant_logFC_3 <- ifelse(
  limma_res$P.Value < 0.05 & abs(limma_res$logFC) >= 3,
  "yes",
  "no"
)

summary_logfc <- data.frame(
  logFC_threshold = c(1, 2, 3),
  significant_genes = c(
    sum(limma_res$significant_logFC_1 == "yes"),
    sum(limma_res$significant_logFC_2 == "yes"),
    sum(limma_res$significant_logFC_3 == "yes")
  )
)

cat("\nКоличество значимых генов при разных порогах logFC:\n")
print(summary_logfc)

cat("\nПервые 10 строк LIMMA результата:\n")
print(head(limma_res, 10))


# =========================
# 6. Сохранение таблиц
# =========================

write.table(
  limma_res,
  file = "results/limma_results_all.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  summary_logfc,
  file = "results/limma_summary_logFC_thresholds.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  limma_res[limma_res$significant_logFC_1 == "yes", ],
  file = "results/limma_significant_logFC_1.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  limma_res[limma_res$significant_logFC_2 == "yes", ],
  file = "results/limma_significant_logFC_2.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  limma_res[limma_res$significant_logFC_3 == "yes", ],
  file = "results/limma_significant_logFC_3.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


# =========================
# 7. Volcano plots
# Строим по обычному P.Value
# =========================

make_volcano <- function(df, threshold) {
  plot_df <- df
  
  plot_df$minus_log10_p <- -log10(plot_df$P.Value)
  
  plot_df$category <- "not significant"
  
  plot_df$category[
    plot_df$P.Value < 0.05 &
      plot_df$logFC >= threshold
  ] <- "up in pCR"
  
  plot_df$category[
    plot_df$P.Value < 0.05 &
      plot_df$logFC <= -threshold
  ] <- "down in pCR"
  
  top_genes <- plot_df %>%
    filter(P.Value < 0.05) %>%
    arrange(P.Value) %>%
    head(10)
  
  p <- ggplot(plot_df, aes(x = logFC, y = minus_log10_p, color = category)) +
    geom_point(alpha = 0.6, size = 1.4) +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = top_genes,
      aes(label = Gene.Symbol),
      size = 3,
      max.overlaps = 20
    ) +
    theme_minimal() +
    labs(
      title = paste0("LIMMA volcano plot, logFC threshold = ", threshold),
      x = "logFC",
      y = "-log10(p-value)",
      color = "Category"
    )
  
  return(p)
}

volcano_1 <- make_volcano(limma_res, 1)
volcano_2 <- make_volcano(limma_res, 2)
volcano_3 <- make_volcano(limma_res, 3)

ggsave(
  "figures/04_limma_volcano_logFC1.png",
  volcano_1,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/05_limma_volcano_logFC2.png",
  volcano_2,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  "figures/06_limma_volcano_logFC3.png",
  volcano_3,
  width = 8,
  height = 6,
  dpi = 300
)

print(volcano_1)


# =========================
# 8. График количества значимых генов
# =========================

barplot_logfc <- ggplot(summary_logfc, aes(x = factor(logFC_threshold), y = significant_genes)) +
  geom_col() +
  theme_minimal() +
  labs(
    title = "Number of significant genes for different logFC thresholds",
    x = "logFC threshold",
    y = "Number of significant genes"
  )

ggsave(
  "figures/07_limma_significant_counts_by_logFC.png",
  barplot_logfc,
  width = 7,
  height = 5,
  dpi = 300
)

cat("\nГотово. Сохранены файлы:\n")
cat("results/limma_results_all.tsv\n")
cat("results/limma_summary_logFC_thresholds.tsv\n")
cat("results/limma_significant_logFC_1.tsv\n")
cat("results/limma_significant_logFC_2.tsv\n")
cat("results/limma_significant_logFC_3.tsv\n")
cat("figures/04_limma_volcano_logFC1.png\n")
cat("figures/05_limma_volcano_logFC2.png\n")
cat("figures/06_limma_volcano_logFC3.png\n")
cat("figures/07_limma_significant_counts_by_logFC.png\n")

