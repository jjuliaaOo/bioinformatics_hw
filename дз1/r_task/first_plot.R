genes <- c("BRCA1", "TP53", "EGFR")
expression <- c(12.5, 45.2, 30.1)
condition <- c("Control", "Treatment", "Treatment")

exp_data <- data.frame(
	genes = genes,
	expression = expression,
	condition = condition,
	stringsAsFactors = FALSE
)

str(exp_data)

png(file.path("r_task", "expression_plot.png"))
barplot(
	exp_data$expression,
	names.arg = exp_data$genes,
	xlab = "Genes",
	ylab = "Expression",
	main = "Gene expression"
)
dev.off()
