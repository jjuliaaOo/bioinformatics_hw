# read CSV (можно запустить из терминала)
args <- commandArgs(trailingOnly = FALSE)
this_file <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
this_dir <- dirname(normalizePath(this_file))

csv_path <- file.path(this_dir, "sample_data.csv")
df <- read.csv(csv_path, stringsAsFactors = FALSE)

# mean
mean_score <- mean(df$Score)
cat("Mean Score:", mean_score, "\n")

# score
treatment_df <- df[df$Group == "Treatment", ]
max_treatment <- max(treatment_df$Score)
cat("Max Score in Treatment:", max_treatment, "\n")

# boxplot
out_png <- file.path(this_dir, "score_boxplot.png")
png(out_png)
boxplot(Score ~ Group,
        data = df,
        main = "Score Distribution by Group",
        xlab = "Group",
        ylab = "Score")
dev.off()

cat("Saved plot to:", out_png, "\n")
