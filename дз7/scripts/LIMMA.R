BiocManager::install("limma")


library(limma)

geo_id <- 'GSE63885'
exp = read.csv(file = '/projects/mipt_dbmp_biotechnology/rna_seq_diff_exp/GSE63885/expression_for_limma2.csv', header = TRUE, row.names = 'Gene.Symbol')
ann = read.csv(file = '/projects/mipt_dbmp_biotechnology/rna_seq_diff_exp/GSE63885/annotation_for_limma.csv', header = TRUE, row.names = 'X')


exp_set <- ExpressionSet(assayData = as.matrix(exp), phenoData = AnnotatedDataFrame(ann))

slope <- factor(ann$clinical.status.post.1st.line.chemotherapy..cr...complete.response..pr...partial.response..sd...stable.disease..p...progression..ch1, levels = c("pCR",
                                                                                                                                                                     "pNC"), labels = c(1, 0))

pCR=as.integer(as.vector(slope))

iterse <- rep(1, length(slope))

design <- cbind(npCR=iterse,pCR=pCR)

fit <- lmFit(exp_set, design)
fit <- eBayes(fit)
top <-topTable(fit, coef="pCR", adjust="BH", n = Inf)