#Call all the librery
library(TCGAbiolinks)
library(DESeq2)
library("parallel")
library("BiocParallel")

#set Directory
setwd("~/Documents/Covid")

#Call the datafile
data = read.table(file = "GSE235262_raw_counts.txt", sep = "\t", header = TRUE, row.names = 1)
length(unique(data$gene_name))
data1 = data[,2:267]
sample_name = colnames(data1)
colData = as.data.frame(sample_name)
variants = c(rep(paste0("BA"),80), rep(paste0("Delta"), 31), rep(paste0("OMEX"), 48), rep(paste0("PTR"), 107))
colData$variant = variants
cts <- as.matrix(data1)
colData$condition <- factor(colData$variant)
rownames(colData) = colData$sample_name
all(rownames(colData) %in% colnames(cts))
raw_data = round(data1, digits = 0)
cts <- as.matrix(raw_data)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = colData,
                              design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 30
table(keep)
dds <- dds[keep,]
dds$definition
dds$condition
dds$condition <- relevel(dds$condition, ref = "BA")
register(MulticoreParam(14))
dds1 <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(14))
resultsNames(dds1)
res <- results(dds1, name="condition_Delta_vs_BA")
summary(res)


sum(res$padj < 0.1 & abs(res$log2FoldChange)>=1, na.rm=TRUE)
deg_BA_Delta = res$padj < 0.1 & abs(res$log2FoldChange)>=1
ids_deg_BA_Delta = which(deg_BA_Delta == "TRUE")
write.table(res[ids_deg_BA_Delta,], file = "covid_deg_BA_Delta.txt", sep = "\t", col.names = NA)
write.table(ids_deg_BA_Delta, file = "ids_covid_deg_BA_Delta.txt", sep = "\t", col.names = NA)
degn = rownames(res[ids_deg_BA_Delta,])
write.table(degn, file = "covid_degnames_ba_delta.txt", sep = "\t", col.names = NA)
degn_expr_ba_delta = raw_data[degn,]
t_degn_expr_ba_delta = t(degn_expr_ba_delta)








res <- results(dds1, name="condition_OMEX_vs_BA")
summary(res)

sum(res$padj < 0.1 & abs(res$log2FoldChange)>=1, na.rm=TRUE)
deg_BA_OMEX = res$padj < 0.1 & abs(res$log2FoldChange)>=1
ids_deg_BA_OMEX = which(deg_BA_OMEX == "TRUE")
write.table(res[ids_deg_BA_OMEX,], file = "covid_deg_BA_OMEX.txt", sep = "\t", col.names = NA)
write.table(ids_deg_BA_OMEX, file = "ids_covid_deg_BA_OMEX.txt", sep = "\t", col.names = NA)
degn = rownames(res[ids_deg_BA_OMEX,])
write.table(degn, file = "covid_degnames_ba_omex.txt", sep = "\t", col.names = NA)
degn_expr_ba_omex = raw_data[degn,]
t_degn_expr_ba_omex = t(degn_expr_ba_omex)

res <- results(dds1, name="condition_PTR_vs_BA")
summary(res)

sum(res$padj < 0.1 & abs(res$log2FoldChange)>=1, na.rm=TRUE)
deg_BA_PTR = res$padj < 0.1 & abs(res$log2FoldChange)>=1
ids_deg_BA_PTR = which(deg_BA_Delta == "TRUE")
write.table(res[ids_deg_BA_PTR,], file = "covid_deg_BA_PTR.txt", sep = "\t", col.names = NA)
write.table(ids_deg_BA_PTR, file = "ids_covid_deg_BA_PTR.txt", sep = "\t", col.names = NA)
degn = rownames(res[ids_deg_BA_PTR,])
write.table(degn, file = "covid_degnames_ba_ptr.txt", sep = "\t", col.names = NA)
degn_expr_ba_ptr = raw_data[degn,]
t_degn_expr_ba_ptr = t(degn_expr_ba_ptr)


dds$condition <- relevel(dds$condition, ref = "OMEX")
register(MulticoreParam(14))
dds1 <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(14))
resultsNames(dds1)
res <- results(dds1, name="condition_Delta_vs_OMEX")
summary(res)
sum(res$padj < 0.1 & abs(res$log2FoldChange)>=1, na.rm=TRUE)
deg_OMEX_Delta = res$padj < 0.1 & abs(res$log2FoldChange)>=1
ids_deg_OMEX_Delta = which(deg_OMEX_Delta == "TRUE")
write.table(res[ids_deg_OMEX_Delta,], file = "covid_deg_OMEX_Delta.txt", sep = "\t", col.names = NA)
write.table(ids_deg_OMEX_Delta, file = "ids_covid_deg_OMEX_Delta.txt", sep = "\t", col.names = NA)
degn = rownames(res[ids_deg_OMEX_Delta,])
write.table(degn, file = "covid_degnames_omex_delta.txt", sep = "\t", col.names = NA)
degn_expr_omex_delta = raw_data[degn,]
t_degn_expr_omex_delta = t(degn_expr_omex_delta)

res <- results(dds1, name="condition_PTR_vs_OMEX")
summary(res)
sum(res$padj < 0.1 & abs(res$log2FoldChange)>=1, na.rm=TRUE)
deg_OMEX_PTR = res$padj < 0.1 & abs(res$log2FoldChange)>=1
ids_deg_OMEX_PTR = which(deg_OMEX_PTR == "TRUE")
write.table(res[ids_deg_OMEX_PTR,], file = "covid_deg_OMEX_PTR.txt", sep = "\t", col.names = NA)
write.table(ids_deg_OMEX_PTR, file = "ids_covid_deg_OMEX_PTR.txt", sep = "\t", col.names = NA)
degn = rownames(res[ids_deg_OMEX_PTR,])
write.table(degn, file = "covid_degnames_omex_ptr.txt", sep = "\t", col.names = NA)
degn_expr_omex_ptr = raw_data[degn,]
t_degn_expr_omex_ptr = t(degn_expr_omex_ptr)

dds$condition <- relevel(dds$condition, ref = "PTR")
dds1 <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(14))
resultsNames(dds1)

res <- results(dds1, name="condition_Delta_vs_PTR")
sum(res$padj < 0.1 & abs(res$log2FoldChange)>=1, na.rm=TRUE)
deg_PTR_Delta = res$padj < 0.1 & abs(res$log2FoldChange)>=1
ids_deg_PTR_Delta = which(deg_PTR_Delta == "TRUE")
write.table(res[ids_deg_PTR_Delta,], file = "covid_deg_PTR_Delta.txt", sep = "\t", col.names = NA)
write.table(ids_deg_PTR_Delta, file = "ids_covid_deg_PTR_Delta.txt", sep = "\t", col.names = NA)
degn = rownames(res[ids_deg_PTR_Delta,])
write.table(degn, file = "covid_degnames_ptr_delta.txt", sep = "\t", col.names = NA)
degn_expr_ptr_delta = raw_data[degn,]
t_degn_expr_ptr_delta = t(degn_expr_ptr_delta)


degn_p1=rownames(res[ids_deg_BA_Delta,])
degn_p2=rownames(res[ids_deg_BA_OMEX,])
degn_p3=rownames(res[ids_deg_BA_PTR,])
degn_p4=rownames(res[ids_deg_OMEX_Delta,])
degn_p5=rownames(res[ids_deg_OMEX_PTR,])
degn_p6=rownames(res[ids_deg_PTR_Delta,])

un_degn=Reduce(union, list(degn_p1, degn_p2, degn_p3, degn_p4, degn_p5, degn_p6))
class(un_degn)
degn_expr_un_degn = raw_data[un_degn, ]
class(degn_expr_un_degn)
View(degn_expr_un_degn)

dds <- DESeqDataSetFromMatrix(countData = raw_data, 
                              colData = colData, 
                              design = ~condition)


vst_data <- vst(dds, blind=FALSE)
View(vst_data)
degn_expr_un_degn_vst = vst_data[un_degn,]