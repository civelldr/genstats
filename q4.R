source("http://www.bioconductor.org/biocLite.R")
biocLite(c("goseq"))
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
colramp = colorRampPalette(c(3,"white",2))(9)

library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
library(goseq)
#Q1
# All reads were realigned to the NCBI m37 version of the mouse genome
# corresponds to NCBI Build 37 -> mm9 mouse

#Q2
library(Biobase)
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata_bot = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata_bot) > 5]
edata_bot = edata_bot[rowMeans(edata_bot) > 5, ]
mod = model.matrix(~ pdata_bot$strain)
fit_limma = lmFit(edata_bot,mod)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma$t)
sum(ebayes_limma$p.value[,2] < 0.05)
#How many genes are differentially expressed at the 5% FDR level using Benjamini-Hochberg correction? 
bh_adj_limma <- p.adjust(ebayes_limma$p.value[,2], method="BH")
sum(bh_adj_limma <= 0.05) 

limma_pvals_adj = topTable(ebayes_limma,number=dim(edata_bot)[1])$adj.P.Val
sum(limma_pvals_adj <= 0.05) 

hist(limma_pvals_adj,col=2)
quantile(limma_pvals_adj)
