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
edata_bot <- log2(edata_bot + 1) # key
mod = model.matrix(~ pdata_bot$strain)
fit_limma = lmFit(edata_bot,mod)
ebayes_limma = eBayes(fit_limma)

#How many genes are differentially expressed at the 5% FDR level using Benjamini-Hochberg correction? 
bh_adj_limma <- p.adjust(ebayes_limma$p.value[,2], method="fdr")
sum(bh_adj_limma < 0.05) #223
sig <- data.frame(bh_adj_limma[bh_adj_limma < 0.05])
head(sig, 1) # ENSMUSG00000000402

# same thing btw, easier to use topTable
# limma_pvals_adj = topTable(ebayes_limma,number=dim(edata_bot)[1])$adj.P.Val
# sum(limma_pvals_adj <= 0.05)

#Q3 perform gene enrichment on the above significantly DE genes
library(goseq)
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
library("org.Mm.eg.db")
# get a list of genes that are below 5% FDR as a T/F array
genes <- as.integer(bh_adj_limma < 0.05)
names(genes) <- rownames(fdata_bot)
pwf <- nullp(genes,"mm9","ensGene")
GO.wall <- goseq(pwf,"mm9","ensGene")
head(GO.wall) # GO:0004888 

#Q4
GO.wall[1,6] # "transmembrane signaling receptor activity"

#Q5
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
edata_bot <- log2(edata_bot + 1) # key
mod <- model.matrix(~ pdata_bot$strain + pdata_bot$lane.number)
fit_limma = lmFit(edata_bot,mod)
ebayes_limma = eBayes(fit_limma)

### multiple testing correction 5% FDR level using Benjamini-Hochberg correction
###
bh_adj_limma <- p.adjust(ebayes_limma$p.value[,2], method="BH")
sum(bh_adj_limma < 0.05) #167

  # enrichment
genes_adj <- as.integer(bh_adj_limma < 0.05)
table(genes_adj)
length(genes_adj)
names(genes_adj) <- rownames(fdata_bot)
pwf_adj <- nullp(genes_adj,"mm9","ensGene")
GO.wall_adj <- goseq(pwf_adj,"mm9","ensGene")
head(GO.wall_adj ,10)
topTenCorrected <- GO.wall_adj [1:10,1]

### unadjusted for multiple correction
###
unadj_limma <- ebayes_limma$p.value[,2]
sum(unadj_limma < 0.05) # 639

  # enrichment
genes_raw <- as.integer(unadj_limma < 0.05)
table(genes_raw)
names(genes_raw) <- rownames(fdata_bot)
length(genes_raw)
pwf_raw <- nullp(genes_raw,"mm9","ensGene")
GO.wall_raw <- goseq(pwf_raw,"mm9","ensGene")
head(GO.wall_raw,10)
topTenUncorrected <- GO.wall_raw[1:10,1]

# how many of the top 10 GO categories are in common between gene enrichment performed on 
# adjusted and unadjusted pvalues
intersect(topTenCorrected, topTenUncorrected) # I get 3 but apparently it's 2. Hm. 
