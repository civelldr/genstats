#(3)
library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
data(sample.ExpressionSet, package = "Biobase")
se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)
assay <- assay(se)
pheno <- colData(se)
rr <- rowRanges(se)


#(Q5)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
bottomly_p=pData(bot)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
bodymap_p=pData(bm)

#(Q7)
edata = exprs(bm)

