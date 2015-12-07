source("http://www.bioconductor.org/biocLite.R")
biocLite(c("SummarizedExperiment"))
tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
colramp = colorRampPalette(c(3,"white",2))(9)

#(3)
library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
data(sample.ExpressionSet, package = "Biobase")
se = makeSummarizedExperimentFromExpressionSet(sample.ExpressionSet)
assay <- assay(se)
pheno <- colData(se)
f <- rowData(se) # defunct!
rr <- rowRanges(se)

# The covariates in the Bottomly data set (experiment number, lane number) are balanced with respect to strain. The covariates 
# in the Bodymap data set (gender, age, number of technical replicates) are not balanced with respect to tissue.

#(5)
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

#(7)
edata = exprs(bm)
row_sums = rowSums(edata)

index = which(rank(-row_sums) <= 500 ) 

# or..
# edata = edata[order(-row_sums),]
# index = 1:500

heatmap(edata[index,],Rowv=NA,Colv=NA)

# yes the samples are next to each other

#(8)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
pdata = pData(bm)
edata = exprs(bm)

library(DESeq2)
library(limma)
s1s2 <- edata[,1:2]

log2_edata <- log2(edata + 1)
rlog_edata <- rlog(edata)

plot(x = (log2_edata[,1] + log2_edata[,2]) / 2, y = log2_edata[,1] - log2_edata[,2], cex = .2)
plot(x = (rlog_edata[,1] + rlog_edata[,2]) / 2, y = rlog_edata[,1] - rlog_edata[,2], col="red", cex = 0.2)
points(x = (rlog_edata[,1] + rlog_edata[,2]) / 2, y = rlog_edata[,1] - rlog_edata[,2], col="red", cex = 0.2)
# rlog on it's own
plot(x = (rlog_edata[,1] + rlog_edata[,2]) / 2, y = rlog_edata[,1] - rlog_edata[,2], col="red", cex = 0.2)

#(9)
library(dendextend) # biocLite(c("dendextend"))
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

table(pdata$study) 

groupCodes <- pdata$study
colorCodes <- c(Montgomery="red", Pickrell="blue")

pdata$sid <- ifelse(pdata$study == "Montgomery", "M", "P")
pdata$sample.id <- paste(pdata$sid, pdata$sample.id, sep="_")
colnames(edata) <- pdata$sample.id

# c1 with no changes to the data 
dist1 = dist(t(edata))
dim(edata)
hclust1 <- hclust(dist1)
plot(hclust1, hang=-1)
dend <- as.dendrogram(hclust1)

labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
plot(dend)

# c2 after filtering all the genes with rowMeans <
edata2 = edata[rowMeans(edata) > 100,]
dim(edata2)
dist2 = dist(t(edata2))
hclust2 <- hclust(dist2)
plot(hclust2, hang=-1)
dend2 <- as.dendrogram(hclust2)
labels_colors(dend2) <- colorCodes[groupCodes][order.dendrogram(dend2)]
plot(dend2) # identical

#c3 after taking the log2 transform of the data without filtering
edata3 <- log2(edata + 1)
dist3 = dist(t(edata3))
hclust3 <- hclust(dist3)
plot(hclust3, hang=-1)
dend3 <- as.dendrogram(hclust3)
labels_colors(dend3) <- colorCodes[groupCodes][order.dendrogram(dend3)]
plot(dend3) # identical

#(10)
library(cluster)
library(fpc)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

groupCodes <- pdata$study
colorCodes <- c(Montgomery="red", Pickrell="blue")

pdata$sid <- ifelse(pdata$study == "Montgomery", "M", "P")
pdata$sample.id <- paste(pdata$sid, pdata$sample.id, sep="_")
colnames(edata) <- pdata$sample.id


edata = edata[rowMeans(edata) > 100,]
edata <- log2(edata + 1)
dist1 = dist(t(edata))

# set.seed(1235)
kmeans1 = kmeans(t(edata),centers=2, set.seed(1235), nstart = 10, iter.max = 100)

df_kmeans <- data.frame(kmeans1$cluster)
df_kmeans$sample.id <- rownames(df_kmeans)
rownames(df_kmeans) <- NULL
km <- merge(pdata, df_kmeans)
table(km$study, km$kmeans1.cluster)

hclust1 <- hclust(dist1)
ct <- cutree(hclust1, 2)

test <- merge(df_kmeans, df_ct)

df_ct <- data.frame(ct)
df_ct$sample.id <- rownames(df_ct)
rownames(df_ct) <- NULL
ct <- merge(pdata, df_ct)
table(ct$study, ct$ct)




# kmeans1 = kmeans(edata,centers=2, set.seed(1235), nstart = 100, iter.max = 100)
# plotcluster(t(edata), kmeans1$cluster)

names(kmeans1)

## ------------------------------------------------------------------------
matplot(t(kmeans1$centers),col=1:2,type="l",lwd=3)

## ------------------------------------------------------------------------
kmeans1$cluster[1:10]
table(kmeans1$cluster)

## ------------------------------------------------------------------------
heatmap(as.matrix(dist1)[order(kmeans1$cluster),],col=colramp,Colv=NA,Rowv=NA)

## ------------------------------------------------------------------------
kmeans2 = kmeans(edata,centers=3, nstart = 40)
table(t(kmeans1$cluster),kmeans2$cluster)

