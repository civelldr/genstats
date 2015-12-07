# Q1
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
# What percentage of variation is explained by the 1st 
# principal component in the data set if you:
# a. do no transformations
# b. do log2(data + 1) transformation
# c. do log2(data + 1) and subtract row means

pc1 = prcomp(edata)
summary(pc1)  # looks like ~ 89% from PC1

edata2 <- log2(edata + 1)
pc2 <- prcomp(edata2)  
summary(pc2) # 97%

edata3 <-  log2(edata + 1)
edata_centered = edata3 - rowMeans(edata) # there must be something wrong with edata_centered
pc3 <- prcomp(edata_centered)
summary(pc3) # weird, I'm seeing 1 for this when the answer looks like it should be 0.35

#Q2
set.seed(333)
kmeans1 = kmeans(t(edata_centered),centers=2)  # forgot to transpose...?
svd1 = svd(edata_centered)
names(svd1) # d u v, d diagnal, v
cor(kmeans1$cluster, svd1$v[,1])  # who knows, I'm getting 0.81
## 0.87 is right

#Q3
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
edata = as.matrix(edata)
lm1 = lm(edata[1,] ~ as.factor(pdata_bm$num.tech.reps))
plot(edata[1,], as.factor(pdata_bm$num.tech.reps))
# There are very few samples with more than 2 replicates so the estimates for those values will not be very good.

#Q4 

lm2 <- lm(edata[1,] ~ pdata_bm$age + pdata_bm$gender)
summary(lm2) 
# -23.91. This coefficient means that there is an average decrease of 23.91 in the count variable per year within each gender.


#Q5
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata = log2(as.matrix(edata) + 1)
mod = model.matrix(~ pdata$population)
fit = lm.fit(mod,t(edata))
names(fit)
dim(fit$residuals)
dim(fit$effects)
dim(fit$coefficients)

# Residual matrix: 129 x 52580
# Effects matrix: 129 x 52580
# Coefficients matrix: 2 x 52580

#Q6 what is "effects"
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata <- log2(edata + 1)
#The estimated fitted values for all samples for each gene, with the values for each gene stored in the columns of the matrix.			

#Q7
library(limma)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)

pdata_bm$age  # 11-13 samples don't have age
edata_eval <- edata[,c(1:10, 14:19)]

mod_adj <- model.matrix(~ pdata_bm$age)
fit_limma = lmFit(edata_eval,mod_adj)
fit_limma$coefficients[1000] # 2469.87

# wrong: -27.61. The model doesn't fit well since there appears to be a non-linear trend in the data.
# 27.61. The model doesn't fit well since there appears to be a non-linear trend in the data.			
# -27.61. The model doesn't fit well since there are two large outlying values and the rest of the values are near zero.			
# wrong: 2469.87. The model doesn't fit well since there are two large outlying values and the rest of the values are near zero.	Inorrect	0.00	
# wrong: -23.25. The model doesn't fit well since there are two large outlying values and the rest of the values are near zero.
# -27.61. The model fits well since there seems to be a flat trend in the counts.

#Q8
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
mod_adj <- model.matrix(~ pdata_bm$age + pdata_bm$tissue.type)
fit_limma = lmFit(edata_eval,mod_adj)

# Since tissue.type is a factor variable with many levels, this model has more coefficients to estimate per gene (18) than data points per gene (16)

#Q9 
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
table(pdata$population, pdata$study)
# > The effects are difficult to distinguish because the study variable and population variable are perfectly correlated.

#Q10
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
set.seed(33353)

edata <- log2(edata + 1)
pdata_bm$age  # 11-13 samples don't have age
edata <- edata[,c(1:10, 14:19)]
edata <- edata[rowMeans(edata) >= 1,]
pdata_bm <- pdata_bm[!is.na(pdata_bm$age),]

mod = model.matrix(~age,data=pdata_bm)
mod0 = model.matrix(~1, data=pdata_bm)
sva1 = sva(edata,mod,mod0,n.sv=1) # comp with the null model, how many surrogate batch effects

names(sva1) 
dim(sva1$sv) # new covariates that are created that are potentially 
# the batch effects.  there's just 1 here 

cor(sva1$sv,pdata_bm$age)
# [1] 0.2028679
summary(lm(sva1$sv ~ pdata_bm$age))
summary(lm(sva1$sv ~ pdata_bm$gender)) # has the smallest pvalue, but not sig..
summary(lm(sva1$sv ~ pdata_bm$race))

boxplot(sva1$sv ~ pdata_bm$gender)
# > Correlation with age: 0.20 More highly correlated with gender.
