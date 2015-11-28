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
edata_centered = edata3 - rowMeans(edata)
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
## wrong: The data are right skewed is wrong
## wrong: The difference between 2 and 5 technical replicates is not the same as the difference between 5 and 6 technical replicates.
# There may be different numbers of counts for different numbers of technical replicates.			
# There is only one data point with a value of 6 so it is likely that the estimated value for that number of technical replicates is highly variable.			
# wrong: The variable num.tech.reps is a continuous variable. 2x

#Q4 

lm2 <- lm(edata[1,] ~ pdata$age + pdata$gender)
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
# wrong: The model coefficients for all samples for each gene, with the values for each gene stored in the columns of the matrix 
# wrong: The model coefficients for all samples for each gene, with the values for each gene stored in the rows of the matrix.
# 2x The estimated fitted values for all samples for each gene, with the values for each gene stored in the columns of the matrix.			
# The model coefficients for all samples for each gene, with the values for each gene stored in the rows of the matrix.
# The shrunken estimated fitted values for all samples for each gene, with the values for each gene stored in the columns of the matrix.

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

#Q8
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata = exprs(bm)
pdata_bm=pData(bm)
mod_adj <- model.matrix(~ pdata_bm$age + pdata_bm$tissue.type)
fit_limma = lmFit(edata_eval,mod_adj)

# either this: tissue.type has 18 levels but there are only 16 data points per gene, so this model can't fit a unique solution.
# or this: Since tissue.type is a factor variable with many levels, this model has more coefficients to estimate per gene (18) than data points per gene (16)
# wrong: Normally this model wouldn't fit well since we have more coefficients (18) than data points per gene (16). But since we have so many genes to estimate with, the model fits well.
# The model doesn't fit well since age should be treated as a factor variable.
# The model doesn't fit well because there are a large number of outliers for the white blood cell tissue.

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
