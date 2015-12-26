# source("https://bioconductor.org/biocLite.R")
# biocLite(c("snpStats", "broom", "sva"))

#Q1
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

snp3 = as.numeric(snpdata[,3])
table(snp3)
snp3[snp3==0] = NA
table(snp3)
# Fit a linear model and a logistic regression model to the data for the 3rd SNP. 
# What are the coefficients for the SNP variable? How are they interpreted? 
# (Hint: Don't forget to recode the 0 values to NA for the SNP data)

## logistic
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3) # Logistic Model = -0.16

## linear
lm3 <- lm(status ~ snp3)
tidy(lm3) #  Linear Model = -0.04 
# Both models are fit on the additive scale. 
# So in the linear model case, the coefficient is the decrease in probability associated with 
# each additional copy of the minor allele, in the logistic regression case, it is the 
# decrease in the log odds ratio associated with each additional copy of the minor allele.

#Q2
#The log odds is always more interpretable than a change in probability on the additive scale.
# >>> If you included more variables it would be possible to get negative estimates for the probability of being a case from the linear model, but this would be prevented with the logistic regression model.
# >> In this case the linear model fits the data much worse than the logistic regression model does. The coefficient is much smaller in absolute value so the effect is clearly not modeled well.
# > It is customary to use logistic regression for case-control data like those obtained from genome-wide association studies.

#Q3
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
status10 <- status[!is.na(snp10)]

snp10 = as.numeric(snpdata[,10])
table(snp10)
snp10[snp10==0] = NA
table(snp10)

# additive model
glm10 <- glm(status ~ snp10, family="binomial")
tidy(glm10)

# recessive model
snp10_rec = (snp10 == 3) 
glm10_rec = glm(status ~ snp10_rec, family="binomial")
tidy(glm10_rec)

table(glm10$fitted.values, status10)
table(glm10_rec$fitted.values, status10)

# > No, in all cases, the fitted values are near 0.5 and there are about an equal number of cases and controls in each group. This is true regardless of whether you fit a recessive or additive model.

#Q4
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
results_e <- list()
results_s <- list()
results_glms <- list()

for (i in 1:ncol(snpdata)){        
  
  snpi = as.numeric(snpdata[,i])
  snpi[snpi==0] = NA
  # additive model
  glmi <- glm(status ~ snpi, family="binomial")
  results_e[i] = tidy(glmi)[[2]][[2]]   
  results_s[i] = tidy(glmi)[[4]][[2]]  
  results_glms[[i]] <- glmi
}    

mean(unlist(results_s))
# 0.007155377
min(unlist(results_s))
#  -4.251469
max(unlist(results_s))
# 3.900891

#Q5
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc

# Fit an additive logistic regression model to each SNP and square the coefficients.
coeffsq <- unlist(results_s)**2

# What is the correlation with the results from using snp.rhs.tests and chi.squared ?
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10) # unadjusted model
slotNames(glm_all)
cs <- chi.squared(glm_all)
cor(cs, coeffsq)

# > 0.99. They are both testing for the same association using the same additive regression model 
# on the logistic scale but using slightly different tests.

#Q6
library(Biobase)
library(GenomicRanges)
library(SummarizedExperiment)
library(devtools)
library(Biobase)
library(limma)
library(edge)
library(genefilter)

con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)

edata_log2 <- log2(edata + 1)
# library(edgeR)
# y <- calcNormFactors(edata)
# plotMDS(edata)

tstats_obj = rowttests(as.matrix(edata_log2),pdata$population)
names(tstats_obj)
hist(tstats_obj$statistic,col=2)

## ------------------------------------------------------------------------
fstats_obj = rowFtests(as.matrix(edata_log2),pdata$population)
names(fstats_obj)
hist(fstats_obj$statistic,col=2)

summary(tstats_obj)
summary(fstats_obj)

# same pvalues, different stats

#Q7
library(DESeq2)
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
edata = edata[rowMeans(edata) > 100,]
fdata = fData(mp)

de = DESeqDataSetFromMatrix(edata, pdata, ~study) 
glm_all_nb = DESeq(de)
result_nb = results(glm_all_nb)

edata_log2 <- log2(edata + 1)
mod = model.matrix(~ pdata$study)
fit_limma = lmFit(edata_log2,mod)
ebayes_limma = eBayes(fit_limma)
head(ebayes_limma$t)

deseq_stats <- result_nb$stat
limma_stats <- ebayes_limma$t[,2]
cor(deseq_stats, limma_stats) # 0.9278701

par(mfrow=c(1,2))
plot(x = (deseq_stats + limma_stats) / 2, y = deseq_stats - limma_stats, cex = .2)
plot(deseq_stats, limma_stats, cex=0.2)
par(mfrow=c(1,1))  # more diff in small stats

#Q8 BH correction for multiple testing 
# using previous deseq and limma results
bh_adj_deseq <- p.adjust(result_nb$pvalue, method="BH")
sum(bh_adj_deseq < 0.05) # 1995

bh_adj_limma <- p.adjust(ebayes_limma$p.value[,2], method="BH")
sum(bh_adj_limma < 0.05) # 2807

#Q9
# >> Yes and no. It is surprising because there are a large fraction of the genes that are significantly different, but it isn't that surprising because we would expect that when comparing measurements from very different batches.

#Q10
# >> The p-values should have a spike near zero (the significant results) and be flat to the right hand side (the null results) so the distribution pushed toward one suggests conservative p-value calculation.
