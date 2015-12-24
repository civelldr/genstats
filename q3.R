source("https://bioconductor.org/biocLite.R")
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
tidy(glm3)

## linear
lm3 <- lm(status ~ snp3)
tidy(lm3)
# Linear Model = -0.04 Logistic Model = -0.16 Both models are fit on the additive scale. 
# So in the linear model case, the coefficient is the decrease in probability associated with 
# each additional copy of the minor allele, in the logistic regression case, it is the 
# decrease in the log odds ratio associated with each additional copy of the minor allele.

#Q2
#The log odds is always more interpretable than a change in probability on the additive scale.
#If you included more variables it would be possible to get negative estimates for the probability of being a case from the linear model, but this would be prevented with the logistic regression model.
# > In this case the linear model fits the data much worse than the logistic regression model does. The coefficient is much smaller in absolute value so the effect is clearly not modeled well.
#...It is customary to use logistic regression for case-control data like those obtained from genome-wide association studies.

#Q3
library(snpStats)
library(broom)
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
status10 <- status[snp10 > 0]

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

table(glm10$fitted.values, status)
table(glm10_rec$fitted.values, status)

# The recessive model shows a strong effect, but the additive model shows no difference so the recessive model is better.
# The recessive model fits much better since there are more parameters to fit and the effect size is so large.
# The recessive model fits much better since it appears that once you aggregate the heterozygotes and homozygous minor alleles, there is a bigger difference in the proportion of cases and controls.
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
coeffsq <- unlist(results_e)**2

# What is the correlation with the results from using snp.rhs.tests and chi.squared ?
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10) # unadjusted model
slotNames(glm_all)
cs <- chi.squared(glm_all)






