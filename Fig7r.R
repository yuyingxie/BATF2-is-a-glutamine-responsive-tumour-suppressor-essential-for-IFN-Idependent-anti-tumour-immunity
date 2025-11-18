library(gee)
library(geepack)
library(survival)
library(survminer)

#################
# Fig 7 r, GEE test
#################
dat = read.csv('dat_BATF2.csv', header = TRUE)


model1 = geeglm(Total ~ as.factor(Genotype)  + Week, id = ID, data = dat,  corstr = "ar1")
summary(model1)
anova(model1)
