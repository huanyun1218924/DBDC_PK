#load packages
library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)
library(ggplot2)
library(ggvenn)
library(ggalt)
library(gridExtra)
library(patchwork)
library(lme4)       #pacakge used for linear mixed model!
library(emmeans)

#read data
load("plasma_sample.RData")

#check data
dim(plasma_sample)
#1125 995

#check covariates
covar

#check metabolites
metvar

#check group and time
table(plasma_sample$group)

table(plasma_sample$time)

#run linear mixed model to identify in response to food intake, adjusted for age, sex,race, and BMI
results <- list()  
results_l <- list()
results_b <- list()
results_p <- list()
results_e <- list()
results_c <- list()

for (i in 1:length(metvar) {
  
  metabolite <- metvar[i]
  data <- plasma_sample[covar]
  data$Concentration <- as.numeric(plasma_sample[[metabolite]])  
  data$group <- as.factor(data$group)
  data$time <- as.factor(data$time)
  
  #run linear mixed model
  fit1 <- lmer(log(Concentration) ~ group*time + ips_age + ps_sex + race + ips_bmi + (1 | ID), data = data)
  fit2 <- lmer(log(Concentration) ~ group + time + ips_age + ps_sex + race + ips_bmi + (1 | ID), data = data)
  
  #save the results
  results[[i]] <- summary(fit1)
  results_l[[i]] <- anova(fit1, fit2)[2,"Pr(>Chisq)"]
  results_b[[i]] <- summary(emmeans(fit1, pairwise ~ group | time))$emmeans
  results_p[[i]] <- summary(emmeans(fit1, pairwise ~ group | time))$contrasts
  results_e[[i]] <- summary(emmeans(fit1, pairwise ~ time | group))$emmeans
  results_c[[i]] <- summary(emmeans(fit1, pairwise ~ time | group))$contrasts
  
}

#save results
save(results,results_b,results_p,results_e,results_c,results_l,file="Res_LMM_Plasma.RData")
