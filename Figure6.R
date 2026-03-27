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
library(showtext)
library(gridExtra)
library(patchwork)
library(nlme)
library(emmeans)
library(cowplot)
library(openxlsx)
library(MESS)
library(PKNCA)
library(pracma)
library(lme4)
library(lmerTest)
library(UpSetR)
library(ggsci)

#define function to calculate NCA parameters
ncafunc <- function(df,food) {
  
  results = data.frame(group=NA,name=NA,Cmax=NA,Tmax=NA,AUC=NA,AUCinf=NA,Half=NA) #define a blank data frame to save results
  data <- df[which(df$group==food),]
  data$time <- as.numeric(data$time)
  data <- arrange(data,time)
  
  for (i in var_group_time) {
    
    tryCatch({
      
      data$Concentration = as.numeric(data[,i])
      data$Dose = 2
      
      d_conc = data[,c("ID","time","Dose","Concentration")]
      d_dose <- d_conc[d_conc$time == 0,]
      d_dose$time <- 0
      id <- unique(intersect(d_conc$ID,d_dose$ID)) #get ID existed in both con and dose dataset
      d_conc <- d_conc[which(d_conc$ID %in% id),]
      d_dose <- d_dose[which(d_dose$ID %in% id),]
      
      conci <- PKNCAconc(d_conc, Concentration~time|ID)
      dosei <- PKNCAdose(d_dose, Concentration~time|ID)
      obji  <- PKNCAdata(conci, dosei, options = list(auc.method = "linear"))
      
      res = pk.nca(obji)$result
      
      cmax_values <- res %>% filter((end == Inf) &(PPTESTCD == "cmax")) %>% select(PPORRES)
      tmax_values <- res %>% filter((end == Inf) &(PPTESTCD == "tmax")) %>% select(PPORRES) + 1 ### tmax produced here is only relative to the start time
      auc_values <- res %>% filter((end == 24) &(PPTESTCD == "auclast")) %>% select(PPORRES)
      aucinf_values <- res %>% filter((end == Inf) &(PPTESTCD == "aucinf.obs")) %>% select(PPORRES)
      halflife_values <- res %>% filter((end == Inf) &(PPTESTCD == "half.life")) %>% select(PPORRES)
      
      results[i,"group"] = food
      results[i,"name"] = i
      results[i,"Cmax"] = mean(cmax_values$PPORRES,na.rm = T)
      results[i,"Tmax"] = mean(tmax_values$PPORRES,na.rm = T)
      results[i,"AUC"] = mean(auc_values$PPORRES,na.rm=T)
      results[i,"AUCinf"] = mean(aucinf_values$PPORRES,na.rm=T)
      results[i,"Half"] = mean(halflife_values$PPORRES,na.rm=T)
      
    }, error = function(e) {
      
      # Handle error: skip this variable and move on
      print(paste("Error running model for", i, ":", e$message))
      
    })
  }
  
  results <- results[-1,] %>% as.data.frame()
  return(results)
}

#define function to calculate NCA parameters for urinary metabolites
ncafunc <- function(df,food) {
  
  results = data.frame(group=NA,name=NA,Cmax=NA,Tmax=NA,AUC=NA,AUCinf=NA,Half=NA) #define a blank data frame to save results
  data <- df[which(df$group==food),]
  data$time <- as.numeric(data$time)
  data <- arrange(data,time)
  
  for (i in var_group_time) {
    
    tryCatch({
      
      data$Concentration = as.numeric(data[,i])
      data$Dose = 2
      
      d_conc = data[,c("ID","time","Dose","Concentration")]
      d_dose <- d_conc[d_conc$time == 2,]
      d_dose$time <- 2
      id <- unique(intersect(d_conc$ID,d_dose$ID)) #get ID existed in both con and dose dataset
      d_conc <- d_conc[which(d_conc$ID %in% id),]
      d_dose <- d_dose[which(d_dose$ID %in% id),]
      
      conci <- PKNCAconc(d_conc, Concentration~time|ID)
      dosei <- PKNCAdose(d_dose, Concentration~time|ID)
      obji  <- PKNCAdata(conci, dosei, options = list(auc.method = "linear"))
      
      res = pk.nca(obji)$result
      
      cmax_values <- res %>% filter((end == Inf) &(PPTESTCD == "cmax")) %>% select(PPORRES)
      tmax_values <- res %>% filter((end == Inf) &(PPTESTCD == "tmax")) %>% select(PPORRES) + 1 ### tmax produced here is only relative to the start time
      auc_values <- res %>% filter((end == 24) &(PPTESTCD == "auclast")) %>% select(PPORRES)
      aucinf_values <- res %>% filter((end == Inf) &(PPTESTCD == "aucinf.obs")) %>% select(PPORRES)
      halflife_values <- res %>% filter((end == Inf) &(PPTESTCD == "half.life")) %>% select(PPORRES)
      
      results[i,"group"] = food
      results[i,"name"] = i
      results[i,"Cmax"] = mean(cmax_values$PPORRES,na.rm = T)
      results[i,"Tmax"] = mean(tmax_values$PPORRES,na.rm = T)
      results[i,"AUC"] = mean(auc_values$PPORRES,na.rm=T)
      results[i,"AUCinf"] = mean(aucinf_values$PPORRES,na.rm=T)
      results[i,"Half"] = mean(halflife_values$PPORRES,na.rm=T)
      
    }, error = function(e) {
      
      # Handle error: skip this variable and move on
      print(paste("Error running model for", i, ":", e$message))
      
    })
  }
  
  results <- results[-1,] %>% as.data.frame()
  return(results)
}

#read data
load("Res_no_filtering.Rdata")

use <- unique(c(bf_all_new,ck_all_new,ch_all_new,yg_all_new,cr_all_new,ww_all_new,ot_all_new,pt_all_new))

var_group_time <- use

df$time <- as.numeric(df$time)

table(df$time)
#1   2   3   4   5   6   7   8   9  10 
#116 118 113 117 110 114 108 112 103 114 

df$time <- df$time-1
table(df$time)
#0   1   2   3   4   5   6   7   8   9 
#116 118 113 117 110 114 108 112 103 114

df$time <- gsub("9","24",df$time)
table(df$time)

df$time <- as.numeric(df$time)

#calculate half, cmax, tmax across different food groups and metabolites
pk_bf <- ncafunc(df,food="Beef") 
pk_cs <- ncafunc(df,food="Cheese")
pk_ck <- ncafunc(df,food="Chicken")
pk_yg <- ncafunc(df,food="Yogurt")
pk_cn <- ncafunc(df,food="Corn")
pk_ot <- ncafunc(df,food="Oats")
pk_pt <- ncafunc(df,food="Potato")
pk_bd <- ncafunc(df,food="Whole Wheat Bread")

res_pk <- rbind(pk_bf[which(pk_bf$name %in% bf_all_new),],pk_cs[which(pk_bf$name %in% ch_all_new),],
                pk_ck[which(pk_bf$name %in% ck_all_new),],pk_cn[which(pk_bf$name %in% cr_all_new),],
                pk_ot[which(pk_bf$name %in% ot_all_new),],pk_pt[which(pk_bf$name %in% pt_all_new),],
                pk_bd[which(pk_bf$name %in% ww_all_new),],pk_yg[which(pk_bf$name %in% yg_all_new),]) %>% as.data.frame()

write.csv(res_pk, file="Res_PK_Plasma_Unfiltered.csv")