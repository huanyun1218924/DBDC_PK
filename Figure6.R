#load packages
library(data.table)
library(dplyr)
library(PKNCA)       #non-compartmental model kinetic analysis

#define function to calculate kinetic parameters for plasma metabolites
ncafunc <- function(plasma_sample,food) {
  
  results = data.frame(group=NA,name=NA,Cmax=NA,Tmax=NA,AUC=NA,AUCinf=NA,Half=NA) #define a blank data frame to save results
  data <- plasma_sample[which(plasma_sample$group==food),]
  data$time <- as.numeric(data$time)
  data <- arrange(data,time)
  
  for (i in metvar) {
    
    tryCatch({
      
      data$Concentration = as.numeric(data[,i])
      data$Dose = 2
      
      d_conc = data[,c("ID","time","Dose","Concentration")]
      d_dose <- d_conc[d_conc$time == 0,]
      d_dose$time <- 0
      id <- unique(intersect(d_conc$ID,d_dose$ID))                                              #get ID existed in both con and dose dataset
      d_conc <- d_conc[which(d_conc$ID %in% id),]
      d_dose <- d_dose[which(d_dose$ID %in% id),]
      
      conci <- PKNCAconc(d_conc, Concentration~time|ID)
      dosei <- PKNCAdose(d_dose, Concentration~time|ID)
      obji  <- PKNCAdata(conci, dosei, options = list(auc.method = "linear"))
      
      res = pk.nca(obji)$result
      
      cmax_values <- res %>% filter((end == Inf) &(PPTESTCD == "cmax")) %>% select(PPORRES)
      tmax_values <- res %>% filter((end == Inf) &(PPTESTCD == "tmax")) %>% select(PPORRES) + 1 #note:tmax produced here is only relative to the start time
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

#define function to calculate kinetic parameters for urinary metabolites
ncafunc2 <- function(urine_sample,food) {
  
  results = data.frame(group=NA,name=NA,Cmax=NA,Tmax=NA,AUC=NA,AUCinf=NA,Half=NA) #define a blank data frame to save results
  data <- urine_sample[which(urine_sample$group==food),]
  data$time <- as.numeric(data$time)
  data <- arrange(data,time)
  
  for (i in metvar) {
    
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

#--------------------------------------------------------------------------------------------
#
#           part 1: calculate kinetic parameters of annotated plasma metabolites 
#
#--------------------------------------------------------------------------------------------
#read data
load("plasma_sample.RData")

#transform the timepoint data into numeric type for analysis
plasma_sample$time <- as.numeric(plasma_sample$time)
table(plasma_sample$time)
#1   2   3   4   5   6   7   8   9  10 
#116 118 113 117 110 114 108 112 103 114 

plasma_sample$time <- plasma_sample$time-1
table(plasma_sample$time)
#0   1   2   3   4   5   6   7   8   9 
#116 118 113 117 110 114 108 112 103 114

plasma_sample$time <- gsub("9","24",plasma_sample$time)
plasma_sample$time <- as.numeric(plasma_sample$time) 
table(plasma_sample$time)
#0   1   2   3   4   5   6   7   8  24 
#116 118 113 117 110 114 108 112 103 114

#calculate half, cmax, tmax across each metabolites by food
pk_bf <- ncafunc(plasma_sample,food="Beef") 
pk_cs <- ncafunc(plasma_sample,food="Cheese")
pk_ck <- ncafunc(plasma_sample,food="Chicken")
pk_yg <- ncafunc(plasma_sample,food="Yogurt")
pk_cn <- ncafunc(plasma_sample,food="Corn")
pk_ot <- ncafunc(plasma_sample,food="Oats")
pk_pt <- ncafunc(plasma_sample,food="Potato")
pk_bd <- ncafunc(plasma_sample,food="Whole Wheat Bread")

#save results
save(pk_bf,pk_ck,pk_ch,pk_yg,pk_cn,pk_ot,pk_pt,pk_bd,file="Results_Plasma_NCA.RData")

#--------------------------------------------------------------------------------------------
#
#           part 2: calculate kinetic parameters of annotated urinary metabolites 
#
#--------------------------------------------------------------------------------------------
#read data
load("urine_sample.RData")

#transform the timepoint data into numeric type for analysis
urine_sample$time <- as.numeric(urine_sample$time)
table(urine_sample$time)
#1   2   3   4   5   6 
#114 114 107 110 105 112 

urine_sample$time <- urine_sample$time-1
table(urine_sample$time)
#0   1   2   3   4   5 
#114 114 107 110 105 112 

urine_sample$time <- gsub("4","8",urine_sample$time)
urine_sample$time <- gsub("3","6",urine_sample$time)
urine_sample$time <- gsub("2","4",urine_sample$time)
urine_sample$time <- gsub("1","2",urine_sample$time)
urine_sample$time <- gsub("5","24",urine_sample$time)

urine_sample$time <- as.numeric(urine_sample$time)
table(urine_sample$time)
#0   2   4   6   8  24 
#114 114 107 110 105 112 

#calculate half, cmax, tmax across each metabolites by food
pk_bf <- ncafunc2(urine_sample,food="Beef") 
pk_cs <- ncafunc2(urine_sample,food="Cheese")
pk_ck <- ncafunc2(urine_sample,food="Chicken")
pk_yg <- ncafunc2(urine_sample,food="Yogurt")
pk_cn <- ncafunc2(urine_sample,food="Corn")
pk_ot <- ncafunc2(urine_sample,food="Oats")
pk_pt <- ncafunc2(urine_sample,food="Potato")
pk_bd <- ncafunc2(urine_sample,food="Whole Wheat Bread")

#save results
save(pk_bf,pk_ck,pk_ch,pk_yg,pk_cn,pk_ot,pk_pt,pk_bd,file="Results_Urine_NCA.RData")
