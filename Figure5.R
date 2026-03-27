#load packages
library(data.table)
library(dplyr)
library(tidyr)
library(cvTools) 
library(glmnet)
library(pROC)
library(PredictABEL)
library(DescTools)
library(tidyverse)

#define function for use
elasfunc <- function(DS, x1, x2, y, z){
  pen1 <- rep(0,length(x1)); pen2 <- rep(1,length(x2)); 
  pen <- c(pen1, pen2)
  DS$predictScore = 0
  elastic_matrix = as.matrix(DS[,c(x1,x2,y,z)])
  set.seed(123456)
  k = 10 
  folds = cvFolds(dim(elastic_matrix)[1], K=k)
  for (i in 1:10) {
    fit = cv.glmnet(as(elastic_matrix[folds$subsets[folds$which != i],c(x2)],"dgCMatrix"),as.vector(elastic_matrix[folds$subsets[folds$which != i],y]),family="binomial",type.measure="auc",alpha=1,penalty.factor=pen2)
    prediction = predict(fit,as(elastic_matrix[folds$subsets[folds$which == i],c(x2)],"dgCMatrix"),type="response",s="lambda.min")
    elastic_matrix[folds$subsets[folds$which == i],"predictScore"] = prediction[,1]
    lambda_min_10F = fit$lambda.min
    cmin = coef(fit, s=lambda_min_10F)
    metabsmin = data.frame(cmin[which(cmin[,1]!=0),])
    write.csv(metabsmin, file = paste(y_bi,"coeff",i,".csv",sep=''))
    write.csv(elastic_matrix, file=paste(y_bi,".csv",sep=''))
  }
}

cvnpred <- function(trainDS, testDS, x1, x2, y, k = 10, r = 100, p = 12){
  pred1 <- NULL; pred2 <- NULL; outcome <- NULL; ROCx1 <- NULL; ROCx2 <- NULL; ROCy <- NULL
  for(i in 1:r) { 
    set.seed(01)
    folds <- sample(cut(seq(1, nrow(trainDS)), breaks = k, labels=F))                                                   
    for(j in 1:k) {                                                              
      trainset <- trainDS[which(folds != j),]
      testset <- trainDS[which(folds == j),]
      #Model1                                                                    
      m1 <- paste(x1, collapse = "+"); fml1 <- paste(as.factor(y), "~", m1, sep = "") 
      fit1 <- glm(fml1, data = trainset, family = "binomial"('logit'), x = TRUE)        
      predictor1 <- predict(fit1, newdata = testset, type= "response") %>% as.data.frame
      #Model2
      m2 <- paste(x2, collapse = "+"); fml2 <- paste(as.factor(y), "~", m2, sep = "") 
      fit2 <- glm(fml2, data = trainset, family = "binomial"('logit'), x = TRUE)
      predictor2 <- predict(fit2, newdata = testset, type= "response") %>% as.data.frame()
      pred1 <- rbind(pred1, predictor1)
      pred2 <- rbind(pred2, predictor2)
      outcome <- rbind(outcome, testset[, p:(p+1)])
      ROCy <- outcome[, 1]; ROCx1 <- pred1[, 1]; ROCx2 <- pred2[, 1] 
    }
  }
  ROC_1 <- roc(ROCy, ROCx1); ROC_2 <- roc(ROCy, ROCx2)
  #NRIIDI <- reclassification(data = outcome, cOutcome = 1, predrisk1 = pred1[,1], predrisk2 = pred2[,1], cutoff = c(0, 0.30, 0.50, 0.80, 1))
  list(ROC_1, ROC_2)
}

putauc <- function(roc, format="%.3f"){
  a <- sprintf(format, ci.auc(roc)[1:3])
  b <- paste(a[2], " (", a[1],"-", a[3],")", sep="")
  return (b)
}

PLOTROC <- function(ROC_1, ROC_2, ROC_3, ROC_4){
  
  par(family = "Calibri") 
  
  plot(ROC_1,grid=F, legacy.axes=T,col="#006EBE", lty=1,family="Calibri")
  plot(ROC_2,grid=F, legacy.axes=T,col="#E69F00", lty=1,family="Calibri", add = TRUE)
  plot(ROC_3,grid=F, legacy.axes=T,col="#009E73", lty=1,family="Calibri", add = TRUE)
  plot(ROC_4,grid=F, legacy.axes=T,col="#FA5555", lty=1,family="Calibri", add = TRUE)
  
  AUC_1 <- sprintf("%.2f",ci.auc(ROC_1)[1:3])
  AUC_2 <- sprintf("%.2f",ci.auc(ROC_2)[1:3])
  AUC_3 <- sprintf("%.2f",ci.auc(ROC_3)[1:3])
  AUC_4 <- sprintf("%.2f",ci.auc(ROC_4)[1:3])
  
  Model1 <- paste("Covariates: ",AUC_1[2], " (", AUC_1[1],"-", AUC_1[3],")",sep="")
  Model2 <- paste("Plasma metabolites: ",AUC_2[2], " (", AUC_2[1],"-", AUC_2[3],")",sep="")
  Model3 <- paste("Urine metabolites: ",AUC_3[2], " (", AUC_3[1],"-", AUC_3[3],")",sep="")
  Model4 <- paste("Plasma and urine metabolites: ",AUC_4[2], " (", AUC_4[1],"-", AUC_4[3],")",sep="")
  
  legend(x=0.90,y=0.20,c(Model1,Model2,Model3,Model4),lty=1,col=c("#006EBE","#E69F00","#009E73","#FA5555"),lwd=3,cex=0.7,bty = "n")
}

run_cv_metabolites <- function(data, y_bi, metabolite_vars, r) {
  
  for(i in 1:r) { 
    set.seed(123)
    k <- 10
    folds <- sample(rep(1:k, length.out = nrow(data)))
    
    auc_results <- data.frame(metabolite = character(), AUC = numeric(), LCI = numeric(), UCI = numeric(), stringsAsFactors = FALSE)
    
    for (met in metabolite_vars) {
      
      cv_pred <- rep(NA, nrow(data))
      
      
      for (k in 1:k) {
        test_idx  <- which(folds == k)
        train_idx <- which(folds != k)
        
        m <- paste(met, collapse = "+"); fml <- paste(as.factor(y_bi), "~", m, sep = "") 
        
        # training logistic regression on a single metabolite
        fit <- glm(fml, data[train_idx, ], family = binomial)
        
        # predicted probabilities for test fold
        cv_pred[test_idx] <- predict(fit, newdata = data[test_idx, ], type = "response")
        
      }
      # average CV AUC
      roc_obj <- roc(response = data[,y_bi], predictor = cv_pred, quiet = TRUE)
      auc_value <- sprintf("%.2f",ci.auc(roc_obj)[1:3])
      
      auc_results <- rbind(auc_results, data.frame(name = met, AUC = auc_value[2], LCI = auc_value[1], UCI = auc_value[3]))
    }    
  }
  auc_results
}

#read data
load("Data_For_ElasticNet.RData")

#run elastic net using plasma data
for (i in out){
  y_bi=i
  fit <- elasfunc(FULLDS, x1, x2, y_bi, z_bi)
}

#run elastic net using urine data
for (i in out){
  y_bi=i
  fit <- elasfunc(FULLDS, x1, x3, y_bi, z_bi)
}

#run elastic net using both plasma and urine data
for (i in out){
  y_bi=i
  fit <- elasfunc(FULLDS, x1, x4, y_bi, z_bi)
}

#read beef results
bf1_1 <- fread("Ifbfcoeff1.csv") %>% as.data.frame()
bf1_2 <- fread("Ifbfcoeff2.csv") %>% as.data.frame()
bf1_3 <- fread("Ifbfcoeff3.csv") %>% as.data.frame()
bf1_4 <- fread("Ifbfcoeff4.csv") %>% as.data.frame()
bf1_5 <- fread("Ifbfcoeff5.csv") %>% as.data.frame()
bf1_6 <- fread("Ifbfcoeff6.csv") %>% as.data.frame()
bf1_7 <- fread("Ifbfcoeff7.csv") %>% as.data.frame()
bf1_8 <- fread("Ifbfcoeff8.csv") %>% as.data.frame()
bf1_9 <- fread("Ifbfcoeff9.csv") %>% as.data.frame()
bf1_10 <- fread("Ifbfcoeff10.csv") %>% as.data.frame()

bf1 <- bind_rows(bf1_1,bf1_2,bf1_3,bf1_4,bf1_5,bf1_6,bf1_7,bf1_8,bf1_9,bf1_10); names(bf1) <- c("Var","Coeff")
bf1 <- bf1 %>% count(Var) %>% arrange(n) 

bf2_1 <- fread("Ifbfcoeff1.csv") %>% as.data.frame()
bf2_2 <- fread("Ifbfcoeff2.csv") %>% as.data.frame()
bf2_3 <- fread("Ifbfcoeff3.csv") %>% as.data.frame()
bf2_4 <- fread("Ifbfcoeff4.csv") %>% as.data.frame()
bf2_5 <- fread("Ifbfcoeff5.csv") %>% as.data.frame()
bf2_6 <- fread("Ifbfcoeff6.csv") %>% as.data.frame()
bf2_7 <- fread("Ifbfcoeff7.csv") %>% as.data.frame()
bf2_8 <- fread("Ifbfcoeff8.csv") %>% as.data.frame()
bf2_9 <- fread("Ifbfcoeff9.csv") %>% as.data.frame()
bf2_10 <- fread("Ifbfcoeff10.csv") %>% as.data.frame()

bf2 <- bind_rows(bf2_1,bf2_2,bf2_3,bf2_4,bf2_5,bf2_6,bf2_7,bf2_8,bf2_9,bf2_10); names(bf2) <- c("Var","Coeff")
bf2 <- bf2 %>% count(Var) %>% arrange(n) 

bf3_1 <- fread("Ifbfcoeff1.csv") %>% as.data.frame()
bf3_2 <- fread("Ifbfcoeff2.csv") %>% as.data.frame()
bf3_3 <- fread("Ifbfcoeff3.csv") %>% as.data.frame()
bf3_4 <- fread("Ifbfcoeff4.csv") %>% as.data.frame()
bf3_5 <- fread("Ifbfcoeff5.csv") %>% as.data.frame()
bf3_6 <- fread("Ifbfcoeff6.csv") %>% as.data.frame()
bf3_7 <- fread("Ifbfcoeff7.csv") %>% as.data.frame()
bf3_8 <- fread("Ifbfcoeff8.csv") %>% as.data.frame()
bf3_9 <- fread("Ifbfcoeff9.csv") %>% as.data.frame()
bf3_10 <- fread("Ifbfcoeff10.csv") %>% as.data.frame()

bf3 <- bind_rows(bf3_1,bf3_2,bf3_3,bf3_4,bf3_5,bf3_6,bf3_7,bf3_8,bf3_9,bf3_10); names(bf3) <- c("Var","Coeff")
bf3 <- bf3 %>% count(Var) %>% arrange(n) 

#read chicken results
ck1_1 <- fread("Ifckcoeff1.csv") %>% as.data.frame()
ck1_2 <- fread("Ifckcoeff2.csv") %>% as.data.frame()
ck1_3 <- fread("Ifckcoeff3.csv") %>% as.data.frame()
ck1_4 <- fread("Ifckcoeff4.csv") %>% as.data.frame()
ck1_5 <- fread("Ifckcoeff5.csv") %>% as.data.frame()
ck1_6 <- fread("Ifckcoeff6.csv") %>% as.data.frame()
ck1_7 <- fread("Ifckcoeff7.csv") %>% as.data.frame()
ck1_8 <- fread("Ifckcoeff8.csv") %>% as.data.frame()
ck1_9 <- fread("Ifckcoeff9.csv") %>% as.data.frame()
ck1_10 <- fread("Ifckcoeff10.csv") %>% as.data.frame()

ck1 <- bind_rows(ck1_1,ck1_2,ck1_3,ck1_4,ck1_5,ck1_6,ck1_7,ck1_8,ck1_9,ck1_10); names(ck1) <- c("Var","Coeff")
ck1 <- ck1 %>% count(Var) %>% arrange(n) 

ck2_1 <- fread("Ifckcoeff1.csv") %>% as.data.frame()
ck2_2 <- fread("Ifckcoeff2.csv") %>% as.data.frame()
ck2_3 <- fread("Ifckcoeff3.csv") %>% as.data.frame()
ck2_4 <- fread("Ifckcoeff4.csv") %>% as.data.frame()
ck2_5 <- fread("Ifckcoeff5.csv") %>% as.data.frame()
ck2_6 <- fread("Ifckcoeff6.csv") %>% as.data.frame()
ck2_7 <- fread("Ifckcoeff7.csv") %>% as.data.frame()
ck2_8 <- fread("Ifckcoeff8.csv") %>% as.data.frame()
ck2_9 <- fread("Ifckcoeff9.csv") %>% as.data.frame()
ck2_10 <- fread("Ifckcoeff10.csv") %>% as.data.frame()

ck2 <- bind_rows(ck2_1,ck2_2,ck2_3,ck2_4,ck2_5,ck2_6,ck2_7,ck2_8,ck2_9,ck2_10); names(ck2) <- c("Var","Coeff")
ck2 <- ck2 %>% count(Var) %>% arrange(n) 

ck3_1 <- fread("Ifckcoeff1.csv") %>% as.data.frame()
ck3_2 <- fread("Ifckcoeff2.csv") %>% as.data.frame()
ck3_3 <- fread("Ifckcoeff3.csv") %>% as.data.frame()
ck3_4 <- fread("Ifckcoeff4.csv") %>% as.data.frame()
ck3_5 <- fread("Ifckcoeff5.csv") %>% as.data.frame()
ck3_6 <- fread("Ifckcoeff6.csv") %>% as.data.frame()
ck3_7 <- fread("Ifckcoeff7.csv") %>% as.data.frame()
ck3_8 <- fread("Ifckcoeff8.csv") %>% as.data.frame()
ck3_9 <- fread("Ifckcoeff9.csv") %>% as.data.frame()
ck3_10 <- fread("Ifckcoeff10.csv") %>% as.data.frame()

ck3 <- bind_rows(ck3_1,ck3_2,ck3_3,ck3_4,ck3_5,ck3_6,ck3_7,ck3_8,ck3_9,ck3_10); names(ck3) <- c("Var","Coeff")
ck3 <- ck3 %>% count(Var) %>% arrange(n) 

#read cheese results
ch1_1 <- fread("Ifchcoeff1.csv") %>% as.data.frame()
ch1_2 <- fread("Ifchcoeff2.csv") %>% as.data.frame()
ch1_3 <- fread("Ifchcoeff3.csv") %>% as.data.frame()
ch1_4 <- fread("Ifchcoeff4.csv") %>% as.data.frame()
ch1_5 <- fread("Ifchcoeff5.csv") %>% as.data.frame()
ch1_6 <- fread("Ifchcoeff6.csv") %>% as.data.frame()
ch1_7 <- fread("Ifchcoeff7.csv") %>% as.data.frame()
ch1_8 <- fread("Ifchcoeff8.csv") %>% as.data.frame()
ch1_9 <- fread("Ifchcoeff9.csv") %>% as.data.frame()
ch1_10 <- fread("Ifchcoeff10.csv") %>% as.data.frame()

ch1 <- bind_rows(ch1_1,ch1_2,ch1_3,ch1_4,ch1_5,ch1_6,ch1_7,ch1_8,ch1_9,ch1_10); names(ch1) <- c("Var","Coeff")
ch1 <- ch1 %>% count(Var) %>% arrange(n) 

ch2_1 <- fread("Ifchcoeff1.csv") %>% as.data.frame()
ch2_2 <- fread("Ifchcoeff2.csv") %>% as.data.frame()
ch2_3 <- fread("Ifchcoeff3.csv") %>% as.data.frame()
ch2_4 <- fread("Ifchcoeff4.csv") %>% as.data.frame()
ch2_5 <- fread("Ifchcoeff5.csv") %>% as.data.frame()
ch2_6 <- fread("Ifchcoeff6.csv") %>% as.data.frame()
ch2_7 <- fread("Ifchcoeff7.csv") %>% as.data.frame()
ch2_8 <- fread("Ifchcoeff8.csv") %>% as.data.frame()
ch2_9 <- fread("Ifchcoeff9.csv") %>% as.data.frame()
ch2_10 <- fread("Ifchcoeff10.csv") %>% as.data.frame()

ch2 <- bind_rows(ch2_1,ch2_2,ch2_3,ch2_4,ch2_5,ch2_6,ch2_7,ch2_8,ch2_9,ch2_10); names(ch2) <- c("Var","Coeff")
ch2 <- ch2 %>% count(Var) %>% arrange(n) 

ch3_1 <- fread("Ifchcoeff1.csv") %>% as.data.frame()
ch3_2 <- fread("Ifchcoeff2.csv") %>% as.data.frame()
ch3_3 <- fread("Ifchcoeff3.csv") %>% as.data.frame()
ch3_4 <- fread("Ifchcoeff4.csv") %>% as.data.frame()
ch3_5 <- fread("Ifchcoeff5.csv") %>% as.data.frame()
ch3_6 <- fread("Ifchcoeff6.csv") %>% as.data.frame()
ch3_7 <- fread("Ifchcoeff7.csv") %>% as.data.frame()
ch3_8 <- fread("Ifchcoeff8.csv") %>% as.data.frame()
ch3_9 <- fread("Ifchcoeff9.csv") %>% as.data.frame()
ch3_10 <- fread("Ifchcoeff10.csv") %>% as.data.frame()

ch3 <- bind_rows(ch3_1,ch3_2,ch3_3,ch3_4,ch3_5,ch3_6,ch3_7,ch3_8,ch3_9,ch3_10); names(ch3) <- c("Var","Coeff")
ch3 <- ch3 %>% count(Var) %>% arrange(n) 

#read yogurt results
yg1_1 <- fread("Ifygcoeff1.csv") %>% as.data.frame()
yg1_2 <- fread("Ifygcoeff2.csv") %>% as.data.frame()
yg1_3 <- fread("Ifygcoeff3.csv") %>% as.data.frame()
yg1_4 <- fread("Ifygcoeff4.csv") %>% as.data.frame()
yg1_5 <- fread("Ifygcoeff5.csv") %>% as.data.frame()
yg1_6 <- fread("Ifygcoeff6.csv") %>% as.data.frame()
yg1_7 <- fread("Ifygcoeff7.csv") %>% as.data.frame()
yg1_8 <- fread("Ifygcoeff8.csv") %>% as.data.frame()
yg1_9 <- fread("Ifygcoeff9.csv") %>% as.data.frame()
yg1_10 <- fread("Ifygcoeff10.csv") %>% as.data.frame()

yg1 <- bind_rows(yg1_1,yg1_2,yg1_3,yg1_4,yg1_5,yg1_6,yg1_7,yg1_8,yg1_9,yg1_10); names(yg1) <- c("Var","Coeff")
yg1 <- yg1 %>% count(Var) %>% arrange(n) 

yg2_1 <- fread("Ifygcoeff1.csv") %>% as.data.frame()
yg2_2 <- fread("Ifygcoeff2.csv") %>% as.data.frame()
yg2_3 <- fread("Ifygcoeff3.csv") %>% as.data.frame()
yg2_4 <- fread("Ifygcoeff4.csv") %>% as.data.frame()
yg2_5 <- fread("Ifygcoeff5.csv") %>% as.data.frame()
yg2_6 <- fread("Ifygcoeff6.csv") %>% as.data.frame()
yg2_7 <- fread("Ifygcoeff7.csv") %>% as.data.frame()
yg2_8 <- fread("Ifygcoeff8.csv") %>% as.data.frame()
yg2_9 <- fread("Ifygcoeff9.csv") %>% as.data.frame()
yg2_10 <- fread("Ifygcoeff10.csv") %>% as.data.frame()

yg2 <- bind_rows(yg2_1,yg2_2,yg2_3,yg2_4,yg2_5,yg2_6,yg2_7,yg2_8,yg2_9,yg2_10); names(yg2) <- c("Var","Coeff")
yg2 <- yg2 %>% count(Var) %>% arrange(n) 

yg3_1 <- fread("Ifygcoeff1.csv") %>% as.data.frame()
yg3_2 <- fread("Ifygcoeff2.csv") %>% as.data.frame()
yg3_3 <- fread("Ifygcoeff3.csv") %>% as.data.frame()
yg3_4 <- fread("Ifygcoeff4.csv") %>% as.data.frame()
yg3_5 <- fread("Ifygcoeff5.csv") %>% as.data.frame()
yg3_6 <- fread("Ifygcoeff6.csv") %>% as.data.frame()
yg3_7 <- fread("Ifygcoeff7.csv") %>% as.data.frame()
yg3_8 <- fread("Ifygcoeff8.csv") %>% as.data.frame()
yg3_9 <- fread("Ifygcoeff9.csv") %>% as.data.frame()
yg3_10 <- fread("Ifygcoeff10.csv") %>% as.data.frame()

yg3 <- bind_rows(yg3_1,yg3_2,yg3_3,yg3_4,yg3_5,yg3_6,yg3_7,yg3_8,yg3_9,yg3_10); names(yg3) <- c("Var","Coeff")
yg3 <- yg3 %>% count(Var) %>% arrange(n) 

#read corn results
cr1_1 <- fread("Ifcrcoeff1.csv") %>% as.data.frame()
cr1_2 <- fread("Ifcrcoeff2.csv") %>% as.data.frame()
cr1_3 <- fread("Ifcrcoeff3.csv") %>% as.data.frame()
cr1_4 <- fread("Ifcrcoeff4.csv") %>% as.data.frame()
cr1_5 <- fread("Ifcrcoeff5.csv") %>% as.data.frame()
cr1_6 <- fread("Ifcrcoeff6.csv") %>% as.data.frame()
cr1_7 <- fread("Ifcrcoeff7.csv") %>% as.data.frame()
cr1_8 <- fread("Ifcrcoeff8.csv") %>% as.data.frame()
cr1_9 <- fread("Ifcrcoeff9.csv") %>% as.data.frame()
cr1_10 <- fread("Ifcrcoeff10.csv") %>% as.data.frame()

cr1 <- bind_rows(cr1_1,cr1_2,cr1_3,cr1_4,cr1_5,cr1_6,cr1_7,cr1_8,cr1_9,cr1_10); names(cr1) <- c("Var","Coeff")
cr1 <- cr1 %>% count(Var) %>% arrange(n) 

cr2_1 <- fread("Ifcrcoeff1.csv") %>% as.data.frame()
cr2_2 <- fread("Ifcrcoeff2.csv") %>% as.data.frame()
cr2_3 <- fread("Ifcrcoeff3.csv") %>% as.data.frame()
cr2_4 <- fread("Ifcrcoeff4.csv") %>% as.data.frame()
cr2_5 <- fread("Ifcrcoeff5.csv") %>% as.data.frame()
cr2_6 <- fread("Ifcrcoeff6.csv") %>% as.data.frame()
cr2_7 <- fread("Ifcrcoeff7.csv") %>% as.data.frame()
cr2_8 <- fread("Ifcrcoeff8.csv") %>% as.data.frame()
cr2_9 <- fread("Ifcrcoeff9.csv") %>% as.data.frame()
cr2_10 <- fread("Ifcrcoeff10.csv") %>% as.data.frame()

cr2 <- bind_rows(cr2_1,cr2_2,cr2_3,cr2_4,cr2_5,cr2_6,cr2_7,cr2_8,cr2_9,cr2_10); names(cr2) <- c("Var","Coeff")
cr2 <- cr2 %>% count(Var) %>% arrange(n) 

cr3_1 <- fread("Ifcrcoeff1.csv") %>% as.data.frame()
cr3_2 <- fread("Ifcrcoeff2.csv") %>% as.data.frame()
cr3_3 <- fread("Ifcrcoeff3.csv") %>% as.data.frame()
cr3_4 <- fread("Ifcrcoeff4.csv") %>% as.data.frame()
cr3_5 <- fread("Ifcrcoeff5.csv") %>% as.data.frame()
cr3_6 <- fread("Ifcrcoeff6.csv") %>% as.data.frame()
cr3_7 <- fread("Ifcrcoeff7.csv") %>% as.data.frame()
cr3_8 <- fread("Ifcrcoeff8.csv") %>% as.data.frame()
cr3_9 <- fread("Ifcrcoeff9.csv") %>% as.data.frame()
cr3_10 <- fread("Ifcrcoeff10.csv") %>% as.data.frame()

cr3 <- bind_rows(cr3_1,cr3_2,cr3_3,cr3_4,cr3_5,cr3_6,cr3_7,cr3_8,cr3_9,cr3_10); names(cr3) <- c("Var","Coeff")
cr3 <- cr3 %>% count(Var) %>% arrange(n) 

#read oats results
ot1_1 <- fread("Ifotcoeff1.csv") %>% as.data.frame()
ot1_2 <- fread("Ifotcoeff2.csv") %>% as.data.frame()
ot1_3 <- fread("Ifotcoeff3.csv") %>% as.data.frame()
ot1_4 <- fread("Ifotcoeff4.csv") %>% as.data.frame()
ot1_5 <- fread("Ifotcoeff5.csv") %>% as.data.frame()
ot1_6 <- fread("Ifotcoeff6.csv") %>% as.data.frame()
ot1_7 <- fread("Ifotcoeff7.csv") %>% as.data.frame()
ot1_8 <- fread("Ifotcoeff8.csv") %>% as.data.frame()
ot1_9 <- fread("Ifotcoeff9.csv") %>% as.data.frame()
ot1_10 <- fread("Ifotcoeff10.csv") %>% as.data.frame()

ot1 <- bind_rows(ot1_1,ot1_2,ot1_3,ot1_4,ot1_5,ot1_6,ot1_7,ot1_8,ot1_9,ot1_10); names(ot1) <- c("Var","Coeff")
ot1 <- ot1 %>% count(Var) %>% arrange(n) 

ot2_1 <- fread("Ifotcoeff1.csv") %>% as.data.frame()
ot2_2 <- fread("Ifotcoeff2.csv") %>% as.data.frame()
ot2_3 <- fread("Ifotcoeff3.csv") %>% as.data.frame()
ot2_4 <- fread("Ifotcoeff4.csv") %>% as.data.frame()
ot2_5 <- fread("Ifotcoeff5.csv") %>% as.data.frame()
ot2_6 <- fread("Ifotcoeff6.csv") %>% as.data.frame()
ot2_7 <- fread("Ifotcoeff7.csv") %>% as.data.frame()
ot2_8 <- fread("Ifotcoeff8.csv") %>% as.data.frame()
ot2_9 <- fread("Ifotcoeff9.csv") %>% as.data.frame()
ot2_10 <- fread("Ifotcoeff10.csv") %>% as.data.frame()

ot2 <- bind_rows(ot2_1,ot2_2,ot2_3,ot2_4,ot2_5,ot2_6,ot2_7,ot2_8,ot2_9,ot2_10); names(ot2) <- c("Var","Coeff")
ot2 <- ot2 %>% count(Var) %>% arrange(n) 

ot3_1 <- fread("Ifotcoeff1.csv") %>% as.data.frame()
ot3_2 <- fread("Ifotcoeff2.csv") %>% as.data.frame()
ot3_3 <- fread("Ifotcoeff3.csv") %>% as.data.frame()
ot3_4 <- fread("Ifotcoeff4.csv") %>% as.data.frame()
ot3_5 <- fread("Ifotcoeff5.csv") %>% as.data.frame()
ot3_6 <- fread("Ifotcoeff6.csv") %>% as.data.frame()
ot3_7 <- fread("Ifotcoeff7.csv") %>% as.data.frame()
ot3_8 <- fread("Ifotcoeff8.csv") %>% as.data.frame()
ot3_9 <- fread("Ifotcoeff9.csv") %>% as.data.frame()
ot3_10 <- fread("Ifotcoeff10.csv") %>% as.data.frame()

ot3 <- bind_rows(ot3_1,ot3_2,ot3_3,ot3_4,ot3_5,ot3_6,ot3_7,ot3_8,ot3_9,ot3_10); names(ot3) <- c("Var","Coeff")
ot3 <- ot3 %>% count(Var) %>% arrange(n) 

#read bread results
ww1_1 <- fread("Ifwwcoeff1.csv") %>% as.data.frame()
ww1_2 <- fread("Ifwwcoeff2.csv") %>% as.data.frame()
ww1_3 <- fread("Ifwwcoeff3.csv") %>% as.data.frame()
ww1_4 <- fread("Ifwwcoeff4.csv") %>% as.data.frame()
ww1_5 <- fread("Ifwwcoeff5.csv") %>% as.data.frame()
ww1_6 <- fread("Ifwwcoeff6.csv") %>% as.data.frame()
ww1_7 <- fread("Ifwwcoeff7.csv") %>% as.data.frame()
ww1_8 <- fread("Ifwwcoeff8.csv") %>% as.data.frame()
ww1_9 <- fread("Ifwwcoeff9.csv") %>% as.data.frame()
ww1_10 <- fread("Ifwwcoeff10.csv") %>% as.data.frame()

ww1 <- bind_rows(ww1_1,ww1_2,ww1_3,ww1_4,ww1_5,ww1_6,ww1_7,ww1_8,ww1_9,ww1_10); names(ww1) <- c("Var","Coeff")
ww1 <- ww1 %>% count(Var) %>% arrange(n) 

ww2_1 <- fread("Ifwwcoeff1.csv") %>% as.data.frame()
ww2_2 <- fread("Ifwwcoeff2.csv") %>% as.data.frame()
ww2_3 <- fread("Ifwwcoeff3.csv") %>% as.data.frame()
ww2_4 <- fread("Ifwwcoeff4.csv") %>% as.data.frame()
ww2_5 <- fread("Ifwwcoeff5.csv") %>% as.data.frame()
ww2_6 <- fread("Ifwwcoeff6.csv") %>% as.data.frame()
ww2_7 <- fread("Ifwwcoeff7.csv") %>% as.data.frame()
ww2_8 <- fread("Ifwwcoeff8.csv") %>% as.data.frame()
ww2_9 <- fread("Ifwwcoeff9.csv") %>% as.data.frame()
ww2_10 <- fread("Ifwwcoeff10.csv") %>% as.data.frame()

ww2 <- bind_rows(ww2_1,ww2_2,ww2_3,ww2_4,ww2_5,ww2_6,ww2_7,ww2_8,ww2_9,ww2_10); names(ww2) <- c("Var","Coeff")
ww2 <- ww2 %>% count(Var) %>% arrange(n) 

ww3_1 <- fread("Ifwwcoeff1.csv") %>% as.data.frame()
ww3_2 <- fread("Ifwwcoeff2.csv") %>% as.data.frame()
ww3_3 <- fread("Ifwwcoeff3.csv") %>% as.data.frame()
ww3_4 <- fread("Ifwwcoeff4.csv") %>% as.data.frame()
ww3_5 <- fread("Ifwwcoeff5.csv") %>% as.data.frame()
ww3_6 <- fread("Ifwwcoeff6.csv") %>% as.data.frame()
ww3_7 <- fread("Ifwwcoeff7.csv") %>% as.data.frame()
ww3_8 <- fread("Ifwwcoeff8.csv") %>% as.data.frame()
ww3_9 <- fread("Ifwwcoeff9.csv") %>% as.data.frame()
ww3_10 <- fread("Ifwwcoeff10.csv") %>% as.data.frame()

ww3 <- bind_rows(ww3_1,ww3_2,ww3_3,ww3_4,ww3_5,ww3_6,ww3_7,ww3_8,ww3_9,ww3_10); names(ww3) <- c("Var","Coeff")
ww3 <- ww3 %>% count(Var) %>% arrange(n) 

#read potato results
pt1_1 <- fread("Ifptcoeff1.csv") %>% as.data.frame()
pt1_2 <- fread("Ifptcoeff2.csv") %>% as.data.frame()
pt1_3 <- fread("Ifptcoeff3.csv") %>% as.data.frame()
pt1_4 <- fread("Ifptcoeff4.csv") %>% as.data.frame()
pt1_5 <- fread("Ifptcoeff5.csv") %>% as.data.frame()
pt1_6 <- fread("Ifptcoeff6.csv") %>% as.data.frame()
pt1_7 <- fread("Ifptcoeff7.csv") %>% as.data.frame()
pt1_8 <- fread("Ifptcoeff8.csv") %>% as.data.frame()
pt1_9 <- fread("Ifptcoeff9.csv") %>% as.data.frame()
pt1_10 <- fread("Ifptcoeff10.csv") %>% as.data.frame()

pt1 <- bind_rows(pt1_1,pt1_2,pt1_3,pt1_4,pt1_5,pt1_6,pt1_7,pt1_8,pt1_9,pt1_10); names(pt1) <- c("Var","Coeff")
pt1 <- pt1 %>% count(Var) %>% arrange(n) 

pt2_1 <- fread("Ifptcoeff1.csv") %>% as.data.frame()
pt2_2 <- fread("Ifptcoeff2.csv") %>% as.data.frame()
pt2_3 <- fread("Ifptcoeff3.csv") %>% as.data.frame()
pt2_4 <- fread("Ifptcoeff4.csv") %>% as.data.frame()
pt2_5 <- fread("Ifptcoeff5.csv") %>% as.data.frame()
pt2_6 <- fread("Ifptcoeff6.csv") %>% as.data.frame()
pt2_7 <- fread("Ifptcoeff7.csv") %>% as.data.frame()
pt2_8 <- fread("Ifptcoeff8.csv") %>% as.data.frame()
pt2_9 <- fread("Ifptcoeff9.csv") %>% as.data.frame()
pt2_10 <- fread("Ifptcoeff10.csv") %>% as.data.frame()

pt2 <- bind_rows(pt2_1,pt2_2,pt2_3,pt2_4,pt2_6,pt2_7,pt2_8,pt2_9,pt2_10); names(pt2) <- c("Var","Coeff")
pt2 <- pt2 %>% count(Var) %>% arrange(n) 

pt3_1 <- fread("Ifptcoeff1.csv") %>% as.data.frame()
pt3_2 <- fread("Ifptcoeff2.csv") %>% as.data.frame()
pt3_3 <- fread("Ifptcoeff3.csv") %>% as.data.frame()
pt3_4 <- fread("Ifptcoeff4.csv") %>% as.data.frame()
pt3_5 <- fread("Ifptcoeff5.csv") %>% as.data.frame()
pt3_6 <- fread("Ifptcoeff6.csv") %>% as.data.frame()
pt3_7 <- fread("Ifptcoeff7.csv") %>% as.data.frame()
pt3_8 <- fread("Ifptcoeff8.csv") %>% as.data.frame()
pt3_9 <- fread("Ifptcoeff9.csv") %>% as.data.frame()
pt3_10 <- fread("Ifptcoeff10.csv") %>% as.data.frame()

pt3 <- bind_rows(pt3_1,pt3_2,pt3_3,pt3_4,pt3_5,pt3_6,pt3_7,pt3_8,pt3_9,pt3_10); names(pt3) <- c("Var","Coeff")
pt3 <- pt3 %>% count(Var) %>% arrange(n) 

#organize elastic net results from plasma data
bf1_1$group <- "Beef"
ck1_1$group <- "Chicken"
ch1_1$group <- "Cheese"
yg1_1$group <- "Yogurt"
cr1_1$group <- "Corn"
ot1_1$group <- "Oats"
ww1_1$group <- "Bread"
pt1_1$group <- "Potato"

use <- rbind(bf1_1,ck1_1,ch1_1,yg1_1,cr1_1,ot1_1,ww1_1,pt1_1)

names(use) <- c("name","coeff","group")

use$uID <- sub("[_].*$", "", use$name)
use <- left_join(use,pla.anno[,c("uID","Metabolite")],by="uID") %>% left_join(uri.anno[,c("uID","Metabolite")],by="uID")

use[grepl("_u", use$name),"Metabolite.x"] <- use[grepl("_u", use$name),"Metabolite.y"]
use[which(is.na(use$Metabolite.x)),"Metabolite.x"] <- use[which(is.na(use$Metabolite.x)),"name"]

names(use)[5] <- "Metabolite"
use <- use[,-6]

use1 <- use[-which(use$name %in% c("(Intercept)","ips_age","ips_bmi","ps_sex","race")),]

use1 <- use1[which(use1$coeff>0),]

use1$type <- ifelse(grepl("_u", use1$name), "Urine", "Plasma")

use1_1 <- use1[grepl("t1", use1$name), ]
use1_2 <- use1[grepl("t2", use1$name), ]
use1_3 <- use1[grepl("t3", use1$name), ]
use1_4 <- use1[grepl("t4", use1$name), ]
use1_5 <- use1[grepl("t5", use1$name), ]
use1_6 <- use1[grepl("t6", use1$name), ]
use1_7 <- use1[grepl("t7", use1$name), ]
use1_8 <- use1[grepl("t8", use1$name), ]
use1_24 <- use1[grepl("t24", use1$name), ]

use1_2 <- use1_2 %>% filter(!grepl("t24", name))

use1_1$time <- "1-hr"
use1_2$time <- "2-hr"
use1_3$time <- "3-hr"
use1_4$time <- "4-hr"
use1_5$time <- "5-hr"
use1_6$time <- "6-hr"
use1_7$time <- "7-hr"
use1_8$time <- "8-hr"
use1_24$time <- "24-hr"

use1 <- rbind(use1_1,use1_2,use1_3,use1_4,use1_5,use1_6,use1_7,use1_8,use1_24)

use1$group <- factor(use1$group,levels=c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Bread","Potato"))

#organize elastic net results from urine data
bf2_1$group <- "Beef"
ck2_1$group <- "Chicken"
ch2_1$group <- "Cheese"
yg2_1$group <- "Yogurt"
cr2_1$group <- "Corn"
ot2_1$group <- "Oats"
ww2_1$group <- "Bread"
pt2_1$group <- "Potato"

use <- rbind(bf2_1,ck2_1,ch2_1,yg2_1,cr2_1,ot2_1,ww2_1,pt2_1)

names(use) <- c("name","coeff","group")

use$uID <- sub("[_].*$", "", use$name)
use <- left_join(use,pla.anno[,c("uID","Metabolite")],by="uID") %>% left_join(uri.anno[,c("uID","Metabolite")],by="uID")

use[grepl("_u", use$name),"Metabolite.x"] <- use[grepl("_u", use$name),"Metabolite.y"]
use[which(is.na(use$Metabolite.x)),"Metabolite.x"] <- use[which(is.na(use$Metabolite.x)),"name"]

names(use)[5] <- "Metabolite"
use <- use[,-6]

use2 <- use[-which(use$name %in% c("(Intercept)","ips_age","ips_bmi","ps_sex","race")),]

use2 <- use2[which(use2$coeff>0),]

use2$type <- ifelse(grepl("_u", use2$name), "Urine", "Plasma")

use2_1 <- use2[grepl("t1", use2$name), ]
use2_2 <- use2[grepl("t2", use2$name), ]
use2_3 <- use2[grepl("t3", use2$name), ]
use2_4 <- use2[grepl("t4", use2$name), ]
use2_5 <- use2[grepl("t5", use2$name), ]
use2_6 <- use2[grepl("t6", use2$name), ]
use2_7 <- use2[grepl("t7", use2$name), ]
use2_8 <- use2[grepl("t8", use2$name), ]
use2_24 <- use2[grepl("t24", use2$name), ]

use2_2 <- use2_2 %>% filter(!grepl("t24", name))

use2_1$time <- "1-hr"
use2_2$time <- "2-hr"
use2_3$time <- "3-hr"
use2_4$time <- "4-hr"
use2_5$time <- "5-hr"
use2_6$time <- "6-hr"
use2_7$time <- "7-hr"
use2_8$time <- "8-hr"
use2_24$time <- "24-hr"

use2 <- rbind(use2_1,use2_2,use2_3,use2_4,use2_5,use2_6,use2_7,use2_8,use2_24)

#organize elastic net results from plasma and urine data
bf3_1$group <- "Beef"
ck3_1$group <- "Chicken"
ch3_1$group <- "Cheese"
yg3_1$group <- "Yogurt"
cr3_1$group <- "Corn"
ot3_1$group <- "Oats"
ww3_1$group <- "Bread"
pt3_1$group <- "Potato"

use <- rbind(bf3_1,ck3_1,ch3_1,yg3_1,cr3_1,ot3_1,ww3_1,pt3_1)

names(use) <- c("name","coeff","group")

use$uID <- sub("[_].*$", "", use$name)
use <- left_join(use,pla.anno[,c("uID","Metabolite")],by="uID") %>% left_join(uri.anno[,c("uID","Metabolite")],by="uID")

use[grepl("_u", use$name),"Metabolite.x"] <- use[grepl("_u", use$name),"Metabolite.y"]
use[which(is.na(use$Metabolite.x)),"Metabolite.x"] <- use[which(is.na(use$Metabolite.x)),"name"]

names(use)[5] <- "Metabolite"
use <- use[,-6]

use3 <- use[-which(use$name %in% c("(Intercept)","ips_age","ips_bmi","ps_sex","race")),]

use3 <- use3[which(use3$coeff>0),]

use3$type <- ifelse(grepl("_u", use3$name), "Urine", "Plasma")

use3_1 <- use3[grepl("t1", use3$name), ]
use3_2 <- use3[grepl("t2", use3$name), ]
use3_3 <- use3[grepl("t3", use3$name), ]
use3_4 <- use3[grepl("t4", use3$name), ]
use3_5 <- use3[grepl("t5", use3$name), ]
use3_6 <- use3[grepl("t6", use3$name), ]
use3_7 <- use3[grepl("t7", use3$name), ]
use3_8 <- use3[grepl("t8", use3$name), ]
use3_24 <- use3[grepl("t24", use3$name), ]

use3_2 <- use3_2 %>% filter(!grepl("t24", name))

use3_1$time <- "1-hr"
use3_2$time <- "2-hr"
use3_3$time <- "3-hr"
use3_4$time <- "4-hr"
use3_5$time <- "5-hr"
use3_6$time <- "6-hr"
use3_7$time <- "7-hr"
use3_8$time <- "8-hr"
use3_24$time <- "24-hr"

use3 <- rbind(use3_1,use3_2,use3_3,use3_4,use3_5,use3_6,use3_7,use3_8,use3_24)

#read predicted score data
setwd("/n/holylfs05/LABS/liang_lab/Lab/huanyun/DBDC/PK/results/temp/ElasticNet/Plasma")
bf1 <- fread("Ifbf.csv") %>% as.data.frame()
ck1 <- fread("Ifck.csv") %>% as.data.frame()
ch1 <- fread("Ifch.csv") %>% as.data.frame()
yg1 <- fread("Ifyg.csv") %>% as.data.frame()
cr1 <- fread("Ifcr.csv") %>% as.data.frame()
ot1 <- fread("Ifot.csv") %>% as.data.frame()
ww1 <- fread("Ifww.csv") %>% as.data.frame()
pt1 <- fread("Ifpt.csv") %>% as.data.frame()

FULLDS$bf_mrs_pla <- bf1$predictScore
FULLDS$ck_mrs_pla <- ck1$predictScore
FULLDS$ch_mrs_pla <- ch1$predictScore
FULLDS$yg_mrs_pla <- yg1$predictScore
FULLDS$cr_mrs_pla <- cr1$predictScore
FULLDS$ot_mrs_pla <- ot1$predictScore
FULLDS$ww_mrs_pla <- ww1$predictScore
FULLDS$pt_mrs_pla <- pt1$predictScore

#read predicted score data
setwd("/n/holylfs05/LABS/liang_lab/Lab/huanyun/DBDC/PK/results/temp/ElasticNet/Urine")
bf2 <- fread("Ifbf.csv") %>% as.data.frame()
ck2 <- fread("Ifck.csv") %>% as.data.frame()
ch2 <- fread("Ifch.csv") %>% as.data.frame()
yg2 <- fread("Ifyg.csv") %>% as.data.frame()
cr2 <- fread("Ifcr.csv") %>% as.data.frame()
ot2 <- fread("Ifot.csv") %>% as.data.frame()
ww2 <- fread("Ifww.csv") %>% as.data.frame()
pt2 <- fread("Ifpt.csv") %>% as.data.frame()

FULLDS$bf_mrs_uri <- bf2$predictScore
FULLDS$ck_mrs_uri <- ck2$predictScore
FULLDS$ch_mrs_uri <- ch2$predictScore
FULLDS$yg_mrs_uri <- yg2$predictScore
FULLDS$cr_mrs_uri <- cr2$predictScore
FULLDS$ot_mrs_uri <- ot2$predictScore
FULLDS$ww_mrs_uri <- ww2$predictScore
FULLDS$pt_mrs_uri <- pt2$predictScore

#read predicted score data
setwd("/n/holylfs05/LABS/liang_lab/Lab/huanyun/DBDC/PK/results/temp/ElasticNet/Plasma_Urine")
bf3 <- fread("Ifbf.csv") %>% as.data.frame()
ck3 <- fread("Ifck.csv") %>% as.data.frame()
ch3 <- fread("Ifch.csv") %>% as.data.frame()
yg3 <- fread("Ifyg.csv") %>% as.data.frame()
cr3 <- fread("Ifcr.csv") %>% as.data.frame()
ot3 <- fread("Ifot.csv") %>% as.data.frame()
ww3 <- fread("Ifww.csv") %>% as.data.frame()
pt3 <- fread("Ifpt.csv") %>% as.data.frame()

FULLDS$bf_mrs_com <- bf3$predictScore
FULLDS$ck_mrs_com <- ck3$predictScore
FULLDS$ch_mrs_com <- ch3$predictScore
FULLDS$yg_mrs_com <- yg3$predictScore
FULLDS$cr_mrs_com <- cr3$predictScore
FULLDS$ot_mrs_com <- ot3$predictScore
FULLDS$ww_mrs_com <- ww3$predictScore
FULLDS$pt_mrs_com <- pt3$predictScore

load("BFI_For_Use_Urine_Unfilerted.RData")

data <- df

data$Ifbf <- ifelse(data$group=="Beef","1","0")
data$Ifck <- ifelse(data$group=="Chicken","1","0")
data$Ifyg <- ifelse(data$group=="Yogurt","1","0")
data$Ifch <- ifelse(data$group=="Cheese","1","0")
data$Ifcr <- ifelse(data$group=="Corn","1","0")
data$Ifot <- ifelse(data$group=="Oats","1","0")
data$Ifww <- ifelse(data$group=="Whole Wheat Bread","1","0")
data$Ifpt <- ifelse(data$group=="Potato","1","0")

var <- c("Ifbf","Ifck","Ifyg","Ifch","Ifcr","Ifot","Ifww","Ifpt")

for (i in var){
  data[,i] <- as.factor(data[,i])
}

fm1 <- cvnpred(data, data, x1, x1, "Ifbf", k = 10, r = 5, p = match("Ifbf",names(data)))
fm2 <- cvnpred(data, data, x1, x1, "Ifck", k = 10, r = 5, p = match("Ifck",names(data)))
fm3 <- cvnpred(data, data, x1, x1, "Ifch", k = 10, r = 5, p = match("Ifch",names(data)))
fm4 <- cvnpred(data, data, x1, x1, "Ifyg", k = 10, r = 5, p = match("Ifyg",names(data)))
fm5 <- cvnpred(data, data, x1, x1, "Ifcr", k = 10, r = 5, p = match("Ifcr",names(data)))
fm6 <- cvnpred(data, data, x1, x1, "Ifot", k = 10, r = 5, p = match("Ifot",names(data)))
fm7 <- cvnpred(data, data, x1, x1, "Ifww", k = 10, r = 5, p = match("Ifww",names(data)))
fm8 <- cvnpred(data, data, x1, x1, "Ifpt", k = 10, r = 5, p = match("Ifpt",names(data)))

#predictive ability by covariates
a0 <- fm1[[1]]
b0 <- fm2[[1]]
c0 <- fm3[[1]]
d0 <- fm4[[1]]
e0 <- fm5[[1]]
f0 <- fm6[[1]]
g0 <- fm7[[1]]
h0 <- fm8[[1]]

#combine covariates and mrs results together
a1 <- roc(FULLDS$Ifbf,FULLDS$bf_mrs_pla)
b1 <- roc(FULLDS$Ifck,FULLDS$ck_mrs_pla)
c1 <- roc(FULLDS$Ifch,FULLDS$ch_mrs_pla)
d1 <- roc(FULLDS$Ifyg,FULLDS$yg_mrs_pla)
e1 <- roc(FULLDS$Ifcr,FULLDS$cr_mrs_pla)
f1 <- roc(FULLDS$Ifot,FULLDS$ot_mrs_pla)
g1 <- roc(FULLDS$Ifww,FULLDS$ww_mrs_pla)

res_ot <- use1[which(use1$group=="Oats"),]

f1 <- cvnpred(FULLDS, FULLDS, x1, res_ot$name, "Ifot", k = 2, r = 5, p = match("Ifot",names(FULLDS)))[[2]]

res_ww <- use1[which(use1$group=="Bread"),]

g1 <- cvnpred(FULLDS, FULLDS, x1, res_ww$name, "Ifww", k = 10, r = 5, p = match("Ifww",names(FULLDS)))[[2]]

h1 <- roc(FULLDS$Ifpt,FULLDS$pt_mrs_pla)

a2 <- roc(FULLDS$Ifbf,FULLDS$bf_mrs_uri)
b2 <- roc(FULLDS$Ifck,FULLDS$ck_mrs_uri)
c2 <- roc(FULLDS$Ifch,FULLDS$ch_mrs_uri)
d2 <- roc(FULLDS$Ifyg,FULLDS$yg_mrs_uri)
e2 <- roc(FULLDS$Ifcr,FULLDS$cr_mrs_uri)
f2 <- roc(FULLDS$Ifot,FULLDS$ot_mrs_uri)
g2 <- roc(FULLDS$Ifww,FULLDS$ww_mrs_uri)

res_yg <- use2[which(use2$group=="Yogurt"),]

d2 <- cvnpred(FULLDS, FULLDS, x1, res_yg$name, "Ifyg", k = 2, r = 5, p = match("Ifyg",names(FULLDS)))[[2]]

res_cr <- use2[which(use2$group=="Corn"),]

e2 <- cvnpred(FULLDS, FULLDS, x1, res_cr$name, "Ifcr", k = 10, r = 5, p = match("Ifcr",names(FULLDS)))[[2]]

res_ot <- use2[which(use2$group=="Oats"),]

f2 <- cvnpred(FULLDS, FULLDS, x1, res_ot$name, "Ifot", k = 10, r = 5, p = match("Ifot",names(FULLDS)))[[2]]

res_ww <- use2[which(use2$group=="Bread"),]

g2 <- cvnpred(FULLDS, FULLDS, x1, res_ww$name, "Ifww", k = 10, r = 5, p = match("Ifww",names(FULLDS)))[[2]]

res_pt <- use2[which(use2$group=="Potato"),]

h2 <- cvnpred(FULLDS, FULLDS, x1, res_pt$name, "Ifpt", k = 10, r = 5, p = match("Ifpt",names(FULLDS)))[[2]]

a3 <- roc(FULLDS$Ifbf,FULLDS$bf_mrs_com)
b3 <- roc(FULLDS$Ifck,FULLDS$ck_mrs_com)
c3 <- roc(FULLDS$Ifch,FULLDS$ch_mrs_com)
d3 <- roc(FULLDS$Ifyg,FULLDS$yg_mrs_com)
e3 <- roc(FULLDS$Ifcr,FULLDS$cr_mrs_com)
f3 <- roc(FULLDS$Ifot,FULLDS$ot_mrs_com)
g3 <- roc(FULLDS$Ifww,FULLDS$ww_mrs_com)

res_ot <- use3[which(use3$group=="Oats"),]
f3 <- cvnpred(FULLDS, FULLDS, x1, res_ot$name, "Ifot", k = 10, r = 5, p = match("Ifot",names(FULLDS)))[[2]]

res_ww <- use3[which(use3$group=="Bread"),]

g3 <- cvnpred(FULLDS, FULLDS, x1, res_ww$name, "Ifww", k = 10, r = 5, p = match("Ifww",names(FULLDS)))[[2]]

h3 <- roc(FULLDS$Ifpt,FULLDS$pt_mrs_com)

#--------------------------------------------------------------------------------------------
#
#            Plot AUC curve
#
#--------------------------------------------------------------------------------------------
#plot plasmaROC curve
png("ROC_Plasma.png",width=450,height=400)
par(family = "Calibri") 

plot(a0,grid=F, legacy.axes=T,col="#FA5555", lty=2,family="Calibri")
plot(b0,grid=F, legacy.axes=T,col="#4F81BD", lty=2,family="Calibri", add = TRUE)
plot(c0,grid=F, legacy.axes=T,col="#cc4bc4", lty=2,family="Calibri", add = TRUE)
plot(d0,grid=F, legacy.axes=T,col="#9f81f7", lty=2,family="Calibri", add = TRUE)
plot(e0,grid=F, legacy.axes=T,col="#orange", lty=2,family="Calibri", add = TRUE)
plot(f0,grid=F, legacy.axes=T,col="#B15928", lty=2,family="Calibri", add = TRUE)
plot(g0,grid=F, legacy.axes=T,col="grey", lty=2,family="Calibri", add = TRUE)
plot(h0,grid=F, legacy.axes=T,col="pink", lty=2,family="Calibri", add = TRUE)

plot(a1,grid=F, legacy.axes=T,col="#FA5555", lty=1,family="Calibri", add = TRUE)
plot(b1,grid=F, legacy.axes=T,col="#4F81BD", lty=1,family="Calibri", add = TRUE)
plot(c1,grid=F, legacy.axes=T,col="#cc4bc4", lty=1,family="Calibri", add = TRUE)
plot(d1,grid=F, legacy.axes=T,col="#9f81f7", lty=1,family="Calibri", add = TRUE)
plot(e1,grid=F, legacy.axes=T,col="#orange", lty=1,family="Calibri", add = TRUE)
plot(f1,grid=F, legacy.axes=T,col="#B15928", lty=1,family="Calibri", add = TRUE)
plot(g1,grid=F, legacy.axes=T,col="grey", lty=1,family="Calibri", add = TRUE)
plot(h1,grid=F, legacy.axes=T,col="pink", lty=1,family="Calibri", add = TRUE)

AUC_1 <- sprintf("%.2f",ci.auc(a0)[1:3])
AUC_2 <- sprintf("%.2f",ci.auc(b0)[1:3])
AUC_3 <- sprintf("%.2f",ci.auc(c0)[1:3])
AUC_4 <- sprintf("%.2f",ci.auc(d0)[1:3])
AUC_5 <- sprintf("%.2f",ci.auc(e0)[1:3])
AUC_6 <- sprintf("%.2f",ci.auc(f0)[1:3])
AUC_7 <- sprintf("%.2f",ci.auc(g0)[1:3])
AUC_8 <- sprintf("%.2f",ci.auc(h0)[1:3])

Model1 <- paste("Beef: ",AUC_1[2], " (", AUC_1[1],"-", AUC_1[3],")",sep="")
Model2 <- paste("Chicken: ",AUC_2[2], " (", AUC_2[1],"-", AUC_2[3],")",sep="")
Model3 <- paste("Cheese: ",AUC_3[2], " (", AUC_3[1],"-", AUC_3[3],")",sep="")
Model4 <- paste("Yogurt: ",AUC_4[2], " (", AUC_4[1],"-", AUC_4[3],")",sep="")
Model5 <- paste("Corn: ",AUC_5[2], " (", AUC_5[1],"-", AUC_5[3],")",sep="")
Model6 <- paste("Oats: ",AUC_6[2], " (", AUC_6[1],"-", AUC_6[3],")",sep="")
Model7 <- paste("Bread: ",AUC_7[2], " (", AUC_7[1],"-", AUC_7[3],")",sep="")
Model8 <- paste("Potato: ",AUC_8[2], " (", AUC_8[1],"-", AUC_8[3],")",sep="")

legend(x=0.90,y=0.40,c(Model1,Model2,Model3,Model4,Model5,Model6,Model7,Model8),lty=2,col=c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"),lwd=3,cex=0.75,bty = "n")

AUC_1 <- sprintf("%.2f",ci.auc(a1)[1:3])
AUC_2 <- sprintf("%.2f",ci.auc(b1)[1:3])
AUC_3 <- sprintf("%.2f",ci.auc(c1)[1:3])
AUC_4 <- sprintf("%.2f",ci.auc(d1)[1:3])
AUC_5 <- sprintf("%.2f",ci.auc(e1)[1:3])
AUC_6 <- sprintf("%.2f",ci.auc(f1)[1:3])
AUC_7 <- sprintf("%.2f",ci.auc(g1)[1:3])
AUC_8 <- sprintf("%.2f",ci.auc(h1)[1:3])

Model1 <- paste("Beef: ",AUC_1[2], " (", AUC_1[1],"-", AUC_1[3],")",sep="")
Model2 <- paste("Chicken: ",AUC_2[2], " (", AUC_2[1],"-", AUC_2[3],")",sep="")
Model3 <- paste("Cheese: ",AUC_3[2], " (", AUC_3[1],"-", AUC_3[3],")",sep="")
Model4 <- paste("Yogurt: ",AUC_4[2], " (", AUC_4[1],"-", AUC_4[3],")",sep="")
Model5 <- paste("Corn: ",AUC_5[2], " (", AUC_5[1],"-", AUC_5[3],")",sep="")
Model6 <- paste("Oats: ",AUC_6[2], " (", AUC_6[1],"-", AUC_6[3],")",sep="")
Model7 <- paste("Bread: ",AUC_7[2], " (", AUC_7[1],"-", AUC_7[3],")",sep="")
Model8 <- paste("Potato: ",AUC_8[2], " (", AUC_8[1],"-", AUC_8[3],")",sep="")

legend(x=0.45,y=0.40,c(Model1,Model2,Model3,Model4,Model5,Model6,Model7,Model8),lty=1,col=c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"),lwd=3,cex=0.75,bty = "n")
dev.off()

#plot urine ROC curve
png("ROC_Urine.png",width=400,height=400)
par(family = "Calibri") 

plot(a0,grid=F, legacy.axes=T,col="#FA5555", lty=2,family="Calibri")
plot(b0,grid=F, legacy.axes=T,col="#4F81BD", lty=2,family="Calibri", add = TRUE)
plot(c0,grid=F, legacy.axes=T,col="#cc4bc4", lty=2,family="Calibri", add = TRUE)
plot(d0,grid=F, legacy.axes=T,col="#9f81f7", lty=2,family="Calibri", add = TRUE)
plot(e0,grid=F, legacy.axes=T,col="#orange", lty=2,family="Calibri", add = TRUE)
plot(f0,grid=F, legacy.axes=T,col="#B15928", lty=2,family="Calibri", add = TRUE)
plot(g0,grid=F, legacy.axes=T,col="grey", lty=2,family="Calibri", add = TRUE)
plot(h0,grid=F, legacy.axes=T,col="pink", lty=2,family="Calibri", add = TRUE)

plot(a2,grid=F, legacy.axes=T,col="#FA5555", lty=1,family="Calibri", add = TRUE)
plot(b2,grid=F, legacy.axes=T,col="#4F81BD", lty=1,family="Calibri", add = TRUE)
plot(c2,grid=F, legacy.axes=T,col="#cc4bc4", lty=1,family="Calibri", add = TRUE)
plot(d2,grid=F, legacy.axes=T,col="#9f81f7", lty=1,family="Calibri", add = TRUE)
plot(e2,grid=F, legacy.axes=T,col="#orange", lty=1,family="Calibri", add = TRUE)
plot(f2,grid=F, legacy.axes=T,col="#B15928", lty=1,family="Calibri", add = TRUE)
plot(g2,grid=F, legacy.axes=T,col="grey", lty=1,family="Calibri", add = TRUE)
plot(h2,grid=F, legacy.axes=T,col="pink", lty=1,family="Calibri", add = TRUE)

AUC_1 <- sprintf("%.2f",ci.auc(a0)[1:3])
AUC_2 <- sprintf("%.2f",ci.auc(b0)[1:3])
AUC_3 <- sprintf("%.2f",ci.auc(c0)[1:3])
AUC_4 <- sprintf("%.2f",ci.auc(d0)[1:3])
AUC_5 <- sprintf("%.2f",ci.auc(e0)[1:3])
AUC_6 <- sprintf("%.2f",ci.auc(f0)[1:3])
AUC_7 <- sprintf("%.2f",ci.auc(g0)[1:3])
AUC_8 <- sprintf("%.2f",ci.auc(h0)[1:3])

Model1 <- paste("Beef: ",AUC_1[2], " (", AUC_1[1],"-", AUC_1[3],")",sep="")
Model2 <- paste("Chicken: ",AUC_2[2], " (", AUC_2[1],"-", AUC_2[3],")",sep="")
Model3 <- paste("Cheese: ",AUC_3[2], " (", AUC_3[1],"-", AUC_3[3],")",sep="")
Model4 <- paste("Yogurt: ",AUC_4[2], " (", AUC_4[1],"-", AUC_4[3],")",sep="")
Model5 <- paste("Corn: ",AUC_5[2], " (", AUC_5[1],"-", AUC_5[3],")",sep="")
Model6 <- paste("Oats: ",AUC_6[2], " (", AUC_6[1],"-", AUC_6[3],")",sep="")
Model7 <- paste("Bread: ",AUC_7[2], " (", AUC_7[1],"-", AUC_7[3],")",sep="")
Model8 <- paste("Potato: ",AUC_8[2], " (", AUC_8[1],"-", AUC_8[3],")",sep="")

legend(x=0.90,y=0.30,c(Model1,Model2,Model3,Model4,Model5,Model6,Model7,Model8),lty=2,col=c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"),lwd=3,cex=0.75,bty = "n")

AUC_1 <- sprintf("%.2f",ci.auc(a2)[1:3])
AUC_2 <- sprintf("%.2f",ci.auc(b2)[1:3])
AUC_3 <- sprintf("%.2f",ci.auc(c2)[1:3])
AUC_4 <- sprintf("%.2f",ci.auc(d2)[1:3])
AUC_5 <- sprintf("%.2f",ci.auc(e2)[1:3])
AUC_6 <- sprintf("%.2f",ci.auc(f2)[1:3])
AUC_7 <- sprintf("%.2f",ci.auc(g2)[1:3])
AUC_8 <- sprintf("%.2f",ci.auc(h2)[1:3])

Model1 <- paste("Beef: ",AUC_1[2], " (", AUC_1[1],"-", AUC_1[3],")",sep="")
Model2 <- paste("Chicken: ",AUC_2[2], " (", AUC_2[1],"-", AUC_2[3],")",sep="")
Model3 <- paste("Cheese: ",AUC_3[2], " (", AUC_3[1],"-", AUC_3[3],")",sep="")
Model4 <- paste("Yogurt: ",AUC_4[2], " (", AUC_4[1],"-", AUC_4[3],")",sep="")
Model5 <- paste("Corn: ",AUC_5[2], " (", AUC_5[1],"-", AUC_5[3],")",sep="")
Model6 <- paste("Oats: ",AUC_6[2], " (", AUC_6[1],"-", AUC_6[3],")",sep="")
Model7 <- paste("Bread: ",AUC_7[2], " (", AUC_7[1],"-", AUC_7[3],")",sep="")
Model8 <- paste("Potato: ",AUC_8[2], " (", AUC_8[1],"-", AUC_8[3],")",sep="")

legend(x=0.45,y=0.30,c(Model1,Model2,Model3,Model4,Model5,Model6,Model7,Model8),lty=1,col=c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"),lwd=3,cex=0.75,bty = "n")
dev.off()

#plot plasma and urine ROC curve
png("ROC_Combined.png",width=450,height=450)
par(family = "Calibri") 

plot(a0,grid=F, legacy.axes=T,col="#FA5555", lty=2,family="Calibri")
plot(b0,grid=F, legacy.axes=T,col="#4F81BD", lty=2,family="Calibri", add = TRUE)
plot(c0,grid=F, legacy.axes=T,col="#cc4bc4", lty=2,family="Calibri", add = TRUE)
plot(d0,grid=F, legacy.axes=T,col="#9f81f7", lty=2,family="Calibri", add = TRUE)
plot(e0,grid=F, legacy.axes=T,col="#orange", lty=2,family="Calibri", add = TRUE)
plot(f0,grid=F, legacy.axes=T,col="#B15928", lty=2,family="Calibri", add = TRUE)
plot(g0,grid=F, legacy.axes=T,col="grey", lty=2,family="Calibri", add = TRUE)
plot(h0,grid=F, legacy.axes=T,col="pink", lty=2,family="Calibri", add = TRUE)

plot(a3,grid=F, legacy.axes=T,col="#FA5555", lty=1,family="Calibri", add = TRUE)
plot(b3,grid=F, legacy.axes=T,col="#4F81BD", lty=1,family="Calibri", add = TRUE)
plot(c3,grid=F, legacy.axes=T,col="#cc4bc4", lty=1,family="Calibri", add = TRUE)
plot(d3,grid=F, legacy.axes=T,col="#9f81f7", lty=1,family="Calibri", add = TRUE)
plot(e3,grid=F, legacy.axes=T,col="#orange", lty=1,family="Calibri", add = TRUE)
plot(f3,grid=F, legacy.axes=T,col="#B15928", lty=1,family="Calibri", add = TRUE)
plot(g3,grid=F, legacy.axes=T,col="grey", lty=1,family="Calibri", add = TRUE)
plot(h3,grid=F, legacy.axes=T,col="pink", lty=1,family="Calibri", add = TRUE)

AUC_1 <- sprintf("%.2f",ci.auc(a0)[1:3])
AUC_2 <- sprintf("%.2f",ci.auc(b0)[1:3])
AUC_3 <- sprintf("%.2f",ci.auc(c0)[1:3])
AUC_4 <- sprintf("%.2f",ci.auc(d0)[1:3])
AUC_5 <- sprintf("%.2f",ci.auc(e0)[1:3])
AUC_6 <- sprintf("%.2f",ci.auc(f0)[1:3])
AUC_7 <- sprintf("%.2f",ci.auc(g0)[1:3])
AUC_8 <- sprintf("%.2f",ci.auc(h0)[1:3])

Model1 <- paste("Beef: ",AUC_1[2], " (", AUC_1[1],"-", AUC_1[3],")",sep="")
Model2 <- paste("Chicken: ",AUC_2[2], " (", AUC_2[1],"-", AUC_2[3],")",sep="")
Model3 <- paste("Cheese: ",AUC_3[2], " (", AUC_3[1],"-", AUC_3[3],")",sep="")
Model4 <- paste("Yogurt: ",AUC_4[2], " (", AUC_4[1],"-", AUC_4[3],")",sep="")
Model5 <- paste("Corn: ",AUC_5[2], " (", AUC_5[1],"-", AUC_5[3],")",sep="")
Model6 <- paste("Oats: ",AUC_6[2], " (", AUC_6[1],"-", AUC_6[3],")",sep="")
Model7 <- paste("Bread: ",AUC_7[2], " (", AUC_7[1],"-", AUC_7[3],")",sep="")
Model8 <- paste("Potato: ",AUC_8[2], " (", AUC_8[1],"-", AUC_8[3],")",sep="")

legend(x=0.90,y=0.30,c(Model1,Model2,Model3,Model4,Model5,Model6,Model7,Model8),lty=2,col=c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"),lwd=3,cex=0.75,bty = "n")

AUC_1 <- sprintf("%.2f",ci.auc(a3)[1:3])
AUC_2 <- sprintf("%.2f",ci.auc(b3)[1:3])
AUC_3 <- sprintf("%.2f",ci.auc(c3)[1:3])
AUC_4 <- sprintf("%.2f",ci.auc(d3)[1:3])
AUC_5 <- sprintf("%.2f",ci.auc(e3)[1:3])
AUC_6 <- sprintf("%.2f",ci.auc(f3)[1:3])
AUC_7 <- sprintf("%.2f",ci.auc(g3)[1:3])
AUC_8 <- sprintf("%.2f",ci.auc(h3)[1:3])

Model1 <- paste("Beef: ",AUC_1[2], " (", AUC_1[1],"-", AUC_1[3],")",sep="")
Model2 <- paste("Chicken: ",AUC_2[2], " (", AUC_2[1],"-", AUC_2[3],")",sep="")
Model3 <- paste("Cheese: ",AUC_3[2], " (", AUC_3[1],"-", AUC_3[3],")",sep="")
Model4 <- paste("Yogurt: ",AUC_4[2], " (", AUC_4[1],"-", AUC_4[3],")",sep="")
Model5 <- paste("Corn: ",AUC_5[2], " (", AUC_5[1],"-", AUC_5[3],")",sep="")
Model6 <- paste("Oats: ",AUC_6[2], " (", AUC_6[1],"-", AUC_6[3],")",sep="")
Model7 <- paste("Bread: ",AUC_7[2], " (", AUC_7[1],"-", AUC_7[3],")",sep="")
Model8 <- paste("Potato: ",AUC_8[2], " (", AUC_8[1],"-", AUC_8[3],")",sep="")

legend(x=0.45,y=0.30,c(Model1,Model2,Model3,Model4,Model5,Model6,Model7,Model8),lty=1,col=c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"),lwd=3,cex=0.75,bty = "n")
dev.off()

#calculate AUC based on plasma individual metabolites
data <- data_use
var <- use1

res_bf <- run_cv_metabolites(data = data,"Ifbf",metabolite_vars = var[which(var$group=="Beef"),]$name,r=5)
res_ck <- run_cv_metabolites(data = data,"Ifck",metabolite_vars = var[which(var$group=="Chicken"),]$name,r=5)
res_ch <- run_cv_metabolites(data = data,"Ifch",metabolite_vars = var[which(var$group=="Cheese"),]$name,r=5)
res_yg <- run_cv_metabolites(data = data,"Ifyg",metabolite_vars = var[which(var$group=="Yogurt"),]$name,r=5)
res_cr <- run_cv_metabolites(data = data,"Ifcr",metabolite_vars = var[which(var$group=="Corn"),]$name,r=5)
res_ot <- run_cv_metabolites(data = data,"Ifot",metabolite_vars = var[which(var$group=="Oats"),]$name,r=5)
res_ww <- run_cv_metabolites(data = data,"Ifww",metabolite_vars = var[which(var$group=="Bread"),]$name,r=5)
res_pt <- run_cv_metabolites(data = data,"Ifpt",metabolite_vars = var[which(var$group=="Potato"),]$name,r=5)

res_bf$group <- "Beef"
res_ck$group <- "Chicken"
res_ch$group <- "Cheese"
res_yg$group <- "Yogurt"
res_cr$group <- "Corn"
res_ot$group <- "Oats"
res_ww$group <- "Bread"
res_pt$group <- "Potato"

res_bf <- left_join(res_bf,var,by=c("group","name"))
res_ck <- left_join(res_ck,var,by=c("group","name"))
res_ch <- left_join(res_ch,var,by=c("group","name"))
res_yg <- left_join(res_yg,var,by=c("group","name"))
res_cr <- left_join(res_cr,var,by=c("group","name"))
res_ot <- left_join(res_ot,var,by=c("group","name"))
res_ww <- left_join(res_ww,var,by=c("group","name"))
res_pt <- left_join(res_pt,var,by=c("group","name"))

res <- rbind(res_bf,res_ck,res_ch,res_yg,res_cr,res_ot,res_ww,res_pt) 

res$Label <- paste(res$Metabolite,"_",res$time,"_",res$type,sep='')

res$group <- factor(res$group,levels=rev(c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Bread","Potato")))

res$id <- seq(1, nrow(res))
res$id <- as.factor(res$id)

res <- res %>% arrange(group,AUC)
res$id <- factor(res$id,levels = res$id)

res1 <- res
res1$AUC <- as.numeric(res1$AUC)

res1$group <- factor(res1$group,levels=rev(c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Bread","Potato")))

p1 <- ggplot(res1, aes(x=id,y=AUC,fill=group)) +       
  geom_col(position = "dodge", width = 0.5) + 
  scale_fill_manual(values=rev(c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"))) +
  geom_text(data=res1, aes(x=id, y=AUC+0.02, label=round(AUC,2)), color="black", family="calibri",size=2.5) +
  labs(title = "", x = "", y = "AUC") +
  scale_x_discrete(breaks = res1$id,labels = res1$Label) +
  theme_classic() +
  theme(strip.text=element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.text.y = element_text(size=8,color="black",family="calibri"),
        axis.text.x = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.title.x = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.title = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        panel.spacing.x = unit(0.01, "cm"),
        panel.spacing.y = unit(0.01, "cm"),
        legend.position = "none",
        legend.title = element_text(color="black",family="calibri",size=5),
        legend.text = element_text(color="black",family="calibri",size=5))+
  coord_flip()

#calculate AUC based on urine individual metabolites
data <- FULLDS
var <- use2

res_bf <- run_cv_metabolites(data = data,"Ifbf",metabolite_vars = var[which(var$group=="Beef"),]$name,r=5)
res_ck <- run_cv_metabolites(data = data,"Ifck",metabolite_vars = var[which(var$group=="Chicken"),]$name,r=5)
res_ch <- run_cv_metabolites(data = data,"Ifch",metabolite_vars = var[which(var$group=="Cheese"),]$name,r=5)
res_yg <- run_cv_metabolites(data = data,"Ifyg",metabolite_vars = var[which(var$group=="Yogurt"),]$name,r=5)
res_cr <- run_cv_metabolites(data = data,"Ifcr",metabolite_vars = var[which(var$group=="Corn"),]$name,r=5)
res_ot <- run_cv_metabolites(data = data,"Ifot",metabolite_vars = var[which(var$group=="Oats"),]$name,r=5)
res_ww <- run_cv_metabolites(data = data,"Ifww",metabolite_vars = var[which(var$group=="Bread"),]$name,r=5)
res_pt <- run_cv_metabolites(data = data,"Ifpt",metabolite_vars = var[which(var$group=="Potato"),]$name,r=5)

res_bf$group <- "Beef"
res_ck$group <- "Chicken"
res_ch$group <- "Cheese"
res_yg$group <- "Yogurt"
res_cr$group <- "Corn"
res_ot$group <- "Oats"
res_ww$group <- "Bread"
res_pt$group <- "Potato"

res_bf <- left_join(res_bf,var,by=c("group","name"))
res_ck <- left_join(res_ck,var,by=c("group","name"))
res_ch <- left_join(res_ch,var,by=c("group","name"))
res_yg <- left_join(res_yg,var,by=c("group","name"))
res_cr <- left_join(res_cr,var,by=c("group","name"))
res_ot <- left_join(res_ot,var,by=c("group","name"))
res_ww <- left_join(res_ww,var,by=c("group","name"))
res_pt <- left_join(res_pt,var,by=c("group","name"))

res <- rbind(res_bf,res_ck,res_ch,res_yg,res_cr,res_ot,res_ww,res_pt) 

res$Label <- paste(res$Metabolite,"_",res$time,"_",res$type,sep='')

res$group <- factor(res$group,levels=rev(c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Bread","Potato")))

res$id <- seq(1, nrow(res))
res$id <- as.factor(res$id)

res <- res %>% arrange(group,AUC)
res$id <- factor(res$id,levels = res$id)

res2 <- res
res2$AUC <- as.numeric(res2$AUC)

p2 <- ggplot(res2, aes(x=id, y=AUC,fill=group)) +       
  geom_col(position = "dodge", width = 0.5) + 
  scale_fill_manual(values=rev(c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"))) +
  geom_text(data=res2, aes(x=id, y=AUC+0.02, label=round(AUC,2)), color="black", family="calibri",size=2.5) +
  labs(title = "", x = "", y = "AUC") +
  scale_x_discrete(breaks = res2$id,labels = res2$Label) +
  theme_classic() +
  theme(strip.text=element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.text.y = element_text(size=8,color="black",family="calibri"),
        axis.text.x = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.title.x = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.title = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        panel.spacing.x = unit(0.01, "cm"),
        panel.spacing.y = unit(0.01, "cm"),
        legend.position = "none",
        legend.title = element_text(color="black",family="calibri",size=5),
        legend.text = element_text(color="black",family="calibri",size=5))+
  coord_flip()

#calculate AUC based on plasma and urine individual metabolites
data <- data_use
var <- use3

res_bf <- run_cv_metabolites(data = data,"Ifbf",metabolite_vars = var[which(var$group=="Beef"),]$name,r=5)
res_ck <- run_cv_metabolites(data = data,"Ifck",metabolite_vars = var[which(var$group=="Chicken"),]$name,r=5)
res_ch <- run_cv_metabolites(data = data,"Ifch",metabolite_vars = var[which(var$group=="Cheese"),]$name,r=5)
res_yg <- run_cv_metabolites(data = data,"Ifyg",metabolite_vars = var[which(var$group=="Yogurt"),]$name,r=5)
res_cr <- run_cv_metabolites(data = data,"Ifcr",metabolite_vars = var[which(var$group=="Corn"),]$name,r=5)
res_ot <- run_cv_metabolites(data = data,"Ifot",metabolite_vars = var[which(var$group=="Oats"),]$name,r=5)
res_ww <- run_cv_metabolites(data = data,"Ifww",metabolite_vars = var[which(var$group=="Bread"),]$name,r=5)
res_pt <- run_cv_metabolites(data = data,"Ifpt",metabolite_vars = var[which(var$group=="Potato"),]$name,r=5)

res_bf$group <- "Beef"
res_ck$group <- "Chicken"
res_ch$group <- "Cheese"
res_yg$group <- "Yogurt"
res_cr$group <- "Corn"
res_ot$group <- "Oats"
res_ww$group <- "Bread"
res_pt$group <- "Potato"

res_bf <- left_join(res_bf,var,by=c("group","name"))
res_ck <- left_join(res_ck,var,by=c("group","name"))
res_ch <- left_join(res_ch,var,by=c("group","name"))
res_yg <- left_join(res_yg,var,by=c("group","name"))
res_cr <- left_join(res_cr,var,by=c("group","name"))
res_ot <- left_join(res_ot,var,by=c("group","name"))
res_ww <- left_join(res_ww,var,by=c("group","name"))
res_pt <- left_join(res_pt,var,by=c("group","name"))

res <- rbind(res_bf,res_ck,res_ch,res_yg,res_cr,res_ot,res_ww,res_pt) 

res$Label <- paste(res$Metabolite,"_",res$time,"_",res$type,sep='')

res$group <- factor(res$group,levels=rev(c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Bread","Potato")))

res$id <- seq(1, nrow(res))
res$id <- as.factor(res$id)

res <- res %>% arrange(group,AUC)
res$id <- factor(res$id,levels = res$id)

res3 <- res
res3$AUC <- as.numeric(res3$AUC)

p3 <- ggplot(res3, aes(x=id, y=AUC,fill=group)) +       
  geom_col(position = "dodge", width = 0.5) + 
  scale_fill_manual(values=rev(c("#FA5555","#4F81BD","#cc4bc4","#9f81f7","orange","#B15928","grey","pink"))) +
  geom_text(data=res3, aes(x=id, y=AUC+0.02, label=round(AUC,2)), color="black", family="calibri",size=2.5) +
  labs(title = "", x = "", y = "AUC") +
  scale_x_discrete(breaks = res3$id,labels = res3$Label) +
  theme_classic() +
  theme(strip.text=element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.text.y = element_text(size=8,color="black",family="calibri"),
        axis.text.x = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.title.x = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        axis.title = element_text(size=8,hjust=0.5,color="black",family="calibri"),
        panel.spacing.x = unit(0.01, "cm"),
        panel.spacing.y = unit(0.01, "cm"),
        legend.position = "none",
        legend.title = element_text(color="black",family="calibri",size=5),
        legend.text = element_text(color="black",family="calibri",size=5))+
  coord_flip()

p <- plot_grid(p1,p2,p3,nrow=1)

ggsave("AUC_Ind.png", plot = p, width = 13, height = 9, dpi = 600)