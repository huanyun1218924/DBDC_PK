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

#--------------------------------------------------------------------------------------------
#
#           Step 1: re-organize data
#
#--------------------------------------------------------------------------------------------
#read plasma data
load("plasma_sample.RData")

plasma_sample1 <- plasma_sample[which(plasma_sample$time=="0"),c("ID","group",metvar)]
plasma_sample2 <- plasma_sample[which(plasma_sample$time=="1"),c("ID","group",metvar)]
plasma_sample3 <- plasma_sample[which(plasma_sample$time=="2"),c("ID","group",metvar)]
plasma_sample4 <- plasma_sample[which(plasma_sample$time=="3"),c("ID","group",metvar)]
plasma_sample5 <- plasma_sample[which(plasma_sample$time=="4"),c("ID","group",metvar)]
plasma_sample6 <- plasma_sample[which(plasma_sample$time=="5"),c("ID","group",metvar)]
plasma_sample7 <- plasma_sample[which(plasma_sample$time=="6"),c("ID","group",metvar)]
plasma_sample8 <- plasma_sample[which(plasma_sample$time=="7"),c("ID","group",metvar)]
plasma_sample9 <- plasma_sample[which(plasma_sample$time=="8"),c("ID","group",metvar)]
plasma_sample10 <- plasma_sample[which(plasma_sample$time=="24"),c("ID","group",metvar)]

names(plasma_sample1)[3:ncol(plasma_sample1)] <- paste(metvar,"_t0",sep='')
names(plasma_sample2)[3:ncol(plasma_sample2)] <- paste(metvar,"_t1",sep='')
names(plasma_sample3)[3:ncol(plasma_sample3)] <- paste(metvar,"_t2",sep='')
names(plasma_sample4)[3:ncol(plasma_sample4)] <- paste(metvar,"_t3",sep='')
names(plasma_sample5)[3:ncol(plasma_sample5)] <- paste(metvar,"_t4",sep='')
names(plasma_sample6)[3:ncol(plasma_sample6)] <- paste(metvar,"_t5",sep='')
names(plasma_sample7)[3:ncol(plasma_sample7)] <- paste(metvar,"_t6",sep='')
names(plasma_sample8)[3:ncol(plasma_sample8)] <- paste(metvar,"_t7",sep='')
names(plasma_sample9)[3:ncol(plasma_sample9)] <- paste(metvar,"_t8",sep='')
names(plasma_sample10)[3:ncol(plasma_sample10)] <- paste(metvar,"_t24",sep='')

plasma_sample_new <- full_join(plasma_sample1,plasma_sample2,by=c("group","ID")) %>%
  full_join(plasma_sample3,by=c("group","ID")) %>%
  full_join(plasma_sample4,by=c("group","ID")) %>%
  full_join(plasma_sample5,by=c("group","ID")) %>%
  full_join(plasma_sample6,by=c("group","ID")) %>%
  full_join(plasma_sample7,by=c("group","ID")) %>%
  full_join(plasma_sample8,by=c("group","ID")) %>%
  full_join(plasma_sample9,by=c("group","ID")) %>%
  full_join(plasma_sample10,by=c("group","ID"))

#read urine data
load("urine_sample.RData")

urine_sample1 <- urine_sample[which(urine_sample$time=="0"),c("ID","group",metvar)]
urine_sample2 <- urine_sample[which(urine_sample$time=="2"),c("ID","group",metvar)]
urine_sample3 <- urine_sample[which(urine_sample$time=="4"),c("ID","group",metvar)]
urine_sample4 <- urine_sample[which(urine_sample$time=="6"),c("ID","group",metvar)]
urine_sample5 <- urine_sample[which(urine_sample$time=="8"),c("ID","group",metvar)]
urine_sample6 <- urine_sample[which(urine_sample$time=="24"),c("ID","group",metvar)]

names(urine_sample1)[3:ncol(urine_sample1)] <- paste(metvar,"_t0_u",sep='')
names(urine_sample2)[3:ncol(urine_sample2)] <- paste(metvar,"_t2_u",sep='')
names(urine_sample3)[3:ncol(urine_sample3)] <- paste(metvar,"_t4_u",sep='')
names(urine_sample4)[3:ncol(urine_sample4)] <- paste(metvar,"_t6_u",sep='')
names(urine_sample5)[3:ncol(urine_sample5)] <- paste(metvar,"_t8_u",sep='')
names(urine_sample6)[3:ncol(urine_sample6)] <- paste(metvar,"_t24_u",sep='')

urine_sample_new <- full_join(urine_sample1,urine_sample2,by=c("group","ID")) %>%
  full_join(urine_sample3,by=c("group","ID")) %>%
  full_join(urine_sample4,by=c("group","ID")) %>%
  full_join(urine_sample5,by=c("group","ID")) %>%
  full_join(urine_sample6,by=c("group","ID"))

#combine plasma and urine data together
com_dat <-  full_join(plasma_sample_new,urine_sample_new,by=c("group","ID"))

#define metabolite list for use
var_plasma <- names(plasma_sample_new)[3:ncol(plasma_sample_new)] #plasma
var_urine <- names(urine_sample_new)[3:ncol(urine_sample_new)]    #urine
var_all <- c(var_plasma,var_urine)                                #combined

#define outcome 
data_use <- com_dat %>% left_join(plasma_sample[!duplicated(plasma_sample$ID),c("ID","ips_age", "ps_sex", "race", "ips_bmi")],by="ID")

data_use$Ifbf <- ifelse(data_use$group=="Beef","1","0")
data_use$Ifck <- ifelse(data_use$group=="Chicken","1","0")
data_use$Ifyg <- ifelse(data_use$group=="Yogurt","1","0")
data_use$Ifch <- ifelse(data_use$group=="Cheese","1","0")
data_use$Ifcr <- ifelse(data_use$group=="Corn","1","0")
data_use$Ifot <- ifelse(data_use$group=="Oats","1","0")
data_use$Ifww <- ifelse(data_use$group=="Whole Wheat Bread","1","0")
data_use$Ifpt <- ifelse(data_use$group=="Potato","1","0")

out <- c("Ifbf","Ifck","Ifyg","Ifch","Ifcr","Ifot","Ifww","Ifpt")

for (i in out){
  data_use[,i] <- as.factor(data_use[,i])
}

#impute covariates: median value for continuous variable
cov <- c("ips_age", "ps_sex", "race", "ips_bmi")

data_use[which(is.na(data_use$ips_age)),"ips_age"] <- 31
data_use[which(is.na(data_use$ips_bmi)),"ips_bmi"] <- 25.4
data_use[which(is.na(data_use$ps_sex)),"ps_sex"] <- 2
data_use[which(is.na(data_use$race)),"race"] <- 1

#standarize each metabolites for elastic net regression
data_use = data_use[,c(cov,var_all,out)] 
data_use = na.omit(data_use)
data_use[var_all] <- apply(data_use[var_all],2,scale)

#--------------------------------------------------------------------------------------------
#
#           Step 2: run elastic net model to identify metablites differentiating foods
#
#--------------------------------------------------------------------------------------------
z_bi = "predictScore"

#run elastic net using plasma data
for (i in out){
  y_bi=i
  fit <- elasfunc(data_use, cov, var_plasma, y_bi, z_bi)
}

#run elastic net using urine data
for (i in out){
  y_bi=i
  fit <- elasfunc(data_use, cov, var_urine, y_bi, z_bi)
}

#run elastic net using both plasma and urine data
for (i in out){
  y_bi=i
  fit <- elasfunc(data_use, cov, var_all, y_bi, z_bi)
}

#--------------------------------------------------------------------------------------------
#
#            Step 3: calculate AUC based on individual metabolites
#
#--------------------------------------------------------------------------------------------
#calculate AUC based on plasma individual metabolites
res_bf <- run_cv_metabolites(data = data_use,"Ifbf",metabolite_vars = var_plasma,r=5)
res_ck <- run_cv_metabolites(data = data_use,"Ifck",metabolite_vars = var_plasma,r=5)
res_ch <- run_cv_metabolites(data = data_use,"Ifch",metabolite_vars = var_plasma,r=5)
res_yg <- run_cv_metabolites(data = data_use,"Ifyg",metabolite_vars = var_plasma,r=5)
res_cr <- run_cv_metabolites(data = data_use,"Ifcr",metabolite_vars = var_plasma,r=5)
res_ot <- run_cv_metabolites(data = data_use,"Ifot",metabolite_vars = var_plasma,r=5)
res_ww <- run_cv_metabolites(data = data_use,"Ifww",metabolite_vars = var_plasma,r=5)
res_pt <- run_cv_metabolites(data = data_use,"Ifpt",metabolite_vars = var_plasma,r=5)

res_bf$group <- "Beef"
res_ck$group <- "Chicken"
res_ch$group <- "Cheese"
res_yg$group <- "Yogurt"
res_cr$group <- "Corn"
res_ot$group <- "Oats"
res_ww$group <- "Bread"
res_pt$group <- "Potato"

res_bf <- left_join(res_bf,var_plasma,by=c("group","name"))
res_ck <- left_join(res_ck,var_plasma,by=c("group","name"))
res_ch <- left_join(res_ch,var_plasma,by=c("group","name"))
res_yg <- left_join(res_yg,var_plasma,by=c("group","name"))
res_cr <- left_join(res_cr,var_plasma,by=c("group","name"))
res_ot <- left_join(res_ot,var_plasma,by=c("group","name"))
res_ww <- left_join(res_ww,var_plasma,by=c("group","name"))
res_pt <- left_join(res_pt,var_plasma,by=c("group","name"))

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
res_bf <- run_cv_metabolites(data = data_use,"Ifbf",metabolite_vars = var_urine,r=5)
res_ck <- run_cv_metabolites(data = data_use,"Ifck",metabolite_vars = var_urine,r=5)
res_ch <- run_cv_metabolites(data = data_use,"Ifch",metabolite_vars = var_urine,r=5)
res_yg <- run_cv_metabolites(data = data_use,"Ifyg",metabolite_vars = var_urine,r=5)
res_cr <- run_cv_metabolites(data = data_use,"Ifcr",metabolite_vars = var_urine,r=5)
res_ot <- run_cv_metabolites(data = data_use,"Ifot",metabolite_vars = var_urine,r=5)
res_ww <- run_cv_metabolites(data = data_use,"Ifww",metabolite_vars = var_urine,r=5)
res_pt <- run_cv_metabolites(data = data_use,"Ifpt",metabolite_vars = var_urine,r=5)

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
res_bf <- run_cv_metabolites(data = data_use,"Ifbf",metabolite_vars = var_all,r=5)
res_ck <- run_cv_metabolites(data = data_use,"Ifck",metabolite_vars = var_all,r=5)
res_ch <- run_cv_metabolites(data = data_use,"Ifch",metabolite_vars = var_all,r=5)
res_yg <- run_cv_metabolites(data = data_use,"Ifyg",metabolite_vars = var_all,r=5)
res_cr <- run_cv_metabolites(data = data_use,"Ifcr",metabolite_vars = var_all,r=5)
res_ot <- run_cv_metabolites(data = data_use,"Ifot",metabolite_vars = var_all,r=5)
res_ww <- run_cv_metabolites(data = data_use,"Ifww",metabolite_vars = var_all,r=5)
res_pt <- run_cv_metabolites(data = data_use,"Ifpt",metabolite_vars = var_all,r=5)

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

ggsave("AUC_Ind.png", plot = p, width = 13, height = 9, dpi = 600)
