#install packages used for time-series clustering analysis
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Mfuzz") 

#load packages based on R version 4.3.3
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(grid)
library(gridExtra)
library(Mfuzz)
library(circlize)
library(ComplexHeatmap)

#define function for use
geo_mean <- function(x) {
  if (any(x <= 0)) {
    return(NA)       
  }
  exp(mean(log(x)))
}  #note: this is used for calcualting geometrix mean

#--------------------------------------------------------------------------------------------
#
#           part 1: run cluster analysis for plasma metabolites
#
#--------------------------------------------------------------------------------------------
#read plasma data
load("plasma_sample.RData") #note: it contains the simulated plasma data same with raw data and its annotation file

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

#calculate geometric mean across each metabolite by group and by time
results <- plasma_sample %>% dplyr::group_by(group, time) %>% dplyr::summarise(dplyr::across(dplyr::all_of(metvar),geo_mean,.names = "geo_mean_{.col}"),.groups = "drop")
names(results)[3:ncol(plasma_sample)] <- metvar

results[3:ncol(results)] <- apply(results[3:ncol(results)],2,function(x)as.numeric(x)) #make sure each metabolite is numeric type

#define food name
var <- c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Whole Wheat Bread","Potato")

#run cluster analysis across each food 
for (i in var){

  #define dataset for use
  use <- results[results$group==i,c("time",metvar)] %>% t() %>% as.data.frame()
  names(use) <- use[1,]
  use <- use[-1,]
  names(use) <- c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t24") #10 timepoints
  
  for (j in 1:10){
    use[,j] <- as.numeric(use[,j])
  }
  
  pd <- data.frame(time = c("0","1","2","3","4","5","6","7","8","24"))
  rownames(pd) <- colnames(use)
  
  #standarize dataset for clustering analysis
  eset <- new("ExpressionSet", exprs = as.matrix(use))
  keep <- rowSums(is.na(exprs(eset))) < ncol(exprs(eset)) - 1
  eset_filtered <- eset[keep, ]
  eset_std <- standardise(eset_filtered)
  
  keep2 <- rowSums(is.na(exprs(eset_std))) < ncol(exprs(eset_std)) - 1
  eset_filtered <- eset_std[keep2, ]
  
  #run model to get the optimized cluster number
  m <- mestimate(eset_std)                                      # Estimate fuzziness parameter
  c  <- which.min(Dmin(eset_filtered, m = m, crange = 2:10))    # Plot to help choose cluster number
  cl <- mfuzz(eset_filtered, c = c, m=m)                        # c = number of clusters
  
  #create membership object
  membership_df <- cl$membership %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "membership") %>%
    group_by(gene) %>%
    slice_max(membership, with_ties = FALSE) 
  
  #extract sample size info
  size <- as.data.frame(cl$size)
  size$cluster <- seq_len(nrow(size))
  size$cluster <- as.numeric(size$cluster)
  names(size)[1] <-  "Num"
  
  #extract cluster assignments info
  exprs_long <- exprs(eset_filtered) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression")
  
  exprs_long <- exprs_long %>%
    left_join(data.frame(gene = names(cl$cluster), cluster = cl$cluster), by = "gene") %>%
    left_join(pd %>% rownames_to_column("sample"), by = "sample") %>%
    left_join(membership_df[,-2], by = "gene") %>%
    left_join(size,by="cluster")
  
  #plot cluster curve and save it
  exprs_long$time <- factor(exprs_long$time, levels=c("0","1","2","3","4","5","6","7","8","24"))
  exprs_long$cluster <- paste("Cluster",exprs_long$cluster,sep='')
  
  exprs_long$cluster <- factor(exprs_long$cluster,levels = unique(exprs_long$cluster[order(as.numeric(gsub("\\D", "", exprs_long$cluster)))]))
  
  p1 <-ggplot(exprs_long, aes(x = time, y = expression, group = gene)) +
    geom_line(alpha=0.3,color="grey40",size=0.7) +
    geom_smooth(aes(group = cluster, color = as.factor(cluster)), method = "loess", se = FALSE, size = 1.2) +
    scale_color_manual(values=c("#1f77b4", "#ff7f0e", "#2ca02c", "#FA5555", "#9467bd","#8c564b", "#bcbd22", "#17becf","#a6cee3", "#fdbf6f","pink","yellow","blue","orange")) +
    facet_wrap(~ cluster, scales = "free_y",nrow=5) +
    labs(x = "Time", y = "Standardized Intensity")+
    theme_minimal() +
    theme(strip.text=element_text(size=12,color="black",family="calibri"),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          axis.text.y = element_text(size=12,color="black",family="calibri"),
          axis.text.x = element_text(size=12,color="black",family="calibri"),
          axis.title = element_text(size=12,color="black",family="calibri"),
          axis.ticks = element_line(color="black"),
          panel.spacing.x = unit(0.01, "cm"),
          panel.spacing.y = unit(0.01, "cm"),
          legend.position = "none",
          legend.title = element_text(color="black",family="calibri",size=12),
          legend.text = element_text(color="black",family="calibri",size=12))
  
  ggsave(paste(i,"Cluster_Plasma.png",sep=''), plot = p1, width = 6, height =12, dpi = 600)
  
  #save results, including standarized dataset, sample size, gene membership
  sz <- cl$size %>% as.data.frame()
  gene_cluster <- cbind(cl$cluster, cl$membership)
  colnames(gene_cluster)[1] <- 'cluster'

  save(eset_filtered,sz,gene_cluster,file=paste("",i,"_Cluster_Plasma.RData",sep=''))
  
  #plot cluster using heatmap
  mem <- as.data.frame(gene_cluster) %>% arrange(cluster)
  mem$name <- rownames(mem)
  mem$cluster  <- factor(mem$cluster,levels=c("1","2","3","4","5","6","7","8","9","10"))
  
  rs <- eset_filtered@assayData[["exprs"]] %>% as.data.frame()
  rs <- rs[match(mem$name,rownames(rs)),]

  breaks <- c(-2,0,2)
  colors <- c("#006EBE", "white", "#FF6666")
  col_fun <- colorRamp2(breaks, colors)
  
  col_ha<-columnAnnotation(foo = names(rs),
                           col = list(foo = c("t0" = "#80CDC1", "t1" = "#87CEFA", "t2"="#D15FEE","t3" = "#EE6A50","t4"="#FFC125","t5"="pink","t6"="#440154", "t7"="#3E4A89", "t8"="#26828E","t24"="#B4DE2C")),gp = gpar(col = "black",family="Calibri",size=9),
                           annotation_name_side = "right",annotation_name_rot=0,
                           annotation_legend_param = list(foo = list(title = "Time",labels = c("0hr", "1hr", "2hr","3hr","4hr","5hr","6hr","7hr","8hr","24hr"))),
                           annotation_label = " ")
  
  row_ha<-rowAnnotation(foo = mem$cluster,
                        col = list(foo = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c", "4" = "#FA5555", "5" = "#9467bd","6" = "#8c564b", "7" = "#bcbd22", "8" = "#17becf","9" = "#a6cee3", "10" = "#fdbf6f")),gp = gpar(col = "black",family="Calibri",size=9),
                        annotation_name_side = "bottom",annotation_name_rot=0,
                        annotation_legend_param = list(foo = list(title = "Cluster")),
                        annotation_label = " ") #trait classification
  
  p2 <- Heatmap(rs, 
                col=col_fun,
                rect_gp = gpar(col = NA, lwd = 1),
                show_row_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                column_dend_reorder = TRUE,
                column_gap = unit(4, "mm"),
                row_gap = unit(2, "mm"),
                show_column_dend = TRUE,
                cluster_rows = F,
                cluster_columns = F,
                column_title = NULL,
                left_annotation = row_ha,
                top_annotation = col_ha,
                show_heatmap_legend = FALSE)
  
  ggsave(paste(i,"_Heatmap_Plasma.png",sep=''), plot = wrap_elements(grid.grabExpr(draw(p2))), width = 6, height = 12, dpi = 300) 
  
}

#check beef results
load("Beef_Cluster_Plasma.RData")

cl1 <- gene_cluster[which(gene_cluster$cluster==1),]
cl1$name <- rownames(cl1)
cl1 <- left_join(cl1,anno_plasma,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

cl2 <- gene_cluster[which(gene_cluster$cluster==2),]
cl2$name <- rownames(cl2)
cl2 <- left_join(cl2,anno_plasma,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

cl3 <- gene_cluster[which(gene_cluster$cluster==3),]
cl3$name <- rownames(cl3)
cl3 <- left_join(cl3,anno_plasma,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

cl8 <- gene_cluster[which(gene_cluster$cluster==8),]
cl8$name <- rownames(cl8)
cl8 <- left_join(cl8,anno_plasma,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

#clear all datasets for next-step urine analysis
rm(list = ls())
                                  
#--------------------------------------------------------------------------------------------
#
#           part 2: run cluster analysis for urinary metabolites
#
#--------------------------------------------------------------------------------------------
#read urine data
load("urine_sample.RData") #note: it contains the simulated urine data same with raw data and its annotation file

#check data
dim(urine_sample)
#662 631

#check covariates
covar

#check metabolites
metvar

#check group and time
table(urine_sample$group)

table(urine_sample$time)

#calculate geometric mean across each metabolite by group and by time
results <- urine_sample %>% dplyr::group_by(group, time) %>% dplyr::summarise(dplyr::across(dplyr::all_of(metvar),geo_mean,.names = "geo_mean_{.col}"),.groups = "drop")
names(results)[3:ncol(urine_sample)] <- metvar

results[3:ncol(results)] <- apply(results[3:ncol(results)],2,function(x)as.numeric(x)) #make sure each metabolite is numeric type

#define food name
var <- c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Whole Wheat Bread","Potato")

#run cluster analysis across each food 
for (i in var){

  #define dataset for use
  use <- results[results$group==i,c("time",metvar)] %>% t() %>% as.data.frame()
  names(use) <- use[1,]
  use <- use[-1,]
  names(use) <- c("t0","t2","t4","t6","t8","t24") #6 timepoints
  
  for (j in 1:10){
    use[,j] <- as.numeric(use[,j])
  }
  
  pd <- data.frame(time = c("0","2","4","6","8","24"))
  rownames(pd) <- colnames(use)
  
  #standarize dataset for clustering analysis
  eset <- new("ExpressionSet", exprs = as.matrix(use))
  keep <- rowSums(is.na(exprs(eset))) < ncol(exprs(eset)) - 1
  eset_filtered <- eset[keep, ]
  eset_std <- standardise(eset_filtered)
  
  keep2 <- rowSums(is.na(exprs(eset_std))) < ncol(exprs(eset_std)) - 1
  eset_filtered <- eset_std[keep2, ]
  
  #run model to get the optimized cluster number
  m <- mestimate(eset_std)                                      # Estimate fuzziness parameter
  c  <- which.min(Dmin(eset_filtered, m = m, crange = 2:10))    # Plot to help choose cluster number
  cl <- mfuzz(eset_filtered, c = c, m=m)                        # c = number of clusters
  
  #create membership object
  membership_df <- cl$membership %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "membership") %>%
    group_by(gene) %>%
    slice_max(membership, with_ties = FALSE) 
  
  #extract sample size info
  size <- as.data.frame(cl$size)
  size$cluster <- seq_len(nrow(size))
  size$cluster <- as.numeric(size$cluster)
  names(size)[1] <-  "Num"
  
  #extract cluster assignments info
  exprs_long <- exprs(eset_filtered) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression")
  
  exprs_long <- exprs_long %>%
    left_join(data.frame(gene = names(cl$cluster), cluster = cl$cluster), by = "gene") %>%
    left_join(pd %>% rownames_to_column("sample"), by = "sample") %>%
    left_join(membership_df[,-2], by = "gene") %>%
    left_join(size,by="cluster")
  
  #plot cluster curve and save it
  exprs_long$time <- factor(exprs_long$time, levels=c("0","1","2","3","4","5","6","7","8","24"))
  exprs_long$cluster <- paste("Cluster",exprs_long$cluster,sep='')
  
  exprs_long$cluster <- factor(exprs_long$cluster,levels = unique(exprs_long$cluster[order(as.numeric(gsub("\\D", "", exprs_long$cluster)))]))
  
  p1 <-ggplot(exprs_long, aes(x = time, y = expression, group = gene)) +
    geom_line(alpha=0.3,color="grey40",size=0.7) +
    geom_smooth(aes(group = cluster, color = as.factor(cluster)), method = "loess", se = FALSE, size = 1.2) +
    scale_color_manual(values=c("#1f77b4", "#ff7f0e", "#2ca02c", "#FA5555", "#9467bd","#8c564b", "#bcbd22", "#17becf","#a6cee3", "#fdbf6f","pink","yellow","blue","orange")) +
    facet_wrap(~ cluster, scales = "free_y",nrow=5) +
    labs(x = "Time", y = "Standardized Intensity")+
    theme_minimal() +
    theme(strip.text=element_text(size=12,color="black",family="calibri"),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          axis.text.y = element_text(size=12,color="black",family="calibri"),
          axis.text.x = element_text(size=12,color="black",family="calibri"),
          axis.title = element_text(size=12,color="black",family="calibri"),
          axis.ticks = element_line(color="black"),
          panel.spacing.x = unit(0.01, "cm"),
          panel.spacing.y = unit(0.01, "cm"),
          legend.position = "none",
          legend.title = element_text(color="black",family="calibri",size=12),
          legend.text = element_text(color="black",family="calibri",size=12))
  
  ggsave(paste(i,"Cluster_Urine.png",sep=''), plot = p1, width = 6, height =12, dpi = 600)
  
  #save results, including standarized dataset, sample size, gene membership
  sz <- cl$size %>% as.data.frame()
  gene_cluster <- cbind(cl$cluster, cl$membership)
  colnames(gene_cluster)[1] <- 'cluster'

  save(eset_filtered,sz,gene_cluster,file=paste("",i,"_Cluster_Urine.RData",sep=''))
  
  #plot cluster using heatmap
  mem <- as.data.frame(gene_cluster) %>% arrange(cluster)
  mem$name <- rownames(mem)
  mem$cluster  <- factor(mem$cluster,levels=sort(unique(mem$cluster)))
  
  rs <- eset_filtered@assayData[["exprs"]] %>% as.data.frame()
  rs <- rs[match(mem$name,rownames(rs)),]

  breaks <- c(-2,0,2)
  colors <- c("#006EBE", "white", "#FF6666")
  col_fun <- colorRamp2(breaks, colors)
  
  col_ha<-columnAnnotation(foo = names(rs),
                           col = list(foo = c("t0" = "#80CDC1", "t2"="#D15FEE","t4"="#FFC125","t6"="#440154", "t8"="#26828E","t24"="#B4DE2C")),gp = gpar(col = "black",family="Calibri",size=9),
                           annotation_name_side = "right",annotation_name_rot=0,
                           annotation_legend_param = list(foo = list(title = "Time",labels = c("0hr", "2hr","4hr","6hr","8hr","24hr"))),
                           annotation_label = " ")
  
  row_ha<-rowAnnotation(foo = mem$cluster,
                        col = list(foo = c("1" = "#1f77b4", "2" = "#ff7f0e", "3" = "#2ca02c", "4" = "#FA5555", "5" = "#9467bd","6" = "#8c564b", "7" = "#bcbd22", "8" = "#17becf","9" = "#a6cee3", "10" = "#fdbf6f")),gp = gpar(col = "black",family="Calibri",size=9),
                        annotation_name_side = "bottom",annotation_name_rot=0,
                        annotation_legend_param = list(foo = list(title = "Cluster")),
                        annotation_label = " ") #trait classification
  
  p2 <- Heatmap(rs, 
                col=col_fun,
                rect_gp = gpar(col = NA, lwd = 1),
                show_row_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                column_dend_reorder = TRUE,
                column_gap = unit(4, "mm"),
                row_gap = unit(2, "mm"),
                show_column_dend = TRUE,
                cluster_rows = F,
                cluster_columns = F,
                column_title = NULL,
                left_annotation = row_ha,
                top_annotation = col_ha,
                show_heatmap_legend = FALSE)
  
  ggsave(paste(i,"_Heatmap_Urine.png",sep=''), plot = wrap_elements(grid.grabExpr(draw(p2))), width = 6, height = 12, dpi = 300) 
  
}
