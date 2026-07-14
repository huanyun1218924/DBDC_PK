#load packages
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(grid)
library(gridExtra)
library(Mfuzz)
library(circlize)
library(pheatmap)

#read data
load("plasma_sample.RData")

#check data
dim(plasma_sample)
#1125 995

#check annotation file
anno_1 <- anno[which(anno$Class %in% c("Azoles","Azolidines","Diazines","Pyridines and derivatives","Isoquinolines and derivatives","Indoles and derivatives","Benzopyrans","Quinolines and derivatives","Benzodioxoles","Benzothiazoles","Cinchona alkaloids","Naphthalenes")),]
anno_2 <- anno[which(anno$Class %in% c("Purine nucleotides","Pyrimidine nucleosides","Purine nucleosides","5'-Deoxyribonucleosides","Nucleoside and nucleotide analogues","Imidazopyrimidines","Pyrazolopyrimidines","Pteridines and derivatives")),]
anno_3 <- anno[which(anno$Class %in% c("Flavonoids","Tropane alkaloids","Cinchona alkaloids","Phenylpropanoic acids","Biotin and derivatives")),]
anno_4 <- anno[which(anno$Class %in% c("Peptidomimetics","Hydroxy acids and derivatives","Keto acids and derivatives","Carboxylic acids and derivatives","Carboximidic acids and derivatives")),]
anno_5 <- anno[which(anno$Class %in% c("Organic sulfonic acids and derivatives","Organic sulfuric acids and derivatives","Sulfinic acids and derivatives","Organooxygen compounds","Organonitrogen compounds","Phenols","Phenol ethers")),]
anno_6 <- anno[which(anno$Class %in% c("Benzene and substituted derivatives","Naphthalenes","Tetrahydroisoquinolines")),]
anno_7 <- anno[which(anno$Class %in% c("Benzocycloheptapyridines","Lactones Oxepanes","Tetrapyrroles and derivatives","Pyrrolidines","Others")),]
anno_8 <- anno[which(anno$Class %in% c("Glycerophospholipids","Glycerolipids","Fatty Acyls","Sphingolipids","Prenol lipids","Steroids and steroid derivatives")),]

anno_1$Label <- "Heterocyclic Compounds"
anno_2$Label <- "Nucleotides & Analogues"
anno_3$Label <- "Plant Metabolites"
anno_4$Label <- "Amino Acid and Peptide-Related Compounds"
anno_5$Label <- "Sulfur- and Oxygen-Containing Organics"
anno_6$Label <- "Aromatic Hydrocarbons and Derivatives"
anno_7$Label <- "Others"
anno_8$Label <- anno_8$Class

anno <- rbind(anno_1,anno_2,anno_3,anno_4,anno_5,anno_6,anno_7,anno_8)

#get geometric mean
geo_mean <- function(x) {
  if (any(x <= 0)) {
    return(NA) # Geometric mean is undefined for non-positive values
  }
  exp(mean(log(x)))
}

results <- df %>% group_by(group, time) %>% summarise(across(metvar, geo_mean, .names = "geo_mean_{col}")) %>% ungroup()
names(results)[3:ncol(df] <- metvar

results[3:35641] <- apply(results[3:35641],2,function(x)as.numeric(x))

metvar <- pla.anno[which(pla.anno$`HMDB_specificity (1=match; 2=representative)` %in% c(1,2)),]$name

#define parameters used for heatmap
breaks <- c(-2,0,2)
colors <- c("#006EBE", "white", "#FF6666")
col_fun <- colorRamp2(breaks, colors)

#make dataset using metabolites as rows and timepoints as columns
var <- c("Beef","Chicken","Cheese","Yogurt","Corn","Oats","Whole Wheat Bread","Potato")

for (i in var){
  
  use <- results[results$group==i,c("time",metvar)] %>% t() %>% as.data.frame()
  names(use) <- use[1,]
  use <- use[-1,]
  names(use) <- c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t24")
  
  for (j in 1:10){
    use[,j] <- as.numeric(use[,j])
  }
  
  pd <- data.frame(time = c("0","1","2","3","4","5","6","7","8","24"))
  rownames(pd) <- colnames(use)
  
  #transform it into gene set and standarize data
  eset <- new("ExpressionSet", exprs = as.matrix(use))
  keep <- rowSums(is.na(exprs(eset))) < ncol(exprs(eset)) - 1
  eset_filtered <- eset[keep, ]
  eset_std <- standardise(eset_filtered)
  
  keep2 <- rowSums(is.na(exprs(eset_std))) < ncol(exprs(eset_std)) - 1
  eset_filtered <- eset_std[keep2, ]
  
  #run model
  m <- mestimate(eset_std)                    # Estimate fuzziness parameter
  c  <- which.min(Dmin(eset_filtered, m = m, crange = 2:10))    # Plot to help choose cluster number
  
  cl <- mfuzz(eset_filtered, c = c, m=m)      # c = number of clusters
  
  # Create membership object
  membership_df <- cl$membership %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "membership") %>%
    group_by(gene) %>%
    slice_max(membership, with_ties = FALSE) 
  
  # Extract sample size info
  size <- as.data.frame(cl$size)
  size$cluster <- seq_len(nrow(size))
  size$cluster <- as.numeric(size$cluster)
  names(size)[1] <-  "Num"
  
  # Extract cluster assignments and tidy data
  exprs_long <- exprs(eset_filtered) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression")
  
  exprs_long <- exprs_long %>%
    left_join(data.frame(gene = names(cl$cluster), cluster = cl$cluster), by = "gene") %>%
    left_join(pd %>% rownames_to_column("sample"), by = "sample") %>%
    left_join(membership_df[,-2], by = "gene") %>%
    left_join(size,by="cluster")
  
  #plot cluster
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
  
  #save results
  sz <- cl$size %>% as.data.frame()
  gene_cluster <- cbind(cl$cluster, cl$membership)
  colnames(gene_cluster)[1] <- 'cluster'
  
  #save results, including standarized dataset, sample size, gene membership
  save(eset_filtered,sz,gene_cluster,file=paste("",i,"_Cluster_Plasma.RData",sep=''))
  
  #plot cluster using heatmap
  mem <- as.data.frame(gene_cluster) %>% arrange(cluster)
  mem$name <- rownames(mem)
  mem$cluster  <- factor(mem$cluster,levels=c("1","2","3","4","5","6","7","8","9","10"))
  
  rs <- eset_filtered@assayData[["exprs"]] %>% as.data.frame()
  rs <- rs[match(mem$name,rownames(rs)),]
  
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
  
  #plot subclass in each cluster
  mem2 <- left_join(mem,pla.anno[,c("name","HMDB_ID","HMDB_specificity (1=match; 2=representative)"),by="name"])
  mem3 <- mem2[which(mem2$`HMDB_specificity (1=match; 2=representative)`%in% c(1,2)),]                   
  mem4 <- left_join(mem3,anno,by="HMDB_ID")
  
  mem4[which(is.na(mem4$Label)),"Label"] <- "Others"
  
  cl_num <- mem4 %>%
    group_by(Label,cluster) %>%
    summarise(count=n())
  
  cl_num$cluster <- paste("Cluster",cl_num$cluster,sep='')
  cl_num$cluster <- factor(cl_num$cluster,levels = unique(cl_num$cluster))
  
  p3 <-  ggplot(cl_num, aes(x = reorder(Label,count), y = count,fill=Label)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values=c("#1f77b4", "#ff7f0e", "#2ca02c", "#FA5555", "#9467bd","#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf","#a6cee3", "#fb9a99", "#fdbf6f")) +
    facet_wrap(~ cluster,ncol=2) +
    scale_x_discrete(drop = TRUE)+
    labs(title = "",x = "", y = "Count") +
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
          legend.text = element_text(color="black",family="calibri",size=12))+
    coord_flip()
  
  ggsave(paste(i,"Cluster_Bar_Plasma.png",sep=''), plot = p3, width = 9, height =12, dpi = 600)
  
}

load("Beef_Cluster_Plasma.RData")

cl1 <- gene_cluster[which(gene_cluster$cluster==1),]
cl1$name <- rownames(cl1)
cl1 <- left_join(cl1,pla.anno,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

cl2 <- gene_cluster[which(gene_cluster$cluster==2),]
cl2$name <- rownames(cl2)
cl2 <- left_join(cl2,pla.anno,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

cl3 <- gene_cluster[which(gene_cluster$cluster==3),]
cl3$name <- rownames(cl3)
cl3 <- left_join(cl3,pla.anno,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

cl8 <- gene_cluster[which(gene_cluster$cluster==8),]
cl8$name <- rownames(cl8)
cl8 <- left_join(cl8,pla.anno,by="name") %>% left_join(anno[,c("HMDB_ID","KEGGs")],by="HMDB_ID")

#supplementary figure 4
var <- c("Chicken","Cheese","Yogurt","Corn","Oats","Whole Wheat Bread","Potato")

plots <- list()

for (i in var){
  
  use <- results[results$group==i,c("time",metvar)] %>% t() %>% as.data.frame()
  names(use) <- use[1,]
  use <- use[-1,]
  names(use) <- c("t0","t1","t2","t3","t4","t5","t6","t7","t8","t24")
  
  for (j in 1:10){
    use[,j] <- as.numeric(use[,j])
  }
  
  pd <- data.frame(time = c("0","1","2","3","4","5","6","7","8","24"))
  rownames(pd) <- colnames(use)
  
  #transform it into gene set and standarize data
  eset <- new("ExpressionSet", exprs = as.matrix(use))
  keep <- rowSums(is.na(exprs(eset))) < ncol(exprs(eset)) - 1
  eset_filtered <- eset[keep, ]
  eset_std <- standardise(eset_filtered)
  
  keep2 <- rowSums(is.na(exprs(eset_std))) < ncol(exprs(eset_std)) - 1
  eset_filtered <- eset_std[keep2, ]
  
  #run model
  m <- mestimate(eset_std)                    # Estimate fuzziness parameter
  c  <- which.min(Dmin(eset_filtered, m = m, crange = 2:10))    # Plot to help choose cluster number
  
  cl <- mfuzz(eset_filtered, c = c, m=m)      # c = number of clusters
  
  # Create membership object
  membership_df <- cl$membership %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cluster", values_to = "membership") %>%
    group_by(gene) %>%
    slice_max(membership, with_ties = FALSE) 
  
  # Extract sample size info
  size <- as.data.frame(cl$size)
  size$cluster <- seq_len(nrow(size))
  size$cluster <- as.numeric(size$cluster)
  names(size)[1] <-  "Num"
  
  # Extract cluster assignments and tidy data
  exprs_long <- exprs(eset_filtered) %>%
    as.data.frame() %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "sample", values_to = "expression")
  
  exprs_long <- exprs_long %>%
    left_join(data.frame(gene = names(cl$cluster), cluster = cl$cluster), by = "gene") %>%
    left_join(pd %>% rownames_to_column("sample"), by = "sample") %>%
    left_join(membership_df[,-2], by = "gene") %>%
    left_join(size,by="cluster")
  
  #plot cluster
  exprs_long$time <- factor(exprs_long$time, levels=c("0","1","2","3","4","5","6","7","8","24"))
  exprs_long$cluster <- paste("Cluster",exprs_long$cluster,sep='')
  
  exprs_long$cluster <- factor(exprs_long$cluster,levels = unique(exprs_long$cluster[order(as.numeric(gsub("\\D", "", exprs_long$cluster)))]))
  
  p <- ggplot(exprs_long, aes(x = time, y = expression, group = gene)) +
    geom_line(alpha=0.3,color="grey70",size=0.2) +
    geom_smooth(aes(group = cluster, color = as.factor(cluster)), method = "loess", se = FALSE, size = 0.5) +
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
  
  plots[[i]] <- p + ggtitle(i)
}

p <- grid.arrange(grobs = plots, nrow = 2)
ggsave("Cluster_Plasma.png", plot = p, width = 12, height =15, dpi = 600)
