
/usr/local/R4.2/bin/R
library(dplyr)
library(trqwe)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(dplyr)
library(rjson)
library(parallel)
library(future)
library(future.apply)
library(data.table)
library(reticulate)
library(ReductionWrappers)
library(s2a)
library(SCP)
options(future.globals.maxSize = 300 * 1024^3)
plan("multiprocess", workers = 10)
plan()
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")
library(trqwe)


*************************内皮空间蛋白组分析********************************************************
*************************内皮空间蛋白组分析********************************************************
*************************内皮空间蛋白组分析********************************************************
*************************内皮空间蛋白组分析********************************************************

sgIGFBP5_VS_sgScr <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Protomics/Processed_Data/SCLC_sgIGFBP5_VS_sgScr_Protein_Matrix_annotaion.csv")
Normal_vs_NSCLC <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Protomics/Processed_Data/Normal_VS_NSCLC_Protein_Matrix_annotaion.csv")

sgIGFBP5_VS_sgScr <- sgIGFBP5_VS_sgScr[,c("UNIPROT","Gene","TRM.Scr.L5","TRM.Scr.1","TRM.I5.1","TRM.I5.4.L1")]
Normal_vs_NSCLC <- Normal_vs_NSCLC[,c("UNIPROT","Gene","Conlung1","Conlung2","Conlung3", "Nsclc.1","Nsclc.2", "Nsclc.3")]
names(sgIGFBP5_VS_sgScr) <- c("UNIPROT","Gene","SCLC_sgScr.1","SCLC_sgScr.2","SCLC_sgIgfbp5.1","SCLC_sgIgfbp5.2")
names(Normal_vs_NSCLC) <-  c("UNIPROT","Gene","Normal.1","Normal.2","Normal.3","NSCLC.1","NSCLC.2","NSCLC.3")

all_protein <- merge(Normal_vs_NSCLC,sgIGFBP5_VS_sgScr,by="UNIPROT",all.x=TRUE,all.y=TRUE)
write.csv(all_protein,file="Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")


all_protein <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Protomics/Processed_Data/Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")

score1 <- sum(all_protein$Normal.1)/10000000
Normalized_Normal.1 <- all_protein$Normal.1/score1
score2 <- sum(all_protein$Normal.2)/10000000
Normalized_Normal.2 <- all_protein$Normal.2/score2
score3 <- sum(all_protein$Normal.3)/10000000
Normalized_Normal.3 <- all_protein$Normal.3/score3

score4 <- sum(all_protein$NSCLC.1)/10000000
Normalized_NSCLC.1 <- all_protein$NSCLC.1/score4
score5 <- sum(all_protein$NSCLC.2)/10000000
Normalized_NSCLC.2 <- all_protein$NSCLC.2/score5
score6 <- sum(all_protein$NSCLC.3)/10000000
Normalized_NSCLC.3 <- all_protein$NSCLC.3/score6
score7 <- sum(all_protein$SCLC_sgScr.1)/10000000
Normalized_SCLC_sgScr.1 <- all_protein$SCLC_sgScr.1/score7
score8 <- sum(all_protein$SCLC_sgScr.2)/10000000
Normalized_SCLC_sgScr.2 <- all_protein$SCLC_sgScr.2/score8
score9 <- sum(all_protein$SCLC_sgIgfbp5.1)/10000000
Normalized_SCLC_sgIgfbp5.1 <- all_protein$SCLC_sgIgfbp5.1/score9
score10 <- sum(all_protein$SCLC_sgIgfbp5.2)/10000000
Normalized_SCLC_sgIgfbp5.2 <- all_protein$SCLC_sgIgfbp5.2/score10


tmp <- data.frame(Normalized_Normal.1,Normalized_Normal.2,Normalized_Normal.3,Normalized_NSCLC.1,Normalized_NSCLC.2,Normalized_NSCLC.3,Normalized_SCLC_sgScr.1,Normalized_SCLC_sgScr.2,Normalized_SCLC_sgIgfbp5.1,Normalized_SCLC_sgIgfbp5.2)
tmp$UNIPROT <- all_protein$UNIPROT
tmp$Gene <- all_protein$Gene
write.csv(tmp,file="Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")




suppressPackageStartupMessages({
	library(Rsamtools)
	library(GenomicFeatures)
	library(GenomicAlignments)
	library(BiocParallel)
	library(pheatmap)
	library(RColorBrewer)
	library(PoiClaClu)
	library(org.Mm.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(pathview)
	library(org.Hs.eg.db)
	library(AnnotationDbi)
	library(DOSE)
	library(clusterProfiler)
	library(topGO)
	library(ggplot2)
})


Protein <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Protomics/Processed_Data/Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")
rownames(Protein) <- Protein$UNIPROT
Protein <- Protein[,c(-1,-2,-3,-4)]
NSCLC <- Protein[,c(1:3)]
SCLC <- Protein[,c(4:5)]
names(SCLC) <- c("SCLC_1","SCLC_2")
names(NSCLC) <- c("NSCLC_1","NSCLC_2","NSCLC_3")
SCLC_VS_NSCLC_logFC <- log2((rowMeans(SCLC)+1) / (rowMeans(NSCLC)+1))
library("future.apply")
p_values <- future_lapply(seq(1,nrow(SCLC)), function(x){
 res <- t.test(x = as.numeric(SCLC[x,][1,]), y = as.numeric(NSCLC[x,][1,]))
 res$p.value
})
p <- unlist(p_values)
p.adj <- p.adjust(p, method = "fdr")
genelist <- as.data.frame(SCLC_VS_NSCLC_logFC)
genelist$SCLC_VS_NSCLC_pvalues <- p
genelist$SCLC_VS_NSCLC_padj <- p.adj
all_exp <- cbind(SCLC,NSCLC)
Protein_basemean <- apply(all_exp,1,mean)
all <- cbind(SCLC,NSCLC,Protein_basemean,genelist, Protein[,c("UNIPROT","Gene","ENSEMBL")] )
all$ENTREZID <- mapIds(x = org.Mm.eg.db,
                        keys = as.character(all$Gene),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
all <- data.frame(all)
all <- apply(all,2,as.character)
write.csv(all, file ="1_SCLC_VS_NSCLC_Protein.csv")





Protein <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Protomics/Processed_Data/Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")
rownames(Protein) <- Protein$UNIPROT
Protein <- Protein[,c(-1,-2,-3,-4)]
sgScr <- Protein[,c(4:5)]
sgIgfbp5 <- Protein[,c(6:7)]
names(sgIgfbp5) <- c("sgIgfbp5_1","sgIgfbp5_2")
names(sgScr) <- c("sgScr_1","sgScr_2")
sgIgfbp5_VS_sgScr_logFC <- log2((rowMeans(sgIgfbp5)+1) / (rowMeans(sgScr)+1))

library("future.apply")
p_values <- future_lapply(seq(1,nrow(sgIgfbp5)), function(x){
 res <- t.test(x = as.numeric(sgIgfbp5[x,][1,]), y = as.numeric(sgScr[x,][1,]))
 res$p.value
})
p <- unlist(p_values)
p.adj <- p.adjust(p, method = "fdr")
genelist <- as.data.frame(sgIgfbp5_VS_sgScr_logFC)
genelist$sgIgfbp5_VS_sgScr_pvalues <- p
genelist$sgIgfbp5_VS_sgScr_padj <- p.adj
all_exp <- cbind(sgIgfbp5,sgScr)
Protein_basemean <- apply(all_exp,1,mean)
all <- cbind(sgIgfbp5,sgScr,Protein_basemean,genelist, Protein[,c("UNIPROT","Gene","ENSEMBL")] )

all$ENTREZID <- mapIds(x = org.Mm.eg.db,
                        keys = as.character(all$Gene),
            keytype ="SYMBOL",
            column ="ENTREZID",
            multiVals="first")
all <- data.frame(all)
all <- apply(all,2,as.character)
write.csv(all, file ="1_sgIgfbp5_VS_sgScr_Protein.csv")




***************************火山图****************************************************************
***************************火山图****************************************************************
***************************火山图****************************************************************

library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")

SCLC_VS_NSCLC_Protein <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/decenber_following_cc/1_SCLC_VS_NSCLC_Protein.csv")
genes_plant <- c("Tjp1", "Tjp2", "Col4a1", "Col4a2", "Col18a1", "Lamc1", "Lama4", "Lamb1", "Cgn", "Epha2", "Ctnnb1", "Afdn")

library(ggplot2)
logFC <- SCLC_VS_NSCLC_Protein$SCLC_VS_NSCLC_logFC
basemean <- SCLC_VS_NSCLC_Protein$Protein_basemean
SCLC_VS_NSCLC_Protein$cc_score <- logFC*basemean/500

SCLC_VS_NSCLC_Protein$group <- ifelse(SCLC_VS_NSCLC_Protein$cc_score > 0,"SCLC_hi","SCLC_low")
SCLC_hi <- subset(SCLC_VS_NSCLC_Protein,group=="SCLC_hi")
SCLC_low <- subset(SCLC_VS_NSCLC_Protein,group=="SCLC_low")

SCLC_hi$final_score <- log(as.numeric(SCLC_hi$cc_score+1,2))
SCLC_low$final_score <- -1*log(as.numeric(SCLC_low$cc_score*-1+1,2))

tmp <- rbind(SCLC_hi,SCLC_low)
tmp <- tmp[!duplicated(tmp$Gene),]
range(tmp$final_score)

data <- data.frame(final_score=tmp$final_score, pvalue=tmp$SCLC_VS_NSCLC_pvalues)
data$sig[(data$pvalue > 0.1 | data$pvalue=="NA" | -1 < data$pvalue |  data$pvalue < -1 )] <- "no"
data$sig[data$pvalue <= 0.1 & data$final_score >= 0.5] <- "up"
data$sig[data$pvalue <= 0.1 & data$final_score <= -0.5] <- "down"
data$symbol <- tmp$Gene
data$log10_pvalue <- -log10(data$pvalue)

range(data$log10_pvalue)

data$log10_pvalue[data$log10_pvalue > 4] <- 4
data$final_score[data$final_score > 5] <- 5
data$final_score[data$final_score < -5] <- -5

sel_genes <- function(gene){
  tmp <- subset(data,symbol==gene)
  return(tmp)
}
aa <- future_lapply(as.list(as.character(genes_plant)),sel_genes)
all_res_plant <- do.call(rbind,aa)

pdf(file ="Protein_volcanomap.pdf", width = 6, height = 6)
ggplot(data, aes(x = final_score, y = log10_pvalue, color = sig)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_colour_manual(values = c('#FF3333', '#0066CC', 'grey'), limits = c('up', 'down', 'no')) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = 'black', fill = 'transparent'),
    plot.title = element_text(hjust = 0.5)
  ) +
  theme(
    legend.key = element_rect(fill = 'transparent'),
    legend.background = element_rect(fill = 'transparent'),
    legend.position = c(0.9, 0.93)
  ) +
  geom_vline(xintercept = c(-0.5, 0.5), color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.1, 10), color = 'gray', size = 0.3) +
  xlim(-5, 5) + ylim(0, 4) +
  labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)\n', color = '', title = '') +
  theme(legend.position = 'right') +
  
  # 添加基因标注
  geom_text_repel(
    data = all_res_plant,
    aes(x = final_score, y = log10_pvalue, label = symbol),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = unit(0.5, "lines"),  # 增大字符与边界的间距
    point.padding = unit(0.3, "lines"), # 控制字符与点的距离
    segment.color = "black",
    segment.size = 0.5,
    max.overlaps = Inf, # 允许文本重叠排布
    force = 5, # 强化文本与点分离的力
    min.segment.length = 0.1 # 避免过短的线段
  ) +
  
  # 添加指示线，避免重叠
  geom_segment(
    data = all_res_plant,
    aes(
      x = final_score, 
      y = log10_pvalue - 1.5, # 指线起点向下偏移
      xend = final_score, 
      yend = log10_pvalue
    ),
    color = "black", # 线的颜色
    size = 0.3,      # 线的粗细
    linetype = "solid" # 线的样式
  )
dev.off()




matrix <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/decenber_following_cc/1_SCLC_VS_NSCLC_Protein.csv")
SCLC_hi <- subset(matrix, SCLC_VS_NSCLC_logFC > 2)
dd <- as.vector(na.omit(SCLC_hi$ENTREZID))
GOupres_1_all <- enrichGO(gene = dd, 
             OrgDb = "org.Mm.eg.db",
             ont = "all", 
                 pvalueCutoff = 1, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 1,
                     minGSSize = 1, 
                     maxGSSize = 1000, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupres_1_all,file="1_SCLC_VS_NSCLC_Protein_SCLC_hi_GOres_log2.csv")






SCLC_hi <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/decenber_following_cc/1_SCLC_VS_NSCLC_Protein_SCLC_hi_GOres_log2.csv")
sel_path <- subset(SCLC_hi, Description=="extracellular matrix" | Description=="establishment of endothelial intestinal barrier" | Description=="adherens junction organization" | 
	Description=="establishment of endothelial barrier" | 
	Description=="regulation of basement membrane organization" | 
	Description=="extracellular matrix organization" | Description=="tight junction assembly" | 
	Description=="tight junction organization" | Description=="basement membrane organization" | 
	Description=="maintenance of blood-brain barrier" | Description=="tight junction" | 
	Description=="basement membrane" | Description=="gap junction")

sel_path <- data.frame(sel_path)
sel_path <- sel_path[order(sel_path$pvalue,decreasing=FALSE),]
sel_path$Description <- factor(sel_path$Description, levels=unique((as.character(sel_path$Description))))

sel_path$pvalue_score <- -log10(sel_path$pvalue)
library(ggplot2)
library(ggpubr)
p1 <- ggbarplot(sel_path, 
  x = "Description", 
  y = "pvalue_score",
  color = "#5B8FCF",            # Set bar border colors to white
  fill ="#5B8FCF",
  sort.val = "asc",          # Sort the value in dscending order
  x.text.angle = 90,           # Rotate vertically x axis texts
  rotate = TRUE,
  title="")+ylim(0,10)
ggsave(p1,file="SCLC_hi_GO_res_log2.svg",width =8, height = 4,dpi=1080)
