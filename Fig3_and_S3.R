$$$$$$$$$$$$$$$$$$Fig3_S3$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$Fig3_S3$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$Fig3_S3$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$Fig3_S3$$$$$$$$$$$$$$$$$$$$

library(dplyr)
library(data.table)
library(tidyverse)

library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(SCP)
library(scde)

library(future)
library(future.apply)
library(parallel)

library(ggplot2)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(patchwork)
library(scales)
library(RColorBrewer)
library(BuenColors)
library(pheatmap)

library(clusterProfiler)
library(DOSE)
library(topGO)
library(pathview)
library(GSEABase)
library(Biobase)
library(genefilter)
library(limma)
library(GSVA)

library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)

library(Rsamtools)
library(GenomicFeatures)
library(GenomicAlignments)
library(BiocParallel)

library(iTALK)
library(nichenetr)

library(reticulate)
library(rjson)

library(trqwe)
library(s2a)
library(tidyr)

library(future)
library(future.apply)

library(ggplot2)
library(ggrepel)


source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")


#Fig 3A

SCLC_VS_NSCLC <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_RNAseq/DEseq2normalized_SCLC_VS_NSCLC_Endo.csv")
SCLC_VS_NSCLC_hi <- subset(SCLC_VS_NSCLC,pvalue < 0.05 & log2FoldChange > 0 )


#Genes highly expressed in tumor cells, immune cells, Fibroblast cells were filtered out using single-cell RNA-sequencing data from SCLC patients (Chan et al.)
Rudin_All_Cells <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/CancerCell_FindAllMarkers.csv")
Rudin_All_Cells <- Rudin_All_Cells %>% mutate(mouse_gene = convert_human_to_mouse_symbols(as.character(Rudin_All_Cells$gene))) %>% drop_na()
Myeloid_Markers <- subset(Rudin_All_Cells,cluster=="Macrophage" | cluster=="Neutrophil" | cluster=="DC" | cluster=="Mast" )
Myeloid_Markers <- subset(Myeloid_Markers,p_val_adj < 0.05 &  avg_log2FC > 0.1)
Fibroblast_Markers <- subset(Rudin_All_Cells,cluster=="Fibroblast")
Fibroblast_Markers <- subset(Fibroblast_Markers,p_val_adj < 0.05 &  avg_log2FC > 0.1)
Tcell_Markers <- subset(Rudin_All_Cells,cluster=="T cell")
Tcell_Markers <- subset(Tcell_Markers,p_val_adj < 0.05 &  avg_log2FC > 0.1)
LUAD_Markers <- subset(Rudin_All_Cells,cluster=="NSCLC")
LUAD_Markers <- subset(LUAD_Markers,p_val_adj < 0.05 &  avg_log2FC > 0.1)
Epi_Markers <- subset(Rudin_All_Cells,cluster=="AE1" | cluster=="AEP" | cluster=="Ciliated" | cluster=="Basal" | cluster=="Club"| cluster=="Mucinous" | cluster=="Neuroendocrine" )
Epi_Markers <- subset(Epi_Markers,p_val_adj < 0.05 &  avg_log2FC > 0.1)
SCLC_Markers <- subset(Rudin_All_Cells,cluster=="SCLC-A" | cluster=="SCLC-N" | cluster=="SCLC-P")
SCLC_Markers <- subset(SCLC_Markers,p_val_adj < 0.05 &  avg_log2FC > 0.1)
Bcell_Markers <- subset(Rudin_All_Cells,cluster=="B cell" | cluster=="Plasma cell")
Bcell_Markers <- subset(Bcell_Markers,p_val_adj < 0.05 &  avg_log2FC > 0.1)


Pure_SCLC_Endo_hi <- setdiff(SCLC_VS_NSCLC_hi$X,SCLC_Markers$mouse_gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"
Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,LUAD_Markers$mouse_gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"
Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,Epi_Markers$mouse_gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"
Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,Fibroblast_Markers$mouse_gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"
Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,Myeloid_Markers$mouse_gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"

Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,Tcell_Markers$mouse_gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"
Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,Bcell_Markers$mouse_gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"
tmp <- subset(SCLC_VS_NSCLC, SCLC_VS_NSCLC$X %in% Pure_SCLC_Endo_hi$Symbol)
tmp$cc_score <- tmp$log2FoldChange* tmp$baseMean * tmp$baseMean
tmp <- tmp[order(-tmp$cc_score),]
write.csv(tmp,file="SCLC_VS_NSCLC_Endothelium_highly_expressed_genes_clean.csv")

GOupall_summry_all <- enrichGO(gene = as.character(na.omit(tmp$ENTREZID)), 
	           OrgDb = "org.Mm.eg.db",
				ont = "all", 
		             pvalueCutoff = 1, 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 1,
                     minGSSize = 1, 
                     maxGSSize = 1000, 
                     readable = T, 
                     pool = FALSE)
write.csv(GOupall_summry_all,"GOres_SCLC_VS_NSCLC_Endothelium_highly_expressed_genes.csv")

bulkRNA_SCLC_hi <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_RNAseq/GOres_SCLC_VS_NSCLC_Endothelium_highly_expressed_genes.csv")
sel_path <- subset(bulkRNA_SCLC_hi, Description=="adherens junction organization" | 
	Description=="tight junction assembly" | 
	Description=="tight junction organization" | 
	Description=="maintenance of blood-brain barrier" | 
	Description=="tight junction" | 
	Description=="basement membrane" | 
	Description=="gap junction")
sel_path <- data.frame(sel_path)
sel_path <- sel_path[order(sel_path$pvalue,decreasing=FALSE),]
sel_path$Description <- factor(sel_path$Description, levels=unique((as.character(sel_path$Description))))
sel_path$pvalue_score <- -log10(sel_path$pvalue)
library(ggplot2)
library(ggpubr)
p1 <- ggbarplot(sel_path, 
  x = "Description", 
  y = "pvalue_score",
  color = "#5B8FCF",            
  fill ="#5B8FCF",
  sort.val = "asc",         
  x.text.angle = 90,           
  rotate = TRUE,
  title="")+ylim(0,8)
ggsave(p1,file="Fig3A.png",width =8, height = 4,dpi=1080)





#Fig S3B

SCLC_VS_NSCLC_RNA <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_RNAseq/DEseq2normalized_SCLC_VS_NSCLC_Endo.csv")
genes_plant <- c("Cldn15", "Col4a5", "Col4a6", "Col8a2", "Des", "Gjb3", "Gjb4", "Lamc3")
library(ggplot2)
logFC <- SCLC_VS_NSCLC_RNA$log2FoldChange
basemean <- SCLC_VS_NSCLC_RNA$baseMean
SCLC_VS_NSCLC_RNA$cc_score <- logFC*basemean*basemean/10000

SCLC_VS_NSCLC_RNA$group <- ifelse(SCLC_VS_NSCLC_RNA$cc_score > 0,"SCLC_hi","SCLC_low")
SCLC_hi <- subset(SCLC_VS_NSCLC_RNA,group=="SCLC_hi")
SCLC_low <- subset(SCLC_VS_NSCLC_RNA,group=="SCLC_low")

SCLC_hi$final_score <- log(as.numeric(SCLC_hi$cc_score+1,10))
SCLC_low$final_score <- -1*log(as.numeric(SCLC_low$cc_score*-1+1,10))

tmp <- rbind(SCLC_hi,SCLC_low)

range(tmp$final_score)

data <- data.frame(final_score=tmp$final_score, pvalue=tmp$pvalue)
data$sig[(data$pvalue > 0.05 | data$pvalue=="NA" | -1 < data$pvalue |  data$pvalue < -1 )] <- "no"
data$sig[data$pvalue <= 0.05 & data$final_score >= 0.5] <- "up"
data$sig[data$pvalue <= 0.05 & data$final_score <= -0.5] <- "down"
data$symbol <- tmp$X
data$log10_pvalue <- -log10(data$pvalue)

range(data$log10_pvalue)
data$log10_pvalue[data$log10_pvalue > 15] <- 15

data$final_score[data$final_score > 15] <- 15
data$final_score[data$final_score < -15] <- -15

sel_genes <- function(gene){
  tmp <- subset(data,symbol==gene)
  return(tmp)
}
aa <- future_lapply(as.list(as.character(genes_plant)),sel_genes)
all_res_plant <- do.call(rbind,aa)


library(ggplot2)
library(ggrepel)

pdf(file ="FigS3B.pdf", width = 6, height = 6)
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
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
  xlim(-15, 15) + ylim(0, 15) +
  labs(x = '\nLog2 Fold Change', y = '-log10(pvalue)\n', color = '', title = '') +
  theme(legend.position = 'right') +
    geom_text_repel(
    data = all_res_plant,
    aes(x = final_score, y = log10_pvalue, label = symbol),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = unit(0.5, "lines"),  
    point.padding = unit(0.3, "lines"), 
    segment.color = "black",
    segment.size = 0.5,
    max.overlaps = Inf, 
    force = 5, 
    min.segment.length = 0.1 
  )
  
dev.off()


#Fig 3B
Cell_endo_marker <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/mouse_endothelial_cell_atlas/Cell_organ_specific_exp_genes.csv")
brain <- subset(Cell_endo_marker, pct.2 < 0.2 & p_val_adj < 0.01 & avg_logFC > 0.5 & cluster=="brain")
brain <- brain[,c("gene","cluster")]
brain <- brain %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(brain$gene))) %>% drop_na()
write.csv(brain,file="Endothelial_barrier_genelist_for_GSEA.csv")


$$$$$$$$$$$$$$$$$$$Proteomics$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$Proteomics$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$Proteomics$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$Proteomics$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


sgIGFBP5_VS_sgScr <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/sgIGFBP5_VS_sgScr_Protein_Matrix_annotaion.csv")
Normal_vs_NSCLC <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/Normal_VS_NSCLC_Protein_Matrix_annotaion.csv")

sgIGFBP5_VS_sgScr <- sgIGFBP5_VS_sgScr[,c("UNIPROT","Gene","TRM.Scr.L5","TRM.Scr.1","TRM.I5.1","TRM.I5.4.L1")]
Normal_vs_NSCLC <- Normal_vs_NSCLC[,c("UNIPROT","Gene","Conlung1","Conlung2","Conlung3", "Nsclc.1","Nsclc.2", "Nsclc.3")]
names(sgIGFBP5_VS_sgScr) <- c("UNIPROT","Gene","SCLC_sgScr.1","SCLC_sgScr.2","SCLC_sgIgfbp5.1","SCLC_sgIgfbp5.2")
names(Normal_vs_NSCLC) <-  c("UNIPROT","Gene","Normal.1","Normal.2","Normal.3","NSCLC.1","NSCLC.2","NSCLC.3")

all_protein <- merge(Normal_vs_NSCLC,sgIGFBP5_VS_sgScr,by="UNIPROT",all.x=TRUE,all.y=TRUE)
write.csv(all_protein,file="Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")


all_protein <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")

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




Protein <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")
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




Protein <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")
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



#Fig 3C
library(future)
library(future.apply)
options(future.globals.maxSize = 3000 * 1024^2)
plan("multiprocess", workers = 8)
plan()
library("ggplot2")
library("ggrepel")

SCLC_VS_NSCLC_Protein <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/1_SCLC_VS_NSCLC_Protein.csv")
genes_plant <- c("Tjp2", "Col4a2", "Col18a1", "Lamc1", "Lamb1")
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

pdf(file ="Fig3C.pdf", width = 6, height = 6)
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
    geom_text_repel(
    data = all_res_plant,
    aes(x = final_score, y = log10_pvalue, label = symbol),
    size = 5,
    fontface = "bold",
    color = "black",
    box.padding = unit(0.5, "lines"),  
    point.padding = unit(0.3, "lines"), 
    segment.color = "black",
    segment.size = 0.5,
    max.overlaps = Inf, 
    force = 5, 
    min.segment.length = 0.1 
  ) 
dev.off()



#Fig 3D
matrix <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/1_SCLC_VS_NSCLC_Protein.csv")
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
write.csv(GOupres_1_all,file="/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/1_SCLC_VS_NSCLC_Protein_GOres_UP.csv")

SCLC_hi <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/1_SCLC_VS_NSCLC_Protein_GOres_UP.csv")
sel_path <- subset(SCLC_hi, 
	Description=="tight junction" | 
	Description=="basement membrane" |
	Description=="adherens junction organization" | 
	Description=="gap junction" |
 	Description=="establishment of endothelial barrier" | 
 	Description=="tight junction assembly" |
	Description=="tight junction organization" | 
	Description=="maintenance of blood-brain barrier")

sel_path <- data.frame(sel_path)
sel_path <- sel_path[order(sel_path$pvalue,decreasing=FALSE),]
sel_path$Description <- factor(sel_path$Description, levels=unique((as.character(sel_path$Description))))

sel_path$pvalue_score <- -log10(sel_path$pvalue)
library(ggplot2)
library(ggpubr)
p1 <- ggbarplot(sel_path, 
  x = "Description", 
  y = "pvalue_score",
  color = "#5B8FCF",           
  fill ="#5B8FCF",
  sort.val = "asc",          
  x.text.angle = 90,         
  rotate = TRUE,
  title="")+ylim(0,5)
ggsave(p1,file="Fig3D.png",width =8, height = 4,dpi=1080)





$$$$$$$$$$$$$$$$$$$$$Define_IEV_signature$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$Define_IEV_signature$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$Define_IEV_signature$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$Define_IEV_signature$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#R3.6
#9022

#FigS3E
library(trqwe)
library(Seurat)
library(ggplot2)
library(pheatmap)
require(RColorBrewer)
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")


UMAP15 <- mcreadRDS("/mnt/data/user_data/abao/1_project/00000000_SCLC_Endothelium/old_submit/All_CancerType_Endo_UMAP15.rds",mc.cores=20)

Pancancer_tmp <- subset(UMAP15,Cancer_Type=="Bladder" | Cancer_Type=="BRCA" | Cancer_Type=="CRC" | Cancer_Type=="ESCC" | Cancer_Type=="OstSRC" | Cancer_Type=="OV" | Cancer_Type=="Prostate" | Cancer_Type=="HCC" | Cancer_Type=="LUSC" | Cancer_Type=="Melanoma" |Cancer_Type=="RCC" | Cancer_Type=="STAD")  
Normal_tmp <- subset(UMAP15, Cancer_Type=="Normal_Breast" | Cancer_Type=="Normal_Lung")        
SCLC_tmp <- subset(UMAP15,Cancer_Type=="SCLC")
LUAD_tmp <- subset(UMAP15,Cancer_Type=="LUAD")

pseudo_bulk_seurat_mean <- function(seurat_obj=seurat_obj,num_split=num_split,seed.use=seed.use,prefix=prefix,slot=slot){
  set.seed(seed.use)
  require(Seurat)
  genes.use <- rownames(seurat_obj)
  cell.sets1 <- split(colnames(seurat_obj), sort(1:length(colnames(seurat_obj))%%num_split))
  profile.set1 = matrix(, nrow = length(genes.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")[genes.use, this.set]
    if (length(this.set) > 1) {
      this.profile <- as.numeric(apply(sub.matrix, 1, function(x) mean(x)))
      profile.set1[, i] <- this.profile
    } else {
      profile.set1[, i] <- sub.matrix
    }
  }
  rownames(profile.set1) <- genes.use
  colnames(profile.set1) <- paste(prefix, 1:length(cell.sets1),sep="_")
  return(profile.set1)
}

SCLC_sudobulk1 <- pseudo_bulk_seurat_mean(seurat_obj=SCLC_tmp,num_split=10,seed.use=1,slot="counts",prefix=paste("SCLC"))
LUAD_sudobulk1 <- pseudo_bulk_seurat_mean(seurat_obj=LUAD_tmp,num_split=10,seed.use=1,slot="counts",prefix=paste("LUAD"))
Normal_sudobulk1 <- pseudo_bulk_seurat_mean(seurat_obj=Normal_tmp,num_split=10,seed.use=1,slot="counts",prefix=paste("Normal"))
PanCancer_sudobulk1 <- pseudo_bulk_seurat_mean(seurat_obj=Pancancer_tmp,num_split=10,seed.use=1,slot="counts",prefix=paste("PanCancer"))

All_sudoEndo <- cbind(Normal_sudobulk1,PanCancer_sudobulk1,LUAD_sudobulk1,SCLC_sudobulk1)
All_sudoEndo <- data.frame(All_sudoEndo)

normalized_matrix <- sweep(All_sudoEndo, 2, colSums(All_sudoEndo), FUN = "/") * 1e6
mcsaveRDS(normalized_matrix,file="sudobulk10_Endo_normalized_matrix.rds" )

expr_mat <- normalized_matrix
group <- c(rep(0, 30), rep(1, 10))
wilcox_pvals <- numeric(nrow(expr_mat))
logFC <- numeric(nrow(expr_mat))
wilcox_pvals <- numeric(nrow(expr_mat))
logFC <- numeric(nrow(expr_mat))

for (i in seq_len(nrow(expr_mat))) {
  gene_expr <- as.numeric(expr_mat[i, ])
  group0_expr <- gene_expr[group == 0]
  group1_expr <- gene_expr[group == 1]

  wt <- tryCatch(
    wilcox.test(group1_expr, group0_expr, exact = FALSE),
    error = function(e) NULL
  )

  wilcox_pvals[i] <- if (!is.null(wt)) wt$p.value else NA
  logFC[i] <- log2(
    mean(group1_expr + 1e-6) /
    mean(group0_expr + 1e-6)
  )
}
wilcox_results <- data.frame(
  Gene = rownames(expr_mat),
  log2FC = logFC,
  P_value = wilcox_pvals,
  FDR = p.adjust(wilcox_pvals, method = "BH")
)
head(wilcox_results[order(wilcox_results$P_value), ], 10)
tmp <- subset(wilcox_results, P_value < 0.01 & log2FC > 0)
write.csv(tmp,file="SCLC_endo_specific_gene_signature.csv")




#Using single-cell RNA-sequencing data from Rudin’s patient cohort, 
#genes highly expressed in SCLC and LUAD tumor cells and immune cells were filtered out.

SCLC_hi <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/SCLC_endo_specific_gene_signature.csv")
Rudin_All_Cells <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/CancerCell_FindAllMarkers.csv")

SCLC_Markers <- subset(Rudin_All_Cells,cluster=="SCLC-A" | cluster=="SCLC-N" | cluster=="SCLC-P")
SCLC_Markers <- subset(SCLC_Markers,p_val_adj < 0.05 &  avg_log2FC > 0)

LUAD_Markers <- subset(Rudin_All_Cells,cluster=="NSCLC")
LUAD_Markers <- subset(LUAD_Markers,p_val_adj < 0.05 &  avg_log2FC > 0)

Pure_SCLC_Endo_hi <- setdiff(SCLC_hi$Gene,SCLC_Markers$gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"
Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,LUAD_Markers$gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"

Tcell_Markers <- subset(Rudin_All_Cells,cluster=="T cell")
Tcell_Markers <- subset(Tcell_Markers,p_val_adj < 0.05 &  avg_log2FC > 0)

Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,Tcell_Markers$gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"


Macrophage_Markers <- subset(Rudin_All_Cells,cluster=="Macrophage")
Macrophage_Markers <- subset(Macrophage_Markers,p_val_adj < 0.05 &  avg_log2FC > 0)

Pure_SCLC_Endo_hi <- setdiff(Pure_SCLC_Endo_hi$Symbol,Macrophage_Markers$gene)
Pure_SCLC_Endo_hi <- data.frame(Pure_SCLC_Endo_hi)
names(Pure_SCLC_Endo_hi) <- "Symbol"

Pure_SCLC_Endo_hi <- Pure_SCLC_Endo_hi %>% mutate(mouse_gene = convert_human_to_mouse_symbols(as.character(Pure_SCLC_Endo_hi$Symbol))) %>% drop_na()
write.csv(Pure_SCLC_Endo_hi,file="IEV_signature.csv")


Pure_SCLC_Endo_hi <- read.csv("./IEV_signature.csv")
normalized_matrix <- mcreadRDS("./sudobulk10_Endo_normalized_matrix.rds",mc.cores=20)
normalized_matrix$Symbol <- rownames(normalized_matrix)

tmp <- subset(normalized_matrix, normalized_matrix$Symbol %in% Pure_SCLC_Endo_hi$Symbol)
tmp <- tmp[,-41]

sel_genes <- c("NID1", "LAMA4", "LAMC1", "LAMB1", "COL18A1","TJP2")

chonglai_zscore_1 <- t(apply(tmp, 1, function(x) (x-mean(x))/sd(x)))
range(chonglai_zscore_1)
chonglai_zscore_1[chonglai_zscore_1 > 1] <- 1
chonglai_zscore_1[chonglai_zscore_1 < -1] <- -1


SeuratObject <- CreateSeuratObject(counts = chonglai_zscore_1, project = "TM")
gene <- rownames(SeuratObject)
SeuratObject@meta.data$group <- rownames(SeuratObject@meta.data)

pdf("FigS3E.pdf",width=8,height=9)
XY_heatmap(seurat_obj=SeuratObject,group="orig.ident",gene = gene,all_num=FALSE,new_names=NULL,labels_rot=90,
  assay_sel="RNA",color=c( "#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac"),
  min_and_max_cut=2,show_row_names=FALSE,mark_gene=sel_genes,label_size=0,scale = FALSE)
dev.off()



#Fig S3H
Mice_DEGs <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/mouse_LUAD_VS_SCLC_FindAllmarker_DEGs_TumorCells.csv")
human_DEGs <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/Tumor/human_SCLC_VS_LUAD_FindAllmarkers.csv")
Mice_DEGs <- Mice_DEGs %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(Mice_DEGs$gene))) %>% drop_na()

Mice_SCLC_hi <- subset(Mice_DEGs,cluster=="SCLC" & p_val_adj < 0.01 & avg_log2FC > 0 )
human_SCLC_hi <- subset(human_DEGs,cluster=="SCLC" & p_val_adj < 0.01 & avg_logFC > 0 )
Mice_SCLC_hi <- Mice_SCLC_hi[,c("human_gene","avg_log2FC")]
human_SCLC_hi <- human_SCLC_hi[,c("gene","avg_logFC")]
names(Mice_SCLC_hi) <- c("Symbol","avg_log2FC")
names(human_SCLC_hi) <- c("Symbol","avg_log2FC")

Mice_LUAD_hi <- subset(Mice_DEGs,cluster=="LUAD" & p_val_adj < 0.01 & avg_log2FC > 0 )
human_LUAD_hi <- subset(human_DEGs,cluster=="LUAD" & p_val_adj < 0.01 & avg_logFC > 0 )

Mice_LUAD_hi <- Mice_LUAD_hi[,c("human_gene","avg_log2FC")]
human_LUAD_hi <- human_LUAD_hi[,c("gene","avg_logFC")]
names(Mice_LUAD_hi) <- c("Symbol","avg_log2FC")
names(human_LUAD_hi) <- c("Symbol","avg_log2FC")
Mice_LUAD_hi$avg_log2FC <- Mice_LUAD_hi$avg_log2FC*-1
human_LUAD_hi$avg_log2FC <- human_LUAD_hi$avg_log2FC*-1
mice_all_genes <- rbind(Mice_SCLC_hi,Mice_LUAD_hi)
human_all_genes <- rbind(human_SCLC_hi,human_LUAD_hi)
names(mice_all_genes) <- c("Symbol","mouse_avg_log2FC")
names(human_all_genes) <- c("Symbol","human_avg_log2FC")
all_tmp <- merge(mice_all_genes,human_all_genes,by="Symbol")
all_tmp <- all_tmp[!duplicated(all_tmp$Symbol),]
rownames(all_tmp) <- all_tmp$Symbol

write.csv(all_tmp,file="mouse_and_human_SCLC_VS_NSCLC_genelist.csv")
library(ggplot2)
library(ggpubr)

ff <- ggplot(all_tmp, aes(x = mouse_avg_log2FC, y = human_avg_log2FC)) +
  geom_point(alpha = 0.5, color = "lightgray") +  
  geom_point(data = all_tmp[all_tmp$Symbol %in% c("ASCL1"), ], aes(x = mouse_avg_log2FC, y = human_avg_log2FC), 
             color = "red", size = 3) +  
  geom_text(data = all_tmp[all_tmp$Symbol %in% c("ASCL1"), ], 
            aes(label = Symbol), vjust = -0.5, color = "red") +  
  stat_smooth(method = "lm", se = FALSE, color = "blue") +  
  stat_cor(data = all_tmp, method = "spearman") +  
  theme_minimal() +  
  labs(title = "Spearman Correlation Scatter Plot", 
       x = "Mouse Avg Log2FC", 
       y = "Human Avg Log2FC")
ggsave(ff,file="FigS3H.png",width=8,height=8)



$$$$$$$$$$$$$$$$$$$$$Impower133_survival$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$Impower133_survival$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$Impower133_survival$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$Impower133_survival$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#Fig3H
library(survival)
library(survminer)
library(dplyr)

setwd("/mnt/data/user_data/ailing/SCLC_IEV_code/IMpower133")
deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IMpower133/IMPOWER133_SCLC_deconv_fraction.csv")
IEV_sig <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/IEV_signature.csv")

IMPOWER133_Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IMpower133/IMPOWER133_Clinical_file.csv")
IMPOWER133_Clinical <- subset(IMPOWER133_Clinical,ACTARM.2=="atezo")
IMPOWER133_Clinical <- IMPOWER133_Clinical[,c("trunc_anonymized_sample_ids","Gay_MDACC_subtypes","PFS_MONTHS", "PFS_CENSOR","BOR")]
IMPOWER133_Clinical$X <- gsub("\\-","\\.",IMPOWER133_Clinical$trunc_anonymized_sample_ids)

ICI_Normalized_data <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IMpower133/IMPOWER133_SCLC_TPM_file.csv")

rownames(ICI_Normalized_data) <- ICI_Normalized_data$UNNAMED..0
ICI_Normalized_data <- ICI_Normalized_data[,c(-1,-2)]
ICI_Normalized_data <- ICI_Normalized_data[,IMPOWER133_Clinical$X]


IEV_marker <- intersect(rownames(ICI_Normalized_data),IEV_sig$Symbol)

IEV_gene <- data.frame(ICI_Normalized_data[IEV_marker,1:132])
IEV_mean <- data.frame(apply(IEV_gene,2,mean))
names(IEV_mean) <- "IEV_sig"
IEV_mean$X <- rownames(IEV_mean)   
tmp_1 <- merge(IMPOWER133_Clinical,IEV_mean,by="X")

#Given the substantial CAF infiltration in NMF4 shown in Figure 3B of Nabet et al., 2024, 
#the fibroblast-derived contribution to the IEV signature was removed.

tmp <- merge(tmp_1,deconv_fraction,by="X")
tmp$IEV_rmFibro_score <- tmp$IEV_sig*(1-tmp$Fibroblast)
tmp$dead_date <- as.numeric(tmp$PFS_MONTHS)*30
tmp$status <- tmp$PFS_CENSOR
genes <- data.frame(t(ICI_Normalized_data[c("IGFBP5","ASCL1"),]))
genes$ID <- rownames(genes)   
final_matrix <- merge(tmp,genes,by.x="X", by.y="ID")
final_matrix$IGFBP5_rmFibro_score <- final_matrix$IGFBP5*(1-final_matrix$Fibroblast)

mtx <- subset(final_matrix, Gay_MDACC_subtypes=="SCLC-A" | Gay_MDACC_subtypes=="SCLC-I")
my_comparisons <- list(c("SCLC-A","SCLC-I"))
ff <- ggboxplot(
  mtx, x = "Gay_MDACC_subtypes", y = "IEV_rmFibro_score",
  color = "Gay_MDACC_subtypes",
  add = "jitter") +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    legend.position = "none"
  )
ggsave(file="Fig3H.png", ff, width = 4, height = 5)




library("survival")
library("survminer")

#To quantify the endothelial contribution to the IEV signature, 
#the total IEV signature score of each sample was weighted by the proportion of endothelial cell infiltration in each patient.

final_matrix$IEV_with_Endo <- final_matrix$IEV_sig * final_matrix$Endothelial
case_info_ALL <- final_matrix
case_info_ALL <- case_info_ALL[order(case_info_ALL$dead_date,decreasing=TRUE),]
event <- gsub("alive",0,case_info_ALL$status)
event <- gsub("dead",1,event)
case_info_ALL$event <- as.numeric(event)
all_case_and_info <- na.omit(case_info_ALL)
case_info_ALL <- all_case_and_info
case_info_ALL.cut <- surv_cutpoint(
   case_info_ALL,
   time = "dead_date",
   event = "event",
   variables = c("IEV_with_Endo"),
   progressbar=TRUE,
   minprop=0.1
)
summary(case_info_ALL.cut)
IGF_cut <- surv_categorize(case_info_ALL.cut) 
library(survival)
fit <- survfit(Surv(dead_date, event) ~ IEV_with_Endo, data = IGF_cut)
pdf("Fig3F_1.pdf",height=8,width=8)
ggsurvplot(fit, data = IGF_cut,
surv.median.line = "hv",
pval = TRUE,
ggtheme = theme_bw(),
risk.table=TRUE)
dev.off()




#Placebo group

deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IMpower133/IMPOWER133_SCLC_deconv_fraction.csv")
IEV_sig <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/IEV_signature.csv")

IMPOWER133_Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IMpower133/IMPOWER133_Clinical_file.csv")
IMPOWER133_Clinical <- subset(IMPOWER133_Clinical,ACTARM.2=="placebo")
IMPOWER133_Clinical <- IMPOWER133_Clinical[,c("trunc_anonymized_sample_ids","Gay_MDACC_subtypes","PFS_MONTHS", "PFS_CENSOR","OS_MONTHS", "OS_CENSOR","BOR")]
IMPOWER133_Clinical$X <- gsub("\\-","\\.",IMPOWER133_Clinical$trunc_anonymized_sample_ids)

ICI_Normalized_data <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IMpower133/IMPOWER133_SCLC_TPM_file.csv")

rownames(ICI_Normalized_data) <- ICI_Normalized_data$UNNAMED..0
ICI_Normalized_data <- ICI_Normalized_data[,c(-1,-2)]
ICI_Normalized_data <- ICI_Normalized_data[,IMPOWER133_Clinical$X]

IEV_marker <- intersect(rownames(ICI_Normalized_data),IEV_sig$Symbol)
IEV_gene <- data.frame(ICI_Normalized_data[IEV_marker,1:139])
IEV_mean <- data.frame(apply(IEV_gene,2,mean))
names(IEV_mean) <- "IEV_sig"
IEV_mean$X <- rownames(IEV_mean)   
tmp_1 <- merge(IMPOWER133_Clinical,IEV_mean,by="X")

tmp <- merge(tmp_1,deconv_fraction,by="X")
tmp$IEV_rmFibro_score <- tmp$IEV_sig*(1-tmp$Fibroblast)
tmp$dead_date <- as.numeric(tmp$OS_MONTHS)*30
tmp$status <- tmp$OS_CENSOR
genes <- data.frame(t(ICI_Normalized_data[c("IGFBP5","ASCL1"),]))
genes$ID <- rownames(genes)   
final_matrix <- merge(tmp,genes,by.x="X", by.y="ID")
final_matrix$IGFBP5_rmFibro_score <- final_matrix$IGFBP5*(1-final_matrix$Fibroblast)

#Patients were stratified according to the grouping strategy used in the atezolizumab group, 
#with an IEV-high to IEV-low ratio of 9:1.
final_matrix$IEV_with_Endo <- final_matrix$IEV_sig * final_matrix$Endothelial
final_matrix <- final_matrix[order(final_matrix$IEV_with_Endo,decreasing=TRUE),]
final_matrix$IEV_Group <- "Low"
final_matrix$IEV_Group[1:125] <- "High"
final_matrix$IEV_Group <- factor(final_matrix$IEV_Group, levels = c("Low","High"))

df <- final_matrix
df$PFS_days <- df$PFS_MONTHS * 30
surv_obj <- Surv(time = df$PFS_days,
                 event = df$PFS_CENSOR)

fit <- survfit(surv_obj ~ IEV_Group, data = df)
pdf("Fig3F_Placebo_1000days.pdf",width=6,height=5)
ggsurvplot(
  fit,
  data = df,
  pval = TRUE,               
  risk.table = TRUE,          
  palette = "jco",            
  xlab = "PFS (Months)",
  ylab = "Survival Probability",
  legend.title = "IEV",
  surv.median.line = "h", xlim = c(0, 1000)
)
dev.off()



df <- final_matrix
df$OS_days <- df$OS_MONTHS * 30
surv_obj <- Surv(time = df$OS_days,
                 event = df$OS_CENSOR)

fit <- survfit(surv_obj ~ IEV_Group, data = df)
pdf("Fig3F_OS.pdf",width=6,height=5)
ggsurvplot(
  fit,
  data = df,
  pval = TRUE,               
  risk.table = TRUE,          
  palette = "jco",            
  xlab = "OS (Months)",
  ylab = "Survival Probability",
  legend.title = "IEV",
  surv.median.line = "h"
)
dev.off()


#FigS3F

final_matrix$IEV_rmFibro <- final_matrix$IEV_sig*(1-final_matrix$Fibroblast)
ff1 <- ggplot(final_matrix, aes(x=IEV_rmFibro, y= T.cell))+geom_point() + stat_smooth(method=lm)+stat_cor(data=final_matrix, method = "spearman")
ggsave(ff1,file="FigS3F.png",width=8,height=8)




library(survival)
library(survminer)
library(dplyr)
IEV_sig <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/IEV_signature.csv")
deconv_fraction <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/first_submit/ImmuneCheckpoint_blokage_Clinical_data/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_SCLC_deconv_fraction.csv")
IMPOWER133_Clinical <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/first_submit/ImmuneCheckpoint_blokage_Clinical_data/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_Clinical_file.csv")
IMPOWER133_Clinical <- subset(IMPOWER133_Clinical,ACTARM.2=="atezo")
IMPOWER133_Clinical <- IMPOWER133_Clinical[,c("trunc_anonymized_sample_ids","Gay_MDACC_subtypes","PFS_MONTHS", "PFS_CENSOR","BOR")]
IMPOWER133_Clinical$X <- gsub("\\-","\\.",IMPOWER133_Clinical$trunc_anonymized_sample_ids)

ICI_Normalized_data <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/first_submit/ImmuneCheckpoint_blokage_Clinical_data/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_SCLC_TPM_file.csv")
rownames(ICI_Normalized_data) <- ICI_Normalized_data$UNNAMED..0
ICI_Normalized_data <- ICI_Normalized_data[,c(-1,-2)]
ICI_Normalized_data <- ICI_Normalized_data[,IMPOWER133_Clinical$X]

IEV_marker <- intersect(rownames(ICI_Normalized_data),IEV_sig$Symbol)
IEV_gene <- data.frame(ICI_Normalized_data[IEV_marker,1:132])
IEV_mean <- data.frame(apply(IEV_gene,2,mean))
names(IEV_mean) <- "IEV_sig"
IEV_mean$X <- rownames(IEV_mean)   
tmp_1 <- merge(IMPOWER133_Clinical,IEV_mean,by="X")

tmp <- merge(tmp_1,deconv_fraction,by="X")
genes <- data.frame(t(ICI_Normalized_data[c("IGFBP5","ASCL1"),]))
genes$ID <- rownames(genes)   
final_matrix <- merge(tmp,genes,by.x="X", by.y="ID")

final_matrix$IEV_rmFibro <- final_matrix$IEV_sig*(1-final_matrix$Fibroblast)
ff1 <- ggplot(final_matrix, aes(x=IEV_rmFibro, y= T.cell))+geom_point() + stat_smooth(method=lm)+stat_cor(data=final_matrix, method = "spearman")
ggsave(ff1,file="IEV_Tcell_Spearman.svg",width=8,height=8)





#Fig3I
all <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/CancerCell_patient_SCLC_subtype_ANP.csv")
Endo <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/Endothelium/Pure_Endo_CancerCell.rds",mc.cores=20)
Endo_num <- data.frame(table(Endo$patient))
names(Endo_num) <- c("ID","Endo_num")
all_1 <- merge(all,Endo_num,by="ID")
all_2 <- subset(all_1,Endo_num > 0)

SCLC.A <- subset(all_2,subtype=="SCLC.A")
SCLC.P <- subset(all_2,subtype=="SCLC.P")
SCLC.N <- subset(all_2,subtype=="SCLC.N")
NSCLC <- subset(all_2,subtype=="NSCLC")

Idents(Endo) <- Endo$patient
SCLC_A <- subset(Endo, idents=c(as.character(SCLC.A$ID)))
SCLC_P <- subset(Endo, idents=c(as.character(SCLC.P$ID)))
SCLC_N <- subset(Endo, idents=c(as.character(SCLC.N$ID)))
NSCLC <- subset(Endo, idents=c(as.character(NSCLC$ID)))

SCLC_A$Group <- "SCLC.A"
SCLC_P$Group <- "SCLC.P"
SCLC_N$Group <- "SCLC.N"
NSCLC$Group <- "NSCLC"
Endo_renew <- merge(x=SCLC_A,y=c(SCLC_P,SCLC_N,NSCLC))
Endo_renew$Group <- factor(Endo_renew$Group,levels=c("NSCLC","SCLC.P","SCLC.N","SCLC.A"))

Endo_renew$Group_1 <- ifelse(Endo_renew$Group=="SCLC.A","ASCL1_hi","ASCL1_low")
Endo_renew$Group_1 <- factor(Endo_renew$Group_1,levels=c("ASCL1_low","ASCL1_hi"))
Idents(Endo_renew) <- Endo_renew$Group_1

ff <- VlnPlot(Endo_renew,features = c("COL4A1","COL4A2","NID1","LAMC1","LAMB1","LAMA4"),pt.size=0,group="Group_1",ncol=3)
ggsave(ff,file="/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/Endothelium/Fig3I.png",width=8,height=6)





%%%%%%%%%%%%%%%%%%%%%%%%%%T细胞浸润在各个亚型里面回归%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T细胞浸润在各个亚型里面回归%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T细胞浸润在各个亚型里面回归%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T细胞浸润在各个亚型里面回归%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%T细胞浸润在各个亚型里面回归%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ZhangPeng_RNA <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/Cell_ZhangPeng_RNA_TPM_log2_Matrix.csv")
Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/Cell_ZhangPeng_Patients_statification.csv")
ZhangPeng_RNA <- as.data.frame(ZhangPeng_RNA)
rownames(ZhangPeng_RNA) <- ZhangPeng_RNA$Gene.Sample.ID
ZhangPeng_RNA <- ZhangPeng_RNA[,c(-1)]
ZhangPeng_RNA <- ZhangPeng_RNA %>% mutate_all(~ ifelse(is.na(.), 0, .))
ZhangPeng_RNA <- data.frame(t(ZhangPeng_RNA))
Tumor_ID <- subset(Clinical,Immune.subtype !="NAT-enriched")
ZhangPeng_RNA_1 <- ZhangPeng_RNA[Tumor_ID$Sample.ID,]
ZhangPeng_RNA <- data.frame(t(ZhangPeng_RNA_1))
ZP_SCLC <- scale(ZhangPeng_RNA)
ZP_SCLC <- data.frame(t(ZP_SCLC[c("ASCL1","YAP1","NEUROD1","POU2F3"),]))
ZP_SCLC$subtype <- apply(ZP_SCLC, 1, function(row) {colnames(ZP_SCLC)[which.max(row)]})
ZP_group <- data.frame(ID=rownames(ZP_SCLC),subtype=ZP_SCLC$subtype)


Nature_RNA <- fread("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/SCLC_Nature_RNA.csv")
Nature_RNA <- as.data.frame(Nature_RNA)
Nature_RNA <- Nature_RNA[!duplicated(Nature_RNA$gene),]
rownames(Nature_RNA) <- Nature_RNA$gene
Nature_RNA <- Nature_RNA[,c(-1,-2)]
Nature_RNA_matrix <- data.frame(Nature_RNA)
Nature_SCLC <- scale(Nature_RNA_matrix)
Nature_SCLC <- data.frame(t(Nature_SCLC[c("ASCL1","YAP1","NEUROD1","POU2F3"),]))
Nature_SCLC$subtype <- apply(Nature_SCLC, 1, function(row) {colnames(Nature_SCLC)[which.max(row)]})
Nature_group <- data.frame(ID=rownames(Nature_SCLC),subtype=Nature_SCLC$subtype)
All_Group <- rbind(ZP_group,Nature_group)

ZhangPeng_deconv <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/113_ZhangPeng_Cell_bulkRNA_deconv_SCLC_fraction.csv")
Nature_deconv <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/Nature_deconv_SCLC_fraction.csv")
all_deconv <- rbind(ZhangPeng_deconv, Nature_deconv)
T_ratio <- all_deconv[,c("X","T.cell","Fibroblast")]
tmp <- merge(All_Group,T_ratio,by.x="ID",by.y="X" )
tmp$Group <- ifelse(tmp$subtype=="ASCL1", "ASCL1","non-ASCL1")


library(ggpubr)
library(ggalluvial)
library(ggpubr)
library(ggplot2)


my_comparisons <- list(c("NE","non-NE"))
ff <- ggboxplot(tmp, x = "Group", y = "T.cell",
               color = "Group",
                add = "jitter",palette = c("#08519c","#e31a1c")) 
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="SCLC_all_Tratio_n194_1.pdf",width=4,height=6)

my_comparisons <- list(c("NE","non-NE"))
ff <- ggboxplot(tmp, x = "Group", y = "T.cell",
               color = "subtype",
                add = "jitter",palette = c("#08519c","#e31a1c","#000000","#c51b7d")) 
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="SCLC_all_Tratio_n194_2.pdf",width=4,height=6)
