#FigS5

library(data.table)
library(sva)
library(limma)
library(Seurat)
library(GSEABase)
library(Biobase)
library(genefilter)
library(limma)
library(RColorBrewer)
library(GSVA)
library(dplyr)
library(Seurat)
library(trqwe)
library(tidyverse)
library(nichenetr)


deconv <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/TCGA_deconv_LUAD_fraction.csv")
IEV_sig <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/IEV_signature.csv")
IGF1_IGF1R_sig <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IGF1_IGF1R_sig.csv")

TCGA_LUAD_matrix <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/TCGA_LUAD_FPKM_matrix.rds",mc.cores=20)
batch_2 <- data.frame(ID=c(rownames(TCGA_LUAD_matrix)),batch="LUAD")
batch_2$malignant <- substring(batch_2$ID,14,15)
batch_2 <- subset(batch_2,malignant < 10)
batch_2 <- batch_2[,-3]
TCGA_LUAD_matrix_1 <- data.frame(TCGA_LUAD_matrix)

sel_genes <- TCGA_LUAD_matrix_1[,c("IGFBP5","ASCL1")]
IEV_sig <- intersect(IEV_sig$Symbol,colnames(TCGA_LUAD_matrix_1))
IEV_sig <- TCGA_LUAD_matrix_1[,IEV_sig]
IEV_sig <- data.frame(apply(IEV_sig,1,mean))
names(IEV_sig) <- "IEV_sig"

tmp <- cbind(sel_genes,IEV_sig)
tmp$X <- rownames(tmp)
tmp <- merge(tmp, deconv, by="X")
tmp$IEV_rmFibro <- tmp$IEV_sig*(1-tmp$Fibroblast)

tmp$group <- ifelse(tmp$ASCL1 > 200,"ASCL1_high","ASCL1_high")
tmp$group <- factor(tmp$group,levels = c("ASCL1_high","ASCL1_exp"))

cols <- c("ASCL1_exp" = "#969696","ASCL1_high" = "#a6cee3")
my_comparisons <- list(c("ASCL1_exp","ASCL1_high"))
ff <- ggboxplot(tmp, x = "group", y = "IEV_rmFibro",color = "group",add = "jitter",palette = cols) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),legend.position = "none")
ggsave("FigS5_IEV_sig_expression.svg", ff, width = 4, height = 6)

ff <- ggboxplot(tmp, x = "group", y = "ASCL1",color = "group",add = "jitter",palette = cols) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),legend.position = "none")
ggsave("FigS5_ASCL1_expression.svg", ff, width = 4, height = 6)

ff <- ggboxplot(tmp, x = "group", y = "IGFBP5",color = "group",add = "jitter",palette = cols) +
  stat_compare_means(comparisons = my_comparisons, label = "p.format", method = "wilcox.test") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),legend.position = "none")
ggsave("FigS5_IGFBP5_expression.svg", ff, width = 4, height = 6)

