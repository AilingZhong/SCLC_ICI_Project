&&&&&&&&&&&&&&&&&&&&&&&&&&human_and_mouse_SCLC_hi_genes&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&human_and_mouse_SCLC_hi_genes&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&human_and_mouse_SCLC_hi_genes&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&human_and_mouse_SCLC_hi_genes&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&human_and_mouse_SCLC_hi_genes&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


Mice_DEGs <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/How_to_find_IGFBP5/mouse_LUAD_VS_SCLC_FindAllmarker_DEGs_TumorCells.csv")
human_DEGs <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/CancerCell_Rudin/Tumor/human_SCLC_VS_LUAD_FindAllmarkers.csv")
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

# 绘制散点图并标记 ASCL1 和 IGFBP5 为红色点
ff <- ggplot(all_tmp, aes(x = mouse_avg_log2FC, y = human_avg_log2FC)) +
  geom_point(alpha = 0.5, color = "lightgray") +  # 散点图
  geom_point(data = all_tmp[all_tmp$Symbol %in% c("ASCL1", "IGFBP5","YAP1","NEUROD1","POU2F3"), ], aes(x = mouse_avg_log2FC, y = human_avg_log2FC), 
             color = "red", size = 3) +  # 红色点标注ASCL1和IGFBP5
  geom_text(data = all_tmp[all_tmp$Symbol %in% c("ASCL1", "IGFBP5","YAP1","NEUROD1","POU2F3"), ], 
            aes(label = Symbol), vjust = -0.5, color = "red") +  # 为这些基因添加标签
  stat_smooth(method = "lm", se = FALSE, color = "blue") +  # 线性回归趋势线
  stat_cor(data = all_tmp, method = "spearman") +  # Spearman相关性计算
  theme_minimal() +  # 使用简洁的主题
  labs(title = "Spearman Correlation Scatter Plot", 
       x = "Mouse Avg Log2FC", 
       y = "Human Avg Log2FC")

# 展示图形
ggsave(ff,file="mouse_human_commone_genes_5genes.svg",width=8,height=8)





&&&&&&&&&&&&&&&&&&&&&&&&&&&&如何定义BBB-like的signature&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&如何定义BBB-like的signature&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&如何定义BBB-like的signature&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&如何定义BBB-like的signature&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

cd /mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/BBB_zen_me_lai_de_correction

/usr/local/R4.2/bin/R

library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(future)
library(future.apply)
options(future.globals.maxSize = 400 * 1024^3)
plan("multiprocess", workers = 15)
plan()
library(scales)
library(BuenColors)
library(ggplot2)
library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(data.table)
library(future.apply)
library(trqwe)
library(ggplot2)
library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(data.table)
library(future.apply)
library(trqwe)

library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(trqwe)
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(iTALK)
library(future.apply)
library(scde)
library(ggplot2)
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")



#brain
final_ctrl_ec_cb <- readRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Reference_Data/Human_BBB/Winkler_et_al_Science_2022/final_ctrl_ec_cb.rds")
final_ctrl_ec_cb$celltype <- "Endothelial"
final_ctrl_ec_cb@meta.data <- final_ctrl_ec_cb@meta.data[,c("orig.ident","nCount_SCT","nFeature_SCT", "Sample","clusters")]
names(final_ctrl_ec_cb@meta.data) <- c("orig.ident","nCount_RNA","nFeature_RNA","sample","CellType")
final_ctrl_ec_cb$Cancer_Type <- "Normal_Brain"
mcsaveRDS(final_ctrl_ec_cb,file="/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Brain/Normal_Brain_Science.rds",mc.cores=20)

Nature_Endo <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Brain/Manolis_Kellis_2022_Nature/ROSMAP_Nature_Endo_Harmony_UMAP20.rds",mc.cores=20)
Nature_Endo$orig.ident <- Nature_Endo$SampleID
Nature_Endo$sample <- Nature_Endo$SampleID
Nature_Endo$CellType <- "Endothelial"
Nature_Endo@meta.data <- Nature_Endo@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","sample","CellType")]
Nature_Endo$Cancer_Type <- "Normal_Brain"
mcsaveRDS(Nature_Endo,file="/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Brain/Normal_Brain_Nature.rds",mc.cores=20)

#Lung
Lung <- Read10X(data.dir = "/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Lung/Circulation_GSE164829/", gene.column = 1)
Lung_object <- CreateSeuratObject(counts = Lung, project = "Lung", min.cells = 3, min.features = 200)
Lung_object$sample <- Lung_object$orig.ident
Lung_object$CellType <- "Endo"
names(Lung_object@meta.data) <- c("orig.ident","nCount_RNA","nFeature_RNA","sample","CellType")
Lung_object$Cancer_Type <- "Normal_Lung"
mcsaveRDS(Lung_object,file="/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Lung/Lung_Normal.rds",mc.cores=20)



Brain_Endo_1 <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Brain/Normal_Brain_Nature.rds",mc.cores=20)
Brain_Endo_2 <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Brain/Normal_Brain_Science.rds",mc.cores=20)
Normal_Lung <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/Normal_Endo/Lung/Lung_Normal.rds",mc.cores=20)
Brain_Endo_2@assays$RNA <- Brain_Endo_2@assays$SCT


All_Merge <- merge(x=Brain_Endo_1,y=c(Brain_Endo_2,Normal_Lung))
All_MERGE_DATA <- All_Merge %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
All_MERGE_DATA.pca <- RunPCA(All_MERGE_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(All_MERGE_DATA.pca,file="human_Brain_Lung_MERGE_DATA_pca.rds",mc.cores=20)

All_FindNeighbors <- FindNeighbors(All_MERGE_DATA.pca,reduction = "pca", dims = 1:50)
All_FindClusters <- FindClusters(All_FindNeighbors,resolution = c(0.1,0.2,0.5,1))
All_UMAP20 <- RunUMAP(All_FindClusters, dims = 1:20)
mcsaveRDS(All_UMAP20,file="human_Brain_Lung_AllCells_UMAP20.rds", mc.cores = 20)


Idents(All_UMAP20) <- All_UMAP20$Cancer_Type
all.markers <- FindAllMarkers(All_UMAP20, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"human_Brain_and_Lung_FindAllMarker.csv")



BBB_signature <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/define_BBB_sig/human_Brain_and_Lung_FindAllMarker.csv")
BBB_signature <- subset(BBB_signature, cluster=="Normal_Brain" & p_val_adj < 0.01 & avg_log2FC > 0 )
SCLC_Endo_DEGs <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/CancerCell_Rudin/Endothelium/FindAllmarker_CancerCell_Endothelium_Hiso.csv")
SCLC_Endo_hi <- subset(SCLC_Endo_DEGs,cluster=="SCLC" & p_val_adj < 0.01 & avg_logFC > 0 )
common_genes <- data.frame(intersect(BBB_signature$gene, SCLC_Endo_hi$gene))
names(common_genes) <- "human_gene"
common_genes <- common_genes %>% mutate(mouse_gene = convert_human_to_mouse_symbols(as.character(common_genes$human_gene))) %>% drop_na()
write.csv(common_genes,file="final_human_Brain_hi_and_SCLC_hi_common_genes.csv")




BBB_signature <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/define_BBB_sig/human_Brain_and_Lung_FindAllMarker.csv")
Brain_Endo_hi <- subset(BBB_signature, cluster=="Normal_Brain" & p_val_adj < 0.01 & avg_log2FC > 0 )
SCLC_Endo_DEGs <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/CancerCell_Rudin/Endothelium/FindAllmarker_CancerCell_Endothelium_Hiso.csv")
SCLC_Endo_hi <- subset(SCLC_Endo_DEGs,cluster=="SCLC" & p_val_adj < 0.01 & avg_logFC > 0 )
common_genes <- data.frame(intersect(Brain_Endo_hi$gene, SCLC_Endo_hi$gene))
names(common_genes) <- "human_gene"
common_genes <- common_genes %>% mutate(mouse_gene = convert_human_to_mouse_symbols(as.character(common_genes$human_gene)))
write.csv(common_genes,file="490_genes_BBB_like_signature_genes.csv")


require(gmp)
enrich_pvalue <- function(N, A, B, k)
{
    m <- A + k
    n <- B + k
    i <- k:min(m,n)

    as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
}

aa2 <- enrich_pvalue(24421,8669,657,490)

library(Vennerable)
data<-Venn(list("Brain.Endo"=unique(na.omit((Brain_Endo_hi$gene))),"SCLC.Endo"=unique(na.omit(SCLC_Endo_hi$gene))))
pdf("Brain_SCLC_Endo_hi_Veen.pdf")
plot(data,doWeight=T)
dev.off()



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


#在两个队列中,使用limma
Nature_RNA <- fread("/mnt/data/user_data/xiangyu/workshop/DATABASE/SCLC_Patient_data/SCLC_Nature_RNA.csv")
Nature_RNA <- as.data.frame(Nature_RNA)
Nature_RNA <- Nature_RNA[!duplicated(Nature_RNA$gene),]
rownames(Nature_RNA) <- Nature_RNA$gene
Nature_RNA <- Nature_RNA[,c(-1,-2)]
Nature_RNA_matrix <- data.frame(Nature_RNA)
ZhangPeng_RNA <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_RNA_TPM_log2_Matrix.csv")
ZhangPeng_RNA <- as.data.frame(ZhangPeng_RNA)
rownames(ZhangPeng_RNA) <- ZhangPeng_RNA$Gene.Sample.ID
ZhangPeng_RNA <- ZhangPeng_RNA[,c(-1)]
ZhangPeng_RNA <- ZhangPeng_RNA %>% mutate_all(~ ifelse(is.na(.), 0, .))
ZhangPeng_RNA <- data.frame(t(ZhangPeng_RNA))
ZhangPeng_RNA$Group <- substring(rownames(ZhangPeng_RNA),1,1)
ZhangPeng_RNA <- subset(ZhangPeng_RNA,Group=="T")
ZhangPeng_RNA_matrix <- ZhangPeng_RNA[,-17236]
ZhangPeng_RNA_matrix <- data.frame(t(ZhangPeng_RNA_matrix))
genes <- intersect(rownames(Nature_RNA_matrix),rownames(ZhangPeng_RNA_matrix))
Nature_RNA_matrix <- Nature_RNA_matrix[genes,]
ZhangPeng_RNA_matrix <- ZhangPeng_RNA_matrix[genes,]
all_matrix <- cbind(Nature_RNA_matrix,ZhangPeng_RNA_matrix)


library(limma)
batch <- factor(c(rep(1, ncol(Nature_RNA_matrix)), rep(2, ncol(ZhangPeng_RNA_matrix))))
all_matrix_corrected <- removeBatchEffect(all_matrix, batch = batch)
mcsaveRDS(all_matrix_corrected,file="Nature2014_and_Zhangpeng_RNA_limma.rds",mc.cores=20)



sel_genes <- intersect(rownames(all_matrix_corrected),common_genes$human_gene)
RNAseq_1 <- all_matrix_corrected[sel_genes,]
SCLC_exp_1 <- data.frame(apply(RNAseq_1,2,mean)) 
names(SCLC_exp_1) <- "BBB_sig"
all_matrix_corrected_1 <- data.frame(t(all_matrix_corrected))
all_matrix_corrected_1$BBB_sig <- SCLC_exp_1$BBB_sig

All_patients <- all_matrix_corrected_1[,c("YAP1","NEUROD1" ,"ASCL1" ,"POU2F3","BBB_sig")]
gene_name<-c()
cor_r<-c()
pvalue<-c()
for (i in 1:ncol(All_patients)){
    g1=colnames(All_patients)[i]
    c_r=cor(as.numeric(All_patients$BBB_sig),as.numeric(All_patients[,i]),method="spearman")
    p=cor.test(as.numeric(All_patients$BBB_sig),as.numeric(All_patients[,i]),method ="spearman")[[3]]
    gene_name=c(gene_name,g1)
    cor_r=c(cor_r,c_r)
    pvalue=c(pvalue,p)
       }   
data_cor <- data.frame(gene_name,cor_r,pvalue)
data_cor <- data_cor[order(-data_cor[,2]),] 
data_cor$order <- c(1:length(cor_r))
key_genes_limma <- subset(data_cor,gene_name=="YAP1" | gene_name=="NEUROD1" | gene_name=="ASCL1" | gene_name=="POU2F3")

key_genes_limma$p_value <- -log10(key_genes_limma$pvalue)

library(ggplot2)


ff <- ggplot(key_genes_limma, aes(x = gene_name, y = 1, size = p_value, color = cor_r)) +
  geom_point() +  # 绘制点图
  scale_size_continuous(
    range = c(0.5,5),  # 设置圆圈大小范围
    breaks = c(0.5,1,1.5,2,2.5,3,3.5,4,4.5,5),  # 设置三个梯度的breaks
    labels = c("0.5","1","1.5","2","2.5","3","3.5","4","4.5","5")  # 为这些大小设置标签
  ) + 
  scale_color_gradientn(
    colors = c("#2166ac","#ffffbf","#fee090", "#d6604d","#b2182b"), 
    values = scales::rescale(c(-0.15,0,0.15,0.3,0.35)), 
    limits = c(-0.15, 0.35)  # 设置颜色渐变
  ) +  
  labs(title = "Dotplot of Gene Correlation and P-value", x = "Gene Name", y = "") +  # 设置标题，去掉y轴标签
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  # 使X轴标签倾斜，避免重叠

# 保存为SVG文件
ggsave("final_in_2cohorts_correlation.svg", plot = ff)



#在各种组学中验证BBB-LIKE SIG
BBB.sig <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/final_human_Brain_hi_and_SCLC_hi_common_genes.csv")

****************************在蛋白水平验证*****************************************************************
****************************在蛋白水平验证*****************************************************************
****************************在蛋白水平验证*****************************************************************
****************************在蛋白水平验证*****************************************************************

normalized_protein <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Protomics/Processed_Data/Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")
common_scRNA_Protein <- intersect(BBB.sig$mouse_gene,normalized_protein$Gene)

common_scRNA_Protein_matrix <- subset(normalized_protein, Gene %in% common_scRNA_Protein)
rownames(common_scRNA_Protein_matrix) <- common_scRNA_Protein_matrix$Gene

Normal <- common_scRNA_Protein_matrix[,c("Normalized_Normal.1","Normalized_Normal.2","Normalized_Normal.3")]
LUAD <- common_scRNA_Protein_matrix[,c("Normalized_NSCLC.1","Normalized_NSCLC.2","Normalized_NSCLC.3")]
SCLC <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgScr.1","Normalized_SCLC_sgScr.2")]
sgIGFBP5 <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgIgfbp5.1","Normalized_SCLC_sgIgfbp5.2")]

Normal_basemean <- data.frame(apply(Normal,1,mean))
LUAD_basemean <- data.frame(apply(LUAD,1,mean))
SCLC_basemean <- data.frame(apply(SCLC,1,mean))
sgIGFBP5_basemean <- data.frame(apply(sgIGFBP5,1,mean))

names(Normal_basemean) <- "Normal_value"
names(LUAD_basemean) <- "LUAD_value"
names(SCLC_basemean) <- "SCLC_value"
names(sgIGFBP5_basemean) <- "sgIGFBP5_value"

tmp <- cbind(Normal_basemean,LUAD_basemean,SCLC_basemean,sgIGFBP5_basemean)

Normal_basemean$Group <- "Normal"
LUAD_basemean$Group <- "LUAD"
SCLC_basemean$Group <- "SCLC"
sgIGFBP5_basemean$Group <- "sgIGFBP5"


Normal <- common_scRNA_Protein_matrix[,c("Normalized_Normal.1","Normalized_Normal.2","Normalized_Normal.3")]
LUAD <- common_scRNA_Protein_matrix[,c("Normalized_NSCLC.1","Normalized_NSCLC.2","Normalized_NSCLC.3")]
SCLC <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgScr.1","Normalized_SCLC_sgScr.2")]
sgIGFBP5 <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgIgfbp5.1","Normalized_SCLC_sgIgfbp5.2")]

Normal_basemean <- data.frame(apply(Normal,1,mean))
LUAD_basemean <- data.frame(apply(LUAD,1,mean))
SCLC_basemean <- data.frame(apply(SCLC,1,mean))
sgIGFBP5_basemean <- data.frame(apply(sgIGFBP5,1,mean))

names(Normal_basemean) <- "MS_value"
names(LUAD_basemean) <- "MS_value"
names(SCLC_basemean) <- "MS_value"
names(sgIGFBP5_basemean) <- "MS_value"

Normal_basemean$Group <- "Normal"
LUAD_basemean$Group <- "LUAD"
SCLC_basemean$Group <- "SCLC"
sgIGFBP5_basemean$Group <- "sgIGFBP5"

tmp <- rbind(Normal_basemean,LUAD_basemean,SCLC_basemean,sgIGFBP5_basemean)
my_comparisons <- list(c("SCLC", "Normal"),c("SCLC","LUAD"),c("SCLC","sgIGFBP5"))



# 加载必要的包
library(ggpubr)
library(dplyr)
library(ggplot2)
library(pheatmap)

tmp$MS_value[tmp$MS_value > 5000] <- 5000
median_values <- tmp %>%
  group_by(Group) %>%
  summarise(median_value = median(MS_value))
ff <- ggboxplot(tmp, x = "Group", y = "MS_value",add = "", color = "Group",palette = c("#08519c","#e31a1c","#000000","#1c9099"))+ylim(0,8000)+
stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = TRUE, method="wilcox.test") 
ggsave(ff,file="Protein_BBB_like_1.svg",width=6,height=6)


ff <- ggbarplot(tmp, x = "Group", y = "MS_value", fill="Group",
 add = c("mean_se"), legend = "none", ylim=c(0,2000),title="")+
 rotate_x_text(angle = 45)+scale_fill_manual(values=c("#377eb8","#e41a1c","#1b9e77","#d95f02","#7570b3","#e7298a","#e6ab02","#a6761d","#666666"))+
   geom_line(data = median_values, aes(x = Group, y = median_value, group = 1), color = "gray", size = 1) +
  geom_point(data = median_values, aes(x = Group, y = median_value), color = "gray", size = 3)
ggsave(ff,file="Protein_BBB_like_2.svg",width=6,height=6)









IGF1_IGF1R_sig <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/IGF1_1R_signature_correlation/IGF1_IGF1R_signaling_pathway_renew.csv")

normalized_protein <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Protomics/Processed_Data/Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")
common_scRNA_Protein <- intersect(IGF1_IGF1R_sig$mouse_gene,normalized_protein$Gene)

common_scRNA_Protein_matrix <- subset(normalized_protein, Gene %in% common_scRNA_Protein)
rownames(common_scRNA_Protein_matrix) <- common_scRNA_Protein_matrix$Gene

Normal <- common_scRNA_Protein_matrix[,c("Normalized_Normal.1","Normalized_Normal.2","Normalized_Normal.3")]
LUAD <- common_scRNA_Protein_matrix[,c("Normalized_NSCLC.1","Normalized_NSCLC.2","Normalized_NSCLC.3")]
SCLC <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgScr.1","Normalized_SCLC_sgScr.2")]
sgIGFBP5 <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgIgfbp5.1","Normalized_SCLC_sgIgfbp5.2")]

Normal_basemean <- data.frame(apply(Normal,1,mean))
LUAD_basemean <- data.frame(apply(LUAD,1,mean))
SCLC_basemean <- data.frame(apply(SCLC,1,mean))
sgIGFBP5_basemean <- data.frame(apply(sgIGFBP5,1,mean))

names(Normal_basemean) <- "Normal_value"
names(LUAD_basemean) <- "LUAD_value"
names(SCLC_basemean) <- "SCLC_value"
names(sgIGFBP5_basemean) <- "sgIGFBP5_value"

tmp <- cbind(Normal_basemean,LUAD_basemean,SCLC_basemean,sgIGFBP5_basemean)

Normal_basemean$Group <- "Normal"
LUAD_basemean$Group <- "LUAD"
SCLC_basemean$Group <- "SCLC"
sgIGFBP5_basemean$Group <- "sgIGFBP5"


Normal <- common_scRNA_Protein_matrix[,c("Normalized_Normal.1","Normalized_Normal.2","Normalized_Normal.3")]
LUAD <- common_scRNA_Protein_matrix[,c("Normalized_NSCLC.1","Normalized_NSCLC.2","Normalized_NSCLC.3")]
SCLC <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgScr.1","Normalized_SCLC_sgScr.2")]
sgIGFBP5 <- common_scRNA_Protein_matrix[,c("Normalized_SCLC_sgIgfbp5.1","Normalized_SCLC_sgIgfbp5.2")]

Normal_basemean <- data.frame(apply(Normal,1,mean))
LUAD_basemean <- data.frame(apply(LUAD,1,mean))
SCLC_basemean <- data.frame(apply(SCLC,1,mean))
sgIGFBP5_basemean <- data.frame(apply(sgIGFBP5,1,mean))

names(Normal_basemean) <- "MS_value"
names(LUAD_basemean) <- "MS_value"
names(SCLC_basemean) <- "MS_value"
names(sgIGFBP5_basemean) <- "MS_value"

Normal_basemean$Group <- "Normal"
LUAD_basemean$Group <- "LUAD"
SCLC_basemean$Group <- "SCLC"
sgIGFBP5_basemean$Group <- "sgIGFBP5"

tmp <- rbind(Normal_basemean,LUAD_basemean,SCLC_basemean,sgIGFBP5_basemean)
my_comparisons <- list(c("SCLC", "Normal"),c("SCLC","LUAD"),c("SCLC","sgIGFBP5"))



# 加载必要的包
library(ggpubr)
library(dplyr)
library(ggplot2)
library(pheatmap)

tmp$MS_value[tmp$MS_value > 5000] <- 5000
median_values <- tmp %>%
  group_by(Group) %>%
  summarise(median_value = median(MS_value))
ff <- ggboxplot(tmp, x = "Group", y = "MS_value",add = "", color = "Group",palette = c("#08519c","#e31a1c","#000000","#1c9099"))+ylim(0,8000)+
stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = TRUE, method="wilcox.test") 
ggsave(ff,file="Protein_IGF1_1R_like_1.svg",width=6,height=6)


ff <- ggbarplot(tmp, x = "Group", y = "MS_value", fill="Group",
 add = c("mean_se"), legend = "none", ylim=c(0,3000),title="")+
 rotate_x_text(angle = 45)+scale_fill_manual(values=c("#377eb8","#e41a1c","#1b9e77","#d95f02","#7570b3","#e7298a","#e6ab02","#a6761d","#666666"))+
   geom_line(data = median_values, aes(x = Group, y = median_value, group = 1), color = "gray", size = 1) +
  geom_point(data = median_values, aes(x = Group, y = median_value), color = "gray", size = 3)
ggsave(ff,file="Protein_IGF1_1R_like_2.svg",width=6,height=6)











****************************在单细胞水平验证*****************************************************************
****************************在单细胞水平验证*****************************************************************
****************************在单细胞水平验证*****************************************************************
****************************在单细胞水平验证*****************************************************************

BBB.sig <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/final_human_Brain_hi_and_SCLC_hi_common_genes.csv")
Mice_model <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_SCLC_Final_Figures/Fig1_S1/Final_Mice_only_LUAD_SCLC.rds",mc.cores=20)
Idents(Mice_model) <- Mice_model$sub_CellType
Endo <- subset(Mice_model, idents=c("Endo"),invert=FALSE)

BBB_like.sig <- as.character(BBB.sig$mouse_gene)
BBB_like.sig <- intersect(rownames(GetAssayData(object = Endo)),BBB_like.sig)
speci_raw <- FetchData(object = Endo, vars = BBB_like.sig,slot="data")
Endo[["BBB_like.sig"]] <- (rowSums(speci_raw))/length(BBB_like.sig)

Mice_Endo <- Endo[[]]

library(ggpubr)
library(ggplot2)

my_comparisons <- list(c("LUAD","SCLC"))
ff <- ggboxplot(Mice_Endo, x = "Group", y = "BBB_like.sig",
               color = "Group",palette = c("#08519c","#e31a1c","#000000","#238b45"))+
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="Mice_scRNA_BBB_sig.svg",width=4,height=8)

  
sgScr_vs_sgIgfbp5 <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Mice_scRNAseq/sgIgfbp5_vs_sgScr_scRNA_ProcessedData/Renames_AllCells_UMAP20.rds",mc.cores=20)
Idents(sgScr_vs_sgIgfbp5) <- sgScr_vs_sgIgfbp5$CellType
Pure_Endo <- subset(sgScr_vs_sgIgfbp5,idents=c("Endo"),invert=FALSE)
BBB_like.sig <- as.character(BBB.sig$mouse_gene)
BBB_like.sig <- intersect(rownames(GetAssayData(object = Pure_Endo)),BBB_like.sig)
speci_raw <- FetchData(object = Pure_Endo, vars = BBB_like.sig,slot="data")
Pure_Endo[["BBB_like.sig"]] <- (rowSums(speci_raw))/length(BBB_like.sig)
sgIGFBP5_Pure_Endo <- Pure_Endo[[]]


library(ggpubr)
library(ggplot2)
my_comparisons <- list(c("sgScr","sgIgfbp5"))
ff <- ggboxplot(sgIGFBP5_Pure_Endo, x = "Group", y = "BBB_like.sig",
               color = "Group",palette = c("#08519c","#e31a1c","#000000","#238b45"))+
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="Mice_scRNA_sgIGFBP5_BBB_like.sig.svg",width=4,height=8)











****************************在2个SCLC队列数据中验证ASCL1—IGFBP5排名*****************************************************************
****************************在2个SCLC队列数据中验证ASCL1—IGFBP5排名*****************************************************************
****************************在2个SCLC队列数据中验证ASCL1—IGFBP5排名*****************************************************************
****************************在2个SCLC队列数据中验证ASCL1—IGFBP5排名*****************************************************************


BBB.sig <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/define_BBB_sig/final_human_Brain_hi_and_SCLC_hi_common_genes.csv")
RNAseq <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/Nature2014_and_Zhangpeng_RNA_limma.rds",mc.cores=20)

sel_genes <- intersect(rownames(RNAseq),BBB.sig$human_gene)
RNAseq_1 <- RNAseq[sel_genes,]
SCLC_exp_1 <- data.frame(apply(RNAseq_1,2,mean)) 
names(SCLC_exp_1) <- "BBB.sig"
RNAseq_1 <- data.frame(t(RNAseq))
RNAseq_1$BBB.sig <- SCLC_exp_1$BBB.sig

gene_name<-c()
cor_r<-c()
pvalue<-c()
for (i in 1:ncol(RNAseq_1)){
    g1=colnames(RNAseq_1)[i]
    c_r=cor(as.numeric(RNAseq_1$BBB.sig),as.numeric(RNAseq_1[,i]),method="spearman")
    p=cor.test(as.numeric(RNAseq_1$BBB.sig),as.numeric(RNAseq_1[,i]),method ="spearman")[[3]]
    gene_name=c(gene_name,g1)
    cor_r=c(cor_r,c_r)
    pvalue=c(pvalue,p)
       }   
data_cor <- data.frame(gene_name,cor_r,pvalue)
data_cor <- data_cor[order(-data_cor[,2]),] 
data_cor$order <- c(1:length(cor_r))
write.csv(data_cor,file="SCLC_2_cohorts_BBB_sig_data_cor_spearman.csv")



#这里不能用bulk中的SCLC VS NSCLC，因为在单细胞数据中已经是SCLC中高的了，在bulk中比较没有意义
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
library(limma)

SCLC_2cohorts <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/Nature2014_and_Zhangpeng_RNA_limma.rds",mc.cores=20)

SCLC_matrix <- data.frame(t(SCLC_2cohorts))
expression_abundance <- data.frame(apply(SCLC_matrix,1,sum)) 
names(expression_abundance) <- "expression_abundance"
SCLC_matrix$expression_abundance <- expression_abundance$expression_abundance


data_cor <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_BBB_zen_me_lai_de_correction/SCLC_2_cohorts_BBB_sig_data_cor_spearman.csv")
Mice_DEGs <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/How_to_find_IGFBP5/mouse_LUAD_VS_SCLC_FindAllmarker_DEGs_TumorCells.csv")
human_DEGs <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/CancerCell_Rudin/Tumor/human_SCLC_VS_LUAD_FindAllmarkers.csv")
Mice_DEGs <- Mice_DEGs %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(Mice_DEGs$gene))) %>% drop_na()

Mice_SCLC_hi <- subset(Mice_DEGs,cluster=="SCLC" & p_val_adj < 0.01 & avg_log2FC > 0 )
human_SCLC_hi <- subset(human_DEGs,cluster=="SCLC" & p_val_adj < 0.01 & avg_logFC > 0 )
Mice_SCLC_hi <- Mice_SCLC_hi[,c("human_gene","avg_log2FC")]
human_SCLC_hi <- human_SCLC_hi[,c("gene","avg_logFC")]
names(Mice_SCLC_hi) <- c("Symbol","avg_log2FC")
names(human_SCLC_hi) <- c("Symbol","avg_log2FC")



#这里是人和小鼠SCLC都高的基因
genes_ovelap <- intersect(Mice_SCLC_hi$Symbol,human_SCLC_hi$Symbol)
genes_ovelap <- intersect(genes_ovelap,rownames(SCLC_2cohorts))
DEGs_1 <- data.frame(SCLC_2cohorts[genes_ovelap,])
expression_abundance <- data.frame(apply(DEGs_1,1,sum)) 
names(expression_abundance) <- "expression_abundance"
DEGs_1$expression_abundance <- expression_abundance$expression_abundance

rownames(data_cor) <- data_cor$gene_name
cor_1 <- data_cor[genes_ovelap,]
tmp_rank <- cbind(DEGs_1,cor_1)

tmp_rank <- merge(tmp_rank,human_DEGs,by.x="gene_name",by.y="gene")
all_tmp <- tmp_rank[,c("gene_name","expression_abundance","cor_r","pvalue","avg_logFC")]

all_tmp$cor.R_log2FC_expression_abundance <-  all_tmp$cor_r*all_tmp$expression_abundance*10000*all_tmp$avg_logFC
all_tmp$Group <- ifelse(all_tmp$cor.R_log2FC_expression_abundance < 0, "Negtive","Positive")
Negtive_group <- subset(all_tmp,Group=="Negtive")
Negtive_group$log10_cor.R_log2FC_expression_abundance <- log(abs(Negtive_group$cor.R_log2FC_expression_abundance),10)*-1
Positive_group <- subset(all_tmp,Group=="Positive")
Positive_group$log10_cor.R_log2FC_expression_abundance <- log(abs(Positive_group$cor.R_log2FC_expression_abundance),10)*1
all_tmp <- rbind(Negtive_group, Positive_group)

all_tmp <- all_tmp[order(-all_tmp$log10_cor.R_log2FC_expression_abundance),] 
write.csv(all_tmp,file="correct_IGFBP5_rank_in_ZhangPeng_and_Nature_log10_cor.R_log2FC_expression_abundance.csv")




all_tmp_1 <- subset(all_tmp,pvalue < 0.05 & cor_r > 0)
all_tmp_1$label <- ifelse(all_tmp_1$gene_name=="IGFBP5" | all_tmp_1$gene_name=="ASCL1","famous","non_famous")
library(ggpubr)

ff <- ggdotchart(all_tmp_1, x = "gene_name", y = "log10_cor.R_log2FC_expression_abundance",
           color = "label",                                # Color by groups
           palette = c("#00AFBB", "#E7B800"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = FALSE,                                # Rotate vertically
           dot.size = 2,                                 # Large dot size
           label = round(all_tmp_1$log10_cor.R_log2FC_expression_abundance),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 4,
                             vjust = 0.3),               # Adjust label parameters
           ggtheme = theme_pubr())                 # ggplot2 theme
ggsave(ff,file="ASCL1_rank_log10.svg",height=10,width=10)





all_tmp <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_BBB_zen_me_lai_de_correction/key_files/correct_IGFBP5_rank_in_ZhangPeng_113_and_Nature_81_log10_cor.R_log2FC_expression_abundance.csv")
ISG_pos <- subset(all_tmp,pvalue < 0.05 & cor_r > 0)

#ASCL1调控基因
ASCL1_binding <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_ASCL1_CUTTag/GSE155692_Trudy_Oliver_ASCL1_CUTtag/macs2_ASCL1_binding_promoter.csv")
ASCL1_binding <- ASCL1_binding %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(ASCL1_binding$SYMBOL))) %>% drop_na()
TRMA_VS_TRM <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_ASCL1_CUTTag/GSE155692_Trudy_Oliver_ASCL1_CUTtag/PRMA_VS_PRM_TrudyOliver.csv")
TRMA_VS_TRM <- subset(TRMA_VS_TRM, log2FoldChange < -2.5 & pvalue < 0.05 )


ASCL1_downstream <- data.frame(intersect(unique(na.omit(TRMA_VS_TRM$human_homolog)), unique(na.omit(ASCL1_binding$human_gene))))
names(ASCL1_downstream) <- "Symbol"

library(Vennerable)
data<-Venn(list("RNA_TRMA_DN"=unique(na.omit(TRMA_VS_TRM$human_homolog)),"ASCL1.chip"=unique(na.omit(ASCL1_binding$human_gene))))
pdf("ASCL1_Chip_RNA_overlap.pdf")
plot(data,doWeight=T)
dev.off()


library(Vennerable)
data <- Venn(list("Cor_genes"=unique(na.omit(sig_tmp$gene_name)),
  "RNA_TRMA_DN"=unique(na.omit(TRMA_VS_TRM$human_homolog)),
  "ASCL1.chip"=unique(na.omit(ASCL1_binding$human_gene))))
pdf("ASCL1_Regulation_genes_and_RankC_genes_overlap.pdf")
plot(data,doWeight=T)
dev.off()

tmp <- subset(sig_tmp, sig_tmp$gene_name %in% ASCL1_downstream$Symbol )
tmp$label <- ifelse(tmp$gene_name=="IGFBP5" | tmp$gene_name=="ASCL1","famous","non_famous")
library(ggpubr)
ff <- ggdotchart(tmp, x = "gene_name", y = "log10_cor.R_log2FC_expression_abundance",
           color = "label",                                # Color by groups
           palette = c("#00AFBB", "#E7B800"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = FALSE,                                # Rotate vertically
           dot.size = 5,                                 # Large dot size
           label = round(tmp$log10_cor.R_log2FC_expression_abundance),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 4,
                             vjust = 0.3),               # Adjust label parameters
           ggtheme = theme_pubr())                 # ggplot2 theme
ggsave(ff,file="IGFBP5_rank_log10.svg",height=10,width=10)




#data <- read.table(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_ASCL1_CUTTag/GSE155692_Trudy_Oliver_ASCL1_CUTtag/GSE155691_15688R_raw_counts.txt.gz"), header = TRUE, sep = "\t")  # 读取文件

com_172genes <- intersect(sig_tmp$gene_name,ASCL1_downstream$Symbol)
TRMA_VS_TRM <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_ASCL1_CUTTag/GSE155692_Trudy_Oliver_ASCL1_CUTtag/PRMA_VS_PRM_TrudyOliver.csv")
TRMA_VS_TRM <- subset(TRMA_VS_TRM, log2FoldChange < -2.5 & pvalue < 0.05 )
exp_172genes <- subset(TRMA_VS_TRM, TRMA_VS_TRM $ human_homolog %in% com_172genes)
exp_172genes <- exp_172genes[,c("human_homolog","log2FoldChange")]
rownames(exp_172genes) <- exp_172genes$human_homolog
exp_172genes <- data.frame(exp_172genes[,-1])

library(pheatmap)
pdf("172_common_genes_heatmap_logFC.pdf",width=3,height=8)
pheatmap(exp_172genes,color = colorRampPalette(c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b"))(50),
  show_rownames=TRUE,cluster_row = FALSE,cluster_col= FALSE,border=FALSE)
dev.off()






