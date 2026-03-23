****************************Fig4_and_FigS4*****************************************************************
****************************Fig4_and_FigS4*****************************************************************
****************************Fig4_and_FigS4*****************************************************************
****************************Fig4_and_FigS4*****************************************************************
#Fig 4A
all_tmp <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/1_SCLC_bulk_cohort_log10_cor.R_log2FC_expression_abundance.csv")
sig_tmp <- subset(all_tmp,pvalue < 0.05 & cor_r > 0)

#ASCL1 binding genes
ASCL1_binding <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/TrudyOliver_ASCL1_RNAseq_Chip/macs2_ASCL1_binding_promoter.csv")
ASCL1_binding <- ASCL1_binding %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(ASCL1_binding$SYMBOL))) %>% drop_na()
TRMA_VS_TRM <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/TrudyOliver_ASCL1_RNAseq_Chip/PRMA_VS_PRM_TrudyOliver.csv")
TRMA_VS_TRM <- subset(TRMA_VS_TRM, log2FoldChange < -2.5 & pvalue < 0.05 )

ASCL1_downstream <- data.frame(intersect(unique(na.omit(TRMA_VS_TRM$human_homolog)), unique(na.omit(ASCL1_binding$human_gene))))
names(ASCL1_downstream) <- "Symbol"

library(Vennerable)
data<-Venn(list("RNA_TRMA_DN"=unique(na.omit(TRMA_VS_TRM$human_homolog)),"ASCL1.chip"=unique(na.omit(ASCL1_binding$human_gene))))
pdf("Fig4A_1.pdf")
plot(data,doWeight=T)
dev.off()

library(Vennerable)
data <- Venn(list("Cor_genes"=unique(na.omit(sig_tmp$gene_name)),
  "RNA_TRMA_DN"=unique(na.omit(TRMA_VS_TRM$human_homolog)),
  "ASCL1.chip"=unique(na.omit(ASCL1_binding$human_gene))))
pdf("Fig4A_2.pdf")
plot(data,doWeight=T)
dev.off()

tmp <- subset(sig_tmp, sig_tmp$gene_name %in% ASCL1_downstream$Symbol )
tmp$label <- ifelse(tmp$gene_name=="IGFBP5" | tmp$gene_name=="ASCL1","famous","non_famous")
write.csv(tmp,file="Fig4A.csv")

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
ggsave(ff,file="Fig4A_3.svg",height=10,width=10)




#data <- read.table(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/0_ASCL1_CUTTag/GSE155692_Trudy_Oliver_ASCL1_CUTtag/GSE155691_15688R_raw_counts.txt.gz"), header = TRUE, sep = "\t")  # 读取文件

com_genes <- intersect(sig_tmp$gene_name,ASCL1_downstream$Symbol)
TRMA_VS_TRM <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/TrudyOliver_ASCL1_RNAseq_Chip/PRMA_VS_PRM_TrudyOliver.csv")
TRMA_VS_TRM <- subset(TRMA_VS_TRM, log2FoldChange < -2.5 & pvalue < 0.05 )
exp_genes <- subset(TRMA_VS_TRM, TRMA_VS_TRM $ human_homolog %in% com_genes)
exp_genes <- exp_genes[,c("human_homolog","log2FoldChange")]
rownames(exp_genes) <- exp_genes$human_homolog
exp_genes <- data.frame(exp_genes[,-1])
library(pheatmap)
pdf("FigS4A.pdf",width=3,height=8)
pheatmap(exp_genes,color = colorRampPalette(c("#f7fbff","#deebf7","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#08519c","#08306b"))(50),
  show_rownames=TRUE,cluster_row = FALSE,cluster_col= FALSE,border=FALSE)
dev.off()


IEV.sig <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/IEV_signature.csv")
RNAseq <- mcreadRDS("./Nature2014_and_Zhangpeng_RNA_limma.rds",mc.cores=20)
sel_genes <- intersect(rownames(RNAseq),IEV.sig$Symbol)
RNAseq_1 <- RNAseq[sel_genes,]
SCLC_exp_1 <- data.frame(apply(RNAseq_1,2,mean)) 
names(SCLC_exp_1) <- "IEV.sig"
RNAseq_1 <- data.frame(t(RNAseq))
RNAseq_1$IEV.sig <- SCLC_exp_1$IEV.sig

tmp <- RNAseq_1[,c("ASCL1","IGFBP5","IEV.sig")]
ff1 <- ggplot(tmp, aes(x=ASCL1, y= IGFBP5))+geom_point() + 
stat_smooth(method=lm)+stat_cor(data=tmp, method = "spearman") 
ggsave(ff1,file="FigS4E.svg")

tmp$IGFBP5[tmp$IGFBP5 > 300] <- 300
ff1 <- ggplot(tmp, aes(x=ASCL1, y= IGFBP5))+geom_point() + 
stat_smooth(method=lm)+stat_cor(data=tmp, method = "spearman") 
ggsave(ff1,file="FigS4F.svg")



&&&&&&&&&&&&&&&&&&&&&&&&&&&CancerCell.scRNA&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&CancerCell.scRNA&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&CancerCell.scRNA&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&CancerCell.scRNA&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
AllCells <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/CancerCell_Rudin_All_Combine.rds",mc.cores=20)
Idents(AllCells) <- AllCells$cell_type_fine
Tumor <- subset(AllCells,idents=c("SCLC-N","SCLC-A","SCLC-P","NSCLC"))
pdf("FigS3D.pdf",width=5,height=5)
VlnPlot(Tumor, features = c("IGFBP5"), ncol = 1,pt.size=0,group.by="cell_type_fine")
dev.off()

Tumor <- subset(AllCells,idents=c("SCLC-N","SCLC-A","SCLC-P","NSCLC"))
Idents(Tumor) <- Tumor$histo
pdf("FigS3C.pdf",width=5,height=5)
VlnPlot(Tumor, features = c("IGFBP5"), ncol = 1,pt.size=0,group.by="histo")
dev.off()




















#Fig 4I
#in Protein level
Pure_SCLC_Endo_hi <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/IEV_signature.csv")
normalized_protein <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Endothelium_protein_data/Normalized_Normal_vs_NSCLC_vs_sgScr_sgIGFBP5_all_protein.csv")

common_scRNA_Protein <- intersect(Pure_SCLC_Endo_hi$mouse_gene,normalized_protein$Gene)

common_scRNA_Protein_matrix <- subset(normalized_protein, Gene%in% common_scRNA_Protein)
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
tmp$MS_value[tmp$MS_value > 5000] <- 5000
median_values <- tmp %>%
  group_by(Group) %>%
  summarise(median_value = median(MS_value))
ff <- ggboxplot(tmp, x = "Group", y = "MS_value",add = "", color = "Group",palette = c("#08519c","#e31a1c","#000000","#1c9099"))+ylim(0,8000)+
stat_compare_means(comparisons = my_comparisons, label = "p.format", paired = TRUE, method="wilcox.test") 
ggsave(ff,file="Fig4I_1.png",width=6,height=6)

ff <- ggbarplot(tmp, x = "Group", y = "MS_value", fill="Group",
 add = c("mean_se"), legend = "none", ylim=c(0,2000),title="")+
 rotate_x_text(angle = 45)+scale_fill_manual(values=c("#377eb8","#e41a1c","#1b9e77","#d95f02","#7570b3","#e7298a","#e6ab02","#a6761d","#666666"))+
   geom_line(data = median_values, aes(x = Group, y = median_value, group = 1), color = "gray", size = 1) +
  geom_point(data = median_values, aes(x = Group, y = median_value), color = "gray", size = 3)
ggsave(ff,file="Fig4I_2.png",width=6,height=6)



#Fig 4J
#in scRNAseq level
Pure_SCLC_Endo_signature <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IEV_gene_signature/IEV_signature.csv")
Mice_model <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/NSCLC_VS_SCLC/Final_Mice_only_LUAD_SCLC.rds",mc.cores=20)
Idents(Mice_model) <- Mice_model$sub_CellType
Endo <- subset(Mice_model, idents=c("Endo"),invert=FALSE)
abnormal_vessel_sig <- as.character(Pure_SCLC_Endo_signature$mouse_gene)
abnormal_vessel_sig <- intersect(rownames(GetAssayData(object = Endo)),abnormal_vessel_sig)
speci_raw <- FetchData(object = Endo, vars = abnormal_vessel_sig,slot="data")
Endo[["Abnormal_Vessel"]] <- (rowSums(speci_raw))/length(abnormal_vessel_sig)
Mice_Endo <- Endo[[]]
library(ggpubr)
library(ggplot2)
my_comparisons <- list(c("LUAD","SCLC"))
ff <- ggboxplot(Mice_Endo, x = "Group", y = "Abnormal_Vessel",
               color = "Group",palette = c("#08519c","#e31a1c","#000000","#238b45"))+
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="Fig4I_1.pdf",width=4,height=8)

sgScr_vs_sgIgfbp5 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/sgIgfbp5_vs_sgScr/Renames_AllCells_UMAP20.rds",mc.cores=20)
Idents(sgScr_vs_sgIgfbp5) <- sgScr_vs_sgIgfbp5$RNA_snn_res.0.2
Pure_Endo <- subset(sgScr_vs_sgIgfbp5,idents=c("9"),invert=FALSE)
abnormal_vessel_sig <- as.character(Pure_SCLC_Endo_signature$mouse_gene)
abnormal_vessel_sig <- intersect(rownames(GetAssayData(object = Pure_Endo)),abnormal_vessel_sig)
speci_raw <- FetchData(object = Pure_Endo, vars = abnormal_vessel_sig,slot="data")
Pure_Endo[["Abnormal_Vessel"]] <- (rowSums(speci_raw))/length(abnormal_vessel_sig)
sgIGFBP5_Pure_Endo <- Pure_Endo[[]]
sgIGFBP5_Pure_Endo$Group <- ifelse(sgIGFBP5_Pure_Endo$sample=="Insitu_sgScr","sgScr","sgIgfbp5")
my_comparisons <- list(c("sgScr","sgIgfbp5"))
ff <- ggboxplot(sgIGFBP5_Pure_Endo, x = "Group", y = "Abnormal_Vessel",
               color = "Group",palette = c("#08519c","#e31a1c","#000000","#238b45"))+
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="Fig4I_2.pdf",width=4,height=8)


$$$$$$$$$$$$$$$$$$$$$$$$$$$TU-cohort$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$TU-cohort$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$TU-cohort$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$TU-cohort$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

IGF1_IGF1R_sig <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/IGF1_IGF1R_sig.csv")

Clinical_Group <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/Cell_Patients_statification.csv")
Phosphosite <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/Cell_Phosphosite_Matrix.csv")
Protein <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/Cell_Protein_Matrix.csv")
RNAseq <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/SCLC_human_bulkRNAseq/Cell_RNA_TPM_log2_Matrix.csv")

RNAseq <- as.data.frame(RNAseq)
rownames(RNAseq) <- RNAseq$Gene.Sample.ID
RNAseq <- RNAseq[,c(-1)]
RNAseq <- RNAseq %>% mutate_all(~ ifelse(is.na(.), 0, .))

library(ggplot2)
library(ggpubr)

sel_genes <- intersect(rownames(RNAseq),IGF1_IGF1R_sig$human_gene)
RNAseq_2 <- RNAseq[sel_genes,]
SCLC_exp_1 <- data.frame(apply(RNAseq_2,2,mean)) 
names(SCLC_exp_1) <- "IGF1_IGF1R"
SCLC_exp_1$Sample.ID <- rownames(SCLC_exp_1)

RNAseq_1 <- data.frame(t(RNAseq[c("IGFBP5","ASCL1"),]))
RNAseq_1$Sample.ID <- rownames(RNAseq_3)
SCLC_exp_2 <- RNAseq_1

All_SCLC_tmp <- merge(SCLC_exp_1,SCLC_exp_2,by="Sample.ID")
All_SCLC_tmp <- merge(All_SCLC_tmp,Clinical_Group,by="Sample.ID")
All_SCLC_tmp <- subset(All_SCLC_tmp,Immune.subtype!="NAT-enriched")

my_comparisons <- list(c("Cold-tumor-enriched","Hot-tumor-enriched"))
ff1 <- ggboxplot(All_SCLC_tmp, x = "Immune.subtype", y = "IGF1_IGF1R",add = "",fill = "Immune.subtype") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Fig4O_SCLC_IGF1_IGF1R.svg",height=6,width=4)


my_comparisons <- list(c("Cold-tumor-enriched","Hot-tumor-enriched"))
ff1 <- ggboxplot(All_SCLC_tmp, x = "Immune.subtype", y = "ASCL1",add = "",fill = "Immune.subtype") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Fig4O_SCLC_ASCL1.svg",height=6,width=4)


my_comparisons <- list(c("Cold-tumor-enriched","Hot-tumor-enriched"))
ff1 <- ggboxplot(All_SCLC_tmp, x = "Immune.subtype", y = "IGFBP5",add = "",fill = "Immune.subtype") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Fig4O_SCLC_IGFBP5.svg",height=6,width=4)

ff <- ggplot(All_SCLC_tmp, aes(x=IGFBP5, y= ASCL1 )) +
  geom_point(alpha = 0.5,color = "lightgray") +  # 散点图
  theme_minimal()+ stat_smooth(method=lm)+stat_cor(data=All_SCLC_tmp, method = "spearman")
ggsave(ff,file="Fig4O_ascl1_Igfbp5_spearman_cor.svg")


*******************************磷酸化水平********************************************************
*******************************磷酸化水平********************************************************
*******************************磷酸化水平********************************************************

rownames(Phosphosite) <- Phosphosite$Phosphosite.Sample.ID
Phosphosite <- Phosphosite[,1:113]
Phosphosite <- Phosphosite[,-1]
Phosphosite <- data.frame(t(Phosphosite))

Cold_group <- subset(Clinical_Group,Immune.subtype=="Cold-tumor-enriched")
Hot_group <- subset(Clinical_Group,Immune.subtype=="Hot-tumor-enriched")
Cold_group$Type <- substring(Cold_group$Sample.ID,1,1)
Hot_group$Type <- substring(Hot_group$Sample.ID,1,1)
Cold_group <- subset(Cold_group,Type !="N")
Hot_group <- subset(Hot_group,Type !="N")

Cold_Ph <- Phosphosite[Cold_group$Sample.ID,]
Hot_Ph <- Phosphosite[Hot_group$Sample.ID,]

Cold_Ph <- data.frame(t(Cold_Ph))
Hot_Ph <- data.frame(t(Hot_Ph))


split_data <- strsplit(rownames(Cold_Ph), "\\.")
Phosphosite_split <- data.frame(do.call(rbind, split_data))
Phosphosite_split <- Phosphosite_split[,-3]
names(Phosphosite_split) <- c("Gene","Position")
Cold_Ph_final <- cbind(Cold_Ph, Phosphosite_split)
common <- intersect(IGF1_IGF1R_signature$human_gene,Cold_Ph_final$Gene)
Cold_Ph_IGF1_IR <- subset(Cold_Ph_final, Gene %in% common)
Cold_Ph_IGF1_IR <- Cold_Ph_IGF1_IR %>% mutate_all(~ ifelse(is.na(.), 0, .))
Cold_Ph_IGF1_IR <- Cold_Ph_IGF1_IR[,c(-91,-92)]
Cold_Ph_mean <- data.frame(apply(Cold_Ph_IGF1_IR,1,mean))
names(Cold_Ph_mean) <- "Ph_level"
Cold_Ph_mean$Group <- "Cold"


split_data <- strsplit(rownames(Hot_Ph), "\\.")
Phosphosite_split <- data.frame(do.call(rbind, split_data))
Phosphosite_split <- Phosphosite_split[,-3]
names(Phosphosite_split) <- c("Gene","Position")
Hot_Ph_final <- cbind(Hot_Ph, Phosphosite_split)
common <- intersect(IGF1_IGF1R_signature$human_gene,Hot_Ph_final$Gene)
Hot_Ph_IGF1_IR <- subset(Hot_Ph_final, Gene %in% common)
Hot_Ph_IGF1_IR <- Hot_Ph_IGF1_IR %>% mutate_all(~ ifelse(is.na(.), 0, .))
Hot_Ph_IGF1_IR <- Hot_Ph_IGF1_IR[,c(-17,-18)]
Hot_Ph_mean <- data.frame(apply(Hot_Ph_IGF1_IR,1,mean))
names(Hot_Ph_mean) <- "Ph_level"
Hot_Ph_mean$Group <- "Hot"

Ph_tmp <- rbind(Cold_Ph_mean,Hot_Ph_mean)
Ph_tmp$Ph_ID <- rownames(Ph_tmp)

library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("Cold","Hot"))
ff1 <- ggboxplot(Ph_tmp, x = "Group", y = "Ph_level",add = "",fill = "Group") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Cell_ZhangPeng_IGF1_1R_Ph_level.svg",height=6,width=4)



my_comparisons <- list(c("Cold","Hot"))
ff <- ggpaired(Ph_tmp, x = "Group", y = "Ph_level",
         color = "Group", line.color = "gray", line.size = 0.4,
         palette = "jco")+stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="Cell_ZhangPeng_IGF1_1R_Ph_level.svg",width=3,height=6.4)


ff1 <- ggboxplot(Ph_tmp, x = "Group", y = "Ph_level",add = "",fill = "Group") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Cell_ZhangPeng_IGF1_1R_Ph_level_1.svg",height=6,width=4)





**************************蛋白表达水平*********************************************************
**************************蛋白表达水平*********************************************************
**************************蛋白表达水平*********************************************************
**************************蛋白表达水平*********************************************************

Clinical_Group <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_Patients_statification.csv")
Protein <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_Protein_Matrix.csv")

SCLC_signature <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/PanCancer_Endo_scRNA_level/0_ALL_DATA_INCLUNDE_NOR_AND_PANCANCER_ENDO/Finaaaal_Pure_SCLC_Endo_signature.csv")
IGF1_IGF1R_signature <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/IGF1_1R_signature_correlation/IGF1_IGF1R_signaling_pathway_renew.csv")



library("survival")
library("survminer")

MS <- as.data.frame(Protein)
rownames(MS) <- MS$Protein.Sample.ID
MS <- MS[,c(-1)]
MS <- MS %>% mutate_all(~ ifelse(is.na(.), 0, .))


sel_genes <- intersect(rownames(MS),SCLC_signature$human_gene)
MS_1 <- MS[sel_genes,]
SCLC_exp_1 <- data.frame(apply(MS_1,2,mean)) 
names(SCLC_exp_1) <- "SCLC_sig"
SCLC_exp_1$Sample.ID <- rownames(SCLC_exp_1)

sel_genes <- intersect(rownames(MS),IGF1_IGF1R_signature$human_gene)
MS_2 <- MS[sel_genes,]
SCLC_exp_2 <- data.frame(apply(MS_2,2,mean)) 
names(SCLC_exp_2) <- "IGF1_IGF1R"
SCLC_exp_2$Sample.ID <- rownames(SCLC_exp_2)

MS_3 <- data.frame(t(MS[c("IGFBP5"),]))
MS_3$Sample.ID <- rownames(MS_3)
SCLC_exp_3 <- MS_3


All_SCLC_tmp <- merge(SCLC_exp_1,SCLC_exp_2,by="Sample.ID")
All_SCLC_tmp <- merge(All_SCLC_tmp,SCLC_exp_3,by="Sample.ID")
All_SCLC_tmp <- merge(All_SCLC_tmp,Clinical_Group,by="Sample.ID")
All_SCLC_tmp <- subset(All_SCLC_tmp,Immune.subtype !="NAT-enriched")

library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("Cold-tumor-enriched","Hot-tumor-enriched"))
ff1 <- ggboxplot(All_SCLC_tmp, x = "Immune.subtype", y = "SCLC_sig",add = "",fill = "Immune.subtype") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Protein_Cell_ZhangPeng_SCLC_Endosig.svg",height=6,width=4)


library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("Cold-tumor-enriched","Hot-tumor-enriched"))
ff1 <- ggboxplot(All_SCLC_tmp, x = "Immune.subtype", y = "IGF1_IGF1R",add = "",fill = "Immune.subtype") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Protein_Cell_ZhangPeng_SCLC_IGF1_IGF1R.svg",height=6,width=4)


library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("Cold-tumor-enriched","Hot-tumor-enriched"))
ff1 <- ggboxplot(All_SCLC_tmp, x = "Immune.subtype", y = "IGFBP5",add = "",fill = "Immune.subtype") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#fc8d59","#cb181d"))
ggsave(ff1,file="Protein_Cell_ZhangPeng_SCLC_IGFBP5.svg",height=6,width=4)


ff <- ggplot(All_SCLC_tmp, aes(x=SCLC_sig, y= IGFBP5 )) +
  geom_point(alpha = 0.5,color = "lightgray") +  # 散点图
  theme_minimal()+ stat_smooth(method=lm)+stat_cor(data=All_SCLC_tmp, method = "spearman")
ggsave(ff,file="Protein_Cell_ZhangPeng_SCLCendo_signature_Igfbp5_spearman_cor.svg")

