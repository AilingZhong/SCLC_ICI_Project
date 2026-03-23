#Fig7-S7
library(dplyr)
library(trqwe)
library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(rjson)
library(parallel)
library(future)
library(future.apply)
library(data.table)
library(reticulate)
library(ReductionWrappers)
library(s2a)
library(SCP)
library(ggpubr)
library(sva)
library(limma)
library(GSEABase)
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(GSVA)
library(tidyverse)
library(nichenetr)
library(org.Hs.eg.db)
library(ggplot2)
library(ggpubr)
library(data.table)
library(AnnotationDbi)


normalized_DEseq2 <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/GSE98894_Pancreas_normalized_DEseq2.csv")
rownames(normalized_DEseq2) <- normalized_DEseq2$X
normalized_DEseq2 <- normalized_DEseq2[,-1]
all_tmp <- data.frame(t(normalized_DEseq2))
all_tmp <- log(all_tmp+1,2)
range(all_tmp$ASCL1)
tmp <- all_tmp[,c("ASCL1","IGFBP5")]
data <- tmp
add_mean_difference <- function(data, cutoff_values, gene, outcome) {
  results_list <- list()
  for (cutoff in cutoff_values) {
    data$Group <- ifelse(data[[gene]] > cutoff, "High", "Low")
        if (length(unique(data$Group)) == 2) {
      group_sizes <- table(data$Group)
      if (all(group_sizes >= 2)) {
        test_result <- wilcox.test(data[[outcome]] ~ data$Group)
        high_mean <- mean(data[data$Group == "High", outcome])
        low_mean <- mean(data[data$Group == "Low", outcome])
        diff <- high_mean - low_mean
        results_list[[length(results_list) + 1]] <- c(
          Cutoff = cutoff, 
          P_Value = test_result$p.value, 
          High_Mean = high_mean, 
          Low_Mean = low_mean, 
          Difference = diff
        )
      }
    }
  }
  results <- do.call(rbind, results_list)
  return(as.data.frame(results))
}

cutoff_values <- unique(data$ASCL1)
results_with_means <- add_mean_difference(data, cutoff_values, gene = "ASCL1", outcome = "IGFBP5")
filtered_results <- results_with_means[results_with_means$Difference > 0, ]
filtered_results <- filtered_results[order(filtered_results$P_Value), ]
print(filtered_results)
all_tmp$Group <- ifelse(all_tmp$ASCL1 < 8.631211, "ASCL1_low","ASCL1_hi")
table(all_tmp$Group)
all_tmp$Group <- factor(all_tmp$Group,levels=c("ASCL1_low","ASCL1_hi"))
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("ASCL1_low","ASCL1_hi"))
ff1 <- ggboxplot(all_tmp, x = "Group", y = "IGFBP5",add = "jitter",fill = "Group") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#cb181d"))
ggsave(ff1,file="Fig5D_Pancreas.png",height=6,width=4)





file_path <- "/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/GSE98894_NatureGenetics_PanEndocrine/"
files <- list.files(file_path, pattern = "counts\\.txt\\.gz$", full.names = TRUE)
read_counts <- function(file) {
  data <- read.table(gzfile(file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_name <- gsub("_counts\\.txt\\.gz$", "", basename(file))
  colnames(data)[2] <- sample_name
  return(data)
}

count_list <- lapply(files, read_counts)

gene_names <- count_list[[1]][, 1]
expression_matrix <- do.call(cbind, lapply(count_list, function(df) df[, 2, drop = FALSE]))

merged_counts <- data.frame(Gene = gene_names, expression_matrix, check.names = FALSE)
rownames(merged_counts) <- merged_counts$Gene

merged_counts$SYMBOL <- mapIds(
  x = org.Hs.eg.db,
  keys = as.character(merged_counts$Gene),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

merged_counts <- na.omit(merged_counts)
rownames(merged_counts) <- merged_counts$SYMBOL
expr_mat <- merged_counts[, !(colnames(merged_counts) %in% c("Gene", "SYMBOL"))]
colnames(expr_mat) <- sub("^GSM[0-9]+_", "", colnames(expr_mat))
write.csv(expr_mat, file = "GSE98894_pan_NE_counts.csv")


info_clinical <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/41588_2018_138_MOESM4_ESM.csv")
inte_rect_info <- subset(info_clinical,Origin=="small intestine"| Origin=="rectum")
inte_rect_info$SampleID
inte_rect <- expr_mat[,inte_rect_info$SampleID]
order1 <- data.frame(order=c(1:44),sample="A")
order2 <- data.frame(order=c(45:99),sample="B")
sampleTable <- rbind(order1,order2)
rownames(sampleTable) <- colnames(inte_rect)
library(DESeq2)
dds_1 <- DESeqDataSetFromMatrix(countData = inte_rect,
                              colData = sampleTable,
                              design = ~sample)
DEseq_1 <- DESeq(dds_1)
normalized_DEseq2 <- counts(DEseq_1, normalized=TRUE)
normalized_DEseq2 <- data.frame(normalized_DEseq2)
write.csv(normalized_DEseq2, "GSE98894_Intestine_Rectum_normalized_DEseq2.csv")


normalized_DEseq2 <- read.csv("./GSE98894_Intestine_Rectum_normalized_DEseq2.csv")
rownames(normalized_DEseq2) <- normalized_DEseq2$X
normalized_DEseq2 <- normalized_DEseq2[,-1]
all_tmp <- data.frame(t(normalized_DEseq2))
all_tmp <- log(all_tmp+1,2)
range(all_tmp$ASCL1)
all_tmp$Group <- ifelse(all_tmp$ASCL1 < 0.5, "ASCL1_low","ASCL1_hi")
table(all_tmp$Group)
my_comparisons <- list(c("ASCL1_low","ASCL1_hi"))
ff1 <- ggboxplot(all_tmp, x = "Group", y = "IGFBP5",add = "jitter",fill = "Group") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#cb181d"))
ggsave(ff1,file="Fig6E_SINEC_RENEC.png",height=6,width=4)



nepc_fpkm <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/NatureMedcine_HimishaBeltran_from_cbioportal_nepc_wcm_2016/data_mrna_seq_fpkm.txt",sep="\t",head=TRUE)
nepc_fpkm <- nepc_fpkm[!duplicated(nepc_fpkm$Hugo_Symbol),]
rownames(nepc_fpkm) <- nepc_fpkm$Hugo_Symbol 
nepc_fpkm <- nepc_fpkm[,c(-1,-2)]
nepc_fpkm <- data.frame(t(nepc_fpkm))
sel_genes <- nepc_fpkm[,c("ASCL1","IGFBP5")]
all_tmp <- log(sel_genes+1,2)
range(all_tmp$ASCL1)

all_tmp$Group <- ifelse(all_tmp$ASCL1 < -1.562851, "ASCL1_low","ASCL1_hi")
table(all_tmp$Group)
my_comparisons <- list(c("ASCL1_low","ASCL1_hi"))
ff1 <- ggboxplot(all_tmp, x = "Group", y = "IGFBP5",add = "jitter",fill = "Group") +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")+
scale_fill_manual(values=c("#2171b5","#cb181d"))
ggsave(ff1,file="Fig6F_CRPC.png",height=6,width=4)




tmp <- all_tmp[,c("ASCL1","IGFBP5")]
data <- tmp
add_mean_difference <- function(data, cutoff_values, gene, outcome) {
  results <- data.frame(Cutoff = numeric(), P_Value = numeric(), High_Mean = numeric(), Low_Mean = numeric(), Difference = numeric())
  
  for (cutoff in cutoff_values) {
    data$Group <- ifelse(data[[gene]] > cutoff, "High", "Low")
    
    if (length(unique(data$Group)) == 2) {
      group_sizes <- table(data$Group)
      if (all(group_sizes >= 2)) {
        test_result <- wilcox.test(data[[outcome]] ~ data$Group)
        high_mean <- mean(data[data$Group == "High", outcome])
        low_mean <- mean(data[data$Group == "Low", outcome])
        diff <- high_mean - low_mean
        results <- rbind(results, data.frame(
          Cutoff = cutoff, 
          P_Value = test_result$p.value, 
          High_Mean = high_mean, 
          Low_Mean = low_mean, 
          Difference = diff
        ))
      }
    }
  }
  return(results)
}



file_path <- "/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE178398_pancreatic_NE/FPKM/"

files <- list.files(file_path, pattern = "FPKM\\.txt\\.gz$", full.names = TRUE)

read_fpkm <- function(file) {
  data <- read.table(gzfile(file), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_name <- gsub("_FPKM\\.txt\\.gz$", "", basename(file))
  colnames(data)[2] <- sample_name
  return(data)
}

fpkm_list <- lapply(files, read_fpkm)

gene_names <- fpkm_list[[1]][, 1]
expression_matrix <- do.call(cbind, lapply(fpkm_list, function(df) df[, -1, drop = FALSE]))
merged_fpkm <- data.frame(Gene = gene_names, expression_matrix, check.names = FALSE)

library(data.table)
expr_mat <- as.data.table(merged_fpkm)

merged_fpkm_avg_dt <- expr_mat[, lapply(.SD, mean, na.rm = TRUE), by = Gene]
merged_fpkm_avg <- as.data.frame(merged_fpkm_avg_dt)

rownames(merged_fpkm_avg) <- merged_fpkm_avg$Gene
merged_fpkm_avg <- merged_fpkm_avg[, -1]

samples <- gsub("_FPKM\\.txt\\.gz$", "", basename(files))
colnames(merged_fpkm_avg) <- samples
out_file <- "/mnt/data/user_data/ailing/SCLC_IEV_code/GSE178398_Pancreatic_NE_FPKM.csv"
write.csv(merged_fpkm_avg, file = out_file)









$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$把412个患者分亚型ASCL1_hi和low$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$把412个患者分亚型ASCL1_hi和low$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$把412个患者分亚型ASCL1_hi和low$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$把412个患者分亚型ASCL1_hi和low$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#把消化道的NE合在一起做
normalized_DEseq2 <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/GSE98894_pan_NE_normalized_DEseq2.csv")
rownames(normalized_DEseq2) <- normalized_DEseq2$X
normalized_DEseq2 <- normalized_DEseq2[,-1]
all_tmp <- data.frame(t(normalized_DEseq2))
all_tmp <- log(all_tmp+1,2)
range(all_tmp$ASCL1)
all_tmp$Group <- ifelse(all_tmp$ASCL1 < 0.7776319, "ASCL1_low","ASCL1_hi")
InteRecPancreas <- all_tmp[,c("ASCL1","Group")]
InteRecPancreas <- data.frame(ID=rownames(InteRecPancreas),Group=InteRecPancreas$Group)


#SCLC分亚型
ZhangPeng_RNA <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_RNA_TPM_log2_Matrix.csv")
Clinical <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_Patients_statification.csv")
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
ZP_group <- data.frame(ID=rownames(ZP_SCLC),Group=ZP_SCLC$subtype)

Nature_RNA <- fread("/mnt/data/user_data/xiangyu/workshop/DATABASE/SCLC_Patient_data/SCLC_Nature_RNA.csv")
Nature_RNA <- as.data.frame(Nature_RNA)
Nature_RNA <- Nature_RNA[!duplicated(Nature_RNA$gene),]
rownames(Nature_RNA) <- Nature_RNA$gene
Nature_RNA <- Nature_RNA[,c(-1,-2)]
Nature_RNA_matrix <- data.frame(Nature_RNA)
Nature_SCLC <- scale(Nature_RNA_matrix)
Nature_SCLC <- data.frame(t(Nature_SCLC[c("ASCL1","YAP1","NEUROD1","POU2F3"),]))
Nature_SCLC$subtype <- apply(Nature_SCLC, 1, function(row) {colnames(Nature_SCLC)[which.max(row)]})
Nature_group <- data.frame(ID=rownames(Nature_SCLC),Group=Nature_SCLC$subtype)
SCLC_Group <- rbind(ZP_group,Nature_group)


#PNAS_前列腺癌
exp_data <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/PNAS_CharlesSawyer_Prostate/data_mrna_seq_fpkm_capture_zscores_ref_all_samples.txt",sep="\t",,head=TRUE)
exp_data <- exp_data[!duplicated(exp_data$Hugo_Symbol),]
rownames(exp_data) <- exp_data$Hugo_Symbol
exp_data <- exp_data[,c(-1)]
colnames(exp_data) <- substring(colnames(exp_data),1,7)
clinical <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/PNAS_CharlesSawyer_Prostate/data_clinical_sample.txt",sep="\t",head=TRUE)
clinical$ID <- substring(clinical$SAMPLE_ID,1,7)
NE <- subset(clinical,PATHOLOGY_CLASSIFICATION=="Small cell" | PATHOLOGY_CLASSIFICATION=="Adenocarcinoma with NE features")
NE_ID <- intersect(NE$ID,colnames(exp_data))
NEPC <- data.frame(ID=NE_ID,Group="ASCL1_hi")

Group_ASCL1 <- rbind(InteRecPancreas,SCLC_Group,NEPC)
names(Group_ASCL1) <- c("ID","Type")
Group_ASCL1$Group <- ifelse(Group_ASCL1$Type=="ASCL1" | Group_ASCL1$Type=="ASCL1_hi","ASCL1_hi","ASCL1_low")
write.csv(Group_ASCL1,file="412_NE_ASCL1_SubType.csv")






$$$$$$$$$$$$$$$$$$$$combat_NE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$combat_NE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$combat_NE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$$$$$$combat_NE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

NG_Intestinal_Rec <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/GSE98894_Intestine_Rectum_normalized_DEseq2.csv")
rownames(NG_Intestinal_Rec) <- NG_Intestinal_Rec$X
NG_Intestinal_Rec <- NG_Intestinal_Rec[,-1]
NG_Intestinal_Rec <- data.frame(NG_Intestinal_Rec)

NG_Pancreas <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/GSE98894_Pancreas_normalized_DEseq2.csv")
rownames(NG_Pancreas) <- NG_Pancreas$X
NG_Pancreas <- NG_Pancreas[,-1]
NG_Pancreas <- data.frame(NG_Pancreas)


SCLC <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/PanCancer_Endo_scRNA_level/TCGA_Deconvolution/Deconvolution/SCLC/Nature2014_bulk_SCLC_matrix.rds",mc.cores=20)
Nature2014_SCLC <- data.frame(t(SCLC))

Clinical <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_Patients_statification.csv")
ZhangPeng_RNA <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_RNA_TPM_log2_Matrix.csv")
ZhangPeng_RNA <- as.data.frame(ZhangPeng_RNA)
rownames(ZhangPeng_RNA) <- ZhangPeng_RNA$Gene.Sample.ID
ZhangPeng_RNA <- ZhangPeng_RNA[,c(-1)]
ZhangPeng_RNA <- ZhangPeng_RNA %>% mutate_all(~ ifelse(is.na(.), 0, .))
ZhangPeng_RNA <- data.frame(t(ZhangPeng_RNA))
Tumor_ID <- subset(Clinical,Immune.subtype !="NAT-enriched")
ZhangPeng_RNA_1 <- ZhangPeng_RNA[Tumor_ID$Sample.ID,]
ZhangPeng_SCLC <- data.frame(t(ZhangPeng_RNA_1))


#PNAS_ProstateCancer
exp_data <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/PNAS_CharlesSawyer_Prostate/data_mrna_seq_fpkm_capture.txt",sep="\t",,head=TRUE)
exp_data <- exp_data[!duplicated(exp_data$Hugo_Symbol),]
rownames(exp_data) <- exp_data$Hugo_Symbol
exp_data <- exp_data[,c(-1)]
colnames(exp_data) <- substring(colnames(exp_data),1,7)
exp_data <- data.frame(t(exp_data))
clinical <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/PNAS_CharlesSawyer_Prostate/data_clinical_sample.txt",sep="\t",head=TRUE)
clinical$ID <- substring(clinical$SAMPLE_ID,1,7)
NE_ID <- subset(clinical,PATHOLOGY_CLASSIFICATION=="Small cell" | PATHOLOGY_CLASSIFICATION=="Adenocarcinoma with NE features")
NE_ID <- intersect(NE_ID$ID,rownames(exp_data))
PNAS_NE <- exp_data[NE_ID,]
PNAS_NEPC <- data.frame(t(PNAS_NE))


common_genes <- Reduce(intersect, list(rownames(NG_Intestinal_Rec),rownames(NG_Pancreas),rownames(Nature2014_SCLC),rownames(ZhangPeng_SCLC),rownames(PNAS_NEPC)))

Nature2014_SCLC_1 <- Nature2014_SCLC[common_genes,]
ZhangPeng_RNA_matrix_1 <- ZhangPeng_SCLC[common_genes,]
PNAS_NEPC_1 <- PNAS_NEPC[common_genes,]
NG_Intestinal_Rec_1 <- NG_Intestinal_Rec[common_genes,]
NG_Pancreas_1 <- NG_Pancreas[common_genes,]
All_NE_Matrix <- cbind(Nature2014_SCLC_1,ZhangPeng_RNA_matrix_1,PNAS_NEPC_1,NG_Intestinal_Rec_1,NG_Pancreas_1)
mcsaveRDS(All_NE_Matrix,file="All_NE_412patients_matrix_without_NM_nepc.rds",mc.cores=20)


library(edgeR)
dge <- DGEList(counts = All_NE_Matrix)
dge <- calcNormFactors(dge, method = "TMM")
normalized_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
library(sva)
batch <- factor(c(rep("Nature_SCLC", ncol(Nature2014_SCLC_1)),
rep("CancerCell_SCLC", ncol(ZhangPeng_RNA_matrix_1)),
rep("PNAS_NEPC", ncol(PNAS_NEPC_1)),
rep("NG_InteRec", ncol(NG_Intestinal_Rec_1)),
rep("NG_PanNE", ncol(NG_Pancreas_1))))
combat_matrix <- ComBat(dat = normalized_matrix, batch = batch, mod = NULL)
mcsaveRDS(combat_matrix,file="All_NE_412patients_combat_without_NM_nepc.rds",mc.cores=20)

matrix_1 <- as.matrix(combat_matrix)
c5_geneSets <- getGmt("/mnt/data/user_data/xiangyu/workshop/DATABASE/GSVA_7.1/c5.all.v7.1.symbols.gmt")
c5_GSVA_res <- gsva(matrix_1, c5_geneSets, min.sz=5, max.sz=500, verbose=FALSE, parallel.sz=10)
saveRDS(c5_GSVA_res,file="All_NE_412patients_combat_without_NM_nepc_GSVA_sudobulk_c5_res.rds")










ZP_SCLC_fraction <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/0_BBB_zen_me_lai_de_correction/113_ZhangPeng_Cell_bulkRNA_deconv_SCLC_fraction.csv")
Nature_SCLC_fraction <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/PanCancer_Endo_scRNA_level/TCGA_Deconvolution/Deconvolution/SCLC/CancerCell_deconv_SCLC_fraction.csv")
PanNEC_fraction <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/Deconvolution/NatureGenetics_deconv_Pancreas_fraction.csv")
SmallInte_fraction <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/Deconvolution/NatureGenetics_deconv_SINE_REC_fraction.csv")

Prostate_fraction <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/Deconvolution/PNAS_Prostate_deconv_Pancreas_fraction.csv")
Prostate_fraction$ID <- substring(Prostate_fraction$X,1,7)
clinical <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/PNAS_CharlesSawyer_Prostate/data_clinical_sample.txt",sep="\t",head=TRUE)
clinical$ID <- substring(clinical$SAMPLE_ID,1,7)
NE <- subset(clinical,PATHOLOGY_CLASSIFICATION=="Small cell" | PATHOLOGY_CLASSIFICATION=="Adenocarcinoma with NE features")
NE$ID <- substring(NE$SAMPLE_ID,1,7)
NE_ID <- intersect(NE$ID,Prostate_fraction$ID)
AD_ID <- setdiff(Prostate_fraction$ID,NE$ID)

NEPC_fraction <- subset(Prostate_fraction,  Prostate_fraction$ID %in% NE_ID  )
PRPC_fraction <- subset(Prostate_fraction, Prostate_fraction$ID %in% AD_ID  )


CC_SCLC_fraction <- ZP_SCLC_fraction[,c("X","T.cell")]
Na_SCLC_fraction <- Nature_SCLC_fraction[,c("X","T.cell")]
PanNEC_fraction <- PanNEC_fraction[,c("X","T.cell")]
SINE_fraction <- SmallInte_fraction[,c("X","T.cell")]
NEPC_fraction <- NEPC_fraction[,c("ID","Tcell")]

names(CC_SCLC_fraction) <- c("ID","T.cell")
names(Na_SCLC_fraction) <- c("ID","T.cell")
names(PanNEC_fraction) <-c("ID","T.cell") 
names(SINE_fraction) <-c("ID","T.cell") 
names(NEPC_fraction) <-c("ID","T.cell") 

CC_SCLC_fraction$TumorType <- "SCLC"
Na_SCLC_fraction$TumorType <- "SCLC"
PanNEC_fraction$TumorType <-  "PanNE"
SINE_fraction$TumorType <-  "SINE"
NEPC_fraction$TumorType <-  "NEPC"

CC_SCLC_fraction$Group <- "NE"
Na_SCLC_fraction$Group <- "NE"
PanNEC_fraction$Group <-  "NE"
SINE_fraction$Group <-  "NE"
NEPC_fraction$Group <-  "NE"

All_Tcell_ratio <- rbind(CC_SCLC_fraction,Na_SCLC_fraction,PanNEC_fraction,SINE_fraction,NEPC_fraction)
mcsaveRDS(All_Tcell_ratio,file="final_412_NE_cohorts_Tcell_ratio.rds")





ASCL1_Group <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/0_BBB_zen_me_lai_de_correction/412_NE_ASCL1_SubType.csv")

Tcell_ratio <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/0_BBB_zen_me_lai_de_correction/final_412_NE_cohorts_Tcell_ratio.rds",mc.cores=20)
combat_matrix <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Vvvvvery_important_BBB_signature_define/BBB_in_panNE/4cohorts_All_NE_combat_without_NM_nepc.rds",mc.cores=20)

scale_exp <- scale(combat_matrix)
ASCL1_IGFBP5 <- data.frame(t(scale_exp[c("ASCL1","IGFBP5"),]))
ASCL1_IGFBP5$ID <- rownames(ASCL1_IGFBP5)

all_tmp <- merge(ASCL1_IGFBP5,Tcell_ratio,by="ID")

write.csv(all_tmp,file="Fig5g_ASCL1_associate_IGFBP5.csv")


final_matrix$TumorType <- factor(final_matrix$TumorType,levels=c("PanNE",'SINE',"NEPC","SCLC"))
final_matrix$TumorType <- as.factor(final_matrix$TumorType)
library(ggpubr)
library(ggplot2)
library(BuenColors)
final_matrix <- final_matrix[order(final_matrix$ASCL1, decreasing = F),]
final_matrix$order <- c(1:412)


final_matrix$Group.y <- as.factor(final_matrix$Group.y)
library(ggpubr)
library(ggplot2)
library(BuenColors)
final_matrix$Group.y <- as.factor(final_matrix$Group.y)
tumor_levels <- levels(final_matrix$Group.y)
y_positions <- seq(from = min(final_matrix$IGFBP5) - 0.1, 
                   by = -0.05, length.out = length(tumor_levels))
tumor_y_mapping <- data.frame(Group.y = tumor_levels, y_pos = y_positions)
final_matrix <- merge(final_matrix, tumor_y_mapping, by = "Group.y")
p1 <- ggplot(final_matrix, aes(order, IGFBP5)) +
  geom_point(alpha = 0.5, size = 0.5, colour = "#969696") +
  scale_color_gradientn(colours = c("#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#d73027")) +
  geom_rug(alpha = 0.8, position = "jitter", sides = "b") +
  geom_smooth(colour = "orange", se = TRUE) +
  theme_classic() +
  xlab("Tcell_low_to_high_order") +
  labs(title = "IGFBP5") + 
  NoLegend() + scale_y_continuous(limits = c(-1, 15), breaks = seq(-1, 15, by = 0.5))
p2 <- p1 + geom_segment(data = final_matrix, aes(x = order, xend = order, y = y_pos, yend = y_pos + 0.02, color = Group.y),inherit.aes = FALSE) + 
  scale_color_manual(values = c("ASCL1_hi" = "#a50f15", "ASCL1_low" ="#000000")) + 
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(), 
        legend.position = "bottom") +  
  geom_text(data = tumor_y_mapping, aes(x = max(final_matrix$order) + 10,  y = y_pos + 0.01, label = Group.y, color = Group.y),inherit.aes = FALSE, hjust = 0, size = 3)
ggsave(p2,file="ASCL1_associate_IGFBP5.svg",width=8,height=8)



&&&&&&&&&&&&&&&&&&&&&&&&&&&&&把GEP_NET拆分开来分ASCL1的组别&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&把GEP_NET拆分开来分ASCL1的组别&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&把GEP_NET拆分开来分ASCL1的组别&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&把GEP_NET拆分开来分ASCL1的组别&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


normalized_DEseq2 <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/GSE98894_Pancreas_normalized_DEseq2.csv")
rownames(normalized_DEseq2) <- normalized_DEseq2$X
normalized_DEseq2 <- normalized_DEseq2[,-1]
PNEC <- data.frame(t(normalized_DEseq2))
PNEC <- log(PNEC+1,2)
range(PNEC$ASCL1)
PNEC_Group <- data.frame(PNEC[,c("ASCL1")])
names(PNEC_Group) <- "ASCL1"
rownames(PNEC_Group) <- rownames(PNEC)
PNEC_Group$Group <- ifelse(PNEC_Group$ASCL1 < 8.631211, "ASCL1_low","ASCL1_hi")
table(PNEC_Group$Group)
PNEC_Group$Group <- factor(PNEC_Group$Group,levels=c("ASCL1_low","ASCL1_hi"))
PNEC_Group$ID <- rownames(PNEC_Group)
PNEC_Group <- PNEC_Group[,c("ID","Group")]


normalized_DEseq2 <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/GSE98894_Pan_cancer_NE_Califo/GSE98894_Intestine_Rectum_normalized_DEseq2.csv")
rownames(normalized_DEseq2) <- normalized_DEseq2$X
normalized_DEseq2 <- normalized_DEseq2[,-1]
SINE <- data.frame(t(normalized_DEseq2))
SINE <- log(SINE+1,2)
range(SINE$ASCL1)
SINE_Group <- data.frame(SINE[,c("ASCL1","IGFBP5")])
names(SINE_Group) <- c("ASCL1","IGFBP5")
rownames(SINE_Group) <- rownames(SINE)
SINE_Group$Group <- ifelse(SINE_Group$ASCL1 < 0.5, "ASCL1_low","ASCL1_hi")
table(SINE_Group$Group)
SINE_Group$Group <- factor(SINE_Group$Group,levels=c("ASCL1_low","ASCL1_hi"))
SINE_Group$ID <- rownames(SINE_Group)

tmp <- SINE_Group[,c("IGFBP5","ASCL1","Group")]
write.csv(tmp,file="Fig5d_NG_GSE98894_SINE_IGFBP5.csv")

SINE_Group <- SINE_Group[,c("ID","Group")]




#SCLC-subtype
ZhangPeng_RNA <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_RNA_TPM_log2_Matrix.csv")
Clinical <- read.csv("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Reference_Data/human_SCLC/Cell_ZhangPeng/Abao_Transfromed_matrix/Cell_ZhangPeng_Patients_statification.csv")
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
ZP_group <- data.frame(ID=rownames(ZP_SCLC),Group=ZP_SCLC$subtype)

Nature_RNA <- fread("/mnt/data/user_data/xiangyu/workshop/DATABASE/SCLC_Patient_data/SCLC_Nature_RNA.csv")
Nature_RNA <- as.data.frame(Nature_RNA)
Nature_RNA <- Nature_RNA[!duplicated(Nature_RNA$gene),]
rownames(Nature_RNA) <- Nature_RNA$gene
Nature_RNA <- Nature_RNA[,c(-1,-2)]
Nature_RNA_matrix <- data.frame(Nature_RNA)
Nature_SCLC <- scale(Nature_RNA_matrix)
Nature_SCLC <- data.frame(t(Nature_SCLC[c("ASCL1","YAP1","NEUROD1","POU2F3"),]))
Nature_SCLC$subtype <- apply(Nature_SCLC, 1, function(row) {colnames(Nature_SCLC)[which.max(row)]})
Nature_group <- data.frame(ID=rownames(Nature_SCLC),Group=Nature_SCLC$subtype)
SCLC_Group <- rbind(ZP_group,Nature_group)


#PNAS_ProstateCancer
exp_data <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/PNAS_CharlesSawyer_Prostate/data_mrna_seq_fpkm_capture_zscores_ref_all_samples.txt",sep="\t",,head=TRUE)
exp_data <- exp_data[!duplicated(exp_data$Hugo_Symbol),]
rownames(exp_data) <- exp_data$Hugo_Symbol
exp_data <- exp_data[,c(-1)]
colnames(exp_data) <- substring(colnames(exp_data),1,7)
clinical <- read.table("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/Neuronendocrine_Cancer/Take_NE_Cohort/bulkRNA/PNAS_CharlesSawyer_Prostate/data_clinical_sample.txt",sep="\t",head=TRUE)
clinical$ID <- substring(clinical$SAMPLE_ID,1,7)
NE <- subset(clinical,PATHOLOGY_CLASSIFICATION=="Small cell" | PATHOLOGY_CLASSIFICATION=="Adenocarcinoma with NE features")
NE_ID <- intersect(NE$ID,colnames(exp_data))
NEPC <- data.frame(ID=NE_ID,Group="ASCL1_hi")

Group_ASCL1 <- rbind(PNEC_Group,SINE_Group, SCLC_Group,NEPC)
names(Group_ASCL1) <- c("ID","Type")
Group_ASCL1$Group <- ifelse(Group_ASCL1$Type=="ASCL1" | Group_ASCL1$Type=="ASCL1_hi","ASCL1_hi","ASCL1_low")
write.csv(Group_ASCL1,file="split_PNEC_SINE_412_NE_ASCL1_SubType.csv")




GSVA <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/PanNE/All_NE_412patients_combat_without_NM_nepc_GSVA_sudobulk_c5_res.rds",mc.cores=20)
ASCL1_Group <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/PanNE/split_PNEC_SINE_412_NE_ASCL1_SubType.csv")
Tcell_ratio <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/PanNE/final_412_NE_cohorts_Tcell_ratio.rds",mc.cores=20)
GSVA <- data.frame(t(GSVA))
GSVA <- GSVA[,c("GO_ESTABLISHMENT_OF_BLOOD_BRAIN_BARRIER","GO_BASEMENT_MEMBRANE_ASSEMBLY")]
GSVA$ID <- rownames(GSVA)
all_tmp <- merge(GSVA,ASCL1_Group,by="ID")
all_tmp <- merge(all_tmp,Tcell_ratio,by="ID")
all_tmp <- all_tmp[order(all_tmp$GO_ESTABLISHMENT_OF_BLOOD_BRAIN_BARRIER, decreasing = F),]
all_tmp$order <- c(1:412)
write.csv(all_tmp,file="Fig5J.csv")

aa_tmp <- all_tmp[,c("T.cell","TumorType","GO_ESTABLISHMENT_OF_BLOOD_BRAIN_BARRIER","order","Group.x")]
aa_tmp$TumorType <- factor(aa_tmp$TumorType,levels=c("PanNE",'SINE',"NEPC","SCLC"))
aa_tmp$Group.x <- as.factor(aa_tmp$Group.x)
tumor_levels <- levels(aa_tmp$Group.x)

group_y_mapping <- data.frame(Group.x = c("ASCL1_low", "ASCL1_hi"),y_pos = c(0.0012, 0.0038))
aa_tmp <- merge(aa_tmp, group_y_mapping, by = "Group.x")
aa_tmp <- aa_tmp[order(aa_tmp$GO_ESTABLISHMENT_OF_BLOOD_BRAIN_BARRIER, decreasing = FALSE), ]
aa_tmp$order <- 1:nrow(aa_tmp)

p2 <- p1 +
  geom_segment(
    data = aa_tmp,
    aes(x = order, xend = order, y = y_pos, yend = y_pos + 0.0012, color = Group.x),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = group_y_mapping,
    aes(
      x = max(aa_tmp$order) + 10,
      y = y_pos + 0.0006,
      label = Group.x,
      color = Group.x
    ),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +
  scale_color_manual(values = c("ASCL1_low" = "#000000", "ASCL1_hi" = "#1f78b4")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(0, 0.03))
ggsave(plot = p2,file = "Fig6J_1.svg",width = 6,height = 6)



all_tmp <- all_tmp[order(all_tmp$GO_ESTABLISHMENT_OF_BLOOD_BRAIN_BARRIER, decreasing = F),]
all_tmp$order <- c(1:412)
aa_tmp <- all_tmp[,c("T.cell","TumorType","GO_ESTABLISHMENT_OF_BLOOD_BRAIN_BARRIER","order","Group.x")]
aa_tmp$TumorType <- factor(aa_tmp$TumorType,levels=c("PanNE",'SINE',"NEPC","SCLC"))
aa_tmp$TumorType <- as.factor(aa_tmp$TumorType)
tumor_levels <- levels(aa_tmp$TumorType)
y_positions <- seq(from = min(aa_tmp$T.cell) - 0.1, 
                   by = -0.05, length.out = length(tumor_levels))
tumor_y_mapping <- data.frame(TumorType = tumor_levels, y_pos = y_positions)
aa_tmp <- merge(aa_tmp, tumor_y_mapping, by = "TumorType")
aa_tmp <- aa_tmp[order(aa_tmp$GO_ESTABLISHMENT_OF_BLOOD_BRAIN_BARRIER, decreasing = F),]
aa_tmp$order <- c(1:412)
p1 <- ggplot(aa_tmp, aes(order, T.cell)) +
  geom_point(alpha = 0.5, size = 0.5, colour = "#969696") +
  scale_color_gradientn(colours = c("#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#d73027")) +
  geom_rug(alpha = 0.8, position = "jitter", sides = "b") +
  geom_smooth(colour = "orange", se = TRUE) +
  theme_classic() +
  xlab("Tcell_hi_to_low_order") +
  labs(title = "T.cell") +
  NoLegend() 
p2 <- p1 + 
  geom_segment(
    data = aa_tmp,
    aes(x = order, xend = order, y = y_pos, yend = y_pos + 0.02, color = TumorType),
    inherit.aes = FALSE
  ) + 
  scale_color_manual(values = c("SCLC" = "#1f78b4", "NEPC" = "#000000", "SINE" = "#e41a1c", "PanNE" = "#8c6bb1")) + 
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  ) +  
  geom_text(
    data = tumor_y_mapping,
    aes(x = max(aa_tmp$order) + 10, y = y_pos + 0.01, label = TumorType, color = TumorType),
    inherit.aes = FALSE,
    hjust = 0,
    size = 3
  ) +
  coord_cartesian(ylim = c(0, 0.03))
ggsave(p2, file = "Fig6J_2.svg", width = 6, height = 6)












***********************CRPC_scRNAseq******************************************
***********************CRPC_scRNAseq******************************************
***********************CRPC_scRNAseq******************************************
***********************CRPC_scRNAseq******************************************
***********************CRPC_scRNAseq******************************************

library(harmony)
library(reticulate)
library(ReductionWrappers)
library(s2a)
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


Adeno_1 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428953_1779_HMP05_IGO_10726_3_dense.csv.gz"))
Adeno_2 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428954_1845_HMP-08_IGO_10837_19_dense.csv.gz"))
Adeno_3 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428955_1968_HMP11_1_IGO_11247_3_dense.csv.gz"))
Adeno_4 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428956_1969_HMP11_2_IGO_11247_4_dense.csv.gz"))
Adeno_5 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428957_1846_JZHP_3_IGO_10837_21_dense.csv.gz"))
Adeno_6 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428958_2284_HMP_14B_IGO_11588_19_dense.csv.gz"))
Adeno_7 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428962_2518_HMP_19_IGO_11874_45_dense.csv.gz"))
Adeno_8 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428963_2624_HMP20_IGO_12065_37_dense.csv.gz"))
Adeno_9 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428964_3226_SZ-1327_HMP25_IGO_12437_254_dense.csv.gz"))
Adeno_10 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Adenocarcinoma/GSM6428965_3309_SZ-1357_HMP26_IGO_12437_I_19_dense.csv.gz"))


NE_1 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/NEPC/GSM6428952_1778_JZ_HMP_04_IGO_10726_2_dense.csv.gz"))
NE_2 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/NEPC/GSM6428960_2016_JZHP_04_IGO_11245_7_dense.csv.gz"))
NE_3 <- read.csv(gzfile("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/NEPC/GSM6428961_2513_HMP_17_IGO_11874_35_dense.csv.gz"))

Adeno_1$X <-  paste0("Adeno_",Adeno_1$X,sep="")
rownames(Adeno_1) <- as.character(Adeno_1$X)
Adeno_1 <- Adeno_1[,c(-1,-2)]
Adeno_1 <- data.frame(t(Adeno_1))


Adeno_2$X <-  paste0("Adeno_",Adeno_2$X,sep="")
rownames(Adeno_2) <- as.character(Adeno_2$X)
Adeno_2 <- Adeno_2[,c(-1,-2)]
Adeno_2 <- data.frame(t(Adeno_2))

Adeno_3$X <-  paste0("Adeno_",Adeno_3$X,sep="")
rownames(Adeno_3) <- as.character(Adeno_3$X)
Adeno_3 <- Adeno_3[,c(-1,-2)]
Adeno_3 <- data.frame(t(Adeno_3))

Adeno_4$X <-  paste0("Adeno_",Adeno_4$X,sep="")
rownames(Adeno_4) <- as.character(Adeno_4$X)
Adeno_4 <- Adeno_4[,c(-1,-2)]
Adeno_4 <- data.frame(t(Adeno_4))

Adeno_5$X <-  paste0("Adeno_",Adeno_5$X,sep="")
rownames(Adeno_5) <- as.character(Adeno_5$X)
Adeno_5 <- Adeno_5[,c(-1,-2)]
Adeno_5 <- data.frame(t(Adeno_5))

Adeno_6$X <-  paste0("Adeno_",Adeno_6$X,sep="")
rownames(Adeno_6) <- as.character(Adeno_6$X)
Adeno_6 <- Adeno_6[,c(-1,-2)]
Adeno_6 <- data.frame(t(Adeno_6))

Adeno_7$X <-  paste0("Adeno_",Adeno_7$X,sep="")
rownames(Adeno_7) <- as.character(Adeno_7$X)
Adeno_7 <- Adeno_7[,c(-1,-2)]
Adeno_7 <- data.frame(t(Adeno_7))

Adeno_8$X <-  paste0("Adeno_",Adeno_8$X,sep="")
rownames(Adeno_8) <- as.character(Adeno_8$X)
Adeno_8 <- Adeno_8[,c(-1,-2)]
Adeno_8 <- data.frame(t(Adeno_8))

Adeno_9$X <-  paste0("Adeno_",Adeno_9$X,sep="")
rownames(Adeno_9) <- as.character(Adeno_9$X)
Adeno_9 <- Adeno_9[,c(-1,-2)]
Adeno_9 <- data.frame(t(Adeno_9))

Adeno_10$X <-  paste0("Adeno_",Adeno_10$X,sep="")
rownames(Adeno_10) <- as.character(Adeno_10$X)
Adeno_10 <- Adeno_10[,c(-1,-2)]
Adeno_10 <- data.frame(t(Adeno_10))


NE_1$X <-  paste0("NE_",NE_1$X,sep="")
rownames(NE_1) <- as.character(NE_1$X)
NE_1 <- NE_1[,c(-1,-2)]
NE_1 <- data.frame(t(NE_1))

NE_2$X <-  paste0("NE_",NE_2$X,sep="")
rownames(NE_2) <- as.character(NE_2$X)
NE_2 <- NE_2[,c(-1,-2)]
NE_2 <- data.frame(t(NE_2))

NE_3$X <-  paste0("NE_",NE_3$X,sep="")
rownames(NE_3) <- as.character(NE_3$X)
NE_3 <- NE_3[,c(-1,-2)]
NE_3 <- data.frame(t(NE_3))




Adeno_1_object <- CreateSeuratObject(counts =Adeno_1, project = "Adeno_1", min.cells = 3, min.features = 200)
Adeno_1_object$sample <- Idents(object = Adeno_1_object)
Adeno_1_object$Group <- "Adeno"

Adeno_2_object <- CreateSeuratObject(counts =Adeno_2, project = "Adeno_2", min.cells = 3, min.features = 200)
Adeno_2_object$sample <- Idents(object = Adeno_2_object)
Adeno_2_object$Group <- "Adeno"

Adeno_3_object <- CreateSeuratObject(counts =Adeno_3, project = "Adeno_3", min.cells = 3, min.features = 200)
Adeno_3_object$sample <- Idents(object = Adeno_3_object)
Adeno_3_object$Group <- "Adeno"

Adeno_4_object <- CreateSeuratObject(counts =Adeno_4, project = "Adeno_4", min.cells = 3, min.features = 200)
Adeno_4_object$sample <- Idents(object = Adeno_4_object)
Adeno_4_object$Group <- "Adeno"

Adeno_5_object <- CreateSeuratObject(counts =Adeno_5, project = "Adeno_5", min.cells = 3, min.features = 200)
Adeno_5_object$sample <- Idents(object = Adeno_5_object)
Adeno_5_object$Group <- "Adeno"

Adeno_6_object <- CreateSeuratObject(counts =Adeno_6, project = "Adeno_6", min.cells = 3, min.features = 200)
Adeno_6_object$sample <- Idents(object = Adeno_6_object)
Adeno_6_object$Group <- "Adeno"

Adeno_7_object <- CreateSeuratObject(counts =Adeno_7, project = "Adeno_7", min.cells = 3, min.features = 200)
Adeno_7_object$sample <- Idents(object = Adeno_7_object)
Adeno_7_object$Group <- "Adeno"

Adeno_8_object <- CreateSeuratObject(counts =Adeno_8, project = "Adeno_8", min.cells = 3, min.features = 200)
Adeno_8_object$sample <- Idents(object = Adeno_8_object)
Adeno_8_object$Group <- "Adeno"

Adeno_9_object <- CreateSeuratObject(counts =Adeno_9, project = "Adeno_9", min.cells = 3, min.features = 200)
Adeno_9_object$sample <- Idents(object = Adeno_9_object)
Adeno_9_object$Group <- "Adeno"

Adeno_10_object <- CreateSeuratObject(counts =Adeno_10, project = "Adeno_10", min.cells = 3, min.features = 200)
Adeno_10_object$sample <- Idents(object = Adeno_10_object)
Adeno_10_object$Group <- "Adeno"



NE_1_object <- CreateSeuratObject(counts =NE_1, project = "NE_1", min.cells = 3, min.features = 200)
NE_1_object$sample <- Idents(object = NE_1_object)
NE_1_object$Group <- "NE"

NE_2_object <- CreateSeuratObject(counts =NE_2, project = "NE_2", min.cells = 3, min.features = 200)
NE_2_object$sample <- Idents(object = NE_2_object)
NE_2_object$Group <- "NE"

NE_3_object <- CreateSeuratObject(counts =NE_3, project = "NE_3", min.cells = 3, min.features = 200)
NE_3_object$sample <- Idents(object = NE_3_object)
NE_3_object$Group <- "NE"


merge.data <- merge(x = Adeno_1_object, y=c(Adeno_2_object,Adeno_3_object,Adeno_4_object,Adeno_5_object,
	Adeno_6_object,Adeno_7_object,Adeno_8_object,Adeno_9_object,Adeno_10_object,NE_1_object,NE_2_object,NE_3_object))
All_obj <- merge.data
All_obj[["percent.MT"]] <- PercentageFeatureSet(All_obj, pattern = "^MT-")



All_MERGE_DATA <- All_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
All_MERGE_DATA.pca <- RunPCA(All_MERGE_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(All_MERGE_DATA.pca,file="All_CRPC_MERGE_DATA_pca.rds",mc.cores=20)
library(harmony)
library(reticulate)
library(ReductionWrappers)
library(s2a)
ALL_MERGE_DATA_harmony  <- RunHarmony(All_MERGE_DATA.pca,"orig.ident")
FindNeighbors <- FindNeighbors(ALL_MERGE_DATA_harmony,reduction = "harmony", dims = 1:20)
FindClusters <- FindClusters(FindNeighbors,resolution = c(0.1,0.2))
mcsaveRDS(FindClusters,file="All_CRPC_Harmony_FindClusters.rds", mc.cores = 20)

Harmony_UMAP30 <- RunUMAP(FindClusters,reduction = "harmony", dims = 1:20)
mcsaveRDS(Harmony_UMAP30,file="NoFilter_All_CRPC_Harmony_UMAP30.rds",mc.cores=20)


Idents(Harmony_UMAP30) <- Harmony_UMAP30$RNA_snn_res.0.2
all.markers <- FindAllMarkers(Harmony_UMAP30, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.markers,"All_CRPC_UMAP20_res02_FindAllmarkers.csv")



p2 <- DimPlot(object = Harmony_UMAP30, reduction = "umap",label=TRUE,
    group.by="RNA_snn_res.0.2") +NoLegend()+labs(title="UMAP")
ggsave(p2,file="All_CRPC_Harmony_UMAP30.png", width=8,height=8)


p2 <- DimPlot(object = Harmony_UMAP15, reduction = "umap",label=TRUE,
    group.by="RNA_snn_res.0.2",split.by="Group") +NoLegend()+labs(title="UMAP")


Idents(Harmony_UMAP30) <- Harmony_UMAP30$RNA_snn_res.0.2
Endo <- subset(Harmony_UMAP30, idents="7", invert=FALSE)




Idents(Harmony_UMAP30) <- Harmony_UMAP30$RNA_snn_res.0.2
new.cluster.ids <- c("Macro","Fibro","Tcell","Epi","Epi","Epi","Epi","Endo","Epi","Epi","Epi","Epi","Epi","Epi","Epi","Tcell","Epi","Tcell","Epi","Fibro","Epi","Bcell","Epi","Macro")
names(new.cluster.ids) <- levels(Harmony_UMAP30)
Harmony_UMAP30 <- RenameIdents(Harmony_UMAP30, new.cluster.ids)
Harmony_UMAP30$CellType <- Idents(Harmony_UMAP30)
Idents(Harmony_UMAP30) <- Harmony_UMAP30$CellType
Harmony_UMAP30$CellType <- factor(Harmony_UMAP30$CellType,levels=c("Epi","Fibro","Endo","Macro","Tcell","Bcell"))

mcsaveRDS(Harmony_UMAP30,file="Renames_CRPC_UMAP30.rds",mc.cores=20)


ff <- DimPlot(object = Harmony_UMAP30, reduction = "umap",pt.size=0.5,
    cols=c("#d73027","#8c510a","#8073ac","#4575b4","#de77ae","#5aae61","#92c5de","#4393c3"),label=TRUE, label.size=5)
ggsave("Renames_CellType.png",width=7.2,height=6.7,dpi=1080)




All_Cells <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Renames_CRPC_UMAP30.rds",mc.cores=20)
meta.data <- All_Cells@meta.data
Adeno <- subset(meta.data,Group=="Adeno")
NE <- subset(meta.data,Group=="NE")
Adeno <- data.frame(table(Adeno$CellType))
NE <- data.frame(table(NE$CellType))
Adeno$group <- "Adeno"
NE$group <- "NE"
Adeno$Ratio <- Adeno$Freq/sum(Adeno$Freq)
NE$Ratio <- NE$Freq/sum(NE$Freq)
tmp1_data <- rbind(Adeno,NE)
tmp1_data$group <- factor(tmp1_data$group,levels=c("Adeno","NE"))
tmp1_data$Var1 <- factor(tmp1_data$Var1,levels=c("Epi", "Fibro",  "Endo", "Macro", "Tcell", "Bcell"))
library(ggpubr)
library(ggalluvial)
ff <- ggplot(tmp1_data, aes( x = group,y=Ratio,fill = Var1, stratum = Var1, alluvium = Var1))+
  geom_col(position = 'stack', width = 0.6)+  scale_y_continuous(breaks=seq(0,1,0.2)) +
  theme_bw() + scale_fill_manual(values=c("#d73027","#8c510a","#8073ac","#4575b4","#de77ae","#5aae61","#92c5de","#4393c3")) + geom_flow(alpha = 0.5) #绘制同类别之间的连接线
ggsave(ff,file="CRPC_CellType_Ratio.pdf",width=4,height=9)














Idents(Harmony_UMAP30) <- Harmony_UMAP30$CellType
Epi <- subset(Harmony_UMAP30,idents=c("Epi"),invert=FALSE)
Epi_object <- CreateSeuratObject(counts =GetAssayData(Epi, slot = "counts",assay="RNA")[,rownames(Epi@meta.data)], meta.data = Epi@meta.data)

Epi_DATA <- Epi_object %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
Epi_DATA.pca <- RunPCA(Epi_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(Epi_DATA.pca,file="Epi_CRPC_MERGE_DATA_pca.rds",mc.cores=20)
library(harmony)
library(reticulate)
library(ReductionWrappers)
library(s2a)
Epi_DATA_harmony  <- RunHarmony(Epi_DATA.pca,"orig.ident")
FindNeighbors <- FindNeighbors(Epi_DATA_harmony,reduction = "harmony", dims = 1:20)
FindClusters <- FindClusters(FindNeighbors,resolution = c(0.1,0.2,0.5))
mcsaveRDS(FindClusters,file="Epi_CRPC_Harmony_FindClusters.rds", mc.cores = 20)

Harmony_UMAP20 <- RunUMAP(FindClusters,reduction = "harmony", dims = 1:20)
mcsaveRDS(Harmony_UMAP20,file="Epi_CRPC_Harmony_UMAP20.rds",mc.cores=20)

Idents(Harmony_UMAP20) <- Harmony_UMAP20$RNA_snn_res.0.5
all.markers <- FindAllMarkers(Harmony_UMAP20, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.markers,"Epi_CRPC_UMAP20_res05_FindAllmarkers.csv")

p2 <- DimPlot(object = Harmony_UMAP20, reduction = "umap",label=TRUE,
    group.by="RNA_snn_res.0.5") +NoLegend()+labs(title="UMAP")
ggsave(p2,file="Epi_CRPC_UMAP20.png", width=8,height=8)




ff <- FeaturePlot(Harmony_UMAP20, features = c("CHGA","CHGB","ASCL1"),order=TRUE,pt.size=0.1,ncol=3,
cols=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#c7e9b4","#fed976","#fd8d3c","#fc4e2a","#e31a1c","#b10026"))
ggsave(ff,file="NE_Marker.png",width=16,height=5)



Idents(Harmony_UMAP20) <- Harmony_UMAP20$RNA_snn_res.0.5
markers.to.plot <- as.character(c("PTPRC","CD3E","CD3G","CD79A","CD79B","CSF1R","CD68","COL1A1","COL1A2","PECAM1","PLVAP"))
ff <- DotPlot(Harmony_UMAP20,features = markers.to.plot,cols = c("#ffffff","#cb181d"),
    assay="RNA",col.min =0,col.max =3,dot.min=0.05,scale=TRUE) + RotatedAxis()
ggsave(ff,file="Marker_Filter_Epi_QC.pdf",width=9.5,height=4.5)




Idents(Harmony_UMAP20) <- Harmony_UMAP20$RNA_snn_res.0.5
Pure_Epi <- subset(Harmony_UMAP20,idents=c("19","20","22","27"),invert=TRUE)
Pure_Epi<- RunUMAP(Pure_Epi,reduction = "harmony", dims = 1:20)

Idents(Pure_Epi) <- Pure_Epi$RNA_snn_res.0.5
new.cluster.ids <- c("Adeno","NE","Adeno","Adeno","Adeno","Adeno","Adeno","NE","Adeno","Adeno",
"Adeno","Adeno","Adeno","NE","NE","Adeno","Adeno","Adeno","Adeno","NE","Adeno","Adeno","NE","Adeno","Adeno")
names(new.cluster.ids) <- levels(Pure_Epi)
Pure_Epi <- RenameIdents(Pure_Epi, new.cluster.ids)
Pure_Epi$sub_CellType <- Idents(Pure_Epi)
Idents(Pure_Epi) <- Pure_Epi$sub_CellType
Pure_Epi$sub_CellType <- factor(Pure_Epi$sub_CellType,levels=c("Adeno","NE"))
mcsaveRDS(Pure_Epi,file="Renames_Pure_NE_UMAP20.rds",mc.cores=20)


p2 <- DimPlot(object = Pure_Epi, reduction = "umap",label=TRUE,
    group.by="sub_CellType") +NoLegend()+labs(title="UMAP")
ggsave(p2,file="Epi_pure_subCellType.png", width=8,height=8)



ff <- FeaturePlot(Pure_Epi, features = c("CHGA","CHGB","ASCL1"),order=TRUE,pt.size=0.1,ncol=3,
cols=c("#0c2c84","#225ea8","#1d91c0","#41b6c4","#c7e9b4","#fed976","#fd8d3c","#fc4e2a","#e31a1c","#b10026"))
ggsave(ff,file="NE_Marker.png",width=16,height=5)
speci_raw <- FetchData(object = Pure_Epi, vars = "IGFBP5")
Pure_Epi[["IGFBP5"]] <- (rowSums(speci_raw))/1
tmp <- Pure_Epi@meta.data
tmp$sub_CellType <- factor(tmp$sub_CellType,levels=c("Adeno","NE"))
library(ggpubr)
library(ggplot2)
my_comparisons <- list(c("Adeno","NE"))
ff <- ggboxplot(tmp, x = "sub_CellType", y = "IGFBP5",
               color = "sub_CellType",
                add = "",palette = c("#08519c","#e31a1c")) +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="wilcox.test")
ggsave(ff,file="IGFBP5.pdf",width=5.5,height=5)




Harmony_UMAP30 <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Neuronendocrine_Prostate_cancer/Chen_Science/Renames_CRPC_UMAP30.rds",mc.cores=20)
Idents(Harmony_UMAP30) <- Harmony_UMAP30$CellType
Endo <- subset(Harmony_UMAP30,idents=c("Endo"),invert=FALSE)
Endo_object <- CreateSeuratObject(counts =GetAssayData(Endo, slot = "counts",assay="RNA")[,rownames(Endo@meta.data)], meta.data = Endo@meta.data)

Endo_DATA <- Endo_object %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
Endo_DATA.pca <- RunPCA(Endo_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(Endo_DATA.pca,file="Endo_CRPC_MERGE_DATA_pca.rds",mc.cores=20)
library(harmony)
library(reticulate)
library(ReductionWrappers)
library(s2a)
Endo_DATA_harmony  <- RunHarmony(Endo_DATA.pca,"orig.ident")
FindNeighbors <- FindNeighbors(Endo_DATA_harmony,reduction = "harmony", dims = 1:20)
FindClusters <- FindClusters(FindNeighbors,resolution = c(0.1,0.2,0.5))
mcsaveRDS(FindClusters,file="Endo_CRPC_Harmony_FindClusters.rds", mc.cores = 20)

Harmony_UMAP20 <- RunUMAP(FindClusters,reduction = "harmony", dims = 1:20)
mcsaveRDS(Harmony_UMAP20,file="Endo_CRPC_Harmony_UMAP20.rds",mc.cores=20)

Idents(Harmony_UMAP20) <- Harmony_UMAP20$RNA_snn_res.0.5
all.markers <- FindAllMarkers(Harmony_UMAP20, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.markers,"Endo_CRPC_UMAP20_res05_FindAllmarkers.csv")