


https://ega-archive.org/studies/EGAS00001005549

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figure1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
library(Seurat)
library(ggplot2)
library(pheatmap)
require(RColorBrewer)
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/MyBestFunction_scRNA.R")
source("/mnt/data/user_data/xiangyu/programme/R_PACKAGES/my_code/Pseudo_CNV_series.R")
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



#Fig 1A
tmp <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_clinical_trail.csv")
tmp_1 <- tmp[,c("ICI_vs_Contrl_ORR","ICI_vs_Control_OS","Tumor_type","Trial.name","pat.No")]
tmp_1 <- na.omit(tmp_1)
tmp_1$ICI_vs_Contrl_ORR <- as.numeric(as.character(tmp_1$ICI_vs_Contrl_ORR))
tmp_1$ICI_vs_Control_OS <- as.numeric(as.character(tmp_1$ICI_vs_Control_OS))

library(ggplot2)
ff <- ggplot(tmp_1, aes(x = ICI_vs_Control_OS,y=ICI_vs_Contrl_ORR,size = pat.No, color = Tumor_type)) +
  geom_point(alpha = 1) +  
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") + 
  scale_size_continuous(range = c(2, 10)) +  
  theme_classic() +  
  labs(x = "ICI_vs_Control_OS", y = "ICI_vs_Contrl_ORR", 
       title = "Correlation of ICI Response with OS",
       size = "Sample Size", color = "Tumor_type") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  )
ggsave(ff, file="Fig1A.png",width=7,height=6)



#Fig 1B
00000000000000000 Bladder ICB 0000000000000000000000000000000000000
00000000000000000 Bladder ICB 0000000000000000000000000000000000000
00000000000000000 Bladder ICB 0000000000000000000000000000000000000
00000000000000000 Bladder ICB 0000000000000000000000000000000000000

Bladder_cancer_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Bladder/renew_Bladder_cancer_deconv_fraction.csv")
Bladder.clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Bladder/GSE176307.sample.info.clean.csv")
PURE01_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Bladder/PURE01_info/PURE01_Bladder_cancer_deconv_fraction.csv")
PURE01_clinical_info <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Bladder/PURE01_info/PURE01_clinical_info.csv")
bladder_GSE176307 <- Bladder_cancer_deconv_fraction[,c("X","T.cell")]
bladder_PURE01 <- PURE01_deconv_fraction[,c("X","T.cell")]
BLCA_T_ratio <- rbind(bladder_GSE176307, bladder_PURE01)
names(BLCA_T_ratio) <- c("ID","T.cell")


00000000000000000 SCLC Nature2024 AND IMpower133 0000000000000000000000000000000000000
00000000000000000 SCLC Nature2024 AND IMpower133 0000000000000000000000000000000000000
00000000000000000 SCLC Nature2024 AND IMpower133 0000000000000000000000000000000000000
00000000000000000 SCLC Nature2024 AND IMpower133 0000000000000000000000000000000000000

IMPOWER133_Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_Clinical_file.csv")
IMPOWER133_Clinical <- subset(IMPOWER133_Clinical,ACTARM.2=="atezo")
IMPOWER133_Clinical <- IMPOWER133_Clinical[,c("trunc_anonymized_sample_ids","PFS_MONTHS","PFS_CENSOR")]
IMPOWER133_Clinical$X <- gsub("\\-","\\.",IMPOWER133_Clinical$trunc_anonymized_sample_ids)
IMPOWER133_SCLC_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_SCLC_deconv_fraction.csv")
SCLC_IMPOWER133_PFS_all <- merge(IMPOWER133_Clinical,IMPOWER133_SCLC_deconv_fraction,by="X")
SCLC_IMPOWER133_PFS_all <- data.frame(SCLC_IMPOWER133_PFS_all[,c("X","T.cell")])
names(SCLC_IMPOWER133_PFS_all) <- c("ID","T.cell")


SCLC_Nature_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/Nature_Evolution_SCLC_deconv_fraction.csv")
#remove PDX,CDX data
SCLC_Nature_deconv_fraction$sample <- substring(SCLC_Nature_deconv_fraction$X,1,6)
library(dplyr)
SCLC_Nature_deconv_fraction <- SCLC_Nature_deconv_fraction %>% dplyr::arrange(sample)
SCLC_Nature_deconv_fraction_mean <- aggregate(x = SCLC_Nature_deconv_fraction[,2:17],
                            by = list(SCLC_Nature_deconv_fraction$sample),
                            FUN = mean)
SCLC_Nature_PFS_all <- data.frame(SCLC_Nature_deconv_fraction_mean[,c("Group.1","T.cell")])
names(SCLC_Nature_PFS_all) <- c("ID","T.cell")
SCLC_T_ratio <- rbind(SCLC_IMPOWER133_PFS_all, SCLC_Nature_PFS_all)


00000000000000000 LUSC AND LUAD 0000000000000000000000000000000000000
00000000000000000 LUSC AND LUAD 0000000000000000000000000000000000000
00000000000000000 LUSC AND LUAD 0000000000000000000000000000000000000
00000000000000000 LUSC AND LUAD 0000000000000000000000000000000000000

NSCLC_NatureGenetics_clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NatureGenetics_clincal.info.csv")
LUAD_id <- subset(NSCLC_NatureGenetics_clinical,Histology_Harmonized=="Adeno")
LUSC_id <- subset(NSCLC_NatureGenetics_clinical,Histology_Harmonized=="Squamous")
NSCLC_NatureGenetics_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NatureGenetics_deconv_fraction.csv")
rownames(NSCLC_NatureGenetics_deconv_fraction) <- NSCLC_NatureGenetics_deconv_fraction$Sample.ID
LUAD_1 <- NSCLC_NatureGenetics_deconv_fraction[na.omit(LUAD_id$Harmonized_SU2C_Participant_ID_v2),]
LUSC_1 <- NSCLC_NatureGenetics_deconv_fraction[na.omit(LUSC_id$Harmonized_SU2C_Participant_ID_v2),]
LUAD_1 <- LUAD_1[,c("Sample.ID", "T.cell")]
LUSC_1 <- LUSC_1[,c("Sample.ID", "T.cell")]
names(LUAD_1) <- c("ID", "T.cell")
names(LUSC_1) <- c("ID", "T.cell")

ICI_CancerType <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NSCLC_ICI_GSE190265.csv")
LUSC_id <- subset(ICI_CancerType,Tumor_Type=="Squa")
LUAD_id <- subset(ICI_CancerType,Tumor_Type=="Non_Squa")
NSCLC_GSE190265_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/GSE190265_deconv_fraction.csv")
rownames(NSCLC_GSE190265_deconv_fraction) <- NSCLC_GSE190265_deconv_fraction$X
LUSC_2 <- NSCLC_GSE190265_deconv_fraction[LUSC_id$Group,]
LUAD_2 <- NSCLC_GSE190265_deconv_fraction[LUAD_id$Group,]
LUAD_2 <- LUAD_2[,c("X", "T.cell")]
LUSC_2 <- LUSC_2[,c("X", "T.cell")]
names(LUAD_2) <- c("ID", "T.cell")
names(LUSC_2) <- c("ID", "T.cell")

LUAD_T_ratio <- rbind(LUAD_1,LUAD_2)
LUSC_T_ratio <- rbind(LUSC_1,LUSC_2)



00000000000000000 GBM 0000000000000000000000000000000000000
00000000000000000 GBM 0000000000000000000000000000000000000
00000000000000000 GBM 0000000000000000000000000000000000000
00000000000000000 GBM 0000000000000000000000000000000000000

GBM_deconvolution <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/GBM/GBM_cancer_deconv_fraction.csv")
GBM_deconvolution$T.cell <- GBM_deconvolution$CD4..T.cells+GBM_deconvolution$CD8..T.cells
GBM_deconvolution <- GBM_deconvolution[,c("X","T.cell")]
names(GBM_deconvolution) <- c("ID", "T.cell")



00000000000000000 Melanoma 0000000000000000000000000000000000000
00000000000000000 Melanoma 0000000000000000000000000000000000000
00000000000000000 Melanoma 0000000000000000000000000000000000000
00000000000000000 Melanoma 0000000000000000000000000000000000000

All_filter_clinical_info <- readRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/already_finished/0_Endo_SCLC/first_submit/ImmuneCheckpoint_blokage_Clinical_data/ICI_Cohort/Melanoma/Melanoma_unknown/All_filter_clinical_info.rds")

Melanoma_GSE78220 <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/GSE78220/Melanoma_GSE78220_deconv_fraction.csv")
Melanoma_GSE91061 <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/GSE91061/GSE91061_Melanoma_deconv_fraction.csv")
Melanoma_GSE100797 <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/GSE100797/Melanoma_GSE100797_deconv_fraction.csv")
#PRJEB23709
Melanoma_unkown <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/Melanoma_unknown/Melanoma_deconv_fraction.csv")
Melanoma_Phs000452 <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/Phs000452/Melanoma_phs000452_deconv_fraction.csv")
Melanoma_GSE78220 <- Melanoma_GSE78220[,c("X","T_Cell")]
Melanoma_GSE91061 <- Melanoma_GSE91061[,c("X","T_Cell")]
Melanoma_GSE100797 <- Melanoma_GSE100797[,c("X","T_Cell")]
Melanoma_unkown <- Melanoma_unkown[,c("X","T_Cell")]
Melanoma_Phs000452 <- Melanoma_Phs000452[,c("X","T_Cell")]

all_MM <- rbind(Melanoma_GSE78220,
Melanoma_GSE91061,
Melanoma_GSE100797,
Melanoma_unkown,
Melanoma_Phs000452)

names(all_MM) <- c("ID","T.cell")



00000000000000000 SRC 0000000000000000000000000000000000000

SRC_ICI <- readRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SRC/SRC_ICI_bulkRNA.rds")

SRC <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SRC/SRC_deconv_fraction.csv")
SRC <- SRC[,c("X","T.cell")]
names(SRC) <- c("ID","T.cell")
BLCA_T_ratio$CancerType <- "BLCA"
SCLC_T_ratio$CancerType <- "SCLC"
LUAD_T_ratio$CancerType <- "LUAD"
LUSC_T_ratio$CancerType <- "LUSC"
GBM_deconvolution$CancerType <- "GBM"
all_MM$CancerType <- "Melanoma"
SRC$CancerType <- "SRC"

All <- rbind(BLCA_T_ratio,SCLC_T_ratio,LUAD_T_ratio,LUSC_T_ratio,GBM_deconvolution,all_MM,SRC)
All_1 <- na.omit(All)
library(ggplot2)
library(dplyr)

All_1 <- All_1 %>% mutate(T.cell.percent = T.cell * 100)
All_1$CancerType <- factor(All_1$CancerType, levels = c("GBM", "SCLC", "LUSC", "SRC", "BLCA", "LUAD", "Melanoma"))
ff <- ggplot(All_1, aes(x = T.cell.percent, y = CancerType)) +
  geom_linerange(aes(xmin = T.cell.percent, xmax = T.cell.percent, ymin = as.numeric(CancerType) - 0.4, ymax = as.numeric(CancerType) + 0.4),
                 color = "brown", alpha = 0.6) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50,60,70,80,90,100)) +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    plot.margin = margin(5, 5, 5, 5)
  ) +
  labs(x = "T cell infiltration")
ggsave(ff, file="Fig1B.png",width=10,height=8)






#Fig 1C

#Melanoma
Melanoma_unkown_deconv <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/Melanoma_unknown/Melanoma_deconv_fraction.csv")
Melanoma_unkown_clinical <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/Melanoma_unknown/All_filter_clinical_info.rds",mc.cores=20)
Melanoma_Tratio <- Melanoma_unkown_deconv[,c("X","T_Cell")]
clinical_1 <- Melanoma_unkown_clinical[,c("sample_name", "Progression.Free.Survival..Days.")]
clinical_1 <- na.omit(clinical_1)
Melanoma_PFS_all <- merge(Melanoma_Tratio,clinical_1,by.x="X",by.y="sample_name")
Melanoma_PFS_all$Status <- "NA"
Melanoma_PFS_all$CancerType <- "Melanoma"
names(Melanoma_PFS_all) <- c("ID","T.cell","PFS","Status","CancerType")
Melanoma_PFS_all$Cohort <- "unknown"

Melanoma_GSE100797 <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/GSE100797/Melanoma_GSE100797_deconv_fraction.csv")
GSE100797_clinical <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/GSE100797/GSE100797.sample.info.clean.rds",mc.cores=20)
GSE100797_clinical$X <- rownames(GSE100797_clinical)
GSE100797_Tratio <- Melanoma_GSE100797[,c("X","T_Cell")]
GSE100797_PFS_all <- merge(GSE100797_Tratio,GSE100797_clinical,by="X")
GSE100797_PFS_all$PFS <- as.numeric(GSE100797_PFS_all$pfs.time)*30
GSE100797_PFS_all$CancerType <- "Melanoma"
GSE100797_PFS_all$Status <- "NA"
GSE100797_PFS_all <- GSE100797_PFS_all[,c("X","PFS","T_Cell","Status","CancerType")]
names(GSE100797_PFS_all) <- c("ID","PFS","T.cell","Status","CancerType")
GSE100797_PFS_all$Cohort <- "GSE100797"




#NSCLC
NatureGenetics_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NatureGenetics_deconv_fraction.csv")
NatureGenetics_clincal.info <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NatureGenetics_clincal.info.csv")
LUAD_id <- subset(NatureGenetics_clincal.info,Histology_Harmonized=="Adeno")
LUSC_id <- subset(NatureGenetics_clincal.info,Histology_Harmonized=="Squamous")
NatureGenetics_clincal.info <- NatureGenetics_clincal.info[,c("Harmonized_SU2C_Participant_ID_v2","Harmonized_PFS_Days","Harmonized_PFS_Event")]
rownames(NatureGenetics_clincal.info) <- NatureGenetics_clincal.info$Harmonized_SU2C_Participant_ID_v2
NatureGenetics_clincal.info <- na.omit(NatureGenetics_clincal.info)
names(NatureGenetics_clincal.info) <- c("ID","PFS","Status")
NSCLC_NatureGenetics_PFS_all <- merge(NatureGenetics_clincal.info,NatureGenetics_deconv_fraction,by.x="ID",by.y="Sample.ID")
NSCLC_NatureGenetics_PFS_all <- data.frame(NSCLC_NatureGenetics_PFS_all[,c("ID","PFS","Status","T.cell")])
names(NSCLC_NatureGenetics_PFS_all) <- c("ID","PFS","Status","T.cell")
rownames(NSCLC_NatureGenetics_PFS_all) <- NSCLC_NatureGenetics_PFS_all$ID
LUAD_NG_PFS_All <- na.omit(NSCLC_NatureGenetics_PFS_all[LUAD_id$Harmonized_SU2C_Participant_ID_v2, ])
LUSC_NG_PFS_All <- na.omit(NSCLC_NatureGenetics_PFS_all[LUSC_id$Harmonized_SU2C_Participant_ID_v2, ])
LUAD_NG_PFS_All$CancerType <- "LUAD"
LUAD_NG_PFS_All$Cohort <- "Nature_SCLC"
LUSC_NG_PFS_All$CancerType <- "LUSC"
LUSC_NG_PFS_All$Cohort <- "Nature_SCLC"


ICI_CancerType <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NSCLC_ICI_GSE190265.csv")
NSCLC_GSE190265 <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/GSE190265_deconv_fraction.csv")
GSE190265_clinical <- read.csv(gzfile("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/GSE190265_samples_Clinical_info_France3.csv.gz"),sep=";",row.names=1)
LUSC_id <- subset(ICI_CancerType,Tumor_Type=="Squa")
LUAD_id <- subset(ICI_CancerType,Tumor_Type=="Non_Squa")
NSCLC_GSE190265_Tratio <- NSCLC_GSE190265[,c("X","T.cell")]
GSE190265_clinical$X <- rownames(GSE190265_clinical)
GSE190265_PFS_all <- merge(NSCLC_GSE190265_Tratio,GSE190265_clinical,by="X")
GSE190265_PFS_all$PFS <- as.numeric(GSE190265_PFS_all$time_PFS)*30
GSE190265_PFS_all <- data.frame(GSE190265_PFS_all[,c("X","PFS","evtPFS","T.cell")])
names(GSE190265_PFS_all) <- c("ID","PFS","Status","T.cell")
rownames(GSE190265_PFS_all) <- GSE190265_PFS_all$ID
GSE190265_PFS_all$ID <- paste0("NSCLC_",GSE190265_PFS_all$ID,sep="")
GSE190265_PFS_all$ID <- gsub("\\-","\\.",GSE190265_PFS_all$ID)
LUAD_GSE190265_PFS_All <- na.omit(GSE190265_PFS_all[LUAD_id$Group, ])
LUSC_GSE190265_PFS_All <- na.omit(GSE190265_PFS_all[LUSC_id$Group, ])
LUAD_GSE190265_PFS_All$CancerType <- "LUAD"
LUAD_GSE190265_PFS_All$Cohort <- "LUAD_GSE190265"
LUSC_GSE190265_PFS_All$CancerType <- "LUSC"
LUSC_GSE190265_PFS_All$Cohort <- "LUSC_GSE190265"


LUAD_PFS_All <- rbind(LUAD_NG_PFS_All,LUAD_GSE190265_PFS_All)
LUSC_PFS_All <- rbind(LUSC_GSE190265_PFS_All,LUSC_NG_PFS_All)





#SCLC
Nature_SCLC_survival_info <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/Nature_SCLC_survival_info.csv")
Nature_SCLC_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/Nature_Evolution_SCLC_deconv_fraction.csv")

Nature_SCLC_sample_info <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/Nature_SCLC_sample_info.csv")
sample.1 <- subset(Nature_SCLC_sample_info, Sequencing..sample.1. == "WES, RNA-seq" |  Sequencing..sample.1. =="WGS, RNA-seq")
sample.1 <- subset(sample.1, Site_acquisition_sample.1 != "CTC (PDX model)")
sample.2 <- subset(Nature_SCLC_sample_info, Sequencing..sample.2. == "WES, RNA-seq" |  Sequencing..sample.2. =="WES, WGS, RNA-seq")
sample.2 <- subset(sample.2, Site_acquisition_sample.2 != "CTC (PDX model)")
sample.3 <- subset(Nature_SCLC_sample_info, Sequencing..sample.3. == "WES, RNA-seq" |  Sequencing..sample.3. =="WES, WGS, RNA-seq")
sample.3 <- subset(sample.3, Site_acquisition_sample.3!= "pleura metastasis (PDX model)"  & Site_acquisition_sample.3!= "CTC (PDX model)" )
sample.4 <- subset(Nature_SCLC_sample_info, Sequencing..sample.4. == "WES, WGS, RNA-seq")
sample.4 <- subset(sample.4, Site_acquisition_sample.4 != "CTC (PDX model)")
useful_id_1 <- data.frame(sample.1$Sample.1.ID)
useful_id_2 <- data.frame(sample.2$Sample.2.ID)
names(useful_id_1) <- "ID"
names(useful_id_2) <- "ID"
useful_ID <- rbind(useful_id_1,useful_id_2)

Nature_SCLC_tmp <- merge(useful_ID,Nature_SCLC_deconv_fraction,by.x="ID",by.y="X")
Nature_SCLC_tmp$ID.2 <- substring(Nature_SCLC_tmp$ID,1,6)
Nature_SCLC_tmp <- Nature_SCLC_tmp[order(Nature_SCLC_tmp$ID.2, decreasing = F),]

library(dplyr)
Nature_SCLC_tmp <- Nature_SCLC_tmp %>% dplyr::arrange(ID.2)
Nature_SCLC_tmp_mean <- aggregate(x = Nature_SCLC_tmp[,2:17],
                            by = list(Nature_SCLC_tmp$ID.2),
                            FUN = mean)
SCLC_PFS_all <- merge(Nature_SCLC_tmp_mean,Nature_SCLC_survival_info,by.x="Group.1",by.y="Sample.ID")
immune_ID <- subset(SCLC_PFS_all, X2nd.line.treatment=="immunotherapy (anti-PD-1)" | X2nd.line.treatment=="immunotherapy (anti-PD-1 + anti-CTLA-4)")
immune_ID[10,39] <- "421"
SCLC_Nature_PFS <- immune_ID[,c("Group.1","T.cell","Status..at.time.of.last.follow.up.","Progression.free_survival..days.")]
SCLC_Nature_PFS$CancerType <- "SCLC"
names(SCLC_Nature_PFS) <- c("ID","T.cell","Status","PFS","CancerType")
SCLC_Nature_PFS$Cohort <- "Nature_SCLC"
#mcsaveRDS(SCLC_Nature_PFS, file="SCLC_Nature_evo_bulkRNA.rds")



#Bladder
Bladder_cancer_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Bladder/renew_Bladder_cancer_deconv_fraction.csv")
Bladder.clinical <- readRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Bladder/GSE176307.sample.info.clean.rds")
Bladder.clinical$X <- rownames(Bladder.clinical)
Bladder.clinical <- Bladder.clinical[,c("X","pfs","alive")]
Bladder_PFS_all <- merge(Bladder.clinical,Bladder_cancer_deconv_fraction,by="X")
Bladder_PFS_all$Status <- ifelse(Bladder_PFS_all$alive=="No","1","0")
Bladder_PFS_all$ID <- Bladder_PFS_all$X
Bladder_PFS_all <- data.frame(Bladder_PFS_all[,c("ID","T.cell","Status","pfs")])
names(Bladder_PFS_all) <- c("ID","T.cell","Status","PFS")
Bladder_PFS_all$CancerType <- "Bladder"
Bladder_PFS_all$Cohort <- "GSE176307_Bladder"



#SRC
SRC_exp <- read.table("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SRC/scRNA_as_reference/immunotherapy/GSE213065_Sarcoma_ecotyper_tpm_bulksamples.txt",sep="\t",head=TRUE)
rownames(SRC_exp) <- SRC_exp$gene
SRC_exp <- SRC_exp[,-1]
SRC_exp <- data.frame(t(SRC_exp))
SRC_exp$ID <- rownames(SRC_exp)
SRC_exp$ID <- substr(SRC_exp$ID, 1, nchar(SRC_exp$ID)- 4)
SRC_exp_mean <- aggregate(x = SRC_exp[,1:18952],
                            by = list(SRC_exp$ID),
                            FUN = mean)
rownames(SRC_exp_mean) <- SRC_exp_mean$Group.1
SRC_exp <- SRC_exp_mean[,-1]
#mcsaveRDS(SRC_exp,file="SRC_ICI_bulkRNA.rds",mc.cores=20)


SRC_Clinical_info <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SRC/scRNA_as_reference/immunotherapy/SRC_Clinical_info_Stanford_immunotherapy_cohort.csv")
SRC_Clinical_info$Status <- SRC_Clinical_info$PFS.Status
SRC_Clinical_info$PFS <- SRC_Clinical_info$PFS.Months*30
SRC_Clinical_info$ID <- SRC_Clinical_info$Patient

SRC_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SRC/SRC_deconv_fraction.csv")
SRC_deconv_fraction$ID <- SRC_deconv_fraction$X
SRC_deconv_fraction$ID <- substr(SRC_deconv_fraction$ID, 1, nchar(SRC_deconv_fraction$ID) - 4)

library(dplyr)
SRC_deconv_fraction <- SRC_deconv_fraction %>% dplyr::arrange(ID)
SRC_deconv_fraction_mean <- aggregate(x = SRC_deconv_fraction[,2:9],
                            by = list(SRC_deconv_fraction$ID),
                            FUN = mean)
rownames(SRC_deconv_fraction_mean) <- SRC_deconv_fraction_mean$Group.1
SRC_deconv_fraction <- SRC_deconv_fraction_mean[,-1]
SRC_deconv_fraction$ID <- rownames(SRC_deconv_fraction)
tmp <- merge(SRC_Clinical_info,SRC_deconv_fraction,by="ID")
SRC_all_info <- data.frame(tmp[,c("ID","PFS","Status","T.cell")])
SRC_all_info$CancerType <- c("SRC")
SRC_all_info$Cohort <- c("SRC_Standford")


ALL <- rbind(
Melanoma_PFS_all,
GSE100797_PFS_all,
LUAD_PFS_All,
LUSC_PFS_All,
SCLC_Nature_PFS,
Bladder_PFS_all,
SRC_all_info)
write.csv(ALL,file="/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Fig1C_PFS_final_ICI_cohort.csv")




SRC_ICI_bulkRNA <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SRC/SRC_ICI_bulkRNA.rds",mc.cores=20)
Bladder_GSE176307_tpm <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Bladder/Bladder_GSE176307_bulkRNA.rds",mc.cores=20)
Melanoma_GSE100797_ProcessedData <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/GSE100797/Melanoma_GSE100797_bulkRNA.rds",mc.cores=20)
Melanoma_Tratio_exp <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Melanoma/Melanoma_unknown/Melanoma_unknown_bulkRNA.rds",mc.cores=20)
NSCLC_GSE190265_exp <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NSCLC_GSE190265_bulkRNA.rds",mc.cores=20)
NSCLC_NatureGenetics_RNA_mean <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/NSCLC/NSCLC_NatureGenetics_bulkRNA.rds",mc.cores=20)

SCLC_Nature_RNA_exp_mean <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/SCLC_Nature_evolution_bulkRNA.rds",mc.cores=20)

gene1 <- colnames(SRC_ICI_bulkRNA)
gene2 <- colnames(SCLC_Nature_RNA_exp_mean)
gene3 <- colnames(NSCLC_NatureGenetics_RNA_mean)
gene4 <- colnames(NSCLC_GSE190265_exp)
gene5 <- colnames(Melanoma_GSE100797_ProcessedData)
gene6 <- colnames(Melanoma_Tratio_exp)
gene7 <- colnames(Bladder_GSE176307_tpm)

common_cols <- Reduce(intersect, list(gene1, gene2, gene3, gene4, gene5, gene6,gene7))

SRC_ICI_bulkRNA <- SRC_ICI_bulkRNA[,common_cols]
SCLC_Nature_data <- SCLC_Nature_RNA_exp_mean[,common_cols]
NSCLC_NG_data <- NSCLC_NatureGenetics_RNA_mean[,common_cols]
NSCLC_GSE190265_data <- NSCLC_GSE190265_exp[,common_cols]
Melanoma_GSE100797_data <- Melanoma_GSE100797_ProcessedData[,common_cols]
Melanoma_unknown_data <- Melanoma_Tratio_exp[,common_cols]
Bladder_GSE176307_data <- Bladder_GSE176307_tpm[,common_cols]


all_matrix <- rbind(SCLC_Nature_data,NSCLC_NG_data,NSCLC_GSE190265_data,Melanoma_GSE100797_data,Melanoma_unknown_data,Bladder_GSE176307_data,SRC_ICI_bulkRNA)
all_dat <- data.frame(t(all_matrix))


batch <- c(rep(1, nrow(SCLC_Nature_data)), rep(2, nrow(NSCLC_NG_data)),
    rep(3, nrow(NSCLC_GSE190265_data)),rep(4, nrow(Melanoma_GSE100797_data)),
    rep(5, nrow(Melanoma_unknown_data)),rep(6, nrow(Bladder_GSE176307_data)), rep(7, nrow(SRC_ICI_bulkRNA)))

library(sva)
combat_all_genes <- ComBat(dat = all_dat, batch = batch, mod = NULL)


######Remove patients without PFS data#########################
######Remove patients without PFS data#########################
######Remove patients without PFS data#########################

ALL <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Fig1C_PFS_final_ICI_cohort.csv")
ALL$ID <- gsub("\\-","\\.",ALL$ID)
patient_with_PFS<- intersect(colnames(combat_all_genes),ALL$ID)
all_exp <- combat_all_genes[,patient_with_PFS]

h_geneSets <- getGmt("/mnt/data/user_data/xiangyu/workshop/DATABASE/GSVA_7.1/h.all.v7.1.symbols.gmt")
c5_geneSets <- getGmt("/mnt/data/user_data/xiangyu/workshop/DATABASE/GSVA_7.1/c5.all.v7.1.symbols.gmt")
c2_geneSets <- getGmt("/mnt/data/user_data/xiangyu/workshop/DATABASE/GSVA_7.1/c2.all.v7.1.symbols.gmt")

matrix <- as.matrix(all_exp)

library(GSVA)

h_gsva_param <- gsvaParam(expr = matrix, geneSets = h_geneSets, minSize = 5, maxSize = 500)
h_GSVA_res <- gsva(h_gsva_param, verbose = FALSE)
c2_gsva_param <- gsvaParam(expr = matrix, geneSets = c2_geneSets, minSize = 5, maxSize = 500)
c2_GSVA_res <- gsva(c2_gsva_param, verbose = FALSE)
c5_gsva_param <- gsvaParam(expr = matrix, geneSets = c5_geneSets, minSize = 5, maxSize = 500)
c5_GSVA_res <- gsva(c5_gsva_param, verbose = FALSE)

sudobulk_GSVA_score <- rbind(h_GSVA_res,c2_GSVA_res,c5_GSVA_res)
write.csv(sudobulk_GSVA_score,"Fig1C_GSVA_score.csv")



sudobulk_GSVA_score <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Fig1C_GSVA_score.csv")
rownames(sudobulk_GSVA_score) <- sudobulk_GSVA_score$X
sudobulk_GSVA_score <- sudobulk_GSVA_score[,-1]
sel_path <- c("GO_ESTABLISHMENT_OF_ENDOTHELIAL_BLOOD_BRAIN_BARRIER",
"GO_REGULATION_OF_NK_T_CELL_ACTIVATION",
"HALLMARK_TNFA_SIGNALING_VIA_NFKB",
"GO_ENDOTHELIAL_CELL_MORPHOGENESIS",
"GO_VASCULAR_SMOOTH_MUSCLE_CELL_DEVELOPMENT")

sel_GSVA <- sudobulk_GSVA_score[sel_path,]
sel_GSVA <- data.frame(t(sel_GSVA))
sel_GSVA$ID <- rownames(sel_GSVA) 
TMP <- merge(sel_GSVA,ALL,by="ID")
TMP <- TMP[order(TMP$T.cell, decreasing = TRUE),]
TMP$order <- c(1:413)

library(BuenColors)
p1 <- ggplot(TMP, aes(order, GO_VASCULAR_SMOOTH_MUSCLE_CELL_DEVELOPMENT,color=order)) +
  geom_point(alpha=0.5,size=0.5,colour="#969696") + scale_color_gradientn(colours = c("#4575b4","#abd9e9","#ffffbf","#fdae61","#d73027"))+
  geom_rug(alpha = 0.8, position = "jitter",sides="b")+ xlab("Tcell_infiltration_high_to_low")+ylim(-1,1)+
  theme_classic() + geom_smooth(colour = "orange",se=TRUE) +NoLegend()+ labs(title="GO_VASCULAR_SMOOTH_MUSCLE_CELL_DEVELOPMENT")
ggsave(p1,file="Fig1C_GO_VASCULAR_SMOOTH_MUSCLE_CELL_DEVELOPMENT.png",height=6,width=5)

p1 <- ggplot(TMP, aes(order, GO_REGULATION_OF_NK_T_CELL_ACTIVATION,color=order)) +
  geom_point(alpha=0.5,size=0.5,colour="#969696") + scale_color_gradientn(colours = c("#4575b4","#abd9e9","#ffffbf","#fdae61","#d73027"))+
  geom_rug(alpha = 0.8, position = "jitter",sides="b")+ xlab("Tcell_infiltration_high_to_low")+ylim(-1,1)+
  theme_classic() + geom_smooth(colour = "orange",se=TRUE) +NoLegend()+ labs(title="GO_REGULATION_OF_NK_T_CELL_ACTIVATION")
ggsave(p1,file="Fig1C_GO_REGULATION_OF_NK_T_CELL_ACTIVATION.png",height=6,width=5)

p1 <- ggplot(TMP, aes(order, GO_ENDOTHELIAL_CELL_MORPHOGENESIS,color=order)) +
  geom_point(alpha=0.5,size=0.5,colour="#969696") + scale_color_gradientn(colours = c("#4575b4","#abd9e9","#ffffbf","#fdae61","#d73027"))+
  geom_rug(alpha = 0.8, position = "jitter",sides="b")+ xlab("Tcell_infiltration_high_to_low")+ylim(-1,1)+
  theme_classic() + geom_smooth(colour = "orange",se=TRUE) +NoLegend()+ labs(title="GO_ENDOTHELIAL_CELL_MORPHOGENESIS")
ggsave(p1,file="Fig1C_GO_ENDOTHELIAL_CELL_MORPHOGENESIS.png",height=6,width=5)

p1 <- ggplot(TMP, aes(order, GO_ESTABLISHMENT_OF_ENDOTHELIAL_BLOOD_BRAIN_BARRIER,color=order)) +
  geom_point(alpha=0.5,size=0.5,colour="#969696") + scale_color_gradientn(colours = c("#4575b4","#abd9e9","#ffffbf","#fdae61","#d73027"))+
  geom_rug(alpha = 0.8, position = "jitter",sides="b")+ xlab("Tcell_infiltration_high_to_low")+ylim(-1,1)+
  theme_classic() + geom_smooth(colour = "orange",se=TRUE) +NoLegend()+ labs(title="GO_ESTABLISHMENT_OF_ENDOTHELIAL_BLOOD_BRAIN_BARRIER")
ggsave(p1,file="Fig1C_GO_ESTABLISHMENT_OF_ENDOTHELIAL_BLOOD_BRAIN_BARRIER.png",height=6,width=5)

p1 <- ggplot(TMP, aes(order, HALLMARK_TNFA_SIGNALING_VIA_NFKB,color=order)) +
  geom_point(alpha=0.5,size=0.5,colour="#969696") + scale_color_gradientn(colours = c("#4575b4","#abd9e9","#ffffbf","#fdae61","#d73027"))+
  geom_rug(alpha = 0.1, position = "jitter",sides="b")+ xlab("Tcell_infiltration_high_to_low")+
  theme_classic() + geom_smooth(colour = "orange",se=TRUE) +NoLegend()+ labs(title="HALLMARK_TNFA_SIGNALING_VIA_NFKB")
ggsave(p1,file="Fig1C_HALLMARK_TNFA_SIGNALING_VIA_NFKB.png",height=6,width=5)



#Fig 1D
IMPOWER133_Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_Clinical_file.csv")
IMPOWER133_Clinical <- subset(IMPOWER133_Clinical,ACTARM.2=="atezo")
IMPOWER133_Clinical <- IMPOWER133_Clinical[,c("trunc_anonymized_sample_ids","PFS_MONTHS","PFS_CENSOR")]
IMPOWER133_Clinical$X <- gsub("\\-","\\.",IMPOWER133_Clinical$trunc_anonymized_sample_ids)
IMPOWER133_SCLC_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_SCLC_deconv_fraction.csv")
SCLC_IMPOWER133_PFS_all <- merge(IMPOWER133_Clinical,IMPOWER133_SCLC_deconv_fraction,by="X")
SCLC_IMPOWER133_PFS_all <- data.frame(SCLC_IMPOWER133_PFS_all[,c("X","T.cell")])
names(SCLC_IMPOWER133_PFS_all) <- c("ID","T.cell")

SCLC_Nature_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/Nature_Evolution_SCLC_deconv_fraction.csv")
#remove PDX,CDX
SCLC_Nature_deconv_fraction$sample <- substring(SCLC_Nature_deconv_fraction$X,1,6)
library(dplyr)
SCLC_Nature_deconv_fraction <- SCLC_Nature_deconv_fraction %>% dplyr::arrange(sample)
SCLC_Nature_deconv_fraction_mean <- aggregate(x = SCLC_Nature_deconv_fraction[,2:17],
                            by = list(SCLC_Nature_deconv_fraction$sample),
                            FUN = mean)
SCLC_Nature_PFS_all <- data.frame(SCLC_Nature_deconv_fraction_mean[,c("Group.1","T.cell")])
names(SCLC_Nature_PFS_all) <- c("ID","T.cell")
SCLC_T_ratio <- rbind(SCLC_IMPOWER133_PFS_all, SCLC_Nature_PFS_all)


IMpower133_TPM <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_SCLC_TPM_file.csv")
Nature_evo_2024 <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/Nature_evo_RNA_exp_symbol_mean.csv")
rownames(IMpower133_TPM) <- IMpower133_TPM$UNNAMED..0
IMpower133_TPM <- IMpower133_TPM[,c(-1,-2)]
IMpower133_TPM_exp <- IMpower133_TPM[,SCLC_IMPOWER133_PFS_all$ID]
rownames(Nature_evo_2024) <- Nature_evo_2024$Group.1
Nature_evo_2024 <- Nature_evo_2024[,c(-1,-2)]
ID <- substring(colnames(Nature_evo_2024),1,6)
names(Nature_evo_2024) <- as.character(ID)
Nature_evo_2024_exp <- Nature_evo_2024[,SCLC_Nature_PFS_all$ID]
mcsaveRDS(IMpower133_TPM_exp,file="IMpower133_TPM_exp.rds",mc.cores=20)
mcsaveRDS(Nature_evo_2024_exp,file="Nature_evo_2024_exp.rds",mc.cores=20)


IMPOWER133_Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_Clinical_file.csv")
IMPOWER133_Clinical <- subset(IMPOWER133_Clinical,ACTARM.2=="atezo")
IMPOWER133_Clinical <- IMPOWER133_Clinical[,c("trunc_anonymized_sample_ids","PFS_MONTHS","PFS_CENSOR")]
IMPOWER133_Clinical$X <- gsub("\\-","\\.",IMPOWER133_Clinical$trunc_anonymized_sample_ids)

IMPOWER133_SCLC_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/IMPOWER133/IMPOWER133_SCLC_deconv_fraction.csv")
SCLC_IMPOWER133_PFS_all <- merge(IMPOWER133_Clinical,IMPOWER133_SCLC_deconv_fraction,by="X")
SCLC_IMPOWER133_PFS_all <- data.frame(SCLC_IMPOWER133_PFS_all[,c("X","T.cell")])
names(SCLC_IMPOWER133_PFS_all) <- c("ID","T.cell")

SCLC_Nature_deconv_fraction <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC/Nature/Nature_Evolution_SCLC_deconv_fraction.csv")
#remove PDX,CDX
SCLC_Nature_deconv_fraction$sample <- substring(SCLC_Nature_deconv_fraction$X,1,6)
library(dplyr)
SCLC_Nature_deconv_fraction <- SCLC_Nature_deconv_fraction %>% dplyr::arrange(sample)
SCLC_Nature_deconv_fraction_mean <- aggregate(x = SCLC_Nature_deconv_fraction[,2:17],
                            by = list(SCLC_Nature_deconv_fraction$sample),
                            FUN = mean)
SCLC_Nature_PFS_all <- data.frame(SCLC_Nature_deconv_fraction_mean[,c("Group.1","T.cell")])
names(SCLC_Nature_PFS_all) <- c("ID","T.cell")
SCLC_T_ratio <- rbind(SCLC_IMPOWER133_PFS_all, SCLC_Nature_PFS_all)



IMpower133_data <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/IMpower133_TPM_exp.rds",mc.cores=20)
Nature_evo_data <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Nature_evo_2024_exp.rds",mc.cores=20)

matrix_1 <- as.matrix(IMpower133_data)
h_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/h.all.v7.1.symbols.gmt")
c5_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/c5.all.v7.1.symbols.gmt")
c2_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/c2.all.v7.1.symbols.gmt")

h_GSVA_res <- gsva(matrix_1, h_geneSets, min.sz=5, max.sz=500, verbose=FALSE, parallel.sz=10)
c2_GSVA_res <- gsva(matrix_1, c2_geneSets, min.sz=5, max.sz=500, verbose=FALSE, parallel.sz=10)
c5_GSVA_res <- gsva(matrix_1, c5_geneSets, min.sz=5, max.sz=500, verbose=FALSE, parallel.sz=10)

#saveRDS(h_GSVA_res,file="IMpower133_All_GSVA_h_GSVA_res.rds")
#saveRDS(c2_GSVA_res,file="IMpower133_All_GSVA_c2_GSVA_res.rds")
#saveRDS(c5_GSVA_res,file="IMpower133_All_GSVA_c5_GSVA_res.rds")




matrix_1 <- as.matrix(Nature_evo_data)
h_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/h.all.v7.1.symbols.gmt")
c5_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/c5.all.v7.1.symbols.gmt")
c2_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/c2.all.v7.1.symbols.gmt")

h_GSVA_res <- gsva(matrix_1, h_geneSets, min.sz=5, max.sz=500, verbose=FALSE, parallel.sz=10)
c2_GSVA_res <- gsva(matrix_1, c2_geneSets, min.sz=5, max.sz=500, verbose=FALSE, parallel.sz=10)
c5_GSVA_res <- gsva(matrix_1, c5_geneSets, min.sz=5, max.sz=500, verbose=FALSE, parallel.sz=10)

#saveRDS(h_GSVA_res,file="Nature_evo_All_GSVA_h_GSVA_res.rds")
#saveRDS(c2_GSVA_res,file="Nature_evo_All_GSVA_c2_GSVA_res.rds")
#saveRDS(c5_GSVA_res,file="Nature_evo_All_GSVA_c5_GSVA_res.rds")

IMpower133_All_GSVA_c2 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/IMpower133_All_GSVA_c2_GSVA_res.rds",mc.cores=20)
IMpower133_All_GSVA_c5 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/IMpower133_All_GSVA_c5_GSVA_res.rds",mc.cores=20)
IMpower133_All_GSVA_h <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/IMpower133_All_GSVA_h_GSVA_res.rds",mc.cores=20)

Nature_evo_All_GSVA_c2 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Nature_evo_All_GSVA_c2_GSVA_res.rds",mc.cores=20)
Nature_evo_All_GSVA_c5 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Nature_evo_All_GSVA_c5_GSVA_res.rds",mc.cores=20)
Nature_evo_All_GSVA_h <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Nature_evo_All_GSVA_h_GSVA_res.rds",mc.cores=20)

sudobulk_GSVA_rds <- rbind(Nature_evo_All_GSVA_c5,Nature_evo_All_GSVA_c2,Nature_evo_All_GSVA_h)
write.csv(sudobulk_GSVA_rds,file="Nature_evo_GSVA.csv")

IMpower133_sudobulk_GSVA_rds <- rbind(IMpower133_All_GSVA_c5,IMpower133_All_GSVA_c2,IMpower133_All_GSVA_h)
write.csv(IMpower133_sudobulk_GSVA_rds,file="IMpower133_GSVA.csv")


Nature_evo_GSVA <- data.frame(t(sudobulk_GSVA_rds))
IMpoewr133_GSVA <- data.frame(t(IMpower133_sudobulk_GSVA_rds))

endo_pathway <- as.character(c(
"KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION",
"REACTOME_TIGHT_JUNCTION_INTERACTIONS",
"GO_POSITIVE_REGULATION_OF_VASCULAR_PERMEABILITY",
"GO_T_CELL_MEDIATED_IMMUNITY"))

endo_pathway <- intersect(colnames(IMpoewr133_GSVA),endo_pathway)
Nature_evo_sel <- Nature_evo_GSVA[,endo_pathway]
IMpoewr133_sel <- IMpoewr133_GSVA[,endo_pathway]
Nature_evo_sel$Group <- "Nature_evo"
IMpoewr133_sel$Group <- "IMpower133"
all_GSVA <- rbind(Nature_evo_sel,IMpoewr133_sel)
all_GSVA$ID <- rownames(all_GSVA)

T.cell_ratio <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/SCLC_T_ratio.csv")
tmp <- merge(all_GSVA,T.cell_ratio,by="ID")
tmp <- tmp[order(-tmp$T.cell, decreasing = FALSE), ]
tmp$order <- c(1:163)

Nature_evo_Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/Nature_SCLC_sample_info.csv")
IMPOWER133_Clinical <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/ICI_Cohort/IMPOWER133_Clinical_file.csv")

subtype <- IMPOWER133_Clinical[,c("trunc_anonymized_sample_ids","Gay_MDACC_subtypes")]
subtype$trunc_anonymized_sample_ids <- gsub("-", ".", subtype$trunc_anonymized_sample_ids)
tmp_1 <- merge(tmp, subtype, by.x= "ID", by.y="trunc_anonymized_sample_ids")

library(ggpubr)
library(ggplot2)
library(BuenColors)
library(Seurat)

tmp_2 <- tmp_1[,c("GO_POSITIVE_REGULATION_OF_VASCULAR_PERMEABILITY","T.cell","Gay_MDACC_subtypes")]
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
tmp_2$order <- c(1:132)

tmp_2$Gay_MDACC_subtypes <- substring(tmp_2$Gay_MDACC_subtypes,5,6)
tmp_2$Gay_MDACC_subtypes <- as.factor(tmp_2$Gay_MDACC_subtypes)
tumor_levels <- levels(tmp_2$Gay_MDACC_subtypes)
y_positions <- seq(from = min(tmp_2$GO_POSITIVE_REGULATION_OF_VASCULAR_PERMEABILITY) - 0.1, 
                   by = -0.05, length.out = length(tumor_levels))
tumor_y_mapping <- data.frame(Gay_MDACC_subtypes = tumor_levels, y_pos = y_positions)
tmp_2 <- merge(tmp_2, tumor_y_mapping, by = "Gay_MDACC_subtypes")
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
p1 <- ggplot(tmp_2, aes(order, GO_POSITIVE_REGULATION_OF_VASCULAR_PERMEABILITY)) +
  geom_point(alpha = 0.5, size = 0.5, colour = "#969696") +
  scale_color_gradientn(colours = c("#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#d73027")) +
  geom_rug(alpha = 0.8, position = "jitter", sides = "b") +
  geom_smooth(colour = "orange", se = TRUE) +
  theme_classic() +
  xlab("Tcell_high_to_low") +
  labs(title = "GO_POSITIVE_REGULATION_OF_VASCULAR_PERMEABILITY") +
  NoLegend()
p2 <- p1 + geom_segment(data = tmp_2, aes(x = order, xend = order, y = y_pos, yend = y_pos + 0.02, color = Gay_MDACC_subtypes),inherit.aes = FALSE) + 
  scale_color_manual(values = c("ASCL1_hi" = "#a50f15", "ASCL1_low" ="#000000")) + 
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(), 
        legend.position = "bottom") + 
  geom_text(data = tumor_y_mapping, aes(x = max(tmp_2$order) + 10,  y = y_pos + 0.01, label = Gay_MDACC_subtypes, color = Gay_MDACC_subtypes),inherit.aes = FALSE, hjust = 0, size = 3)
ggsave(p2,file="Fig1D_1.svg",width=4,height=5)



tmp_2 <- tmp_1[,c("GO_T_CELL_MEDIATED_IMMUNITY","T.cell","Gay_MDACC_subtypes")]
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
tmp_2$order <- c(1:132)
tmp_2$Gay_MDACC_subtypes <- substring(tmp_2$Gay_MDACC_subtypes,5,6)
tmp_2$Gay_MDACC_subtypes <- as.factor(tmp_2$Gay_MDACC_subtypes)
tumor_levels <- levels(tmp_2$Gay_MDACC_subtypes)
y_positions <- seq(from = min(tmp_2$GO_T_CELL_MEDIATED_IMMUNITY) - 0.1, 
                   by = -0.05, length.out = length(tumor_levels))
tumor_y_mapping <- data.frame(Gay_MDACC_subtypes = tumor_levels, y_pos = y_positions)
tmp_2 <- merge(tmp_2, tumor_y_mapping, by = "Gay_MDACC_subtypes")
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
p1 <- ggplot(tmp_2, aes(order, GO_T_CELL_MEDIATED_IMMUNITY)) +
  geom_point(alpha = 0.5, size = 0.5, colour = "#969696") +
  scale_color_gradientn(colours = c("#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#d73027")) +
  geom_rug(alpha = 0.8, position = "jitter", sides = "b") +
  geom_smooth(colour = "orange", se = TRUE) +
  theme_classic() +
  xlab("Tcell_high_to_low") +
  labs(title = "GO_T_CELL_MEDIATED_IMMUNITY") +
  NoLegend()
p2 <- p1 + geom_segment(data = tmp_2, aes(x = order, xend = order, y = y_pos, yend = y_pos + 0.02, color = Gay_MDACC_subtypes),inherit.aes = FALSE) + 
  scale_color_manual(values = c("ASCL1_hi" = "#a50f15", "ASCL1_low" ="#000000")) + 
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(), 
        legend.position = "bottom") + 
  geom_text(data = tumor_y_mapping, aes(x = max(tmp_2$order) + 10,  y = y_pos + 0.01, label = Gay_MDACC_subtypes, color = Gay_MDACC_subtypes),inherit.aes = FALSE, hjust = 0, size = 3)
ggsave(p2,file="Fig1D_2.png",width=4,height=5)


tmp_2 <- tmp_1[,c("KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION","T.cell","Gay_MDACC_subtypes")]
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
tmp_2$order <- c(1:132)
tmp_2$Gay_MDACC_subtypes <- substring(tmp_2$Gay_MDACC_subtypes,5,6)
tmp_2$Gay_MDACC_subtypes <- as.factor(tmp_2$Gay_MDACC_subtypes)
tumor_levels <- levels(tmp_2$Gay_MDACC_subtypes)
y_positions <- seq(from = min(tmp_2$KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION) - 0.1, 
                   by = -0.05, length.out = length(tumor_levels))
tumor_y_mapping <- data.frame(Gay_MDACC_subtypes = tumor_levels, y_pos = y_positions)
tmp_2 <- merge(tmp_2, tumor_y_mapping, by = "Gay_MDACC_subtypes")
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
p1 <- ggplot(tmp_2, aes(order, KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION)) +
  geom_point(alpha = 0.5, size = 0.5, colour = "#969696") +
  scale_color_gradientn(colours = c("#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#d73027")) +
  geom_rug(alpha = 0.8, position = "jitter", sides = "b") +
  geom_smooth(colour = "orange", se = TRUE) +
  theme_classic() +
  xlab("Tcell_high_to_low") +
  labs(title = "KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION") +
  NoLegend()
p2 <- p1 + geom_segment(data = tmp_2, aes(x = order, xend = order, y = y_pos, yend = y_pos + 0.02, color = Gay_MDACC_subtypes),inherit.aes = FALSE) + 
  scale_color_manual(values = c("ASCL1_hi" = "#a50f15", "ASCL1_low" ="#000000")) + 
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(), 
        legend.position = "bottom") + 
  geom_text(data = tumor_y_mapping, aes(x = max(tmp_2$order) + 10,  y = y_pos + 0.01, label = Gay_MDACC_subtypes, color = Gay_MDACC_subtypes),inherit.aes = FALSE, hjust = 0, size = 3)
ggsave(p2,file="Fig1D_3.png",width=4,height=5)



tmp_2 <- tmp_1[,c("REACTOME_TIGHT_JUNCTION_INTERACTIONS","T.cell","Gay_MDACC_subtypes")]
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
tmp_2$order <- c(1:132)
tmp_2$Gay_MDACC_subtypes <- substring(tmp_2$Gay_MDACC_subtypes,5,6)
tmp_2$Gay_MDACC_subtypes <- as.factor(tmp_2$Gay_MDACC_subtypes)
tumor_levels <- levels(tmp_2$Gay_MDACC_subtypes)
y_positions <- seq(from = min(tmp_2$REACTOME_TIGHT_JUNCTION_INTERACTIONS) - 0.1, 
                   by = -0.05, length.out = length(tumor_levels))
tumor_y_mapping <- data.frame(Gay_MDACC_subtypes = tumor_levels, y_pos = y_positions)
tmp_2 <- merge(tmp_2, tumor_y_mapping, by = "Gay_MDACC_subtypes")
tmp_2 <- tmp_2[order(-tmp_2$T.cell, decreasing = F),]
p1 <- ggplot(tmp_2, aes(order, REACTOME_TIGHT_JUNCTION_INTERACTIONS)) +
  geom_point(alpha = 0.5, size = 0.5, colour = "#969696") +
  scale_color_gradientn(colours = c("#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#d73027")) +
  geom_rug(alpha = 0.8, position = "jitter", sides = "b") +
  geom_smooth(colour = "orange", se = TRUE) +
  theme_classic() +
  xlab("Tcell_high_to_low") +
  labs(title = "REACTOME_TIGHT_JUNCTION_INTERACTIONS") +
  NoLegend()
p2 <- p1 + geom_segment(data = tmp_2, aes(x = order, xend = order, y = y_pos, yend = y_pos + 0.02, color = Gay_MDACC_subtypes),inherit.aes = FALSE) + 
  scale_color_manual(values = c("ASCL1_hi" = "#a50f15", "ASCL1_low" ="#000000")) + 
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank(), 
        legend.position = "bottom") + 
  geom_text(data = tumor_y_mapping, aes(x = max(tmp_2$order) + 10,  y = y_pos + 0.01, label = Gay_MDACC_subtypes, color = Gay_MDACC_subtypes),inherit.aes = FALSE, hjust = 0, size = 3)
ggsave(p2,file="Fig1D_4.png",width=4,height=5)



#Fig 1E

#Convert("/mnt/data/user_data/xiangyu/workshop/DATABASE/HTAN/level4/MSK/combined.mnnc/adata.combined.mnnc.010920.h5ad", dest = "h5seurat", overwrite = TRUE)
All_Combine <- LoadH5Seurat("/mnt/data/user_data/xiangyu/workshop/DATABASE/HTAN/level4/MSK/Convert_To_Seurat_Object/adata.combined.mnnc.010920.h5seurat")
Idents(All_Combine) <- All_Combine$cell_type_fine

Bcell <- subset(All_Combine,idents=c("B cell"),invert=FALSE)
Idents(Bcell) <- Bcell$histo
Bcell <- subset(Bcell,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Bcell, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Bcell_DEGs.csv")


Fibroblast <- subset(All_Combine,idents=c("Fibroblast"),invert=FALSE)
Idents(Fibroblast) <- Fibroblast$histo
Fibroblast <- subset(Fibroblast,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Fibroblast, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Fibroblast_DEGs.csv")

Macrophage <- subset(All_Combine,idents=c("Macrophage"),invert=FALSE)
Idents(Macrophage) <- Macrophage$histo
Macrophage <- subset(Macrophage,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Macrophage, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Macrophage_DEGs.csv")
             
Neutrophil <- subset(All_Combine,idents=c("Neutrophil"),invert=FALSE)
Idents(Neutrophil) <- Neutrophil$histo
Neutrophil <- subset(Neutrophil,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Neutrophil, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Neutrophil_DEGs.csv")
     
Plasma <- subset(All_Combine,idents=c("Plasma cell"),invert=FALSE)
Idents(Plasma) <- Plasma$histo
Plasma <- subset(Plasma,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Plasma, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Plasma_DEGs.csv")
        
Tcell <- subset(All_Combine,idents=c("T cell"),invert=FALSE)
Idents(Tcell) <- Tcell$histo
Tcell <- subset(Tcell,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Tcell_DEGs.csv")


Tcell <- subset(All_Combine,idents=c("T cell"),invert=FALSE)
Idents(Tcell) <- Tcell$histo
Tcell <- subset(Tcell,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Tcell, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Tcell_DEGs.csv")


Neutrophil <- subset(All_Combine,idents=c("Neutrophil"),invert=FALSE)
Idents(Neutrophil) <- Neutrophil$histo
Neutrophil <- subset(Neutrophil,idents=c("normal"),invert=TRUE)
all.markers <- FindAllMarkers(Neutrophil, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"Neutrophil_DEGs.csv")



LUAD_Tcell <- subset(Tcell_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Tcell <- subset(Tcell_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Tcell_human <- data.frame(DEG_num=c(length(LUAD_Tcell$gene),length(SCLC_Tcell$gene)),group=c("LUAD.Tcell","SCLC.Tcell"),cell_type="Tcell")

LUAD_Bcell <- subset(Bcell_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Bcell <- subset(Bcell_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Bcell_human <- data.frame(DEG_num=c(length(LUAD_Bcell$gene),length(SCLC_Bcell$gene)),group=c("LUAD.Bcell","SCLC.Bcell"),cell_type="Bcell")


LUAD_Plasma <- subset(Plasma_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Plasma <- subset(Plasma_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Plasma_human <- data.frame(DEG_num=c(length(LUAD_Plasma$gene),length(SCLC_Plasma$gene)),group=c("LUAD.Plasma","SCLC.Plasma"),cell_type="Plasma")


LUAD_Macro <- subset(Macro_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Macro <- subset(Macro_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Macro_human <- data.frame(DEG_num=c(length(LUAD_Macro$gene),length(SCLC_Macro$gene)),group=c("LUAD.Macro","SCLC.Macro"),cell_type="Macro")


LUAD_Neutro <- subset(Neutro_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Neutro <- subset(Neutro_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Neutro_human <- data.frame(DEG_num=c(length(LUAD_Neutro$gene),length(SCLC_Neutro$gene)),group=c("LUAD.Neutro","SCLC.Neutro"),cell_type="Neutro")

LUAD_Tumor <- subset(Tumor_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Tumor <- subset(Tumor_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Tumor_human <- data.frame(DEG_num=c(length(LUAD_Tumor$gene),length(SCLC_Tumor$gene)),group=c("LUAD.Tumor","SCLC.Tumor"),cell_type="Tumor")

LUAD_Fibro <- subset(Fibro_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Fibro <- subset(Fibro_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Fibro_human <- data.frame(DEG_num=c(length(LUAD_Fibro$gene),length(SCLC_Fibro$gene)),group=c("LUAD.Fibro","SCLC.Fibro"),cell_type="Fibro")

LUAD_Endo <- subset(Endo_DEGs,cluster=="LUAD" & avg_logFC > 0.1 & p_val_adj < 0.01)
SCLC_Endo <- subset(Endo_DEGs,cluster=="SCLC" & avg_logFC > 0.1 & p_val_adj < 0.01)
Endo_human <- data.frame(DEG_num=c(length(LUAD_Endo$gene),length(SCLC_Endo$gene)),group=c("LUAD.Endo","SCLC.Endo"),cell_type="Endo")

all_human <- rbind(Tcell_human,Bcell_human,Plasma_human,Macro_human,Neutro_human,Fibro_human,Endo_human)
all_human$histo <- substring(all_human$group,1,4)

library(dplyr)

df_plot <- all_human %>%
  mutate(
    value = ifelse(histo == "LUAD", -DEG_num, DEG_num),
    histo = factor(histo, levels = c("LUAD", "SCLC"))
  )

library(ggplot2)

ff <- ggplot(df_plot, aes(x = value, y = cell_type, fill = histo)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_vline(xintercept = 0, color = "black") +
 scale_x_continuous(
  breaks = seq(-800, 800, by = 400),
  labels = abs
) +
  scale_fill_manual(
    values = c("LUAD" = "#2C5AA0", "SCLC" = "#C83E3A")
  ) +
  labs(
    x = "The number of DEGs",
    y = NULL,
    fill = NULL
  ) +

  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 13, face = "bold")
  )
ggsave(ff,file="Fig1E.png",width=8,height=5)





#Fig1F

Endo_pure <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/Endothelium/Pure_Endo_CancerCell.rds",mc.cores=20)
Idents(Endo_pure) <- Endo_pure$histo

LUAD <- subset(Endo_pure,ident="LUAD")
SCLC <- subset(Endo_pure,ident="SCLC")
normal <- subset(Endo_pure,ident="normal")

pseudo_bulk_seurat_mean <- function(seurat_obj=seurat_obj,num_split=num_split,seed.use=seed.use,prefix=prefix,slot=slot){
  set.seed(seed.use)
  require(Seurat)
  genes.use <- rownames(seurat_obj)
  cell.sets1 <- split(colnames(seurat_obj), sort(1:length(colnames(seurat_obj))%%num_split))
  profile.set1 = matrix(, nrow = length(genes.use), ncol = length(cell.sets1))
  for (i in 1:length(cell.sets1)) {
    this.set <- cell.sets1[[i]]
    sub.matrix <- GetAssayData(seurat_obj, slot = slot, assay = "RNA")[genes.use, this.set]
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

LUAD_pseudobulk <- pseudo_bulk_seurat_mean(seurat_obj=LUAD,num_split=10,seed.use=1,slot="data",prefix="LUAD")
normal_pseudobulk <- pseudo_bulk_seurat_mean(seurat_obj=normal,num_split=10,seed.use=1,slot="data",prefix="normal")
SCLC_pseudobulk <- pseudo_bulk_seurat_mean(seurat_obj=SCLC,num_split=10,seed.use=1,slot="data",prefix="SCLC")
All_pseudobulk <- cbind(normal_pseudobulk,LUAD_pseudobulk,SCLC_pseudobulk)
mcsaveRDS(All_pseudobulk,file="All_pseudobulk.rds",mc.cores=30)


h_geneSets <- getGmt("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/GSVA_7.1/h.all.v7.1.symbols.gmt")
c5_geneSets <- getGmt("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/GSVA_7.1/c5.all.v7.1.symbols.gmt")
c2_geneSets <- getGmt("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/GSVA_7.1/c2.all.v7.1.symbols.gmt")

matrix <- data.frame(All_pseudobulk)
matrix <- as.matrix(matrix)

h_gsva_param <- gsvaParam(expr = matrix, geneSets = h_geneSets, minSize = 5, maxSize = 500)
h_GSVA_res <- gsva(h_gsva_param, verbose = FALSE)
c2_gsva_param <- gsvaParam(expr = matrix, geneSets = c2_geneSets, minSize = 5, maxSize = 500)
c2_GSVA_res <- gsva(c2_gsva_param, verbose = FALSE)
c5_gsva_param <- gsvaParam(expr = matrix, geneSets = c5_geneSets, minSize = 5, maxSize = 500)
c5_GSVA_res <- gsva(c5_gsva_param, verbose = FALSE)

sudobulk_GSVA_score <- rbind(h_GSVA_res,c2_GSVA_res,c5_GSVA_res)
write.csv(sudobulk_GSVA_score,"Endo_AllCells_GSVA_score.csv")


annotation0 <- data.frame(c(1:10))
annotation0$group <- "normal"
annotation1 <- data.frame(c(11:20))
annotation1$group <- "LUAD"
annotation2 <- data.frame(c(21:30))
annotation2$group <- "SCLC"
names(annotation2) <- c("order","group")
names(annotation1) <- c("order","group")
names(annotation0) <- c("order","group")
annotation <- rbind(annotation0,annotation1,annotation2)

group <-as.factor(annotation$group)
design <- model.matrix(~ group + 0)
rownames(design)<-colnames(matrix)
head(design)

contrasts <- makeContrasts(groupLUAD-(groupnormal+groupSCLC)/2,groupnormal-(groupLUAD+groupSCLC)/2,
	groupSCLC-(groupLUAD+groupnormal)/2,
	levels=design)

fiT <- lmFit(h_GSVA_res, design)
fiT2 <- contrasts.fit(fiT, contrasts)
h_fiT3 <- eBayes(fiT2)

fiT <- lmFit(c2_GSVA_res, design)
fiT2 <- contrasts.fit(fiT, contrasts)
c2_fiT3 <- eBayes(fiT2)

fiT <- lmFit(c5_GSVA_res, design)
fiT2 <- contrasts.fit(fiT, contrasts)
c5_fiT3 <- eBayes(fiT2)

GSVA_h <- topTable(h_fiT3, number=1000,p.value=0.05, adjust="BH")
GSVA_c2 <- topTable(c2_fiT3, number=1000,p.value=0.05, adjust="BH")
GSVA_c5 <- topTable(c5_fiT3, number=1000,p.value=0.05, adjust="BH")
GSVA_hc2c5 <- rbind(GSVA_h,GSVA_c2,GSVA_c5)
write.csv(GSVA_hc2c5,"Endothelium_CancerCell_p005_AllCells_GSVA_score.csv")




GSVA_pathway <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/Endothelium/Endo_AllCells_GSVA_score.csv")
rownames(GSVA_pathway) <- GSVA_pathway$X
select_pathways <- subset(GSVA_pathway,
X=="BIOCARTA_IGF1R_PATHWAY" |
X=="BIOCARTA_ECM_PATHWAY" |
X=="REACTOME_SIGNALING_BY_PDGFR_IN_DISEASE" |
X=="PID_IGF1_PATHWAY" |
X=="GO_MAP_KINASE_ACTIVITY" |
X=="BIOCARTA_IGF1MTOR_PATHWAY" |
X=="GO_REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY" |
X=="NABA_BASEMENT_MEMBRANES" )

select_pathways <- data.frame(select_pathways)
select_pathways <- select_pathways[,-1]

select_pathways <- t(apply(select_pathways, 1, function(x) (x-mean(x))/sd(x)))

range(select_pathways)
library(pheatmap)
select_pathways[select_pathways > 1] <- 1
select_pathways[select_pathways < -1] <- -1

library(pheatmap)
library(dplyr)

pdf("Fig1F_heatmap.pdf",width=15,height=4)
pheatmap(select_pathways,clustering_method="ward.D2",
  color = colorRampPalette(c("#2971B1","#6AACD0","#C1DDEB","#F7F7F7","#FACDB5","#E58267","#BB2933"))(50),
  gaps_col = c(10,20),
  fontsize_row=12,show_rownames=TRUE,show_colnames=TRUE,cluster_row =FALSE,cluster_col= FALSE,border=FALSE)
dev.off()



#Fig S1A

#Convert("/mnt/data/user_data/xiangyu/workshop/DATABASE/HTAN/level4/MSK/combined.mnnc/adata.combined.mnnc.010920.h5ad", dest = "h5seurat", overwrite = TRUE)

All_Combine <- LoadH5Seurat("/mnt/data/user_data/xiangyu/workshop/DATABASE/HTAN/level4/MSK/Convert_To_Seurat_Object/adata.combined.mnnc.010920.h5seurat")
tmp <- All_Combine[[]]
LUAD <- subset(tmp,histo=="LUAD")
SCLC <- subset(tmp,histo=="SCLC")

LUAD_Epi <- subset(LUAD,cell_type_fine=="NSCLC" | cell_type_fine=="SCLC-A" | cell_type_fine=="SCLC-N" | cell_type_fine=="SCLC-P" | cell_type_fine=="AE1" |  cell_type_fine=="AEP" | cell_type_fine=="Basal" | cell_type_fine=="Ciliated" | cell_type_fine=="Hepatocyte" | cell_type_fine=="Ionocyte" |cell_type_fine=="Mucinous" |cell_type_fine=="Neuroendocrine" |  cell_type_fine=="Club" |cell_type_fine=="Tuft")          
LUAD_Fibro <- subset(LUAD,cell_type_fine=="Fibroblast" | cell_type_fine=="Pericyte")        
LUAD_Macro <- subset(LUAD,cell_type_fine=="Macrophage")        
LUAD_Mast <- subset(LUAD,cell_type_fine=="Mast")        
LUAD_Neutrophil <- subset(LUAD,cell_type_fine=="Neutrophil")        
LUAD_Plasma <- subset(LUAD,cell_type_fine=="Plasma cell")        
LUAD_Erythrocyte <- subset(LUAD,cell_type_fine=="Erythrocyte")        
LUAD_Endothelial <- subset(LUAD,cell_type_fine=="Endothelial")        
LUAD_DC <- subset(LUAD,cell_type_fine=="DC")        
LUAD_Tcell <- subset(LUAD,cell_type_fine=="T cell")        
LUAD_Bcell <- subset(LUAD,cell_type_fine=="B cell")        

LUAD_Epi$cell_type_fine_1 <- "Epithelial"
LUAD_Fibro$cell_type_fine_1 <- "Fibro"
LUAD_Endothelial$cell_type_fine_1 <- "Endo"
LUAD_Macro$cell_type_fine_1 <- "Macro"
LUAD_Mast$cell_type_fine_1 <- "Mast"
LUAD_Neutrophil$cell_type_fine_1 <- "Neutro"
LUAD_Plasma$cell_type_fine_1 <- "Plasma"
LUAD_Erythrocyte$cell_type_fine_1 <- "Erythroid"
LUAD_DC$cell_type_fine_1 <- "DC"
LUAD_Tcell$cell_type_fine_1 <- "Tcell"
LUAD_Bcell$cell_type_fine_1 <- "Bcell"
LUAD_tmp <- rbind(LUAD_Epi,LUAD_Fibro,LUAD_Endothelial,LUAD_Macro,LUAD_Mast,LUAD_Neutrophil,LUAD_Plasma,LUAD_Erythrocyte,LUAD_DC,LUAD_Tcell,LUAD_Bcell)
LUAD_cellnumber <- data.frame(table(LUAD_tmp$cell_type_fine_1))
LUAD_cellnumber$All_cells <- "57679"
LUAD_cellnumber$Ratio <- as.numeric(LUAD_cellnumber$Freq)/as.numeric(LUAD_cellnumber$All_cells)
LUAD_cellnumber$type <- "LUAD"


SCLC_Epi <- subset(SCLC,cell_type_fine=="NSCLC" | cell_type_fine=="SCLC-A" | cell_type_fine=="SCLC-N" | cell_type_fine=="SCLC-P" | cell_type_fine=="AE1" |  cell_type_fine=="AEP" | cell_type_fine=="Basal" | cell_type_fine=="Ciliated" | cell_type_fine=="Hepatocyte" | cell_type_fine=="Ionocyte" |cell_type_fine=="Mucinous" |cell_type_fine=="Neuroendocrine" |  cell_type_fine=="Club" |cell_type_fine=="Tuft")          
SCLC_Fibro <- subset(SCLC,cell_type_fine=="Fibroblast" | cell_type_fine=="Pericyte")        
SCLC_Macro <- subset(SCLC,cell_type_fine=="Macrophage")        
SCLC_Mast <- subset(SCLC,cell_type_fine=="Mast")        
SCLC_Neutrophil <- subset(SCLC,cell_type_fine=="Neutrophil")        
SCLC_Plasma <- subset(SCLC,cell_type_fine=="Plasma cell")        
SCLC_Endothelial <- subset(SCLC,cell_type_fine=="Endothelial")        
SCLC_DC <- subset(SCLC,cell_type_fine=="DC")        
SCLC_Tcell <- subset(SCLC,cell_type_fine=="T cell")        
SCLC_Bcell <- subset(SCLC,cell_type_fine=="B cell")        
SCLC_Epi$cell_type_fine_1 <- "Epithelial"
SCLC_Fibro$cell_type_fine_1 <- "Fibro"
SCLC_Endothelial$cell_type_fine_1 <- "Endo"
SCLC_Macro$cell_type_fine_1 <- "Macro"
SCLC_Mast$cell_type_fine_1 <- "Mast"
SCLC_Neutrophil$cell_type_fine_1 <- "Neutro"
SCLC_Plasma$cell_type_fine_1 <- "Plasma"
SCLC_DC$cell_type_fine_1 <- "DC"
SCLC_Tcell$cell_type_fine_1 <- "Tcell"
SCLC_Bcell$cell_type_fine_1 <- "Bcell"

SCLC_tmp <- rbind(SCLC_Epi,SCLC_Fibro,SCLC_Endothelial,SCLC_Macro,SCLC_Mast,SCLC_Neutrophil,SCLC_Plasma,SCLC_DC,SCLC_Tcell,SCLC_Bcell)
SCLC_cellnumber <- data.frame(table(SCLC_tmp$cell_type_fine_1))
SCLC_Erythrocyte <- data.frame(Var1="Erythroid",Freq="0")
SCLC_cellnumber <- rbind(SCLC_cellnumber,SCLC_Erythrocyte)
SCLC_cellnumber$All_cells <- "77143"
SCLC_cellnumber$Ratio <- as.numeric(SCLC_cellnumber$Freq)/as.numeric(SCLC_cellnumber$All_cells)
SCLC_cellnumber$type <- "SCLC"
tmp_data <- rbind(LUAD_cellnumber,SCLC_cellnumber)
tmp_data$type <- factor(tmp_data$type,levels=c('LUAD',"SCLC"))
tmp_data$Var1 <- factor(tmp_data$Var1,levels=c("Tcell","Bcell","Plasma","Neutro","Macro","Mast","DC","Erythroid","Endo","Fibro","Epithelial"))
library(ggalluvial)
ff <- ggplot(tmp_data, aes( x = type,y=Ratio,fill = Var1, stratum = Var1, alluvium = Var1))+
  geom_col(position = 'stack', width = 0.6)+
  #geom_bar(position = "stack", stat = "identity", width = 0.6) 
  theme_bw() + scale_fill_manual(values=c("#cb181d","#2166ac","#67a9cf","#78c679","#a6bddb","#99d8c9","#d6604d","#dadaeb","#f1b6da","#f4a582","#737373")) +
  scale_y_continuous(expand = c(0,0)) + geom_flow(alpha = 0.5) #绘制同类别之间的连接线
ggsave(ff,file="FigS1A_2.png",width=5,height=6)




use_python("/usr/bin/python3")
CancerCell_Rudin_All_Combine <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/CancerCell_Rudin_All_Combine.rds",mc.cores=20)
CancerCell_Rudin_OpenTSNE30 <- DoopenTSNE(CancerCell_Rudin_All_Combine, reduction_use = "pca", reduction_save = "openTSNE",dims_use = 1:30)   
Idents(CancerCell_Rudin_OpenTSNE30) <- CancerCell_Rudin_OpenTSNE30$histo
Rudin_NSCLC_SCLC <- subset(CancerCell_Rudin_OpenTSNE30,idents="normal",invert=TRUE)

tmp <- Rudin_NSCLC_SCLC[[]]
tmp_Epi <- subset(tmp,cell_type_fine=="NSCLC" | cell_type_fine=="SCLC-A" | cell_type_fine=="SCLC-N" | cell_type_fine=="SCLC-P" | cell_type_fine=="AE1" |  cell_type_fine=="AEP" | cell_type_fine=="Basal" | cell_type_fine=="Ciliated" | cell_type_fine=="Hepatocyte" | cell_type_fine=="Ionocyte" |cell_type_fine=="Mucinous" |cell_type_fine=="Neuroendocrine" |  cell_type_fine=="Club" |cell_type_fine=="Tuft")          
tmp_Fibro <- subset(tmp,cell_type_fine=="Fibroblast" | cell_type_fine=="Pericyte")        
tmp_Macro <- subset(tmp,cell_type_fine=="Macrophage")        
tmp_Mast <- subset(tmp,cell_type_fine=="Mast")        
tmp_Neutrophil <- subset(tmp,cell_type_fine=="Neutrophil")        
tmp_Plasma <- subset(tmp,cell_type_fine=="Plasma cell")        
tmp_Erythrocyte <- subset(tmp,cell_type_fine=="Erythrocyte")        
tmp_Endothelial <- subset(tmp,cell_type_fine=="Endothelial")        
tmp_DC <- subset(tmp,cell_type_fine=="DC")        
tmp_Tcell <- subset(tmp,cell_type_fine=="T cell")        
tmp_Bcell <- subset(tmp,cell_type_fine=="B cell")      
tmp_Epi$abao_cell_type <- "Epithelial"
tmp_Fibro$abao_cell_type <- "Fibro"
tmp_Endothelial$abao_cell_type <- "Endo"
tmp_Macro$abao_cell_type <- "Macro"
tmp_Mast$abao_cell_type <- "Mast"
tmp_Neutrophil$abao_cell_type <- "Neutro"
tmp_Plasma$abao_cell_type <- "Plasma"
tmp_DC$abao_cell_type <- "DC"
tmp_Tcell$abao_cell_type <- "Tcell"
tmp_Bcell$abao_cell_type <- "Bcell"
tmp_Erythrocyte$abao_cell_type <- "Erythrocyte"
All_tmp <- rbind(tmp_Epi,tmp_Fibro,tmp_Endothelial,tmp_Macro,tmp_Mast,tmp_Neutrophil,tmp_Plasma,tmp_DC,tmp_Tcell,tmp_Bcell,tmp_Erythrocyte)
All_tmp <- All_tmp[rownames(tmp),]
Rudin_NSCLC_SCLC$abao_cell_type <- All_tmp$abao_cell_type
Rudin_NSCLC_SCLC$abao_cell_type <- factor(Rudin_NSCLC_SCLC$abao_cell_type,levels=c("Tcell","Bcell","Plasma","Neutro","Macro","Mast","DC","Erythrocyte","Endo","Fibro","Epithelial"))
p2 <- DimPlot(object = Rudin_NSCLC_SCLC, reduction = "openTSNE",label=TRUE, pt.size = 0.1,label.size=0,raster = FALSE,
    group.by="abao_cell_type",cols=c("#cb181d","#2166ac","#67a9cf","#78c679","#a6bddb","#99d8c9","#d6604d","#dadaeb","#f1b6da","#f4a582","#737373")) +labs(title="openTSNE")
ggsave(p2,file="FigS1A_1.png",width=5.8,height=5.0,dpi=1080)




#Fig S1B

tmp <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/CancerCell_cell_type_fine.csv")
LUAD <- subset(tmp,histo=="LUAD")
SCLC <- subset(tmp,histo=="SCLC")
normal <- subset(tmp,histo=="normal")

LUAD_all_cells <- data.frame(table(LUAD$patient))
LUAD_all_cells <- subset(LUAD_all_cells,Freq>0)
names(LUAD_all_cells) <- c("ID","All_cells")

LUAD_Tcells <- subset(LUAD,cell_type_fine=="T cell")
LUAD_Tcells_id <- data.frame(table(LUAD_Tcells$patient))
names(LUAD_Tcells_id) <- c("ID","Tcells")
rownames(LUAD_Tcells_id) <- LUAD_Tcells_id$ID
LUAD_Tcells_id <- LUAD_Tcells_id[LUAD_all_cells$ID,]

LUAD_Macrophage <- subset(LUAD,cell_type_fine=="Macrophage")
LUAD_Macrophage_id <- data.frame(table(LUAD_Macrophage$patient))
names(LUAD_Macrophage_id) <- c("ID","Macrophage")
rownames(LUAD_Macrophage_id) <- LUAD_Macrophage_id$ID
LUAD_Macrophage_id <- LUAD_Macrophage_id[LUAD_all_cells$ID,]

LUAD_Mast <- subset(LUAD,cell_type_fine=="Mast")
LUAD_Mast_id <- data.frame(table(LUAD_Mast$patient))
names(LUAD_Mast_id) <- c("ID","Mast")
rownames(LUAD_Mast_id) <- LUAD_Mast_id$ID
LUAD_Mast_id <- LUAD_Mast_id[LUAD_all_cells$ID,]

LUAD_Bcell <- subset(LUAD,cell_type_fine=="B cell")
LUAD_Bcell_id <- data.frame(table(LUAD_Bcell$patient))
names(LUAD_Bcell_id) <- c("ID","Bcell")
rownames(LUAD_Bcell_id) <- LUAD_Bcell_id$ID
LUAD_Bcell_id <- LUAD_Bcell_id[LUAD_all_cells$ID,]


SCLC_all_cells <- data.frame(table(SCLC$patient))
SCLC_all_cells <- subset(SCLC_all_cells,Freq>0)
names(SCLC_all_cells) <- c("ID","All_cells")

SCLC_Tcells <- subset(SCLC,cell_type_fine=="T cell")
SCLC_Tcells_id <- data.frame(table(SCLC_Tcells$patient))
names(SCLC_Tcells_id) <- c("ID","Tcells")
rownames(SCLC_Tcells_id) <- SCLC_Tcells_id$ID
SCLC_Tcells_id <- SCLC_Tcells_id[SCLC_all_cells$ID,]

SCLC_Macrophage <- subset(SCLC,cell_type_fine=="Macrophage")
SCLC_Macrophage_id <- data.frame(table(SCLC_Macrophage$patient))
names(SCLC_Macrophage_id) <- c("ID","Macrophage")
rownames(SCLC_Macrophage_id) <- SCLC_Macrophage_id$ID
SCLC_Macrophage_id <- SCLC_Macrophage_id[SCLC_all_cells$ID,]

SCLC_Mast <- subset(SCLC,cell_type_fine=="Mast")
SCLC_Mast_id <- data.frame(table(SCLC_Mast$patient))
names(SCLC_Mast_id) <- c("ID","Mast")
rownames(SCLC_Mast_id) <- SCLC_Mast_id$ID
SCLC_Mast_id <- SCLC_Mast_id[SCLC_all_cells$ID,]

SCLC_Bcell <- subset(SCLC,cell_type_fine=="B cell")
SCLC_Bcell_id <- data.frame(table(SCLC_Bcell$patient))
names(SCLC_Bcell_id) <- c("ID","Bcell")
rownames(SCLC_Bcell_id) <- SCLC_Bcell_id$ID
SCLC_Bcell_id <- SCLC_Bcell_id[SCLC_all_cells$ID,]



LUAD_T_ratio <- data.frame(LUAD_Tcells_id$Tcells/LUAD_all_cells$All_cells)
names(LUAD_T_ratio) <- "T_ratio"
LUAD_T_ratio$type <- "LUAD_Tcell"
SCLC_T_ratio <- data.frame(SCLC_Tcells_id$Tcells/SCLC_all_cells$All_cells)
names(SCLC_T_ratio) <- "T_ratio"
SCLC_T_ratio$type <- "SCLC_Tcell"
Tcells_ratio <- rbind(LUAD_T_ratio,SCLC_T_ratio)
Tcells_ratio$type <- factor(Tcells_ratio$type,levels=c("LUAD_Tcell","SCLC_Tcell"))
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("SCLC_Tcell","LUAD_Tcell"))
ff1 <- ggboxplot(Tcells_ratio, x = "type", y = "T_ratio",add = "jitter",
               color = "type",palette = c("#000000","#a50026")) +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="t.test")


LUAD_Bcell_ratio <- data.frame(LUAD_Bcell_id$Bcell/LUAD_all_cells$All_cells)
names(LUAD_Bcell_ratio) <- "Bcell_ratio"
LUAD_Bcell_ratio$type <- "LUAD_Bcell"
SCLC_Bcell_ratio <- data.frame(SCLC_Bcell_id$Bcell/SCLC_all_cells$All_cells)
names(SCLC_Bcell_ratio) <- "Bcell_ratio"
SCLC_Bcell_ratio$type <- "SCLC_Bcell"
Bcell_ratio <- rbind(LUAD_Bcell_ratio,SCLC_Bcell_ratio)
Bcell_ratio$type <- factor(Bcell_ratio$type,levels=c("LUAD_Bcell","SCLC_Bcell"))
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("SCLC_Bcell","LUAD_Bcell"))
ff2 <- ggboxplot(Bcell_ratio, x = "type", y = "Bcell_ratio",add = "jitter",
               color = "type",palette = c("#000000","#a50026")) +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="t.test")


LUAD_Macrophage_ratio <- data.frame(LUAD_Macrophage_id$Macrophage/LUAD_all_cells$All_cells)
names(LUAD_Macrophage_ratio) <- "Macrophage_ratio"
LUAD_Macrophage_ratio$type <- "LUAD_Macrophage"
SCLC_Macrophage_ratio <- data.frame(SCLC_Macrophage_id$Macrophage/SCLC_all_cells$All_cells)
names(SCLC_Macrophage_ratio) <- "Macrophage_ratio"
SCLC_Macrophage_ratio$type <- "SCLC_Macrophage"
Macrophage_ratio <- rbind(LUAD_Macrophage_ratio,SCLC_Macrophage_ratio)
Macrophage_ratio$type <- factor(Macrophage_ratio$type,levels=c("LUAD_Macrophage","SCLC_Macrophage"))
library(ggplot2)
library(ggpubr)
my_comparisons <- list(c("SCLC_Macrophage","LUAD_Macrophage"))
ff3 <- ggboxplot(Macrophage_ratio, x = "type", y = "Macrophage_ratio",add = "jitter",
               color = "type",palette = c("#000000","#a50026")) +
stat_compare_means(comparisons = my_comparisons, label = "p.format",method="t.test")

ff <- ff1+ff2+ff3
ggsave(ff,file="FigS1B.png",width=12,height=5)


#FigS1C
Endo_pure <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Rudin_CancerCell/Endothelium/Pure_Endo_CancerCell.rds",mc.cores=20)
Endo_pure$histo <- factor(Endo_pure$histo,levels=c("normal","LUAD","SCLC"))
pdf("FigS1C.pdf",width=15,height=10)
VlnPlot(Endo_pure, features = c("MMP2","NID1","LAMC1","COL4A1","LAMA4","LAMB1"), ncol = 3,pt.size=0,group.by="histo")
dev.off()
