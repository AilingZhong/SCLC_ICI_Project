#Fig6_S6
library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)
library(future)
library(future.apply)
library(scales)
library(BuenColors)
library(ggplot2)
library(data.table)
library(trqwe)
library(harmony)
library(reticulate)
library(ReductionWrappers)
library(s2a)


Med_Luyou <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/scRNA/SCLC_LY_NFF/output_v3/ALL_MERGE_DATA_filter.rds", mc.cores = 20)
Idents(Med_Luyou) <- Med_Luyou$group
Ctrl_ICI <- subset(Med_Luyou,idents=c("Ctrl","ICI"),invert=FALSE)
Ctrl_ICI@meta.data <- Ctrl_ICI@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mito","group","sample")]
Ctrl_ICI@meta.data$percent.mt <- Ctrl_ICI@meta.data$percent.mito*100
Ctrl_ICI@meta.data <- Ctrl_ICI@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","group","sample")]
Ctrl_ICI$location <- "subcutaneous"
CCLYlab <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Mice_scRNAseq/OSI_ICI_Processed_Data/Old_without_Luyou/All_MERGE_DATA_pca.rds", mc.cores = 20)
CCLYlab$group <- CCLYlab$sample
CCLYlab@meta.data <- CCLYlab@meta.data[,c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","group","sample")]
CCLYlab$location <- "subcutaneous"
merge.data <- merge(x = Ctrl_ICI, y=CCLYlab)
All_obj <- CreateSeuratObject(counts =GetAssayData(merge.data, 
	slot = "counts",assay="RNA")[,rownames(merge.data@meta.data)], meta.data = merge.data@meta.data)
All_MERGE_DATA <- All_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
All_MERGE_DATA.pca <- RunPCA(All_MERGE_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(All_MERGE_DATA.pca,file="All_MERGE_DATA_pca.rds",mc.cores=20)

ALL_MERGE_DATA_harmony  <- RunHarmony(All_MERGE_DATA.pca,"orig.ident")
FindNeighbors <- FindNeighbors(ALL_MERGE_DATA_harmony,reduction = "harmony", dims = 1:20)
FindClusters <- FindClusters(FindNeighbors,reduction = "harmony",resolution = c(0.1,0.2,0.5,1))
mcsaveRDS(FindClusters,file="ALL_MERGE_DATA_harmony_FindClusters.rds", mc.cores = 20)

All_UMAP20 <- RunUMAP(FindClusters,reduction = "harmony", dims = 1:20)
mcsaveRDS(All_UMAP20,file="Harmony_AllCells_UMAP20.rds", mc.cores = 20)

ff <- DimPlot(object = All_UMAP20, reduction = "umap",pt.size=0.4,label=TRUE, 
    label.size=7,split.by="sample")

Idents(All_UMAP20) <- All_UMAP20$RNA_snn_res.1
all.markers <- FindAllMarkers(All_UMAP20, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(all.markers,"All_UMAP20_FindAllmarkers_res1.csv")



Idents(All_UMAP20) <- All_UMAP20$RNA_snn_res.1
ff <- DimPlot(object = All_UMAP20, reduction = "umap",pt.size=0.4,label=TRUE, label.size=7)+
    cols=c("#cb181d","#2166ac","#67a9cf","#238b45","#f1b6da","#f4a582"),label=TRUE, label.size=7)
ggsave(ff,file="Renames_All_Harmony_UMAP30.png",width=9.3,height=8.5,dpi=1080)


All_UMAP20$RNA_snn_res.1 <- factor(All_UMAP20$RNA_snn_res.1,levels=c(0:28))
Idents(All_UMAP20) <- All_UMAP20$RNA_snn_res.1
new.cluster.ids <- c("Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","T.cell","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Tumor","Macro","Tumor","NK","T.cell","Neutro","Fibro","Plasma","Tumor","DC","Tumor","Endo","Tumor")
names(new.cluster.ids) <- levels(All_UMAP20)
All_UMAP20 <- RenameIdents(All_UMAP20, new.cluster.ids)
All_UMAP20$CellType <- Idents(All_UMAP20)
Idents(All_UMAP20) <- All_UMAP20$CellType
All_UMAP20$CellType <- factor(All_UMAP20$CellType,levels=c("Tumor","T.cell","NK","Macro","Plasma","DC","Neutro","Fibro","Endo"))
mcsaveRDS(All_UMAP20,file="Add_Med_All_UMAP20_Renames.rds",mc.cores=20)





All_UMAP20 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Add_Med_All_UMAP20_Renames.rds",mc.cores=20)
All_UMAP20$CellType <- factor(All_UMAP20$CellType,levels=c("Tumor","T.cell","Plasma","NK","Neutro","Macro","Fibro","Endo","B.cell"))
Idents(All_UMAP20) <- All_UMAP20$CellType
markers.to.plot <- as.character(c("LUCI2LTR","CAS9","V2TC",
    "Cd3e","Cd3d","Cd3g","Cd8a","Cd8b","Cd4",
    "Jchain","Mzb1",
    "Nkg7","Klre1","Gzma",
    "S100a9","S100a8","G0s2",
    "Csf1r","Cd68",  "Col1a1","Col1a2","Dcn","Pecam1","Plvap","Kdr","Cd79a","Cd79b"))
ff <- DotPlot(All_UMAP20,features = markers.to.plot,cols = c("#ffffff","#cb181d"),
    assay="RNA",col.min =0,col.max =3,dot.min=0.05,scale=TRUE) + RotatedAxis()
ggsave(ff,file="Marker_Dotplot.pdf",width=9.5,height=4.5)


meta.data <- All_UMAP20@meta.data
Naive <- subset(meta.data,sample=="Ctrl1" | sample=="Ctrl2")
ICI <- subset(meta.data,sample=="antiPD1" | sample=="ICI1" | sample=="ICI2")
OSI <- subset(meta.data,sample=="OSI")
OSI_ICI <- subset(meta.data,sample=="OSI_antiPD1")

Naive <- data.frame(table(Naive$CellType))
ICI <- data.frame(table(ICI$CellType))
OSI <- data.frame(table(OSI$CellType))
OSI_ICI <- data.frame(table(OSI_ICI$CellType))

Naive$group <- "Naive"
ICI$group <- "ICI"
OSI$group <- "OSI"
OSI_ICI$group <- "OSI_ICI"

Naive$Ratio <- Naive$Freq/sum(Naive$Freq)
ICI$Ratio <- ICI$Freq/sum(ICI$Freq)
OSI$Ratio <- OSI$Freq/sum(OSI$Freq)
OSI_ICI$Ratio <- OSI_ICI$Freq/sum(OSI_ICI$Freq)

tmp1_data <- rbind(Naive,ICI,OSI,OSI_ICI)
tmp1_data$group <- factor(tmp1_data$group,levels=c("Naive","ICI","OSI","OSI_ICI"))

tmp1_data$Var1 <- factor(tmp1_data$Var1,levels=c("Tumor","Fibro","Endo","T.cell","NK","Macro","Plasma","DC","Neutro"))
library(ggpubr)
library(ggalluvial)
ff <- ggplot(tmp1_data, aes( x = group,y=Ratio,fill = Var1, stratum = Var1, alluvium = Var1))+
  geom_col(position = 'stack', width = 0.6)+  scale_y_continuous(breaks=seq(0,1,0.01)) +
  theme_bw() + scale_fill_manual(values=c("#d9d9d9","#6a51a3","#f768a1","#cb181d","#2166ac","#67a9cf","#238b45","#f1b6da","#f4a582")) + geom_flow(alpha = 0.5) 
ggsave(ff,file="Fig6I_CellType_Ratio.pdf",width=6,height=9)





merge.data <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Add_Med_All_UMAP20_Renames.rds",mc.cores=20)
Idents(merge.data) <- merge.data$CellType
Tcell <- subset(merge.data,idents="T.cell",invert=FALSE)
Tcell_object <- CreateSeuratObject(counts =GetAssayData(Tcell, slot = "counts",assay="RNA")[,rownames(Tcell@meta.data)], meta.data = Tcell@meta.data)
Tcell_object[["percent.mt"]] <- PercentageFeatureSet(Tcell_object, pattern = "^mt-")
All_obj <- CreateSeuratObject(counts =GetAssayData(Tcell_object, 
	slot = "counts",assay="RNA")[,rownames(Tcell_object@meta.data)], meta.data = Tcell_object@meta.data)
All_Tcell_DATA <- All_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
All_Tcell_DATA.pca <- RunPCA(All_Tcell_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(All_Tcell_DATA.pca,file="All_Tcell_DATA_pca.rds",mc.cores=20)
FindNeighbors <- FindNeighbors(All_Tcell_DATA.pca, dims = 1:20)
FindClusters <- FindClusters(FindNeighbors,resolution = c(1))
All_UMAP20 <- RunUMAP(FindClusters, dims = 1:20)
mcsaveRDS(All_UMAP20,file="All_UMAP20_Tcell.rds",mc.cores=20)

ff <- DimPlot(object = All_UMAP20, reduction = "umap",pt.size=0.4,label=TRUE, label.size=7)
Idents(All_UMAP20) <- All_UMAP20$RNA_snn_res.1
all.markers <- FindAllMarkers(All_UMAP20, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.markers,"All_UMAP20_FindAllmarkers_res1.csv")

Th1 <- intersect(c("TBX21","CD40LG","IL7R","KLRB1","CD4","RORA","TNFRSF25","ANXA1","LTB","GPR183","CDKN1A","LMNA","S100A11","S100A4","PLIN2"),rownames(only_NKT_T.fil_seurat))
Th2 <- intersect(c("IL1RL1", "GATA3", "IL13", "IL5", "IL17RB", "LTB4R1", "CCR8", "CD200R1", "IL6", "PLAC8", "IGFBP7","CTLA2A", "AHNAK", "RBPJ", "BHLHE40"),rownames(only_NKT_T.fil_seurat))
Tfh <- intersect(c("TOX2","MAF","TBC1D4","TNFRSF4","CD4","UCP2","CORO1B","TIGIT","BATF","NR3C1","ITM2A","TNFRSF18","LIMS1",
    "ICA1","CTSB","SH2D1A","CTLA4","MAGEH1","SESN3","ICOS","NMB","CXCL13","CD200","BTLA","PTPN13"),rownames(only_NKT_T.fil_seurat))
CD8Tm <- intersect(c("EOMES","GZMK","CXCR3"),rownames(only_NKT_T.fil_seurat))
CD8Tem <- intersect(c("GZMA","GZMB","CCR5","HOPX","CXCR6","NKG7","CTSW","CCL5","CST7","PRF1","FASL","ITM2C","CHST12","AOAH","SLAMF7","F2R"),rownames(only_NKT_T.fil_seurat))
Texp <- intersect(c("PDCD1","LAG3","HAVCR2","TIGIT"),rownames(only_NKT_T.fil_seurat))
Tex <- intersect(c("GZMB","ENTPD1","ITGAE"),rownames(only_NKT_T.fil_seurat))
All_genes <- list(Th1,Th2,Tfh,CD8Tm,CD8Tem)
names(All_genes) <- list("Th1","Th2","Tfh","CD8Tm","CD8Tem")

Idents(All_UMAP20) <- All_UMAP20$RNA_snn_res.1
Pure_Tcell<- subset(All_UMAP20,idents=c("5","6","10","12","13"),invert=TRUE)
new.cluster.ids <- c("Naive.like","Early.Active","Cd8.Eff.Mem.","Cd8.Proli","Cd8.Eff.Mem.","Tfh","Treg","Cd8.Eff.Mem.","Th17","Th1")
names(new.cluster.ids) <- levels(Pure_Tcell)
Pure_Tcell <- RenameIdents(Pure_Tcell, new.cluster.ids)
Pure_Tcell$sub_CellType <- Idents(Pure_Tcell)
Idents(Pure_Tcell) <- Pure_Tcell$sub_CellType
Pure_Tcell$sub_CellType <- factor(Pure_Tcell$sub_CellType,levels=c("Naive.like","Early.Active","Cd8.Eff.Mem.","Cd8.Proli","Tfh","Treg","Th17","Th1"))
mcsaveRDS(Pure_Tcell,file="Pure_Tcell_Renames.rds",mc.cores=20)






#FigS6C-D
#9022
antiPD1_sgScr <- Read10X(data.dir="/mnt/data/user_data/abao/1_project/Single_cell_project/7_WangYiyun/Endo_SCLC_ICI/Cellranger/antiPD1_sgScr/outs/filtered_feature_bc_matrix/")
antiPD1_sgIGFBP5 <- Read10X(data.dir="/mnt/data/user_data/abao/1_project/Single_cell_project/7_WangYiyun/Endo_SCLC_ICI/Cellranger/antiPD1_sgIGFBP5/outs/filtered_feature_bc_matrix/")

antiPD1_sgScr_object <- CreateSeuratObject(counts = antiPD1_sgScr, project = "antiPD1_sgScr", min.cells = 3, min.features = 200)
antiPD1_sgScr_object$sample <- Idents(object = antiPD1_sgScr_object)

antiPD1_sgIGFBP5_object <- CreateSeuratObject(counts = antiPD1_sgIGFBP5, project = "antiPD1_sgIGFBP5", min.cells = 3, min.features = 200)
antiPD1_sgIGFBP5_object$sample <- Idents(object = antiPD1_sgIGFBP5_object)

merge.data <- merge(x = antiPD1_sgScr_object,y = c(antiPD1_sgIGFBP5_object))

merge_object <- CreateSeuratObject(counts =GetAssayData(merge.data, slot = "counts",assay="RNA")[,rownames(merge.data@meta.data)], meta.data = merge.data@meta.data)
merge_object[["percent.mt"]] <- PercentageFeatureSet(merge_object, pattern = "^mt-")
pdf("percent_MT_qc_renew.pdf",width=15,height=5)
VlnPlot(merge_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0,group.by="sample")
dev.off()
merge_object <- subset(merge_object, subset = percent.mt < 20)
merge_object <- subset(merge_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7500)
merge_object <- subset(merge_object, subset = nCount_RNA > 500)
All_obj <- CreateSeuratObject(counts =GetAssayData(merge_object, 
	slot = "counts",assay="RNA")[,rownames(merge_object@meta.data)], meta.data = merge_object@meta.data)
All_MERGE_DATA <- All_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
All_MERGE_DATA.pca <- RunPCA(All_MERGE_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(All_MERGE_DATA.pca,file="All_MERGE_DATA_pca.rds",mc.cores=20)
All_FindNeighbors <- FindNeighbors(All_MERGE_DATA.pca,reduction = "pca", dims = 1:50)
All_FindClusters <- FindClusters(All_FindNeighbors,resolution = c(0.1,0.2,0.5,1))
All_UMAP20 <- RunUMAP(All_FindClusters, dims = 1:20)
mcsaveRDS(All_UMAP20,file="AllCells_UMAP20.rds", mc.cores = 20)


All_TSNE10 <- RunTSNE(All_FindClusters, dims = 1:10)

library(ggplot2)
Idents(All_TSNE10) <- All_TSNE10$RNA_snn_res.0.2
new.cluster.ids <- c("Tumor","Tumor","Tumor","Mono&Macro","Tumor","T cell","Tumor","Fibro")
names(new.cluster.ids) <- levels(All_TSNE10)
All_TSNE10 <- RenameIdents(All_TSNE10, new.cluster.ids)
All_TSNE10$CellType <- Idents(All_TSNE10)
Idents(All_TSNE10) <- All_TSNE10$CellType
All_TSNE10$CellType <- factor(All_TSNE10$CellType,levels=c("Tumor","Fibro","Mono&Macro","T cell"))
mcsaveRDS(All_TSNE10,file="AllCells_TSNE10_Renames.rds", mc.cores = 20)


ff <- DimPlot(object = All_TSNE10, reduction = "tsne",pt.size=0.8,
    cols=c("#d6604d","#1f78b4","#b2182b","#33a02c"),label=FALSE, label.size=7)
ggsave(ff,file="Renames_All_TSNE10.png",width=7.3,height=6.7,dpi=1080)

Idents(All_TSNE10) <- All_TSNE10$RNA_snn_res.0.2
markers.to.plot <- as.character(c("LUCI2LTR","V2TC","Chga","Chgb","Ascl1","Epcam","Lgals3","Csf1r","Cd68","Cd3g","Cd3e","Col1a1","Fn1","Jchain",
"Dcn","Pecam1","Vwf","Cd79a","Cd79b","Lyz2","G0s2"))
ff <- DotPlot(All_TSNE10,features = markers.to.plot,cols = c("#ffffff","#cb181d"),
    assay="RNA",col.min =0,col.max =3,dot.min=0.05,scale=TRUE) + RotatedAxis()

Idents(All_TSNE10) <- All_TSNE10$RNA_snn_res.0.2
all.markers <- FindAllMarkers(All_TSNE10, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all.markers,"All_TSNE10_FindAllmarkers_res02.csv")

Idents(All_TSNE10) <- All_TSNE10$CellType
All_TSNE10$CellType <- factor(All_TSNE10$CellType,levels=c("Tumor","Fibro","Mono&Macro","T cell"))
markers.to.plot <- as.character(c("LUCI2LTR","V2TC","Ascl1","Chga","Chgb","Col1a1","Fn1","Csf1r","Cd68","Lyz2","Cd3g","Cd3e"))
ff <- DotPlot(All_TSNE10,features = markers.to.plot,cols = c("#ffffff","#cb181d"),
    assay="RNA",col.min =0,col.max =3,dot.min=0.05,scale=TRUE) + RotatedAxis()
ggsave(ff,file="Marker_scale_renames.pdf",width=7,height=3.5)


meta.data <- All_TSNE10@meta.data
antiPD1_sgScr <- subset(meta.data,sample=="antiPD1_sgScr")
antiPD1_sgIGFBP5 <- subset(meta.data,sample=="antiPD1_sgIGFBP5")
antiPD1_sgScr <- data.frame(table(antiPD1_sgScr$CellType))
antiPD1_sgIGFBP5 <- data.frame(table(antiPD1_sgIGFBP5$CellType))
antiPD1_sgScr$group <- "antiPD1_sgScr"
antiPD1_sgIGFBP5$group <- "antiPD1_sgIGFBP5"
antiPD1_sgScr$Ratio <- antiPD1_sgScr$Freq/sum(antiPD1_sgScr$Freq)
antiPD1_sgIGFBP5$Ratio <- antiPD1_sgIGFBP5$Freq/sum(antiPD1_sgIGFBP5$Freq)
tmp1_data <- rbind(antiPD1_sgScr,antiPD1_sgIGFBP5)
tmp1_data$group <- factor(tmp1_data$group,levels=c("antiPD1_sgScr","antiPD1_sgIGFBP5"))
tmp1_data$Var1 <- factor(tmp1_data$Var1,levels=c("Tumor","Fibro","Mono&Macro","T cell"))
library(ggpubr)
library(ggalluvial)
ff <- ggplot(tmp1_data, aes( x = group,y=Ratio,fill = Var1, stratum = Var1, alluvium = Var1))+
  geom_col(position = 'stack', width = 0.6)+  scale_y_continuous(breaks=seq(0,1,0.01)) +
  theme_bw() + scale_fill_manual(values=c("#878787","#1f78b4","#b2182b","#33a02c")) + geom_flow(alpha = 0.5) #绘制同类别之间的连接线
ggsave(ff,file="FigS6D.pdf",width=4,height=9)



%%%%%%%%%%%%%%%%%%%%%%T_cell_function_marker%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%T_cell_function_marker%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%T_cell_function_marker%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%T_cell_function_marker%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd8_core_genes <- as.character(c(
  "Jak2", "Stat5b", "E2f1", "Ccnd3",
  "Itgb7", "Sell", "Cxcr4", "Cxcr3",
  "Msi2", "Satb1", "Bach2", "Lef1",
  "Fyn", "Fasl", "Nkg7", "Gzmk", "Gzma", "Gzmb"
))

Tcell <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Pure_Tcell_final_update.rds",mc.cores=20)
CD8_Tcell <- subset(Tcell,idents=c('CD8.T_Cycling','CD8.T_TerminalExhausted','CD8.T_IntermediateExhausted','CD8.T_ProgenitorPrecursorExhausted'))
Idents(CD8_Tcell) <- CD8_Tcell$treat
mtx <- GetAssayData(CD8_Tcell, slot = "data",assay="RNA")
mtx <- data.frame(mtx)
#mtx2 <- mtx[rowMeans(mtx > 0) > 0.01, ]
mtx_z <- t(scale(t(mtx)))
mtx_z[is.na(mtx_z)] <- 0
mtx_z <- data.frame(mtx_z)

antiPD1 <- subset(CD8_Tcell,idents=c("antiPD1"))
OSI_antiPD1 <- subset(CD8_Tcell,idents=c("OSI_antiPD1"))
antiPD1_mtx <- GetAssayData(antiPD1, slot = "data",assay="RNA")
ICI_genes <- mtx_z[cd8_core_genes,colnames(data.frame(antiPD1_mtx))]
OSI_antiPD1_mtx <- GetAssayData(OSI_antiPD1, slot = "data",assay="RNA")
OSI_antiPD1_genes <- mtx_z[cd8_core_genes,colnames(data.frame(OSI_antiPD1_mtx))]


ICI_mean <- data.frame(apply(ICI_genes,1,mean)) 
OSI_mean <- data.frame(apply(OSI_antiPD1_genes,1,mean)) 
names(ICI_mean) <- "ICI"
names(OSI_mean) <- "ICI_OSI"
tmp <- cbind(ICI_mean,OSI_mean)
tmp <- na.omit(tmp)
tmp[tmp < -1.2] <- -1.2
tmp[tmp > 1.2] <- 1.2

library(pheatmap)
pdf("FigS6H_1.pdf",width=3,height=7)
pheatmap(tmp,cluster_rows = FALSE,cluster_cols = FALSE,scale = "none",
	color = c("#2166ac", "#67a9cf", "#d1e5f0","#fddbc7", "#ef8a62", "#b2182b"), main = "CD8 Signature")
dev.off()



$$$$$$$$$$$$$$$Tcell_in_sgIGFBP5_combine_aPD1_treatment$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$Tcell_in_sgIGFBP5_combine_aPD1_treatment$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
$$$$$$$$$$$$$$$Tcell_in_sgIGFBP5_combine_aPD1_treatment$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

Tcell_MERGE_UMAP15 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/sgIgfbp5_combine_aPD1/Tcell_MERGE_UMAP15.rds",mc.cores=20)
Idents(Tcell_MERGE_UMAP15) <- Tcell_MERGE_UMAP15$RNA_snn_res.1
Pure_Tcell<- subset(Tcell_MERGE_UMAP15,idents=c("3","4"),invert=TRUE)
new.cluster.ids <- c("CD8T.exh","CD8T.effmem","CD4.Tcell")
names(new.cluster.ids) <- levels(Pure_Tcell)
Pure_Tcell <- RenameIdents(Pure_Tcell, new.cluster.ids)
Pure_Tcell$sub_CellType <- Idents(Pure_Tcell)
Idents(Pure_Tcell) <- Pure_Tcell$sub_CellType
Pure_Tcell$sub_CellType <- factor(Pure_Tcell$sub_CellType,levels=c("CD8T.exh","CD8T.effmem","CD4.Tcell"))

CD8_Tcell <- subset(Pure_Tcell,idents=c("CD8T.exh","CD8T.effmem"))
Idents(CD8_Tcell) <- CD8_Tcell$sample

mtx <- GetAssayData(CD8_Tcell, slot = "data",assay="RNA")
mtx <- data.frame(mtx)
#mtx2 <- mtx[rowMeans(mtx > 0) > 0.01, ]
mtx_z <- t(scale(t(mtx)))
mtx_z[is.na(mtx_z)] <- 0
mtx_z <- data.frame(mtx_z)

antiPD1_sgScr <- subset(CD8_Tcell,idents=c("antiPD1_sgScr"))
antiPD1_sgIGFBP5 <- subset(CD8_Tcell,idents=c("antiPD1_sgIGFBP5"))

antiPD1_sgScr_mtx <- GetAssayData(antiPD1_sgScr, slot = "data",assay="RNA")
sgScr_genes <- mtx_z[cd8_core_genes,colnames(data.frame(antiPD1_sgScr_mtx))]

antiPD1_sgIGFBP5_mtx <- GetAssayData(antiPD1_sgIGFBP5, slot = "data",assay="RNA")
sgIGFBP5_genes <- mtx_z[cd8_core_genes,colnames(data.frame(antiPD1_sgIGFBP5_mtx))]

sgScr_mean <- data.frame(apply(sgScr_genes,1,mean)) 
sgIGFBP5_mean <- data.frame(apply(sgIGFBP5_genes,1,mean)) 
names(sgScr_mean) <- "sgScr"
names(sgIGFBP5_mean) <- "sgIGFBP5"
tmp <- cbind(sgScr_mean,sgIGFBP5_mean)

tmp <- na.omit(tmp)

tmp[tmp < -0.1] <- -0.1
tmp[tmp > 0.1] <- 0.1

library(pheatmap)
pdf("FigS6H_2.pdf",width=4,height=12)
pheatmap(
  tmp,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",   
  color = c("#2166ac", "#67a9cf", "#d1e5f0", "#fddbc7", "#ef8a62", "#b2182b"),
  main = "CD8 Signature (Z-score)"
)
dev.off()



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

All_cells <- mcreadRDS("/mnt/data/user_data/xiangyu/workshop/Zhong_ABao/0_Endo_SCLC/Mice_scRNAseq/OSI_PD1_Processed_Data/Add_Luyou_CRM/Add_CRM_All_UMAP20_Renames.rds",mc.cores=20)
table(All_cells$group)
All_cells$group <- gsub("ICI","antiPD1",All_cells$group)

Idents(All_cells) <- All_cells$CellType

Tumor <- subset(All_cells,idents=c("Tumor"),invert=FALSE)
T.cell <- subset(All_cells,idents=c("T.cell"),invert=FALSE)
Fibro <- subset(All_cells,idents=c("Fibro"),invert=FALSE)
Endo <- subset(All_cells,idents=c("Endo"),invert=FALSE)
Macro <- subset(All_cells,idents=c("Macro"),invert=FALSE)
Neutro <- subset(All_cells,idents=c("Neutro"),invert=FALSE)
DC <- subset(All_cells,idents=c("DC"),invert=FALSE)
Plasma <- subset(All_cells,idents=c("Plasma"),invert=FALSE)
NK <- subset(All_cells,idents=c("NK"),invert=FALSE)

Idents(Tumor) <- Tumor$group
Idents(T.cell) <- T.cell$group
Idents(Endo) <- Endo$group
Idents(Fibro) <- Fibro$group
Idents(Macro) <- Macro$group
Idents(Neutro) <- Neutro$group
Idents(DC) <- DC$group
Idents(Plasma) <- Plasma$group
Idents(NK) <- NK$group


Tumor_markers <- FindAllMarkers(Tumor, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(Tumor_markers,"DEGs_Tumor_markers.csv")

Endo_markers <- FindAllMarkers(Endo, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(Endo_markers,"DEGs_Endo_markers.csv")
Fibro_markers <- FindAllMarkers(Fibro, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(Fibro_markers,"DEGs_Fibro_markers.csv")

DC_markers <- FindAllMarkers(DC, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(DC_markers,"DEGs_DC_markers.csv")
T.cell_markers <- FindAllMarkers(T.cell, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(T.cell_markers,"DEGs_T.cell_markers.csv")
Plasma_markers <- FindAllMarkers(Plasma, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(Plasma_markers,"DEGs_Plasma_markers.csv")
Neutro_markers <- FindAllMarkers(Neutro, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(Neutro_markers,"DEGs_Neutro_markers.csv")
Macro_markers <- FindAllMarkers(Macro, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(Macro_markers,"DEGs_Macro_markers.csv")

NK_markers <- FindAllMarkers(NK, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(NK_markers,"DEGs_NK_markers.csv")



#Fig S6
#cell cell interaction
library(dplyr)
library(Seurat)
library(cowplot)
library(patchwork)
library(iTALK)

All_cells <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Add_Med_All_UMAP20_Renames.rds",mc.cores=20)
table(All_cells$group)
All_cells$group <- gsub("ICI","antiPD1",All_cells$group)

Idents(All_cells) <- All_cells$group
Ctrl <- subset(All_cells, idents=c("Ctrl"),invert=FALSE)
ICI <- subset(All_cells, idents=c("antiPD1"),invert=FALSE)
OSI <- subset(All_cells, idents=c("OSI"),invert=FALSE)
OSI_antiPD1 <- subset(All_cells, idents=c("OSI_antiPD1"),invert=FALSE)

Ctrl_data <- data.frame(GetAssayData(object = Ctrl, slot = 'data'))
Ctrl_data_1 <- Ctrl_data %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(rownames(Ctrl_data)))) %>% drop_na()
Ctrl_data_1 <- Ctrl_data_1[!duplicated(Ctrl_data_1$human_gene),]
rownames(Ctrl_data_1) <- Ctrl_data_1$human_gene
matrix <- Ctrl_data_1[,-21512]
matrix <- data.frame(t(matrix))
meta.data <- data.frame(Ctrl@meta.data)
matrix$cell_type <- meta.data$CellType
matrix$batch <- meta.data$group
highly_exprs_genes <- rawParse(matrix,top_genes=20000,stats='mean')
ICI_Tcell <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='growth factor', database = NULL)
other <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='other', database = NULL)
cytokine <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='cytokine', database = NULL)
checkpoint <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='checkpoint', database = NULL)
aa <- rbind(growth_factor,cytokine,checkpoint,other)
write.csv(aa,"Ctrl_ligand_receptpr.csv")


ICI_data <- data.frame(GetAssayData(object = ICI, slot = 'data'))
ICI_data_1 <- ICI_data %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(rownames(ICI_data)))) %>% drop_na()
ICI_data_1 <- ICI_data_1[!duplicated(ICI_data_1$human_gene),]
rownames(ICI_data_1) <- ICI_data_1$human_gene
matrix <- ICI_data_1[,-32275]
matrix <- data.frame(t(matrix))
meta.data <- data.frame(ICI@meta.data)
matrix$cell_type <- meta.data$CellType
matrix$batch <- meta.data$group
highly_exprs_genes <- rawParse(matrix,top_genes=20000,stats='mean')
growth_factor <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='growth factor', database = NULL)
other <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='other', database = NULL)
cytokine <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='cytokine', database = NULL)
checkpoint <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='checkpoint', database = NULL)
aa <- rbind(growth_factor,cytokine,checkpoint,other)
write.csv(aa,"ICI_ligand_receptpr.csv")


OSI_data <- data.frame(GetAssayData(object = OSI, slot = 'data'))
OSI_data_1 <- OSI_data %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(rownames(OSI_data)))) %>% drop_na()
OSI_data_1 <- OSI_data_1[!duplicated(OSI_data_1$human_gene),]
rownames(OSI_data_1) <- OSI_data_1$human_gene
matrix <- OSI_data_1[,-12988]
matrix <- data.frame(t(matrix))
meta.data <- data.frame(OSI@meta.data)
matrix$cell_type <- meta.data$CellType
matrix$batch <- meta.data$group
highly_exprs_genes <- rawParse(matrix,top_genes=20000,stats='mean')
growth_factor <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='growth factor', database = NULL)
other <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='other', database = NULL)
cytokine <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='cytokine', database = NULL)
checkpoint <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='checkpoint', database = NULL)
aa <- rbind(growth_factor,cytokine,checkpoint,other)
write.csv(aa,"OSI_ligand_receptpr.csv")


OSI_antiPD1_data <- data.frame(GetAssayData(object = OSI_antiPD1, slot = 'data'))
OSI_antiPD1_data_1 <- OSI_antiPD1_data %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(rownames(OSI_antiPD1_data)))) %>% drop_na()
OSI_antiPD1_data_1 <- OSI_antiPD1_data_1[!duplicated(OSI_antiPD1_data_1$human_gene),]
rownames(OSI_antiPD1_data_1) <- OSI_antiPD1_data_1$human_gene
matrix <- OSI_antiPD1_data_1[,-14071]
matrix <- data.frame(t(matrix))
meta.data <- data.frame(OSI_antiPD1@meta.data)
matrix$cell_type <- meta.data$CellType
matrix$batch <- meta.data$group
highly_exprs_genes <- rawParse(matrix,top_genes=20000,stats='mean')
growth_factor <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='growth factor', database = NULL)
other <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='other', database = NULL)
cytokine <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='cytokine', database = NULL)
checkpoint <- FindLR(highly_exprs_genes, data_2 = NULL, datatype="mean count",comm_type='checkpoint', database = NULL)
aa <- rbind(growth_factor,cytokine,checkpoint,other)
write.csv(aa,"OSI_antiPD1_ligand_receptpr.csv")


#Screening based on differentially expressed genes

Ctrl_LR <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/iTalk/Ctrl_ligand_receptpr.csv")
ICI_LR <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/iTalk/ICI_ligand_receptpr.csv")
OSI_LR <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/iTalk/OSI_ligand_receptpr.csv")
ICI_OSI_LR <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/iTalk/OSI_antiPD1_ligand_receptpr.csv")

Ctrl_LR$CellType_Ligand <- paste(Ctrl_LR$cell_from,"_",Ctrl_LR$ligand,sep="")
Ctrl_LR$CellType_Receptor <- paste(Ctrl_LR$cell_to,"_",Ctrl_LR$receptor,sep="")

ICI_LR$CellType_Ligand <- paste(ICI_LR$cell_from,"_",ICI_LR$ligand,sep="")
ICI_LR$CellType_Receptor <- paste(ICI_LR$cell_to,"_",ICI_LR$receptor,sep="")

OSI_LR$CellType_Ligand <- paste(OSI_LR$cell_from,"_",OSI_LR$ligand,sep="")
OSI_LR$CellType_Receptor <- paste(OSI_LR$cell_to,"_",OSI_LR$receptor,sep="")

ICI_OSI_LR$CellType_Ligand <- paste(ICI_OSI_LR$cell_from,"_",ICI_OSI_LR$ligand,sep="")
ICI_OSI_LR$CellType_Receptor <- paste(ICI_OSI_LR$cell_to,"_",ICI_OSI_LR$receptor,sep="")


DEGs_Tumor <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_Tumor_markers.csv")
DEGs_Neutro <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_Neutro_markers.csv")
DEGs_NK <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_NK_markers.csv")
DEGs_Macro <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_Macro_markers.csv")
DEGs_DC <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_DC_markers.csv")
DEGs_Fibro <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_Fibro_markers.csv")
DEGs_Endo <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_Endo_markers.csv")
DEGs_Plasma <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_Plasma_markers.csv")
DEGs_T.cell <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/iTALK_homa_made_method/DEGs_T.cell_markers.csv")

DEGs_Tumor$CellType <- "Tumor"
DEGs_Neutro$CellType <- "Neutro"
DEGs_NK$CellType <- "NK"
DEGs_Macro$CellType <- "Macro"
DEGs_DC$CellType <- "DC"
DEGs_Fibro$CellType <- "Fibro"
DEGs_Endo$CellType <- "Endo"
DEGs_Plasma$CellType <- "Plasma"
DEGs_T.cell$CellType <- "T.cell"

All_DEGs <- rbind(DEGs_Tumor,DEGs_Neutro,DEGs_NK,DEGs_Macro,DEGs_DC,DEGs_Fibro,DEGs_Endo,DEGs_Plasma,DEGs_T.cell)
All_DEGs <- All_DEGs %>% mutate(human_gene = convert_mouse_to_human_symbols(as.character(All_DEGs$gene))) %>% drop_na()
All_DEGs$CellType_Genes <- paste(All_DEGs$CellType,"_",All_DEGs$human_gene,sep="")

Ctrl_hi_genes <- subset(All_DEGs, avg_log2FC > 0 & p_val_adj < 0.01 & cluster=="Ctrl")
ICI_hi_genes <- subset(All_DEGs, avg_log2FC > 0 & p_val_adj < 0.01 & cluster=="antiPD1")
OSI_hi_genes <- subset(All_DEGs, avg_log2FC > 0 & p_val_adj < 0.01 & cluster=="OSI")
ICI_OSI_hi_genes <- subset(All_DEGs, avg_log2FC > 0 & p_val_adj < 0.01 & cluster=="OSI_antiPD1")


common_Ligand <- intersect(Ctrl_LR$CellType_Ligand,Ctrl_hi_genes$CellType_Genes)
Ctrl_LR <- subset(Ctrl_LR, cell_from_mean_exprs > 0 & cell_to_mean_exprs > 0)
Ctrl_LR_filtered <- Ctrl_LR %>% filter(CellType_Ligand %in% common_Ligand)
table(Ctrl_LR_filtered$cell_from)


common_Ligand <- intersect(ICI_LR$CellType_Ligand,ICI_hi_genes$CellType_Genes)
ICI_LR <- subset(ICI_LR, cell_from_mean_exprs > 0 & cell_to_mean_exprs > 0)
ICI_LR_filtered <- ICI_LR %>% filter(CellType_Ligand %in% common_Ligand)
table(ICI_LR_filtered$cell_from)


common_Ligand <- intersect(OSI_LR$CellType_Ligand,OSI_hi_genes$CellType_Genes)
OSI_LR <- subset(OSI_LR, cell_from_mean_exprs > 0 & cell_to_mean_exprs > 0)
OSI_LR_filtered <- OSI_LR %>% filter(CellType_Ligand %in% common_Ligand)
table(OSI_LR_filtered$cell_from)

common_Ligand <- intersect(ICI_OSI_LR$CellType_Ligand,ICI_OSI_hi_genes$CellType_Genes)
ICI_OSI_LR <- subset(ICI_OSI_LR, cell_from_mean_exprs > 0 & cell_to_mean_exprs > 0)
ICI_OSI_LR_filtered <- ICI_OSI_LR %>% filter(CellType_Ligand %in% common_Ligand)
table(ICI_OSI_LR_filtered$cell_from)


Ctrl_Tcell <- subset(Ctrl_LR_filtered,cell_from=="T.cell" & cell_to=="Tumor")
Ctrl_Tcell$ligand <- as.character(Ctrl_Tcell$ligand)
Ctrl_Tcell$receptor <- as.character(Ctrl_Tcell$receptor)
Ctrl_Tcell$cell_from <- as.character(Ctrl_Tcell$cell_from)
Ctrl_Tcell$cell_to <- as.character(Ctrl_Tcell$cell_to)
Ctrl_Tcell$comm_type <- as.character(Ctrl_Tcell$comm_type)
cell_col <- structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b'),names=unique(Ctrl_Tcell$cell_from))
pdf("FigS6C_Ctrl_Tcell_LR.pdf",height=10,width=10)
LRPlot(Ctrl_Tcell,datatype='mean count',link.arr.lwd=Ctrl_Tcell$cell_from_mean_exprs,link.arr.width=Ctrl_Tcell$cell_to_mean_exprs)
dev.off()


ICI_Tcell <- subset(ICI_LR_filtered,cell_from=="T.cell" & cell_to=="Tumor")
ICI_Tcell$ligand <- as.character(ICI_Tcell$ligand)
ICI_Tcell$receptor <- as.character(ICI_Tcell$receptor)
ICI_Tcell$cell_from <- as.character(ICI_Tcell$cell_from)
ICI_Tcell$cell_to <- as.character(ICI_Tcell$cell_to)
ICI_Tcell$comm_type <- as.character(ICI_Tcell$comm_type)
cell_col <- structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b'),names=unique(ICI_Tcell$cell_from))
pdf("FigS6C_ICI_Tcell_LR.pdf",height=10,width=10)
LRPlot(ICI_Tcell,datatype='mean count',link.arr.lwd=ICI_Tcell$cell_from_mean_exprs,link.arr.width=ICI_Tcell$cell_to_mean_exprs)
dev.off()

OSI_Tcell <- subset(OSI_LR_filtered,cell_from=="T.cell" & cell_to=="Tumor")
OSI_Tcell$ligand <- as.character(OSI_Tcell$ligand)
OSI_Tcell$receptor <- as.character(OSI_Tcell$receptor)
OSI_Tcell$cell_from <- as.character(OSI_Tcell$cell_from)
OSI_Tcell$cell_to <- as.character(OSI_Tcell$cell_to)
OSI_Tcell$comm_type <- as.character(OSI_Tcell$comm_type)
cell_col <- structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b'),names=unique(OSI_Tcell$cell_from))
pdf("FigS6C_OSI_Tcell_LR.pdf",height=10,width=10)
LRPlot(OSI_Tcell,datatype='mean count',link.arr.lwd=OSI_Tcell$cell_from_mean_exprs,link.arr.width=OSI_Tcell$cell_to_mean_exprs)
dev.off()


ICI_OSI_Tcell <- subset(ICI_OSI_LR_filtered,cell_from=="T.cell" & cell_to=="Tumor")
ICI_OSI_Tcell$ligand <- as.character(ICI_OSI_Tcell$ligand)
ICI_OSI_Tcell$receptor <- as.character(ICI_OSI_Tcell$receptor)
ICI_OSI_Tcell$cell_from <- as.character(ICI_OSI_Tcell$cell_from)
ICI_OSI_Tcell$cell_to <- as.character(ICI_OSI_Tcell$cell_to)
ICI_OSI_Tcell$comm_type <- as.character(ICI_OSI_Tcell$comm_type)
cell_col <- structure(c('#4a84ad','#4a1dc6','#e874bf','#b79eed', '#ff636b', '#52c63b'),names=unique(ICI_OSI_Tcell$cell_from))
pdf("FigS6C_ICI_OSI_Tcell_LR.pdf",height=10,width=10)
LRPlot(ICI_OSI_Tcell,datatype='mean count',link.arr.lwd=ICI_OSI_Tcell$cell_from_mean_exprs,link.arr.width=ICI_OSI_Tcell$cell_to_mean_exprs)
dev.off()




#FigS6I,S6J,S6K

All_cells <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Add_Med_All_UMAP20_Renames.rds",mc.cores=20)
Idents(All_cells) <- All_cells$CellType
Macro <- subset(All_cells,ident ="Macro")

Macro_obj <- CreateSeuratObject(counts =GetAssayData(Macro, slot = "counts",assay="RNA")[,rownames(Macro@meta.data)], meta.data = Macro@meta.data)
Macro_MERGE_DATA <- Macro_obj %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
    ScaleData(verbose = TRUE,vars.to.regress =NULL)
Macro_MERGE_DATA.pca <- RunPCA(Macro_MERGE_DATA, npcs = 100, verbose = FALSE)
mcsaveRDS(Macro_MERGE_DATA.pca,file="Mice_Macro_MERGE_DATA_pca.rds",mc.cores=20)
JackStraw <- JackStraw(Macro_MERGE_DATA.pca, num.replicate = 100)
ScoreJackStraw <- ScoreJackStraw(JackStraw, dims = 1:20)
FindNeighbors <- FindNeighbors(ScoreJackStraw,dims = 1:20)

FindClusters <- FindClusters(FindNeighbors,resolution = c(0.1,0.2))
mcsaveRDS(FindClusters,file="Mice_Macro_MERGE_DATA_FindClusters.rds", mc.cores = 20)

Macro_MERGE_UMAP10 <- RunUMAP(object = FindClusters, dims = 1:10)
mcsaveRDS(Macro_MERGE_UMAP10,file="Mice_Macro_MERGE_UMAP10.rds", mc.cores = 10)

Idents(Macro_MERGE_UMAP10) <- Macro_MERGE_UMAP10$RNA_snn_res.0.1
ff <- DimPlot(object = Macro_MERGE_UMAP10, reduction = "umap",pt.size=1,
	label=TRUE,ncol=1,label.size=7,group.by="RNA_snn_res.0.1")
ggsave(ff,file="Macro_MERGE_UMAP10_res01.png",width=7.4,height=7)

Idents(Macro_MERGE_UMAP10) <- Macro_MERGE_UMAP10$RNA_snn_res.0.2
all.markers <- FindAllMarkers(Macro_MERGE_UMAP10, only.pos = TRUE, min.pct = 0, logfc.threshold = 0)
write.csv(all.markers,"FindAllmarker_Macro_MERGE_UMAP10_res02.csv")

markers_mouse <- as.character(c(
  "Hla-dra","Hla-drb1","Cd14","Fcgr3","Lyz","Fut4","Cd68","Cd36",
  "Cd163","Marco","Msr1","Junb","Jund","Fos","Tgfb1","Thbs1",
  "Csf1r","Tnfrsf1b","Itgam","Cd33","Adgre1","Vcan","Fcn1",
  "S100a8","S100a9","Nda","Fcer1g","Lst1","Aif1","Ifitm3",
  "Cd1c","Ppa1","Lsp1","Csf2ra","Id2","Gzmb",
  "Jchain","Irf7","Lilra4","Pld4","Tcf4","Il1b","Cxcl8",
  "Ccl3","Ccl4","Nfkbia","Sod2","Gpx1","C1qa","C1qb",
  "C1qc","Apoc1","Apoe","Csf3r","Nampt","Slc25a37",
  "Rnf149","Smchd1","Il1r2","Fam65b","G0s2",
  "Epcam","Chga","Ascl1"
))

Macro_MERGE_UMAP10 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Macro_endo_interaction/Mice_Macro_MERGE_UMAP10.rds",mc.cores=20)
Idents(Macro_MERGE_UMAP10) <- Macro_MERGE_UMAP10$RNA_snn_res.0.2
ff <- DotPlot(Macro_MERGE_UMAP10,features = markers_mouse,cols = c("#ffffff","#cb181d"),
assay="RNA",col.min = 0.1,col.max =2,dot.min=0.05,scale=TRUE) + RotatedAxis()
ggsave(ff,file="Dotplot_qc_marker.pdf",width=18,height=5)


Macro_MERGE_UMAP10 <- mcreadRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Macro_endo_interaction/Mice_Macro_MERGE_UMAP10.rds",mc.cores=20)
Idents(Macro_MERGE_UMAP10) <- Macro_MERGE_UMAP10$RNA_snn_res.0.2
pure_Macro <- subset(Macro_MERGE_UMAP10,idents=c("4","5","6"),invert=TRUE)

FindNeighbors <- FindNeighbors(pure_Macro,dims = 1:20)
FindClusters <- FindClusters(FindNeighbors,resolution = c(0.1,0.2))
Macro_pure_UMAP10 <- RunUMAP(object = FindClusters, dims = 1:10)


Macro_pure_UMAP10$treat <- Macro_pure_UMAP10$group
Macro_pure_UMAP10$treat[Macro_pure_UMAP10$group %in% c("ICI", "antiPD1")] <- "antiPD1"


Idents(Macro_pure_UMAP10) <- Macro_pure_UMAP10$RNA_snn_res.0.2
Macro_pure_UMAP10$treat <- factor(Macro_pure_UMAP10$treat,levels = c("Ctrl","antiPD1","OSI", "OSI_antiPD1"))
ff <- DimPlot(object = Macro_pure_UMAP10, reduction = "umap",pt.size=1,cols=c("#a50f15",
"#045a8d","#c2a5cf","#fb6a4a","#3288bd",
"#525252","#01665e"),label=TRUE,label.size=0)
ggsave(ff, file="Final_SCLC_Macro_pure_UMAP10.png",width=5.2,height=4.6,dpi=1080)


Idents(Macro_pure_UMAP10) <- Macro_pure_UMAP10$RNA_snn_res.0.2
Macro_pure_UMAP10$treat <- factor(Macro_pure_UMAP10$treat,levels = c("Ctrl","antiPD1","OSI", "OSI_antiPD1"))
all.markers <- FindAllMarkers(Macro_pure_UMAP10, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
write.csv(all.markers,"final_Pure_Macro_UMAP15_res02.csv")

Macro_pure_UMAP10$RNA_snn_res.0.2_merge <- Macro_pure_UMAP10$RNA_snn_res.0.2
Macro_pure_UMAP10$RNA_snn_res.0.2_merge[Macro_pure_UMAP10$RNA_snn_res.0.2 %in% c(1, 2)] <- 1
Idents(Macro_pure_UMAP10) <- Macro_pure_UMAP10$RNA_snn_res.0.2_merge
Macro_pure_UMAP10$treat <- factor(Macro_pure_UMAP10$treat,levels = c("Ctrl","antiPD1","OSI", "OSI_antiPD1"))
ff <- DimPlot(object = Macro_pure_UMAP10, reduction = "umap",pt.size=1,
	cols=c("#a50f15","#045a8d","#c2a5cf","#fb6a4a","#3288bd","#525252","#01665e"),label=TRUE,label.size=0)
ggsave(ff, file="FigS6I.png",width=5.2,height=4.6,dpi=1080)



library(dplyr)
library(ggplot2)

tmp <- Macro_pure_UMAP10@meta.data
tmp$treat <- factor(tmp$treat,levels = c("Ctrl","OSI","antiPD1", "OSI_antiPD1"))
subtype_ratio <- tmp %>% group_by(treat, RNA_snn_res.0.2_merge) %>%
  summarise(cell_n = n()) %>%                     
  group_by(treat) %>%
  mutate(total_cells = sum(cell_n),
         ratio = cell_n / total_cells) %>%         
  ungroup()
subtype_ratio <- data.frame(subtype_ratio)
ratio_df <- subtype_ratio[,c("treat","cell_n","RNA_snn_res.0.2_merge")]
library(ggplot2)
ratio_df$treat <- factor(ratio_df$treat,levels = c("Ctrl", "OSI", "antiPD1", "OSI_antiPD1"))
ff <- ggplot(ratio_df,aes(x = treat,stratum   = RNA_snn_res.0.2_merge,alluvium  = RNA_snn_res.0.2_merge,y = cell_n,fill = RNA_snn_res.0.2_merge)) +
  geom_flow(alpha = 0.5, width = 0.7) +   
  geom_stratum(width = 0.7, color = "black") +  
  scale_fill_manual(values = c("#a50f15", "#045a8d", "#c2a5cf", "#fb6a4a")) +
  theme_classic() +
  labs(x = NULL,y = "Number of cells",fill = "Cluster")
ggsave("FigS6J.png",ff,width = 5,height = 5)







Idents(Macro_pure_UMAP10) <- Macro_pure_UMAP10$RNA_snn_res.0.2_merge
Macro_0 <- subset(Macro_pure_UMAP10,idents=c("0"))
Macro_1 <- subset(Macro_pure_UMAP10,idents=c("1"))
Macro_3 <- subset(Macro_pure_UMAP10,idents=c("3"))
Macro_4 <- subset(Macro_pure_UMAP10,idents=c("4"))

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

Macro_0_pseudobulk <- pseudo_bulk_seurat_mean(seurat_obj=Macro_0,num_split=10,seed.use=1,slot="data",prefix="Macro_0")
Macro_1_pseudobulk <- pseudo_bulk_seurat_mean(seurat_obj=Macro_1,num_split=10,seed.use=1,slot="data",prefix="Macro_1")
Macro_3_pseudobulk <- pseudo_bulk_seurat_mean(seurat_obj=Macro_3,num_split=10,seed.use=1,slot="data",prefix="Macro_3")
Macro_4_pseudobulk <- pseudo_bulk_seurat_mean(seurat_obj=Macro_4,num_split=10,seed.use=1,slot="data",prefix="Macro_3")

All_pseudobulk <- cbind(Macro_0_pseudobulk,
Macro_1_pseudobulk,
Macro_3_pseudobulk,Macro_4_pseudobulk)
mcsaveRDS(All_pseudobulk,file="Macrophage_All_pseudobulk_bin10_renew.rds",mc.cores=30)


all_pseudobulk <- readRDS("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Macro_endo_interaction/Macrophage_All_pseudobulk_bin10_renew.rds")
all_pseudobulk <- data.frame(all_pseudobulk)
aa <- all_pseudobulk %>% rownames() %>% convert_mouse_to_human_symbols()
all_pseudobulk$symbol <- as.character(aa)
matrix <- all_pseudobulk[!duplicated(all_pseudobulk$symbol),]
matrix <- na.omit(matrix)
rownames(matrix) <- matrix$symbol
matrix <- matrix[,c(-41)]

matrix <- as.matrix(matrix)
h_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/h.all.v7.1.symbols.gmt")
c5_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/c5.all.v7.1.symbols.gmt")
c2_geneSets <- getGmt("/mnt/data/user_data/abao/2_reference/GSVA_7.1/c2.all.v7.1.symbols.gmt")

h_GSVA_res <- GSVA::gsva(matrix, h_geneSets, min.sz=10, max.sz=500, verbose=FALSE)
c2_GSVA_res <- GSVA::gsva(matrix, c2_geneSets, min.sz=10, max.sz=500, verbose=FALSE)
c5_GSVA_res <- GSVA::gsva(matrix, c5_geneSets, min.sz=10, max.sz=500, verbose=FALSE)
sudobulk_GSVA_score <- rbind(h_GSVA_res,c2_GSVA_res,c5_GSVA_res)
write.csv(sudobulk_GSVA_score,"Macrophage_AllCells_GSVA_score_bin10_renew.csv")



annotation0 <- data.frame(c(1:10))
annotation0$group <- "M0"
annotation1 <- data.frame(c(11:20))
annotation1$group <- "M1"
annotation3 <- data.frame(c(31:40))
annotation3$group <- "M3"
annotation4 <- data.frame(c(41:50))
annotation4$group <- "M4"

names(annotation4) <- c("order","group")
names(annotation3) <- c("order","group")
names(annotation1) <- c("order","group")
names(annotation0) <- c("order","group")
annotation <- rbind(annotation0,annotation1,annotation3,annotation4)
group <-as.factor(annotation$group)
design <- model.matrix(~ group + 0)
rownames(design)<-colnames(matrix)
head(design)

contrasts <- makeContrasts(groupM0-(groupM1+groupM3+groupM4)/3,
	groupM1 -(groupM0+groupM3+groupM4)/3,
    groupM3-(groupM0+groupM1+groupM4)/3,
    groupM4-(groupM0+groupM1+groupM3)/3,
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
write.csv(GSVA_hc2c5,"Macrophage_p005_AllCells_GSVA_score_bin10_renew.csv")

#9122
GSVA_pathway <- read.csv("/mnt/data/user_data/ailing/SCLC_IEV_code/Mouse_Model_scRNAseq/Luyou_Med_as_control_OSI906/Macro_endo_interaction/Macrophage_p005_AllCells_GSVA_score_bin10_renew.csv")
rownames(GSVA_pathway) <- GSVA_pathway$X
select_pathways <- subset(GSVA_pathway,
X=="GO_ACTIVATION_OF_IMMUNE_RESPONSE" | 
X=="GO_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY" |
X=="GO_CYTOKINE_SECRETION_INVOLVED_IN_IMMUNE_RESPONSE" |
X=="HALLMARK_INFLAMMATORY_RESPONSE")

select_pathways <- data.frame(select_pathways[,2:5])

select_pathways <- t(apply(select_pathways, 1, function(x) (x-mean(x))/sd(x)))

range(select_pathways)
library(pheatmap)
select_pathways[select_pathways > 1] <- 1
select_pathways[select_pathways < -1] <- -1

library(pheatmap)
library(dplyr)

pdf("FigS6K.pdf",width=8,height=8)
pheatmap(select_pathways,clustering_method="ward.D2",
  color = colorRampPalette(c("#2971B1","#6AACD0","#C1DDEB","#fed976","#fe9929","#E58267","#BB2933"))(50),
  fontsize_row=12,show_rownames=TRUE,show_colnames=TRUE,cluster_row =TRUE,cluster_col= TRUE,border=FALSE)
dev.off()



















































