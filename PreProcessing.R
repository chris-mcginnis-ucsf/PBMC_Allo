#########################################################
## RNA, MULTI-seq and SCMK Pre-Processing for 8-donor, ##
## 7-donor, Zheng et al, and BD PBMC +/- CD3/CD28 data ##
## Chris McGinnis, Gartner Lab, UCSF, 06/26/2020 ########
#########################################################

library(Seurat)
library(deMULTIplex)
library(DoubletFinder)

##################
## 8-Donor PBMC ##
##################
## Step 1: Define putative cell IDs -------------------------------------------------------------------------------------------------------------------------------------------------------------------
data_l1 <- Read10X('./raw_feature_bc_matrix_L1/')
cell.counts <- Matrix::colSums(data_l1)
cellIDs_l1 <- colnames(data_l1)[which(cell.counts >= 250)]  ## 5,001 cells

data_l2 <- Read10X('./raw_feature_bc_matrix_L2/')
cell.counts <- Matrix::colSums(data_l2)
cellIDs_l2 <- colnames(data_l2)[which(cell.counts >= 250)]  ## 4,326 cells

data_l3 <- Read10X('./raw_feature_bc_matrix_L3/')
cell.counts <- Matrix::colSums(data_l3)
cellIDs_l3 <- colnames(data_l3)[which(cell.counts >= 250)] ## 5,159 cells

data_l4 <- Read10X('./raw_feature_bc_matrix_L4/')
cell.counts <- Matrix::colSums(data_l4)
cellIDs_l4 <- colnames(data_l4)[which(cell.counts >= 250)] ## 5,867 cells


## Step 2: Align MULTI-seq and HashTag FASTQs ---------------------------------------------------------------------------------------------------------------------------------------------------------
## MULTI-seq
readTable_l2_multi <- MULTIseq.preProcess(R1="TGACCA_S2_L001_R1_001.fastq.gz", 
                                          R2="TGACCA_S2_L001_R2_001.fastq.gz",
                                          cellIDs = cellIDs_l2, cell=c(1,16), umi=c(17,28), tag=c(1,8)) 

readTable_l3_multi <- MULTIseq.preProcess(R1="ACAGTG_S3_L001_R1_001.fastq.gz", 
                                          R2="ACAGTG_S3_L001_R2_001.fastq.gz",
                                          cellIDs = cellIDs_l3, cell=c(1,16), umi=c(17,28), tag=c(1,8)) 

bar.ref <- read.csv("LMOlist.csv", header=F, stringsAsFactors=F) 
bar.ref <- bar.ref[9:16,1]
barTable_l2_multi <- MULTIseq.align(readTable_l2_multi, cellIDs_l2, bar.ref)
barTable_l3_multi <- MULTIseq.align(readTable_l3_multi, cellIDs_l3, bar.ref)

## Cell Hashing
readTable_l2_hash <- MULTIseq.preProcess_beta(R1='SL-2792-06_e06_S12_L002_R1_001.fastq.gz', 
                                              R2='SL-2792-06_e06_S12_L002_R2_001.fastq.gz', 
                                              cellIDs = cellIDs_l2, cell=c(1,16), umi=c(17,26), tag=c(26,70))
readTable_l3_hash <- MULTIseq.preProcess_beta(R1='SL-2792-07_e07_S13_L002_R1_001.fastq.gz', 
                                              R2='SL-2792-07_e07_S13_L002_R2_001.fastq.gz', 
                                              cellIDs = cellIDs_l3, cell=c(1,16), umi=c(17,26), tag=c(26,70))

bar.ref_hash <- c("ATTCAAGGGCAGCCGCGTCACGATTGGATACGACTGTTGGACCGG",
                  "TGGATGGGATAAGTGCGTGATGGACCGAAGGGACCTCGTGGCCGG",
                  "CGGCTCGTGCTGCGTCGTCTCAAGTCCAGAAACTCCGTGTATCCT",
                  "ATTGGGAGGCTTTCGTACCGCTGCCGCCACCAGGTGATACCCGCT",
                  "GTGATCCGCGCAGGCACACATACCGACTCAGATGGGTTGTCCAGG",
                  "GCAGCCGGCGTCGTACGAGGCACAGCGGAGACTAGATGAGGCCCC",
                  "CGCGTCCAATTTCCGAAGCCCCGCCCTAGGAGTTCCCCTGCGTGC",
                  "GCCCATTCATTGCACCCGCCAGTGATCGACCCTAGTGGAGCTAAG")
barTable_l2_hash <- MULTIseq.align_BD(readTable_l2_hash, cellIDs_l2, bar.ref_hash)
barTable_l3_hash <- MULTIseq.align_BD(readTable_l3_hash, cellIDs_l3, bar.ref_hash)

## Sample classification workflow continued in "SampleClassification.R"

## Step 3: Remove uninformative genes in read-depth normalized data -----------------------------------------------------------------------------------------------------------------------------------
data_aggr <- Read10X('/Volumes/MULTIseq_1/SLNR_PBMC/aggrSL/outs/raw_feature_bc_matrix/')
cellIDs <- c(cellIDs_l1,cellIDs_l2,cellIDs_l3,cellIDs_l4)
data_aggr <- data_aggr[,cellIDs]
gene.counts <- Matrix::rowSums(data_aggr)
data_aggr <- data_aggr[which(gene.counts >= 3), ]


## Step 4: Pre-process Seurat object ------------------------------------------------------------------------------------------------------------------------------------------------------------------
seu_pbmc_raw <- CreateSeuratObject(data_aggr)
mito.genes <- grep("^MT-", rownames(seu_pbmc_raw@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(seu_pbmc_raw@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(seu_pbmc_raw@assays$RNA@counts)
seu_pbmc_raw@meta.data[,"PercentMito"] <- percent.mito
seu_pbmc_raw <- SCTransform(seu_pbmc_raw)
seu_pbmc_raw <- RunPCA(seu_pbmc_raw)
seu_pbmc_raw <- RunUMAP(seu_pbmc_raw, dims = 1:13)
seu_pbmc_raw <- FindNeighbors(seu_pbmc_raw, dims = 1:13)
seu_pbmc_raw <- FindClusters(seu_pbmc_raw, resolution = 1.0)


## Step 5: Remove low-quality cells (high pMito, low nUMI) --------------------------------------------------------------------------------------------------------------------------------------------
'%ni%' <- Negate('%in%')
seu_pbmc_raw_2 <- SubsetData(seu_pbmc_raw, cells=rownames(seu_pbmc_raw@meta.data)[which(seu_pbmc_raw@active.ident %ni% c(4,10,13,18:20,23,24))]) # Remove 3726 cells
seu_pbmc_raw_2 <- SCTransform(seu_pbmc_raw_2)
seu_pbmc_raw_2 <- RunPCA(seu_pbmc_raw_2)
seu_pbmc_raw_2 <- RunUMAP(seu_pbmc_raw_2, dims = 1:12)
seu_pbmc_raw_2 <- FindNeighbors(seu_pbmc_raw_2, dims = 1:12)
seu_pbmc_raw_2 <- FindClusters(seu_pbmc_raw_2, resolution = 1.0)


## Step 6: Before applying DoubletFinder, subset data by lane -----------------------------------------------------------------------------------------------------------------------------------------
seu_pbmc_raw_2@meta.data[,"LaneID"] <- rep(1, nrow(seu_pbmc_raw_2@meta.data))
seu_pbmc_raw_2@meta.data[grep("-2", rownames(seu_pbmc_raw_2@meta.data)), "LaneID"] <- 2
seu_pbmc_raw_2@meta.data[grep("-3", rownames(seu_pbmc_raw_2@meta.data)), "LaneID"] <- 3
seu_pbmc_raw_2@meta.data[grep("-4", rownames(seu_pbmc_raw_2@meta.data)), "LaneID"] <- 4

## Lane 1
seu_pbmc_raw_l1 <- SubsetData(seu_pbmc_raw_2, cells=rownames(seu_pbmc_raw_2@meta.data)[which(seu_pbmc_raw_2@meta.data$LaneID == 1)]) 
seu_pbmc_raw_l1 <- SCTransform(seu_pbmc_raw_l1)
seu_pbmc_raw_l1 <- RunPCA(seu_pbmc_raw_l1)
seu_pbmc_raw_l1 <- RunUMAP(seu_pbmc_raw_l1, dims = 1:10)
seu_pbmc_raw_l1 <- FindNeighbors(seu_pbmc_raw_l1, dims = 1:10)
seu_pbmc_raw_l1 <- FindClusters(seu_pbmc_raw_l1, resolution = 0.5)

## Lane 2
seu_pbmc_raw_l2 <- SubsetData(seu_pbmc_raw_2, cells=rownames(seu_pbmc_raw_2@meta.data)[which(seu_pbmc_raw_2@meta.data$LaneID == 2)]) 
seu_pbmc_raw_l2 <- SCTransform(seu_pbmc_raw_l2)
seu_pbmc_raw_l2 <- RunPCA(seu_pbmc_raw_l2)
seu_pbmc_raw_l2 <- RunUMAP(seu_pbmc_raw_l2, dims = 1:11)
seu_pbmc_raw_l2 <- FindNeighbors(seu_pbmc_raw_l2, dims = 1:11)
seu_pbmc_raw_l2 <- FindClusters(seu_pbmc_raw_l2, resolution = 0.5)

## Lane 3
seu_pbmc_raw_l3 <- SubsetData(seu_pbmc_raw_2, cells=rownames(seu_pbmc_raw_2@meta.data)[which(seu_pbmc_raw_2@meta.data$LaneID == 3)]) 
seu_pbmc_raw_l3 <- SCTransform(seu_pbmc_raw_l3)
seu_pbmc_raw_l3 <- RunPCA(seu_pbmc_raw_l3)
seu_pbmc_raw_l3 <- RunUMAP(seu_pbmc_raw_l3, dims = 1:13)
seu_pbmc_raw_l3 <- FindNeighbors(seu_pbmc_raw_l3, dims = 1:13)
seu_pbmc_raw_l3 <- FindClusters(seu_pbmc_raw_l3, resolution = 0.5)

## Lane 4
seu_pbmc_raw_l4 <- SubsetData(seu_pbmc_raw_2, cells=rownames(seu_pbmc_raw_2@meta.data)[which(seu_pbmc_raw_2@meta.data$LaneID == 4)]) 
seu_pbmc_raw_l4 <- SCTransform(seu_pbmc_raw_l4)
seu_pbmc_raw_l4 <- RunPCA(seu_pbmc_raw_l4)
seu_pbmc_raw_l4 <- RunUMAP(seu_pbmc_raw_l4, dims = 1:10)
seu_pbmc_raw_l4 <- FindNeighbors(seu_pbmc_raw_l4, dims = 1:10)
seu_pbmc_raw_l4 <- FindClusters(seu_pbmc_raw_l4, resolution = 1.0)


## Step 7: Identify heterotypic doublets in each lane using DoubletFinder -----------------------------------------------------------------------------------------------------------------------------
## Lane 1
sweep.res.list_pbmc_l1 <- paramSweep_v3(seu_pbmc_raw_l1, PCs = 1:10, sct = TRUE)
sweep.stats_pbmc_l1 <- summarizeSweep(sweep.res.list_pbmc_l1, GT = FALSE)
bcmvn_pbmc_l1 <- find.pK(sweep.stats_pbmc_l1)
homotypic.prop.l1 <- modelHomotypic(seu_pbmc_raw_l1@active.ident)           
nExp_poi.l1 <- round(0.04*nrow(seu_pbmc_raw_l1@meta.data))  
nExp_poi.adj.l1 <- round(nExp_poi.l1*(1-homotypic.prop.l1))
seu_pbmc_raw_l1 <- doubletFinder_v3(seu_pbmc_raw_l1, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj.l1, reuse.pANN = FALSE, sct=T)

## Lane 2
sweep.res.list_pbmc_l2 <- paramSweep_v3(seu_pbmc_raw_l2, PCs = 1:11, sct = TRUE)
sweep.stats_pbmc_l2 <- summarizeSweep(sweep.res.list_pbmc_l2, GT = FALSE)
bcmvn_pbmc_l2 <- find.pK(sweep.stats_pbmc_l2)
homotypic.prop.l2 <- modelHomotypic(seu_pbmc_raw_l2@active.ident)           
nExp_poi.l2 <- round(0.04*nrow(seu_pbmc_raw_l2@meta.data))  
nExp_poi.adj.l2 <- round(nExp_poi.l2*(1-homotypic.prop.l2))
seu_pbmc_raw_l2 <- doubletFinder_v3(seu_pbmc_raw_l2, PCs = 1:11, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj.l2, reuse.pANN = FALSE, sct=T)

## Lane 3
sweep.res.list_pbmc_l3 <- paramSweep_v3(seu_pbmc_raw_l3, PCs = 1:13, sct = TRUE)
sweep.stats_pbmc_l3 <- summarizeSweep(sweep.res.list_pbmc_l3, GT = TRUE, GT.calls = final.calls_l3_temp)
bcmvn_pbmc_l3 <- find.pK(sweep.stats_pbmc_l3)
homotypic.prop.l3 <- modelHomotypic(seu_pbmc_raw_l3@active.ident)           
nExp_poi.l3 <- round(0.04*nrow(seu_pbmc_raw_l3@meta.data))  
nExp_poi.adj.l3 <- round(nExp_poi.l3*(1-homotypic.prop.l3))
seu_pbmc_raw_l3 <- doubletFinder_v3(seu_pbmc_raw_l3, PCs = 1:11, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj.l3, reuse.pANN = FALSE, sct=T)

## Lane 4
sweep.res.list_pbmc_l4 <- paramSweep_v3(seu_pbmc_raw_l4, PCs = 1:10, sct = TRUE)
sweep.stats_pbmc_l4 <- summarizeSweep(sweep.res.list_pbmc_l4, GT = FALSE)
bcmvn_pbmc_l4 <- find.pK(sweep.stats_pbmc_l4)
homotypic.prop.l4 <- modelHomotypic(seu_pbmc_raw_l4@active.ident)           
nExp_poi.l4 <- round(0.04*nrow(seu_pbmc_raw_l4@meta.data))  
nExp_poi.adj.l4 <- round(nExp_poi.l4*(1-homotypic.prop.l4))
seu_pbmc_raw_l4 <- doubletFinder_v3(seu_pbmc_raw_l4, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi.adj.l4, reuse.pANN = FALSE, sct=T)


## Step 8: Remove doublets, pre-process cleaned merged Seurat object ----------------------------------------------------------------------------------------------------------------------------------
df.doublets <- c(rownames(seu_pbmc_raw_l1@meta.data)[which(seu_pbmc_raw_l1@meta.data$DF.classifications_0.25_0.01_133 == "Doublet")],
                 rownames(seu_pbmc_raw_l2@meta.data)[which(seu_pbmc_raw_l2@meta.data$DF.classifications_0.25_0.01_132 == "Doublet")],
                 rownames(seu_pbmc_raw_l3@meta.data)[which(seu_pbmc_raw_l3@meta.data$DF.classifications_0.25_0.01_148 == "Doublet")],
                 rownames(seu_pbmc_raw_l4@meta.data)[which(seu_pbmc_raw_l4@meta.data$DF.classifications_0.25_0.01_168 == "Doublet")])
seu_pbmc_raw_2@meta.data[,"DF"] <- rep("Singlet")
seu_pbmc_raw_2@meta.data[df.doublets,"DF"] <- "Doublet"

seu_pbmc_clean <- SubsetData(seu_pbmc_raw_2, cells=rownames(seu_pbmc_raw_2@meta.data)[which(seu_pbmc_raw_2@meta.data$DF == "Singlet")])
seu_pbmc_clean <- SCTransform(seu_pbmc_clean)
seu_pbmc_clean <- RunPCA(seu_pbmc_clean)
seu_pbmc_clean <- RunUMAP(seu_pbmc_clean, dims = 1:12)
seu_pbmc_clean <- FindNeighbors(seu_pbmc_clean, dims = 1:12)
seu_pbmc_clean <- FindClusters(seu_pbmc_clean, resolution = 1.5)


## Step 9: Annotate Cell Types ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seu_pbmc_clean@meta.data[,"CellType"] <- rep("unknown", nrow(seu_pbmc_clean@meta.data))
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(8,21))] <- "B"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(2,14,24))] <- "NK"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(4,6,15,16))] <- "CD14Mono"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident == 20)] <- "DC"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident == 25)] <- "pDC"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(13,22))] <- "CD16Mono"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(0,1,3,5,7,9,12,17,18,23))] <- "CD4T"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(10,11,19))] <- "CD8T"


## Step 10: Subset Ficoll-prepped PBMCs and CD4+ T-cells ----------------------------------------------------------------------------------------------------------------------------------------------
seu_pbmc_ficoll <- SubsetData(seu_pbmc_clean, cells = rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$Donor %in% c("A","B","C"))])
seu_pbmc_ficoll <- SCTransform(seu_pbmc_ficoll)
seu_pbmc_ficoll <- RunPCA(seu_pbmc_ficoll)
seu_pbmc_ficoll <- RunUMAP(seu_pbmc_ficoll, dims = 1:12)
seu_pbmc_ficoll <- FindNeighbors(seu_pbmc_ficoll, dims = 1:12)
seu_pbmc_ficoll <- FindClusters(seu_pbmc_ficoll, resolution = 0.8)

seu_CD4T_ficoll <- SubsetData(seu_pbmc_clean, cells = rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$CellType == "CD4T" & seu_pbmc_clean@meta.data$Donor %in% c("A","B","C"))])
seu_CD4T_ficoll <- SCTransform(seu_CD4T_ficoll)
seu_CD4T_ficoll <- RunPCA(seu_CD4T_ficoll)
seu_CD4T_ficoll <- RunUMAP(seu_CD4T_ficoll, dims = 1:25)
seu_CD4T_ficoll <- FindNeighbors(seu_CD4T_ficoll, dims = 1:25)
seu_CD4T_ficoll <- FindClusters(seu_CD4T_ficoll, resolution = 1.0)


## Step 11: Define CD4+ T-cell subtypes ---------------------------------------------------------------------------------------------------------------------------------------------------------------
seu_CD4T_ficoll@meta.data[,"CD4Tsubset"] <- rep("unknown")
seu_CD4T_ficoll@meta.data$CD4Tsubset[which(seu_CD4T_ficoll@active.ident %in% c(0,2,6,10,11))] <- "Activated"
seu_CD4T_ficoll@meta.data$CD4Tsubset[which(seu_CD4T_ficoll@active.ident %in% c(1,3,5,7:9))] <- "Memory"
seu_CD4T_ficoll@meta.data$CD4Tsubset[which(seu_CD4T_ficoll@active.ident == 4)] <- "Naive"

######################
## Zheng et al PBMC ##
######################
## Step 1: Remove uninformative genes and cell BCs from read-depth normalized Zheng et al (Nat Comm, 2017) dataset ------------------------------------------------------------------------------------
data_zheng <- Read10X('/Users/gartnerlab/Desktop/SLNR_PBMC/REVIEW/zheng/raw_feature_bc_matrix/')
cell.counts <- Matrix::colSums(data_zheng)
data_zheng <- data_zheng[,which(cell.counts >= 250)]
gene.counts <- Matrix::rowSums(data_zheng)
data_zheng <- data_zheng[which(gene.counts > 3), ]


## Step 2: Pre-process Seurat object, remove high %Mito cells -----------------------------------------------------------------------------------------------------------------------------------------
seu_zheng <- CreateSeuratObject(data_zheng)
seu_zheng <- SCTransform(seu_zheng)
seu_zheng <- RunPCA(seu_zheng)
seu_zheng <- RunUMAP(seu_zheng, dims = 1:21)
seu_zheng <- FindNeighbors(seu_zheng, dims = 1:21)
seu_zheng <- FindClusters(seu_zheng, resolution = 0.5)
mito.genes <- grep("MT-", rownames(seu_zheng@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(seu_zheng@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(seu_zheng@assays$RNA@counts)
seu_zheng@meta.data[,"PercentMito"] <- percent.mito

seu_zheng_2 <- SubsetData(seu_zheng, cells=rownames(seu_zheng@meta.data)[which(seu_zheng@active.ident != 9)])
seu_zheng_2 <- SCTransform(seu_zheng_2)
seu_zheng_2 <- RunPCA(seu_zheng_2)
seu_zheng_2 <- RunUMAP(seu_zheng_2, dims = 1:20)
seu_zheng_2 <- FindNeighbors(seu_zheng_2, dims = 1:20)
seu_zheng_2 <- FindClusters(seu_zheng_2, resolution = 0.5)
seu_zheng_2@meta.data[,"LaneID"] <- 1L
seu_zheng_2@meta.data$LaneID[grep("-2",rownames(seu_zheng_2@meta.data))] <- 2
seu_zheng_2@meta.data$LaneID[grep("-3",rownames(seu_zheng_2@meta.data))] <- 3


## Step 3: Before applying DoubletFinder, subset data by lane -----------------------------------------------------------------------------------------------------------------------------------------
seu_zheng_XY <- SubsetData(seu_zheng_2, cells=rownames(seu_zheng_2@meta.data)[which(seu_zheng_2@meta.data$LaneID == 1)])
seu_zheng_XY <- SCTransform(seu_zheng_XY)
seu_zheng_XY <- RunPCA(seu_zheng_XY)
seu_zheng_XY <- RunUMAP(seu_zheng_XY, dims = 1:18)
seu_zheng_XY <- FindNeighbors(seu_zheng_XY, dims = 1:18)
seu_zheng_XY <- FindClusters(seu_zheng_XY, resolution = 0.5)

seu_zheng_X <- SubsetData(seu_zheng_2, cells=rownames(seu_zheng_2@meta.data)[which(seu_zheng_2@meta.data$LaneID == 2)])
seu_zheng_X <- SCTransform(seu_zheng_X)
seu_zheng_X <- RunPCA(seu_zheng_X)
seu_zheng_X <- RunUMAP(seu_zheng_X, dims = 1:16)
seu_zheng_X <- FindNeighbors(seu_zheng_X, dims = 1:16)
seu_zheng_X <- FindClusters(seu_zheng_X, resolution = 0.5)

seu_zheng_Y <- SubsetData(seu_zheng_2, cells=rownames(seu_zheng_2@meta.data)[which(seu_zheng_2@meta.data$LaneID == 3)])
seu_zheng_Y <- SCTransform(seu_zheng_Y)
seu_zheng_Y <- RunPCA(seu_zheng_Y)
seu_zheng_Y <- RunUMAP(seu_zheng_Y, dims = 1:18)
seu_zheng_Y <- FindNeighbors(seu_zheng_Y, dims = 1:18)
seu_zheng_Y <- FindClusters(seu_zheng_Y, resolution = 0.5)


## Step 4: Identify heterotypic doublets in each lane using DoubletFinder -----------------------------------------------------------------------------------------------------------------------------
sweep.res.list_X <- paramSweep_v3(seu_zheng_X, PCs = 1:16, sct = TRUE)
sweep.stats_X <- summarizeSweep(sweep.res.list_X, GT = FALSE)
bcmvn_X <- find.pK(sweep.stats_X) 
nExp_poi <- round(0.07*nrow(seu_zheng_X@meta.data))  
seu_zheng_X <- doubletFinder_v3(seu_zheng_X, PCs = 1:16, pN = 0.25, pK = 0.07, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

sweep.res.list_Y <- paramSweep_v3(seu_zheng_Y, PCs = 1:18, sct = TRUE)
sweep.stats_Y <- summarizeSweep(sweep.res.list_Y, GT = FALSE)
bcmvn_Y <- find.pK(sweep.stats_Y) 
nExp_poi <- round(0.07*nrow(seu_zheng_Y@meta.data)) 
seu_zheng_Y <- doubletFinder_v3(seu_zheng_Y, PCs = 1:18, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)

sweep.res.list_XY <- paramSweep_v3(seu_zheng_XY, PCs = 1:18, sct = TRUE)
sweep.stats_XY <- summarizeSweep(sweep.res.list_XY, GT = FALSE)
bcmvn_XY <- find.pK(sweep.stats_XY)
nExp_poi <- round(0.07*nrow(seu_zheng_XY@meta.data))  
seu_zheng_XY <- doubletFinder_v3(seu_zheng_XY, PCs = 1:18, pN = 0.25, pK = 0.08, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)


## Step 5: Record DoubletFinder results, remove doublets (also souporcell), and low-frequency cell types  ---------------------------------------------------------------------------------------------
seu_zheng_2@meta.data[,"DF"] <- rep("unknown")
seu_zheng_2@meta.data[rownames(seu_zheng_X@meta.data),"DF"] <- seu_zheng_X@meta.data$DF.classifications_0.25_0.07_556
seu_zheng_2@meta.data[rownames(seu_zheng_Y@meta.data),"DF"] <- seu_zheng_Y@meta.data$DF.classifications_0.25_0.09_681
seu_zheng_2@meta.data[rownames(seu_zheng_XY@meta.data),"DF"] <- seu_zheng_XY@meta.data$DF.classifications_0.25_0.08_595

seu_zheng_clean <- SubsetData(seu_zheng_2, cells=rownames(seu_zheng_2@meta.data)[which(seu_zheng_2@meta.data$DF == "Singlet")])
seu_zheng_clean <- SCTransform(seu_zheng_clean)
seu_zheng_clean <- RunPCA(seu_zheng_clean)
seu_zheng_clean <- RunUMAP(seu_zheng_clean, dims = 1:21)
seu_zheng_clean <- FindNeighbors(seu_zheng_clean, dims = 1:21)
seu_zheng_clean <- FindClusters(seu_zheng_clean, resolution = 0.75)

seu_zheng_clean <- SubsetData(seu_zheng_clean, cells=rownames(seu_zheng_clean@meta.data)[which(seu_zheng_clean@meta.data$CellType != "low_frequency" & seu_zheng_clean@meta.data$Donor != "Doublet")])
seu_zheng_clean <- SCTransform(seu_zheng_clean)
seu_zheng_clean <- RunPCA(seu_zheng_clean)
seu_zheng_clean <- RunUMAP(seu_zheng_clean, dims = 1:16)
seu_zheng_clean <- FindNeighbors(seu_zheng_clean, dims = 1:16)
seu_zheng_clean <- FindClusters(seu_zheng_clean, resolution = 0.5)


## Step 6: Annotate cell types ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seu_zheng_clean@meta.data[,"CellType"] <- rep("low_frequency")
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident %in% c(0:2,4,5,11,18))] <- "CD4T"
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident %in% c(3,10,17))] <- "CD8T"
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident == 8)] <- "NK"
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident == 16)] <- "Platelet"
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident %in% c(7,9))] <- "CD14Mono"
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident == 12)] <- "CD16Mono"
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident == 13)] <- "DC"
seu_zheng_clean@meta.data$CellType[which(seu_zheng_clean@active.ident %in% c(6,14))] <- "B"


## Step 7: Subset CD4+ T-cells, define subtypes -------------------------------------------------------------------------------------------------------------------------------------------------------
seu_zheng_CD4T <- SubsetData(seu_zheng_clean, cells=rownames(seu_zheng_clean@meta.data)[which(seu_zheng_clean@meta.data$CellType == "CD4T")])
seu_zheng_CD4T <- SCTransform(seu_zheng_CD4T)
seu_zheng_CD4T <- RunPCA(seu_zheng_CD4T)
seu_zheng_CD4T <- RunUMAP(seu_zheng_CD4T, dims = 1:20)
seu_zheng_CD4T <- FindNeighbors(seu_zheng_CD4T, dims = 1:20)
seu_zheng_CD4T <- FindClusters(seu_zheng_CD4T, resolution = 0.75)

seu_zheng_CD4T@meta.data[,"Subtype"] <- rep("memory")
seu_zheng_CD4T@meta.data$Subtype[which(seu_zheng_CD4T@active.ident %in% c(0,1,3))] <- "activated"
seu_zheng_CD4T@meta.data$Subtype[which(seu_zheng_CD4T@active.ident %in% c(5,7))] <- "naive"


##################
## 7-Donor PBMC ##
##################
## Step 1: Remove uninformative genes and cell BCs ----------------------------------------------------------------------------------------------------------------------------------------------------
rna.raw <- Read10X('.')
cell.counts <- Matrix::colSums(rna.raw)
rna.raw <- rna.raw[, which(cell.counts >= 250)]
gene.counts <- Matrix::rowSums(rna.raw)
rna.raw <- rna.raw[which(gene.counts > 3), ]


## Step 2: Pre-process Seurat object, remove high %Mito cells -----------------------------------------------------------------------------------------------------------------------------------------
seu_gh <- CreateSeuratObject(rna.raw)
seu_gh <- NormalizeData(seu_gh)
seu_gh <- ScaleData(seu_gh)
seu_gh <- FindVariableFeatures(seu_gh, selection.method = "vst", nfeatures = 2000)
seu_gh <- RunPCA(seu_gh)
seu_gh <- RunUMAP(seu_gh, dims=1:18)
seu_gh <- FindNeighbors(seu_gh, dims=1:18)
seu_gh <- FindClusters(seu_gh, resolution = 1.0)
mito.genes <- grep("^MT-", rownames(seu_gh@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(seu_gh@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(seu_gh@assays$RNA@counts)
seu_gh@meta.data[,"PercentMito"] <- percent.mito

seu_gh <- SubsetData(seu_gh, cells=rownames(seu_gh@meta.data)[which(seu_gh@active.ident %ni% c(4,8))])
seu_gh <- NormalizeData(seu_gh)
seu_gh <- ScaleData(seu_gh)
seu_gh <- FindVariableFeatures(seu_gh, selection.method = "vst", nfeatures = 2000)
seu_gh <- RunPCA(seu_gh)
seu_gh <- RunUMAP(seu_gh, dims=1:16)
seu_gh <- FindNeighbors(seu_gh, dims=1:16)
seu_gh <- FindClusters(seu_gh, resolution = 1.0)

seu_gh <- SubsetData(seu_gh, cells=rownames(seu_gh@meta.data)[which(seu_gh@active.ident != 11)])
seu_gh <- NormalizeData(seu_gh)
seu_gh <- ScaleData(seu_gh)
seu_gh <- FindVariableFeatures(seu_gh, selection.method = "vst", nfeatures = 2000)
seu_gh <- RunPCA(seu_gh)
seu_gh <- RunUMAP(seu_gh, dims=1:13)
seu_gh <- FindNeighbors(seu_gh, dims=1:13)
seu_gh <- FindClusters(seu_gh, resolution = 0.5)
    
## Sample classification workflow continued in "SampleClassification.R"

## Step 3: Remove doublets defined by SCMK classifications --------------------------------------------------------------------------------------------------------------------------------------------
seu_gh <- SubsetData(seu_gh, cells=rownames(seu_gh@meta.data)[which(seu_gh@meta.data$class != "Doublet")])
seu_gh <- NormalizeData(seu_gh)
seu_gh <- ScaleData(seu_gh)
seu_gh <- FindVariableFeatures(seu_gh, selection.method = "vst", nfeatures = 2000)
seu_gh <- RunPCA(seu_gh)
seu_gh <- RunUMAP(seu_gh, dims=1:17)
seu_gh <- FindNeighbors(seu_gh, dims=1:17)
seu_gh <- FindClusters(seu_gh, resolution = 0.5)


## Step 4: Annotate PBMC cell types -------------------------------------------------------------------------------------------------------------------------------------------------------------------
seu_gh@meta.data[,'CellType'] <- rep("unknown")
seu_gh@meta.data$CellType[which(seu_gh@active.ident %in% 0:2)] <-"CD4T"
seu_gh@meta.data$CellType[which(seu_gh@active.ident == 3)] <-"CD8T"
seu_gh@meta.data$CellType[which(seu_gh@active.ident == 4)] <-"CMono"
seu_gh@meta.data$CellType[which(seu_gh@active.ident == 5)] <-"B"
seu_gh@meta.data$CellType[which(seu_gh@active.ident %in% c(6,7))] <-"NK"
seu_gh@meta.data$CellType[which(seu_gh@active.ident == 8)] <-"NCMono"
seu_gh@meta.data$CellType[which(seu_gh@active.ident == 9)] <-"DC"
seu_gh@meta.data$CellType[which(seu_gh@active.ident == 11)] <-"Platelet"
seu_gh@meta.data$CellType[which(seu_gh@active.ident == 12)] <-"MKI67+T"

## Step 5: Subset CD4+ T-cells, annotate subtypes -----------------------------------------------------------------------------------------------------------------------------------------------------
seu_gh_cd4t <- SubsetData(seu_gh, cells=rownames(seu_gh@meta.data)[which(seu_gh@active.ident %in% 0:2)])
seu_gh_cd4t <- NormalizeData(seu_gh_cd4t)
seu_gh_cd4t <- ScaleData(seu_gh_cd4t)
seu_gh_cd4t <- FindVariableFeatures(seu_gh_cd4t, selection.method = "vst", nfeatures = 2000)
seu_gh_cd4t <- RunPCA(seu_gh_cd4t)
seu_gh_cd4t <- RunUMAP(seu_gh_cd4t, dims=1:15)
seu_gh_cd4t <- FindNeighbors(seu_gh_cd4t, dims=1:15)
seu_gh_cd4t <- FindClusters(seu_gh_cd4t, resolution = 0.5)

seu_gh_cd4t@meta.data[,'subset'] <- rep("activated")
seu_gh_cd4t@meta.data$subset[which(seu_gh_cd4t@active.ident == 1)] <- "memory"
seu_gh_cd4t@meta.data$subset[which(seu_gh_cd4t@active.ident == 3)] <- "naive"

#########################
## BD PBMC +/- CD/CD28 ##
#########################
## Step 1: Load raw gene expression (long-format) -----------------------------------------------------------------------------------------------------------------------------------------------------
rna.BD.raw <- read.table("Combined_BD-Demo-WTA-AbSeq-SMK_Expression_Data.st", header=T, stringsAsFactors = F)
rna.BD.raw <- rna.BD.raw[,c(1,2,7)]
colnames(rna.BD.raw) <- c("Cell","Gene","Count")
cellIDs <- unique(rna.BD.raw$Cell)
geneIDs <- unique(rna.BD.raw$Gene)


## Step 2: Convert to sparse, wide matrix) ------------------------------------------------------------------------------------------------------------------------------------------------------------
library(Matrix)
raw.mat.BD <- matrix(0L, nrow=length(geneIDs), ncol=length(cellIDs))
raw.mat.BD <- Matrix(raw.mat.BD, sparse = TRUE)
rownames(raw.mat.BD) <- geneIDs; colnames(raw.mat.BD) <- cellIDs
x <- 0
for (i in cellIDs) {
  ind <- which(rna.BD.raw$Cell == i)
  x <- x + 1
  print(paste(x, "out of 5665", sep=" "))
  raw.mat.BD[rna.BD.raw$Gene[ind],x] <- rna.BD.raw$Count[ind]
  rna.BD.raw <- rna.BD.raw[-ind, ]
}


## Step 3: Remove uninformative genes and CITE-seq Ab counts) -----------------------------------------------------------------------------------------------------------------------------------------
gene.counts <- Matrix::rowSums(raw.mat.BD)
raw.mat.BD <- raw.mat.BD[which(gene.counts > 3), ]
ind <- grep("pAb",rownames(raw.mat.BD))
raw.mat.BD <- raw.mat.BD[-ind, ]
seu_bd <- CreateSeuratObject(raw.mat.BD)
seu_bd <- NormalizeData(seu_bd)
seu_bd <- ScaleData(seu_bd)
seu_bd <- FindVariableFeatures(seu_bd, selection.method = "vst", nfeatures = 2000)
seu_bd <- RunPCA(seu_bd)
seu_bd <- RunUMAP(seu_bd, dims=1:17)
seu_bd <- FindNeighbors(seu_bd, dims=1:17)
seu_bd <- FindClusters(seu_bd, resolution = 0.25)
mito.genes <- grep("^MT-", rownames(seu_bd@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(seu_bd@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(seu_bd@assays$RNA@counts)
seu_bd@meta.data[,"PercentMito"] <- percent.mito


## Step 4: Remove low-quality cells) ------------------------------------------------------------------------------------------------------------------------------------------------------------------
cells.to.remove <- rownames(seu_bd@meta.data)[which(seu_bd@active.ident == 10)]
cells.to.remove <- c(cells.to.remove, rownames(seu_bd@meta.data)[which(seu_bd@meta.data$PercentMito >= 0.25)])
ind <- which(rownames(seu_bd@meta.data) %in% cells.to.remove)
cells.to.keep <- rownames(seu_bd@meta.data)[-ind]
seu_bd <- SubsetData(seu_bd, cells=cells.to.keep)
seu_bd <- NormalizeData(seu_bd)
seu_bd <- ScaleData(seu_bd)
seu_bd <- FindVariableFeatures(seu_bd, selection.method = "vst", nfeatures = 2000)
seu_bd <- RunPCA(seu_bd)
seu_bd <- RunUMAP(seu_bd, dims=1:17)
seu_bd <- FindNeighbors(seu_bd, dims=1:17)
seu_bd <- FindClusters(seu_bd, resolution = 0.25)


## Step 5: Read in BD SCMK classifications, remove doublets) ------------------------------------------------------------------------------------------------------------------------------------------
temp <- read.csv('BD-Demo-WTA-AbSeq-SMK_Sample_Tag_Calls.csv')
rownames(temp) <- temp[,1]
seu_bd@meta.data[,"class"] <- temp[rownames(seu_bd@meta.data), "Sample_Tag"]

seu_bd <- SubsetData(seu_bd, cells=rownames(seu_bd@meta.data)[which(seu_bd@meta.data$class != "Multiplet")])
seu_bd <- NormalizeData(seu_bd)
seu_bd <- ScaleData(seu_bd)
seu_bd <- FindVariableFeatures(seu_bd, selection.method = "vst", nfeatures = 2000)
seu_bd <- RunPCA(seu_bd)
seu_bd <- RunUMAP(seu_bd, dims=1:17)
seu_bd <- FindNeighbors(seu_bd, dims=1:17)
seu_bd <- FindClusters(seu_bd, resolution = 1.0)

## Step 6: Annotate cells) ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
seu_bd@meta.data[,'CellType'] <- rep("unknown")
seu_bd@meta.data$CellType[which(seu_bd@active.ident == 0)] <-"CD4T_rest"
seu_bd@meta.data$CellType[which(seu_bd@active.ident %in% c(1,2,9))] <-"CD4T_stim"
seu_bd@meta.data$CellType[which(seu_bd@active.ident == 5)] <-"CD8T"
seu_bd@meta.data$CellType[which(seu_bd@active.ident == 10)] <-"NK_rest"
seu_bd@meta.data$CellType[which(seu_bd@active.ident %in% c(4,12))] <-"NK_stim"
seu_bd@meta.data$CellType[which(seu_bd@active.ident %in% c(3,8,11,13))] <-"RBC"
seu_bd@meta.data$CellType[which(seu_bd@active.ident == 6)] <-"Monocyte"
seu_bd@meta.data$CellType[which(seu_bd@active.ident %in% c(15,17))] <-"Macrophage"
seu_bd@meta.data$CellType[which(seu_bd@active.ident == 7)] <-"B"
seu_bd@meta.data$CellType[which(seu_bd@active.ident == 14)] <-"Granulocyte"
seu_bd@meta.data$CellType[which(seu_bd@active.ident == 16)] <-"Neutrophil"