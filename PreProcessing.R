#############################################################
## Part 1: mRNA, MULTI-seq and Cell Hashing Pre-Processing ##
## Chris McGinnis, Gartner Lab, UCSF, 01/21/2020 ############
#############################################################

library(Seurat)
library(deMULTIplex)

## Step 1: Define putative cell IDs ---------------------------------------------------------------------------------------------------------------------------
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

## Step 2: Align MULTI-seq and HashTag FASTQs -----------------------------------------------------------------------------------------------------------------
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


## Step 3: Remove uninformative genes in read-depth normalized data -------------------------------------------------------------------------------------------
data_aggr <- Read10X('/Volumes/MULTIseq_1/SLNR_PBMC/aggrSL/outs/raw_feature_bc_matrix/')
cellIDs <- c(cellIDs_l1,cellIDs_l2,cellIDs_l3,cellIDs_l4)
data_aggr <- data_aggr[,cellIDs]
gene.counts <- Matrix::rowSums(data_aggr)
data_aggr <- data_aggr[which(gene.counts >= 3), ]


## Step 4: Pre-process Seurat object --------------------------------------------------------------------------------------------------------------------------
seu_pbmc_raw <- CreateSeuratObject(data_aggr)
mito.genes <- grep("^MT-", rownames(seu_pbmc_raw@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(seu_pbmc_raw@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(seu_pbmc_raw@assays$RNA@counts)
seu_pbmc_raw@meta.data[,"PercentMito"] <- percent.mito
seu_pbmc_raw <- SCTransform(seu_pbmc_raw)
seu_pbmc_raw <- RunPCA(seu_pbmc_raw)
seu_pbmc_raw <- RunUMAP(seu_pbmc_raw, dims = 1:13)
seu_pbmc_raw <- FindNeighbors(seu_pbmc_raw, dims = 1:13)
seu_pbmc_raw <- FindClusters(seu_pbmc_raw, resolution = 1.0)


## Step 5: Remove low-quality cells (high pMito, low nUMI) ----------------------------------------------------------------------------------------------------
'%ni%' <- Negate('%in%')
seu_pbmc_raw_2 <- SubsetData(seu_pbmc_raw, cells=rownames(seu_pbmc_raw@meta.data)[which(seu_pbmc_raw@active.ident %ni% c(4,10,13,18:20,23,24))]) # Remove 3726 cells
seu_pbmc_raw_2 <- SCTransform(seu_pbmc_raw_2)
seu_pbmc_raw_2 <- RunPCA(seu_pbmc_raw_2)
seu_pbmc_raw_2 <- RunUMAP(seu_pbmc_raw_2, dims = 1:12)
seu_pbmc_raw_2 <- FindNeighbors(seu_pbmc_raw_2, dims = 1:12)
seu_pbmc_raw_2 <- FindClusters(seu_pbmc_raw_2, resolution = 1.0)


## Step 6: Before applying DoubletFinder, subset data by lane -------------------------------------------------------------------------------------------------
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


## Step 7: Identify heterotypic doublets in each lane using DoubletFinder ------------------------------------------------------------------------------------
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


## Step 8: Remove doublets, pre-process cleaned merged Seurat object ---------------------------------------------------------------------------------------
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


## Step 9: Annotate Cell Types --------------------------------------------------------------------------------------------------------------------------------
seu_pbmc_clean@meta.data[,"CellType"] <- rep("unknown", nrow(seu_pbmc_clean@meta.data))
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(8,21))] <- "B"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(2,14,24))] <- "NK"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(4,6,15,16))] <- "CD14Mono"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident == 20)] <- "DC"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident == 25)] <- "pDC"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(13,22))] <- "CD16Mono"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(0,1,3,5,7,9,12,17,18,23))] <- "CD4T"
seu_pbmc_clean@meta.data$CellType[which(seu_pbmc_clean@active.ident %in% c(10,11,19))] <- "CD8T"


####################
## Visualizations ##
####################

## Fig. S1A: Raw gene expression data
DimPlot(seu_pbmc_raw, label=T) + NoLegend()
FeaturePlot(seu_pbmc_raw, "PercentMito") + NoLegend()
FeaturePlot(seu_pbmc_raw, "nCount_RNA") + NoLegend()

## Fig. S1B: DoubletFinder results
DimPlot(seu_pbmc_raw2, label=T) + NoLegend()
DimPlot(seu_pbmc_raw2, group.by="DF", cols=c("red","black")) + NoLegend()

## Fig. S1C: Cell Type annotations
DimPlot(seu_pbmc_clean, group.by="CellType", cols = c("#e7298a","#1b9e77","#e6ab02","#7570b3","#66a61e","#a6761d")) + NoLegend()
FeaturePlot(seu_pbmc_clean, "IL7R") + NoLegend()
FeaturePlot(seu_pbmc_clean, "CD8A") + NoLegend()
FeaturePlot(seu_pbmc_clean, "SPON2") + NoLegend()
FeaturePlot(seu_pbmc_clean, "MS4A1") + NoLegend()
FeaturePlot(seu_pbmc_clean, "CD14") + NoLegend()
FeaturePlot(seu_pbmc_clean, "CLEC10A") + NoLegend()
FeaturePlot(seu_pbmc_clean, "FCGR3A") + NoLegend()