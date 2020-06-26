#########################################################
## Analysis of Trima-separated PBMCs in 8-donor data ####
## Chris McGinnis, Gartner Lab, UCSF, 06/26/2020 ########
#########################################################

library(Seurat)

## Step 1: Add Ficoll vs Trima prepraton metadata ----------------------------------------------------------------------------------------------------------------------------
seu_pbmc_clean@meta.data[,"Prep"] <- rep("ficoll")
seu_pbmc_clean@meta.data[which(seu_pbmc_clean@meta.data$Donor %in% c("D","E","F","G","H")),"Prep"] <- rep("trima")


## Step 2: Subset by cell type -----------------------------------------------------------------------------------------------------------------------------------------------
seu_NK <- SubsetData(seu_pbmc_clean, cells = rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$CellType == "NK")])
seu_NK <- SCTransform(seu_NK)
seu_NK <- RunPCA(seu_NK)
seu_NK <- RunUMAP(seu_NK, dims = 1:10)
seu_NK <- FindNeighbors(seu_NK, dims = 1:10)
seu_NK <- FindClusters(seu_NK, resolution = 0.8)

seu_CD14Mono <- SubsetData(seu_pbmc_clean, cells = rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$CellType == "CD14Mono")])
seu_CD14Mono <- SCTransform(seu_CD14Mono)
seu_CD14Mono <- RunPCA(seu_CD14Mono)
seu_CD14Mono <- RunUMAP(seu_CD14Mono, dims = 1:17)
seu_CD14Mono <- FindNeighbors(seu_CD14Mono, dims = 1:17)
seu_CD14Mono <- FindClusters(seu_CD14Mono, resolution = 0.8)

## Remove outlier CD14 Mono cluster
seu_CD14Mono <- SubsetData(seu_CD14Mono, cells = rownames(seu_CD14Mono@meta.data)[which(seu_CD14Mono@active.ident != 10)])
seu_CD14Mono <- SCTransform(seu_CD14Mono)
seu_CD14Mono <- RunPCA(seu_CD14Mono)
seu_CD14Mono <- RunUMAP(seu_CD14Mono, dims = 1:18)
seu_CD14Mono <- FindNeighbors(seu_CD14Mono, dims = 1:18)
seu_CD14Mono <- FindClusters(seu_CD14Mono, resolution = 0.8)