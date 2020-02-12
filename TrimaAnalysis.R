###################################################
## Part 3: Analysis of Trima-separated PBMCs ######
## Chris McGinnis, Gartner Lab, UCSF, 01/21/2020 ##
###################################################


## Step 1: Add Ficoll vs Trima prepraton metadata ---------------------------------------------------------------------------------------------------
seu_pbmc_clean@meta.data[,"Prep"] <- rep("ficoll")
seu_pbmc_clean@meta.data[which(seu_pbmc_clean@meta.data$Donor %in% c("D","E","F","G","H")),"Prep"] <- rep("trima")


## Step 2: Subset by cell type ----------------------------------------------------------------------------------------------------------------------
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


## Step 3: Perform sex-corrected marker analysis ----------------------------------------------------------------------------------------------------
## NK
seu_NK@meta.data[,"Sex_Prep"] <- paste(seu_NK@meta.data$Sex, seu_NK@meta.data$Prep, sep="_")
markers_NK_male <- FindMarkers(seu_NK, group.by="Sex_Prep", test.use="bimod", logfc.threshold=log(1.5), ident.1="Male_fresh", ident.2="Male_TRIMA")
markers_NK_female <- FindMarkers(seu_NK, group.by="Sex_Prep", test.use="bimod", logfc.threshold=log(1.5), ident.1="Female_fresh", ident.2="Female_TRIMA")
temp1 <- rownames(markers_NK_male)[which(markers_NK_male$avg_logFC <= 0)]
temp2 <- rownames(markers_NK_female)[which(markers_NK_female$avg_logFC <= 0)]
temp3 <- temp1[which(temp1 %in% temp2)]
markers_NK <- markers_NK_male[temp3, ]

## cMono
seu_CD14Mono@meta.data[,"Sex_Prep"] <- paste(seu_CD14Mono@meta.data$Sex, seu_CD14Mono@meta.data$Prep, sep="_")
markers_CD14Mono_male <- FindMarkers(seu_CD14Mono, group.by="Sex_Prep", test.use="bimod", logfc.threshold=log(1.5), ident.1="Male_fresh", ident.2="Male_TRIMA")
markers_CD14Mono_female <- FindMarkers(seu_CD14Mono, group.by="Sex_Prep", test.use="bimod", logfc.threshold=log(1.5), ident.1="Female_fresh", ident.2="Female_TRIMA")
temp1 <- rownames(markers_CD14Mono_male)[which(markers_CD14Mono_male$avg_logFC <= 0)]
temp2 <- rownames(markers_CD14Mono_female)[which(markers_CD14Mono_female$avg_logFC <= 0)]
temp3 <- temp1[which(temp1 %in% temp2)]
markers_CD14Mono <- markers_CD14Mono_male[temp3, ]


####################
## Visualizations ##
####################

## Fig. S3A: Trima vs Ficoll PBMCs (400 x 400)
DimPlot(seu_pbmc_clean, group.by="CellType", cols = c("#e7298a","#1b9e77","#e6ab02","#7570b3","#66a61e","#a6761d")) + NoLegend()
DimPlot(seu_pbmc_clean, group.by="Prep", cols = c("black","goldenrod"))+NoLegend()

## Fig. S3B: Classical Mono and NK sub-clustering (400 x 400)
DimPlot(seu_CD14Mono, group.by="Prep", cols = c("black","goldenrod"))+NoLegend()
DimPlot(seu_NK, group.by="Prep", cols = c("black","goldenrod"))+NoLegend()

## Fig. S3B: Trima-specific marker analysis (300x300)
library(vioplot)
temp_CD14 <- as.data.frame(t(as.data.frame(seu_CD14Mono@assays$SCT@data[c("CEBPB","HIST1H1C","MNDA"), ])))
temp_CD14[,"Prep"] <- seu_CD14Mono@meta.data$Prep
temp_NK <- as.data.frame(t(as.data.frame(seu_NK@assays$SCT@data[c("JUN","IFNG","GZMA","PRF1"), ])))
temp_NK[,"Prep"] <- seu_NK@meta.data$Prep

par(mar=rep(1,4))
vioplot(CEBPB~Prep, 
        data=temp_CD14[which(temp_CD14$Prep == "ficoll"), ], 
        col="black", plotCentre="line", side="left", drawRect=F,
        ylim=c(min(temp_CD14$CEBPB), max(temp_CD14$CEBPB)))
vioplot(CEBPB~Prep, 
        data=temp_CD14[which(temp_CD14$Prep == "trima"), ], 
        col="goldenrod", plotCentre="line", side="right", drawRect=F, add=T,
        ylim=c(min(temp_CD14$CEBPB), max(temp_CD14$CEBPB)))

par(mar=rep(1,4))
vioplot(HIST1H1C~Prep, 
        data=temp_CD14[which(temp_CD14$Prep == "ficoll"), ], 
        col="black", plotCentre="line", side="left", drawRect=F,
        ylim=c(min(temp_CD14$HIST1H1C), max(temp_CD14$HIST1H1C)))
vioplot(HIST1H1C~Prep, 
        data=temp_CD14[which(temp_CD14$Prep == "trima"), ], 
        col="goldenrod", plotCentre="line", side="right", drawRect=F, add=T,
        ylim=c(min(temp_CD14$HIST1H1C), max(temp_CD14$HIST1H1C)))

par(mar=rep(1,4))
vioplot(MNDA~Prep, 
        data=temp_CD14[which(temp_CD14$Prep == "ficoll"), ], 
        col="black", plotCentre="line", side="left", drawRect=F,
        ylim=c(min(temp_CD14$MNDA), max(temp_CD14$MNDA)))
vioplot(MNDA~Prep, 
        data=temp_CD14[which(temp_CD14$Prep == "trima"), ], 
        col="goldenrod", plotCentre="line", side="right", drawRect=F, add=T,
        ylim=c(min(temp_CD14$MNDA), max(temp_CD14$MNDA)))

par(mar=rep(1,4))
vioplot(JUN~Prep, 
        data=temp_NK[which(temp_NK$Prep == "ficoll"), ], 
        col="black", plotCentre="line", side="left", drawRect=F,
        ylim=c(min(temp_NK$JUN), max(temp_NK$JUN)))
vioplot(JUN~Prep, 
        data=temp_NK[which(temp_NK$Prep == "trima"), ], 
        col="goldenrod", plotCentre="line", side="right", drawRect=F, add=T,
        ylim=c(min(temp_NK$JUN), max(temp_NK$JUN)))

par(mar=rep(1,4))
vioplot(IFNG~Prep, 
        data=temp_NK[which(temp_NK$Prep == "ficoll"), ], 
        col="black", plotCentre="line", side="left", drawRect=F,
        ylim=c(min(temp_NK$IFNG), max(temp_NK$IFNG)))
vioplot(IFNG~Prep, 
        data=temp_NK[which(temp_NK$Prep == "trima"), ], 
        col="goldenrod", plotCentre="line", side="right", drawRect=F, add=T,
        ylim=c(min(temp_NK$IFNG), max(temp_NK$IFNG)))

par(mar=rep(1,4))
vioplot(GZMA~Prep, 
        data=temp_NK[which(temp_NK$Prep == "ficoll"), ], 
        col="black", plotCentre="line", side="left", drawRect=F,
        ylim=c(min(temp_NK$GZMA), max(temp_NK$GZMA)))
vioplot(GZMA~Prep, 
        data=temp_NK[which(temp_NK$Prep == "trima"), ], 
        col="goldenrod", plotCentre="line", side="right", drawRect=F, add=T,
        ylim=c(min(temp_NK$GZMA), max(temp_NK$GZMA)))

par(mar=rep(1,4))
vioplot(PRF1~Prep, 
        data=temp_NK[which(temp_NK$Prep == "ficoll"), ], 
        col="black", plotCentre="line", side="left", drawRect=F,
        ylim=c(min(temp_NK$PRF1), max(temp_NK$PRF1)))
vioplot(PRF1~Prep, 
        data=temp_NK[which(temp_NK$Prep == "trima"), ], 
        col="goldenrod", plotCentre="line", side="right", drawRect=F, add=T,
        ylim=c(min(temp_NK$PRF1), max(temp_NK$PRF1)))

