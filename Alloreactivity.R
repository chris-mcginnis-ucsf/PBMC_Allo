##########################################################
## Part 4: Alloreactivity analysis amongst Ficoll PBMCs ##
## Chris McGinnis, Gartner Lab, UCSF 01/21/2020 ##########
##########################################################

library(Seurat)
library(MASS)
library(philentropy)
library(colorspace)
library(dendextend)

## Step 1: Subset PBMCs, classical monocytes, and CD4+ T-cells by ficoll-separated (provides outgroup for JSD) --------------------------------------
seu_pbmc_ficoll <- SubsetData(seu_pbmc_clean, cells = rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$Donor %in% c("A","B","C"))])
seu_pbmc_ficoll <- SCTransform(seu_pbmc_ficoll)
seu_pbmc_ficoll <- RunPCA(seu_pbmc_ficoll)
seu_pbmc_ficoll <- RunUMAP(seu_pbmc_ficoll, dims = 1:12)
seu_pbmc_ficoll <- FindNeighbors(seu_pbmc_ficoll, dims = 1:12)
seu_pbmc_ficoll <- FindClusters(seu_pbmc_ficoll, resolution = 0.8)

seu_CD14Mono_ficoll <- SubsetData(seu_CD14Mono, cells = rownames(seu_CD14Mono@meta.data)[which(seu_CD14Mono@meta.data$Donor %in% c("A","B","C"))])
seu_CD14Mono_ficoll <- SCTransform(seu_CD14Mono_ficoll)
seu_CD14Mono_ficoll <- RunPCA(seu_CD14Mono_ficoll)
seu_CD14Mono_ficoll <- RunUMAP(seu_CD14Mono_ficoll, dims = 1:17)
seu_CD14Mono_ficoll <- FindNeighbors(seu_CD14Mono_ficoll, dims = 1:17)
seu_CD14Mono_ficoll <- FindClusters(seu_CD14Mono_ficoll, resolution = 0.8)

seu_CD4T_ficoll <- SubsetData(seu_pbmc_clean, cells = rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$CellType == "CD4T" & seu_pbmc_clean@meta.data$Donor %in% c("A","B","C"))])
seu_CD4T_ficoll <- SCTransform(seu_CD4T_ficoll)
seu_CD4T_ficoll <- RunPCA(seu_CD4T_ficoll)
seu_CD4T_ficoll <- RunUMAP(seu_CD4T_ficoll, dims = 1:25)
seu_CD4T_ficoll <- FindNeighbors(seu_CD4T_ficoll, dims = 1:25)
seu_CD4T_ficoll <- FindClusters(seu_CD4T_ficoll, resolution = 0.8)


## Step 2: Subset CD4+ T-cells and classical monocytes to include equal numbers from A-mixed, A-unmixed, B, and C -----------------------------------
## CD4+ T-cells
temp <- paste(seu_CD4T_ficoll@meta.data$Donor_MULTI, seu_CD4T_ficoll@meta.data$LaneID, seu_CD4T_ficoll@meta.data$CD4TSubset, sep="_")
names(temp) <- rownames(seu_CD4T_ficoll@meta.data)
ncell_act <- 111
ncell_mem <- 48
ncell_naive <- 8
cell.vec <- NULL
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_1_Activated",temp)], ncell_act))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_2_Activated",temp)], ncell_act))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_3_Activated",temp)], ncell_act))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_4_Activated",temp)], ncell_act))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_2_Activated",temp)], ncell_act))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_3_Activated",temp)], ncell_act))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_2_Activated",temp)], ncell_act))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_3_Activated",temp)], ncell_act))

cell.vec <- c(cell.vec, sample(names(temp)[grep("A_1_Memory",temp)], ncell_mem))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_2_Memory",temp)], ncell_mem))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_3_Memory",temp)], ncell_mem))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_4_Memory",temp)], ncell_mem))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_2_Memory",temp)], ncell_mem))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_3_Memory",temp)], ncell_mem))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_2_Memory",temp)], ncell_mem))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_3_Memory",temp)], ncell_mem))

cell.vec <- c(cell.vec, sample(names(temp)[grep("A_1_Naive",temp)], ncell_naive))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_2_Naive",temp)], ncell_naive))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_3_Naive",temp)], ncell_naive))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_4_Naive",temp)], ncell_naive))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_2_Naive",temp)], ncell_naive))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_3_Naive",temp)], ncell_naive))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_2_Naive",temp)], ncell_naive))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_3_Naive",temp)], ncell_naive))

cell.vec_CD4T <- cell.vec

## Classical monocytes
temp <- paste(seu_CD14Mono_ficoll@meta.data$Donor_MULTI, seu_CD14Mono_ficoll@meta.data$LaneID, sep="_")
names(temp) <- rownames(seu_CD14Mono_ficoll@meta.data)
ncell_mono <- 48
cell.vec <- NULL
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_1",temp)], ncell_mono))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_2",temp)], ncell_mono))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_3",temp)], ncell_mono))
cell.vec <- c(cell.vec, sample(names(temp)[grep("A_4",temp)], ncell_mono))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_2",temp)], ncell_mono))
cell.vec <- c(cell.vec, sample(names(temp)[grep("B_3",temp)], ncell_mono))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_2",temp)], ncell_mono))
cell.vec <- c(cell.vec, sample(names(temp)[grep("C_3",temp)], ncell_mono))

cell.vec_mono <- cell.vec


## Step 3: Subset and pre-process Seurat objects ---------------------------------------------------------------------------------------------------------------
seu_CD14Mono_jsd <- SubsetData(seu_CD14Mono_ficoll, cells=cell.vec_mono)
seu_CD14Mono_jsd <- SCTransform(seu_CD14Mono_jsd)
seu_CD14Mono_jsd <- RunPCA(seu_CD14Mono_jsd)
seu_CD14Mono_jsd <- RunUMAP(seu_CD14Mono_jsd, dims = 1:10)
seu_CD14Mono_jsd <- FindNeighbors(seu_CD14Mono_jsd, dims = 1:10)
seu_CD14Mono_jsd <- FindClusters(seu_CD14Mono_jsd, resolution = 0.8)

seu_CD4T_jsd <- SubsetData(seu_CD4T_ficoll, cells=cell.vec_CD4T)
seu_CD4T_jsd <- SCTransform(seu_CD4T_jsd)
seu_CD4T_jsd <- RunPCA(seu_CD4T_jsd)
seu_CD4T_jsd <- RunUMAP(seu_CD4T_jsd, dims = 1:14)
seu_CD4T_jsd <- FindNeighbors(seu_CD4T_jsd, dims = 1:14)
seu_CD4T_jsd <- FindClusters(seu_CD4T_jsd, resolution = 0.8)


## Step 4: Extract UMAP and metadata for each cell subset ---------------------------------------------------------------------------------------------------------------
umap_CD4T_jsd <- as.data.frame(seu_CD4T_jsd@reductions$umap@cell.embeddings[,1:2])
umap_CD4T_jsd[,"Donor_Lane"] <- seu_CD4T_jsd@meta.data$Donor_LaneID

umap_CD14Mono_jsd <- as.data.frame(seu_CD14Mono_jsd@reductions$umap@cell.embeddings[,1:2])
umap_CD14Mono_jsd[,"Donor_Lane"] <- seu_CD14Mono_jsd@meta.data$Donor_LaneID


## Step 5: Compute GKDE and JSD in UMAP space for each cell subset ---------------------------------------------------------------------------------------------------------------
groups <- c("A_1","A_2","A_3","A_4","B_2","B_3","C_2","C_3")

p <- ggplot_build(ggplot(umap_CD4T_jsd, aes(x=UMAP_1, y=UMAP_2)))
xrange_CD4T <- p$layout$panel_scales_x[[1]]$range$range
yrange_CD4T <- p$layout$panel_scales_y[[1]]$range$range

p <- ggplot_build(ggplot(umap_CD14Mono_jsd, aes(x=UMAP_1, y=UMAP_2)))
xrange_CD14Mono <- p$layout$panel_scales_x[[1]]$range$range
yrange_CD14Mono <- p$layout$panel_scales_y[[1]]$range$range

kde2d_CD4T <- list()
kde2d_CD14Mono <- list()
JSD.div_CD4T <- data.frame()
JSD.div_CD14Mono <- data.frame()

for (i in 1:8){
  print(i)
  for (j in 1:8){
    kde2d_CD4T[[i]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$Donor_Lane == groups[i])], 1], 
                             umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$Donor_Lane == groups[i])], 2],
                             n=500, lims=c(xrange_CD4T,yrange_CD4T))
    kde2d_CD4T[[j]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$Donor_Lane == groups[j])], 1], 
                             umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$Donor_Lane == groups[j])], 2],
                             n=500, lims=c(xrange_CD4T,yrange_CD4T))
    JSD.div_CD4T[i,j] <- JSD(rbind(as.vector(kde2d_CD4T[[i]]$z), as.vector(kde2d_CD4T[[j]]$z)))
    
    kde2d_CD14Mono[[i]] <- kde2d(umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$Donor_Lane == groups[i])], 1], 
                                 umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$Donor_Lane == groups[i])], 2],
                                 n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
    kde2d_CD14Mono[[j]] <- kde2d(umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$Donor_Lane == groups[j])], 1], 
                                 umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$Donor_Lane == groups[j])], 2],
                                 n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
    JSD.div_CD14Mono[i,j] <- JSD(rbind(as.vector(kde2d_CD14Mono[[i]]$z), as.vector(kde2d_CD14Mono[[j]]$z)))
  }
}

colnames(JSD.div_CD4T) <- groups
rownames(JSD.div_CD4T) <- groups
colnames(JSD.div_CD14Mono) <- groups
rownames(JSD.div_CD14Mono) <- groups


## Step 6: Perform hierarchical clustering ----------------------------------------------------------------------------------------------------------
hclust_CD4T <- as.dendrogram(hclust(as.dist(JSD.div_CD4T), method="ward.D2"))
hclust_CD14Mono <- as.dendrogram(hclust(as.dist(JSD.div_CD14Mono), method="ward.D2"))

## Step 7: Perform donor A permutation test for CD4+ T-cells ----------------------------------------------------------------------------------------
JSD.div_CD4T_permute <- list()
ind_a <- grep("A", umap_CD4T_jsd$Donor_Lane)
umap_CD4T_jsd[,"permute"] <- umap_CD4T_jsd$Donor_Lane)
for (iter in 1:100) {
  print(iter)
  # permute donor A classifications
  ind <- sample(ind_a, length(ind_a)) 
  umap_CD4T_jsd$permute[ind_a] <- umap_CD4T_jsd$permute[ind]
  
  # compute JSD
  kde2d_temp <- list()
  JSD.div_temp <- data.frame()
  for (i in 1:8){
    for (j in 1:8){
      kde2d_temp[[i]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[i])], 1], 
                               umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[i])], 2],
                               n=500, lims=c(xrange_CD4T,yrange_CD4T))
      kde2d_temp[[j]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[j])], 1], 
                               umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[j])], 2],
                               n=500, lims=c(xrange_CD4T,yrange_CD4T))
      JSD.div_temp[i,j] <- JSD(rbind(as.vector(kde2d_temp[[i]]$z), as.vector(kde2d_temp[[j]]$z)))
    }
  }
  
  # store results
  colnames(JSD.div_temp) <- groups
  rownames(JSD.div_temp) <- groups
  JSD.div_temp <- apply(JSD.div_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_temp))))
  JSD.div_CD4T_permute[[iter]] <- JSD.div_temp
  
}

####################
## Visualizations ##
####################

## Fig. 3A: Sample classification density UMAPs
umap_pbmc_ficoll <- as.data.frame(seu_pbmc_ficoll@reductions$umap@cell.embeddings)
umap_pbmc_ficoll[,"Group"] <- paste(seu_pbmc_ficoll@meta.data$Donor, seu_pbmc_ficoll@meta.data$LaneID, sep="_")

DimPlot(seu_pbmc_ficoll, group.by="CellType", cols = c("#e7298a","#1b9e77","#e6ab02","#7570b3","#66a61e","#a6761d")) + NoLegend()

ggplot(umap_pbmc_ficoll[which(umap_pbmc_ficoll$Group %in% c("A_1","A_4")), ], aes(x=UMAP_1, y=UMAP_2)) +
  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=FALSE) +
  scale_fill_gradient2(low="black",mid="#9B2977",high="#FEFCC8",midpoint=0.5) +
  theme(legend.position="none", strip.text=element_blank(), strip.background = element_blank(),
        axis.text=element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

ggplot(umap_pbmc_ficoll[which(umap_pbmc_ficoll$Group %in% c("A_2","A_3")), ], aes(x=UMAP_1, y=UMAP_2)) +
  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=FALSE) +
  scale_fill_gradient2(low="black",mid="#9B2977",high="#FEFCC8",midpoint=0.5) +
  theme(legend.position="none", strip.text=element_blank(), strip.background = element_blank(),
        axis.text=element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

ggplot(umap_pbmc_ficoll[which(umap_pbmc_ficoll$Group %in% c("B_2","B_3","C_2","C_3")), ], aes(x=UMAP_1, y=UMAP_2)) +
  stat_density_2d(geom="raster", aes(fill=stat(ndensity)), contour=FALSE) +
  scale_fill_gradient2(low="black",mid="#9B2977",high="#FEFCC8",midpoint=0.5) +
  theme(legend.position="none", strip.text=element_blank(), strip.background = element_blank(),
        axis.text=element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())

## Fig. 3B: CD4+ T-cell UMAP w/ donor classificaitons
seu_CD4T_jsd@meta.data[,"Donor_LaneID"] <- paste(seu_CD4T_jsd@meta.data$Donor_MULTI, seu_CD4T_jsd@meta.data$LaneID, sep="_")
DimPlot(seu_CD4T_jsd, group.by="Donor_LaneID", cols=c("black","red","red","black","grey","grey","beige","beige"))+NoLegend()

## Fig. 3C: JSD heatmap for CD4+ T-cells (400x400)
jsd.heatmap_CD4T <- as.data.frame(matrix(0L, nrow=64, ncol=3))
colnames(jsd.heatmap_CD4T) <- c("Group_1","Group_2","JSD")
jsd.heatmap_CD4T$Group_1 <- rep(groups, each=8)
jsd.heatmap_CD4T$Group_2 <- rep(groups, 8)
temp <- 0
for (i in 1:8) { temp <- c(temp, as.numeric(JSD.div_CD4T[i, ])) }
jsd.heatmap_CD4T$JSD <- temp[-1]

jsd.heatmap_CD4T$Group_1 <- factor(jsd.heatmap_CD4T$Group_1, levels=rev(c("C_3","C_2","B_3","B_2","A_3","A_1","A_4","A_2")))
jsd.heatmap_CD4T$Group_2 <- factor(jsd.heatmap_CD4T$Group_2, levels=c("C_3","C_2","B_3","B_2","A_3","A_1","A_4","A_2"))
ggplot(jsd.heatmap_CD4T, aes(x=Group_1, y=Group_2, fill=JSD)) + 
  geom_raster() +
  scale_fill_gradient2(low="red", 
                       mid="red", 
                       high="black") +
  theme(legend.position="none", 
        strip.text=element_blank(), 
        strip.background = element_blank(),
        axis.text=element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank())


## Fig. S2B: Classical monocyte UMAP w/ donor classiications
seu_CD14Mono_jsd@meta.data[,"Donor_LaneID"] <- paste(seu_CD14Mono_jsd@meta.data$Donor_MULTI, seu_CD14Mono_jsd@meta.data$LaneID, sep="_")
DimPlot(seu_CD14Mono_jsd, group.by="Donor_LaneID", cols=c("black","red","red","black","grey","grey","beige","beige"))+NoLegend()

## Fig. S2C: JSD heatmap for classical monocytes (400x400)
jsd.heatmap_CD14Mono <- as.data.frame(matrix(0L, nrow=64, ncol=3))
colnames(jsd.heatmap_CD14Mono) <- c("Group_1","Group_2","JSD")
jsd.heatmap_CD14Mono$Group_1 <- rep(groups, each=8)
jsd.heatmap_CD14Mono$Group_2 <- rep(groups, 8)
temp <- 0
for (i in 1:8) { temp <- c(temp, as.numeric(JSD.div_CD14Mono[i, ])) }
jsd.heatmap_CD14Mono$JSD <- temp[-1]

jsd.heatmap_CD14Mono$Group_1 <- factor(jsd.heatmap_CD14Mono$Group_1, levels=rev(c("B_3","B_2","C_3","C_2","A_3","A_2","A_4","A_1")))
jsd.heatmap_CD14Mono$Group_2 <- factor(jsd.heatmap_CD14Mono$Group_2, levels=c("B_3","B_2","C_3","C_2","A_3","A_2","A_4","A_1"))
ggplot(jsd.heatmap_CD14Mono, aes(x=Group_1, y=Group_2, fill=JSD)) + 
  geom_raster() +
  scale_fill_gradient2(low="red", 
                       mid="red", 
                       high="black") +
  theme(legend.position="none", 
        strip.text=element_blank(), 
        strip.background = element_blank(),
        axis.text=element_blank(), 
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank())

## Fig. S2D: Alloreactivity marker analysis in CD4+ T-cells (300x300)
seu_CD4T_ficoll@meta.data[,"TEMP"] <- rep("none")
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Activated" & seu_CD4T_ficoll@meta.data$LaneID %in% c(1,4))] <- "act_unmixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Activated" & seu_CD4T_ficoll@meta.data$LaneID %in% c(2,3))] <- "act_mixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Memory" & seu_CD4T_ficoll@meta.data$LaneID %in% c(1,4))] <- "mem_unmixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Memory" & seu_CD4T_ficoll@meta.data$LaneID %in% c(2,3))] <- "mem_mixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Naive" & seu_CD4T_ficoll@meta.data$LaneID %in% c(1,4))] <- "naive_unmixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Naive" & seu_CD4T_ficoll@meta.data$LaneID %in% c(2,3))] <- "naive_mixed"

seu_CD4T_ficoll@meta.data$TEMP <- factor(seu_CD4T_ficoll@meta.data$TEMP, levels=c("act_unmixed","act_mixed","mem_unmixed","mem_mixed","naive_unmixed","naive_mixed","none"))
seu_CD4T_ficoll <- SetIdent(seu_CD4T_ficoll, value = seu_CD4T_ficoll@meta.data$TEMP)

temp <- as.data.frame(seu_CD4T_jsd@assays$SCT@data[c("IFNG","CD40LG","FOS","DUSP1"), ])
temp <- as.data.frame(t(temp))
temp[,"GROUP"] <- seu_CD4T_ficoll@meta.data$TEMP
ind <- which(temp$GROUP == "none")
temp <- temp[-ind, ]

ggplot(temp[which(temp$GROUP != "none"), ], aes(x=GROUP, y=IFNG, fill=GROUP)) + 
  geom_violin() +
  scale_fill_manual(values = c("black","red","black","red","black","red")) +
  theme_bw() +
  theme(legend.position="none", 
        strip.text=element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank())

ggplot(temp[which(temp$GROUP != "none"), ], aes(x=GROUP, y=CD40LG, fill=GROUP)) + 
  geom_violin() +
  scale_fill_manual(values = c("black","red","black","red","black","red")) +
  theme_bw() +
  theme(legend.position="none", 
        strip.text=element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank())

ggplot(temp[which(temp$GROUP != "none"), ], aes(x=GROUP, y=DUSP1, fill=GROUP)) + 
  geom_violin() +
  scale_fill_manual(values = c("black","red","black","red","black","red")) +
  theme_bw() +
  theme(legend.position="none", 
        strip.text=element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank())

ggplot(temp[which(temp$GROUP != "none"), ], aes(x=GROUP, y=FOS, fill=GROUP)) + 
  geom_violin() +
  scale_fill_manual(values = c("black","red","black","red","black","red")) +
  theme_bw() +
  theme(legend.position="none", 
        strip.text=element_blank(), 
        strip.background = element_blank(),
        axis.title = element_blank(), 
        axis.ticks = element_blank(), 
        axis.line = element_blank())


