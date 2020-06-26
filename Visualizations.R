#####################################################
## Manuscript Visualizations, McGinnis et al, 2020 ##
## Chris McGinnis, Gartner Lab, UCSF, 06/26/2020 ####
#####################################################

library(scales)
library(ggplot2)

## Fig. 2A: Barcode tSNEs colored by (i) deMULTIplex, (ii) demuxEM, and (iii) souporcell ----------------------------------------------------------------------------------------------------
ggplot(barTSNE_multi, aes(x=TSNE1, y=TSNE2, color=MULTI)) + 
  geom_point(size=0.5) +
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","black","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","black")) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

## Fig. 2B: Cell Hashing negative proportion aacross PBMC cell types ------------------------------------------------------------------------------------------------------------------------
temp <- seu_pbmc_ficoll@meta.data[which(seu_pbmc_ficoll@meta.data$LaneID == "3"), c("CellType","HASH")]
temp$HASH[which(temp$HASH %ni% c("Doublet","Negative"))] <- "Singlet"
temp$HASH <- factor(temp$HASH, levels=c("Singlet","Negative","Doublet"))
temp$CellType <- factor(temp$CellType, levels=c("CD4T","CD8T","NK","B","CD14Mono","DC","CD16Mono","pDC"))
ggplot(temp[which(temp$CellType != "pDC" & temp$HASH != "Doublet"), ], aes(x=CellType, fill=HASH)) +
  geom_bar(position = "fill", alpha=0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("black","red")) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

## Fig. 2C: CD4+ T-cell annotations with SCMK negatives -------------------------------------------------------------------------------------------------------------------------------------
seu_CD4T_ficoll@meta.data[,"Fig2C"] <- seu_CD4T_ficoll@meta.data$CD4Tsubset
seu_CD4T_ficoll@meta.data$Fig2C[which(seu_CD4T_ficoll@meta.data$HASH == "Negative")] <- "Negative"
DimPlot(seu_CD4T_ficoll, group.by="Fig2C", order="Negative", cols=c("grey","#DABE93","grey","black","red")) + NoLegend()

## Fig. 3A: Sample classification density UMAPs ---------------------------------------------------------------------------------------------------------------------------------------------
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

## Fig. 3B: Alloreactivity marker analysis, CD4+ T-cells ------------------------------------------------------------------------------------------------------------------------------------
seu_CD4T_ficoll@meta.data[,"TEMP"] <- rep("none")
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Activated" & seu_CD4T_ficoll@meta.data$LaneID %in% c(1,4))] <- "act_unmixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Activated" & seu_CD4T_ficoll@meta.data$LaneID %in% c(2,3))] <- "act_mixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Memory" & seu_CD4T_ficoll@meta.data$LaneID %in% c(1,4))] <- "mem_unmixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Memory" & seu_CD4T_ficoll@meta.data$LaneID %in% c(2,3))] <- "mem_mixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Naive" & seu_CD4T_ficoll@meta.data$LaneID %in% c(1,4))] <- "naive_unmixed"
seu_CD4T_ficoll@meta.data$TEMP[which(seu_CD4T_ficoll@meta.data$Donor_MULTI == "A" & seu_CD4T_ficoll@meta.data$CD4TSubset == "Naive" & seu_CD4T_ficoll@meta.data$LaneID %in% c(2,3))] <- "naive_mixed"

seu_CD4T_ficoll@meta.data$TEMP <- factor(seu_CD4T_ficoll@meta.data$TEMP, levels=c("act_unmixed","act_mixed","mem_unmixed","mem_mixed","naive_unmixed","naive_mixed","none"))
seu_CD4T_ficoll <- SetIdent(seu_CD4T_ficoll, value = seu_CD4T_ficoll@meta.data$TEMP)

VlnPlot(seu_CD4T_ficoll, "IFNG", cols=rep(c("black","red"),3)) + NoLegend()
VlnPlot(seu_CD4T_ficoll, "CD40LG", cols=rep(c("black","red"),3)) + NoLegend()
VlnPlot(seu_CD4T_ficoll, "FOS", cols=rep(c("black","red"),3)) + NoLegend()
VlnPlot(seu_CD4T_ficoll, "DUSP1", cols=rep(c("black","red"),3)) + NoLegend()

## Fig. 3D: JSD analysis bar plots  ---------------------------------------------------------------------------------------------------------------------------------------------------------
ggplot(jsd.summary, aes(x=CellType, y=JSD, fill=Group)) + geom_col(position="dodge", color="black") + ylim(c(0,1)) +
  theme_classic() + scale_fill_manual(values=c("white","#DABE93","grey","black","red")) + theme(legend.position="none") +
  geom_errorbar(aes(ymin=JSD-SD, ymax=JSD+SD), width=.2, position=position_dodge(.9)) 

## Fig. S1A,D: 7-Donor PBMC cell type annotations, classifications, and classification frequency bar-plots ----------------------------------------------------------------------------------
DimPlot(seu_gh, group.by="CellType", cols = c("#e7298a","#7570b3","#66a61e","#1b9e77","#D96002","darkred","#e6ab02","#a6761d","black","grey")) + NoLegend()
DimPlot(seu_gh, cells.highlight = rownames(seu_gh@meta.data)[which(seu_gh@meta.data$class == "Negative")], sizes.highlight = 0.1)+NoLegend()
temp <- seu_gh@meta.data[which(seu_gh@meta.data$CellType %in% c("CD4T","CD8T","NK","B","CMono","DC","NCMono")), c("CellType","class")]
temp$class[which(temp$class != "Negative")] <- "Singlet"
temp$class <- factor(temp$class, levels=c("Singlet","Negative"))
temp$CellType <- factor(temp$CellType, levels=c("NK","CD4T","CD8T","B","DC","CMono","NCMono"))
ggplot(temp, aes(x=CellType, fill=class)) +
  geom_bar(position = "fill", alpha=0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("black","red")) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
# NOTE: Included cell types with at least 100 cells

## Fig. S1B,E: 7-Donor CD4+ T-cell subtype annotations, classifications, and classification frequency bar-plots -----------------------------------------------------------------------------
DimPlot(seu_gh_cd4t, group.by="subtype", cols=c("grey","black","#DABE93")) + NoLegend()
DimPlot(seu_gh_cd4t, cells.highlight = rownames(seu_gh_cd4t@meta.data)[which(seu_gh_cd4t@meta.data$class == "Negative")], sizes.highlight = 0.1)+NoLegend()
temp <- seu_gh_cd4t@meta.data[, c("subset","class")]
temp$class[which(temp$class != "Negative")] <- "Singlet"
temp$class <- factor(temp$class, levels=c("Singlet","Negative"))
ggplot(temp, aes(x=subset, fill=class)) +
  geom_bar(position = "fill", alpha=0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("black","red")) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

## Fig. S1C,F: BD PBMC cell type annotations, classifications, and classification frequency bar-plots ---------------------------------------------------------------------------------------
DimPlot(seu_bd, group.by="CellType", cols = c("#e7298a","#948ee6","#4f4c78","#66a61e","black","#197559","#30e3ac","grey","#e8a527","#6b4c12","darkred")) + NoLegend()
DimPlot(seu_bd, group.by="class", order="Undetermined", cols=c("grey","black","red")) + NoLegend() 

ind <- which(seu_bd@meta.data$CellType != "Granulocyte")
temp <- cbind(as.character(seu_bd@meta.data$CellType[ind]), as.character(seu_bd@meta.data$class[ind]))
colnames(temp) <- c("CellType","class")
temp[which(temp[,2] != "Undetermined"), 2] <- "Classified"
temp <- as.data.frame(temp)
temp$CellType <- factor(temp$CellType, levels=c("CD4T_rest","CD8T","B","NK_rest","NK_stim","CD4T_stim","RBC","Macrophage","Monocyte"))
ggplot(temp, aes(x=CellType, fill=class)) +
  geom_bar(position = "fill", alpha=0.8) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("black","red")) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
# NOTE: Included cell types with at least 100 cells

## Fig. S2A: Trima vs Ficoll PBMC gene expression space -------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_pbmc_clean, group.by="CellType", cols = c("#e7298a","#1b9e77","#e6ab02","#7570b3","#66a61e","#a6761d")) + NoLegend()
DimPlot(seu_pbmc_clean, group.by="Prep", cols = c("black","goldenrod"))+NoLegend()

## Fig. S2B: Classical Mono and NK sub-clustering -------------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_CD14Mono, group.by="Prep", cols = c("black","goldenrod"))+NoLegend()
DimPlot(seu_NK, group.by="Prep", cols = c("black","goldenrod"))+NoLegend()

## Fig. S2B: Trima-specific marker analysis -------------------------------------------------------------------------------------------------------------------------------------------------
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

## Fig. S3A: PBMC cell type frequency plot, 8-donor PBMC ------------------------------------------------------------------------------------------------------------------------------------
ggplot(anno_freq, aes(x=DonorLane, y=Frequency, fill=CellType)) + geom_col(col="black") + theme_bw() + theme(legend.position="none") +
  scale_fill_manual(values = c("#e7298a","#1b9e77","goldenrod","#7570b3","#66a61e","#D96002","#a6761d"))

## Fig. S3B: PBMC cell type pairwise proportion test heatmap, 8-donor PBMC ------------------------------------------------------------------------------------------------------------------
plot.cd4t <- ggplot(proptest_data[which(proptest_data$CellType == "CD4T"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.cd8t <- ggplot(proptest_data[which(proptest_data$CellType == "CD8T"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.cd14mono <- ggplot(proptest_data[which(proptest_data$CellType == "CD14Mono"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.cd16mono <- ggplot(proptest_data[which(proptest_data$CellType == "CD16Mono"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.b <- ggplot(proptest_data[which(proptest_data$CellType == "B"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.nk <- ggplot(proptest_data[which(proptest_data$CellType == "NK"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.dc <- ggplot(proptest_data[which(proptest_data$CellType == "DC"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")

grid.arrange(plot.cd4t, plot.cd8t, plot.cd14mono, plot.cd16mono, plot.b, plot.nk, plot.dc, nrow=2, ncol=4)

## Fig. S3C: CD4+ T-cell subtype frequency plot, 8-donor PBMC -------------------------------------------------------------------------------------------------------------------------------
ggplot(anno_freq_CD4T, aes(x=DonorLane, y=Frequency, fill=Subtype)) + geom_col(col="black") + theme_bw() + theme(legend.position="none") +
  scale_fill_manual(values = c("grey","black","#DABE93"))

## Fig. S3D: CD4+ T-cell subtype pairwise proportion test heatmap, 8-donor PBMC -------------------------------------------------------------------------------------------------------------
plot.act <- ggplot(proptest_CD4T_data[which(proptest_CD4T_data$Subset == "Activated"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.mem <- ggplot(proptest_CD4T_data[which(proptest_CD4T_data$Subset == "Memory"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.naive <- ggplot(proptest_CD4T_data[which(proptest_CD4T_data$Subset == "Naive"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")

grid.arrange(plot.act, plot.mem, plot.naive, nrow=1, ncol=3)

## Fig. S4A: PBMC cell type frequency plot, Zheng et al PBMC --------------------------------------------------------------------------------------------------------------------------------
ggplot(anno_freq_zheng, aes(x=DonorLane, y=Frequency, fill=CellType)) + geom_col(col="black") + theme_bw() + theme(legend.position="none") +
  scale_fill_manual(values = c("#e7298a","#1b9e77","goldenrod","#7570b3","#66a61e","#D96002","#a6761d","black"))

## Fig. S4B: PBMC cell type pairwise proportion test heatmap, Zheng et al PBMC ------------------------------------------------------------------------------------------------------------------
plot.cd4t <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "CD4T"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.cd8t <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "CD8T"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.cd14mono <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "CD14Mono"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.cd16mono <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "CD16Mono"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.b <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "B"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.nk <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "NK"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.dc <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "DC"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.plate <- ggplot(proptest_zheng_data[which(proptest_zheng_data$CellType == "Platelet"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")

grid.arrange(plot.cd4t, plot.cd8t, plot.cd14mono, plot.cd16mono, plot.b, plot.nk, plot.dc, plot.plate, nrow=2, ncol=4)

## Fig. S4C: CD4+ T-cell subtype frequency plot, Zheng et al PBMC ---------------------------------------------------------------------------------------------------------------------------
ggplot(anno_freq_zheng_CD4T, aes(x=Donor_Mix, y=Frequency, fill=Subtype)) + geom_col(col="black") + theme_bw() + theme(legend.position="none") +
  scale_fill_manual(values = c("grey","black","#DABE93"))

## Fig. S4D: CD4+ T-cell subtype pairwise proportion test heatmap, Zheng et al PBMC ---------------------------------------------------------------------------------------------------------
plot.act <- ggplot(proptest_zheng_CD4T_data[which(proptest_zheng_CD4T_data$Subtype == "activated"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.mem <- ggplot(proptest_zheng_CD4T_data[which(proptest_zheng_CD4T_data$Subtype == "memory"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")
plot.naive <- ggplot(proptest_zheng_CD4T_data[which(proptest_zheng_CD4T_data$Subtype == "naive"), ], aes(x=group1, y=group2, fill=Bin)) + 
  geom_tile(color="black") + scale_fill_viridis(limits=c(0,5)) + theme_void() + theme(legend.position="none")

grid.arrange(plot.act, plot.mem, plot.naive, nrow=1, ncol=3)

## Fig. S4F: JSD analysis bar plots, Zheng et al PBMC ---------------------------------------------------------------------------------------------------------------------------------------
ggplot(jsd.summary_zheng, aes(x=CellType, y=JSD, fill=Group)) + geom_col(position="dodge", color="black") + ylim(c(0,1)) +
theme_classic() + scale_fill_manual(values=c("white","#DABE93","red")) + theme(legend.position="none") +
  geom_errorbar(aes(ymin=JSD-SD, ymax=JSD+SD), width=.2, position=position_dodge(.9)) 

## Fig. S5A: Raw gene expression data, 8-donor PBMC -----------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_pbmc_raw, label=T) + NoLegend()
FeaturePlot(seu_pbmc_raw, "PercentMito") + NoLegend()
FeaturePlot(seu_pbmc_raw, "nCount_RNA") + NoLegend()

## Fig. S5B: DoubletFinder results, 8-donor PBMC --------------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_pbmc_raw2, label=T) + NoLegend()
DimPlot(seu_pbmc_raw2, group.by="DF", cols=c("red","black")) + NoLegend()

## Fig. S5C: Cell Type annotations, 8-donor PBMC --------------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_pbmc_clean, group.by="CellType", cols = c("#e7298a","#1b9e77","#e6ab02","#7570b3","#66a61e","#a6761d")) + NoLegend()
FeaturePlot(seu_pbmc_clean, "IL7R") + NoLegend()
FeaturePlot(seu_pbmc_clean, "CD8A") + NoLegend()
FeaturePlot(seu_pbmc_clean, "SPON2") + NoLegend()
FeaturePlot(seu_pbmc_clean, "MS4A1") + NoLegend()
FeaturePlot(seu_pbmc_clean, "CD14") + NoLegend()
FeaturePlot(seu_pbmc_clean, "CLEC10A") + NoLegend()
FeaturePlot(seu_pbmc_clean, "FCGR3A") + NoLegend()

## Fig. S6A: Raw gene expression data, Zheng et al PBMC -------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_zheng, label=T) + NoLegend()
FeaturePlot(seu_zheng, "PercentMito") + NoLegend()
FeaturePlot(seu_zheng, "nCount_RNA") + NoLegend()

## Fig. S6B: DoubletFinder results, Zheng et al PBMC ----------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_zheng_2, label=T) + NoLegend()
DimPlot(seu_zheng_2, group.by="DF", cols=c("red","black")) + NoLegend()

## Fig. S6C: Cell type annotations, Zheng et al PBMC ----------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_zheng_clean, group.by="CellType", cols = c("#e7298a","#1b9e77","goldenrod","#7570b3","#66a61e","#D96002","grey","#a6761d","black")) + NoLegend() + NoAxes()
DimPlot(seu_zheng_clean, group.by="Donor", order="Doublet", cols = c("dodgerblue","goldenrod","black"))+NoLegend()+NoAxes()
FeaturePlot(seu_zheng_clean, "IL7R") + NoLegend()
FeaturePlot(seu_zheng_clean, "CD8A") + NoLegend()
FeaturePlot(seu_zheng_clean, "SPON2") + NoLegend()
FeaturePlot(seu_zheng_clean, "MS4A1") + NoLegend()
FeaturePlot(seu_zheng_clean, "CD14") + NoLegend()
FeaturePlot(seu_zheng_clean, "CLEC10A") + NoLegend()
FeaturePlot(seu_zheng_clean, "FCGR3A") + NoLegend()
FeaturePlot(seu_zheng_clean, "PF4") + NoLegend()
FeaturePlot(seu_zheng_clean, "MKI67") + NoLegend()
FeaturePlot(seu_zheng_clean, "MZB1") + NoLegend()
FeaturePlot(seu_zheng_clean, "LILRA4") + NoLegend()
FeaturePlot(seu_zheng_clean, "GATA2") + NoLegend()

## Fig. S7: CD4+ T-cell annotations, 8-donor PBMC, 7-donor PBMC, and Zheng et al PBMC -------------------------------------------------------------------------------------------------------
DimPlot(seu_CD4T_ficoll, group.by="CD4Tsubset", cols=c("grey","black","#DABE93")) + NoLegend()
DimPlot(seu_zheng_CD4T, group.by="Subtype", cols=c("grey","black","#DABE93")) + NoLegend()
DimPlot(seu_gh_cd4t, group.by="subtype", cols=c("grey","black","#DABE93")) + NoLegend()

## Fig. S8A: Cell type annotations, 7-donor PBMC --------------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_gh, group.by="CellType", cols = c("#e7298a","#7570b3","#66a61e","#1b9e77","#D96002","darkred","#e6ab02","#a6761d","black","grey")) + NoLegend()
FeaturePlot(seu_gh, "IL7R") + NoLegend()
FeaturePlot(seu_gh, "CD8A") + NoLegend()
FeaturePlot(seu_gh, "SPON2") + NoLegend()
FeaturePlot(seu_gh, "MKI67") + NoLegend()
FeaturePlot(seu_gh, "MS4A1") + NoLegend()
FeaturePlot(seu_gh, "CD14") + NoLegend()
FeaturePlot(seu_gh, "CLEC10A") + NoLegend()
FeaturePlot(seu_gh, "FCGR3A") + NoLegend()
FeaturePlot(seu_gh, "PF4") + NoLegend()

## Fig. S8B: Cell type annotations, BD PBMC +/- CD3/CD28 ------------------------------------------------------------------------------------------------------------------------------------
DimPlot(seu_bd, group.by="CellType", cols = c("#e7298a","#948ee6","#4f4c78","#66a61e","black","#197559","#30e3ac","grey","#e8a527","#6b4c12","darkred")) + NoLegend()
DimPlot(seu_bd, group.by="class", order="Undetermined", cols=c("grey","black","red")) + NoLegend() 
FeaturePlot(seu_bd, "IL7R") + NoLegend()
FeaturePlot(seu_bd, "CD8A") + NoLegend()
FeaturePlot(seu_bd, "SPON2") + NoLegend()
FeaturePlot(seu_bd, "MS4A1") + NoLegend()
FeaturePlot(seu_bd, "CD14") + NoLegend()
FeaturePlot(seu_bd, "GBP1") + NoLegend()
FeaturePlot(seu_bd, "ENO1") + NoLegend()
FeaturePlot(seu_bd, "MKI67") + NoLegend()
FeaturePlot(seu_bd, "GNLY") + NoLegend()
FeaturePlot(seu_bd, "HBB") + NoLegend()
FeaturePlot(seu_bd, "LTF") + NoLegend()
FeaturePlot(seu_bd, "GATA2") + NoLegend()