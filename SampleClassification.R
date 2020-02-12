######################################################
## Part 2: MULTI-seq and CellHashing Demultiplexing ##
## Chris McGinnis, Gartner Lab, UCSF, 01/21/2020 #####
######################################################

library(deMULTIplex)

## Step 1: MULTI-seq classification for cells passing QC --------------------------------------------------------------------------------------------
## Lane 2
cells <- unlist(strsplit(rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$LaneID == 2)], split="-2"))
bar.table <- barTable_l2_multi[cells, 1:4]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round1 <- findThresh(call.list=bar.table_sweep.list)
round1.calls <- classifyCells(bar.table, q=findQ(res_round1$res, res_round1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]

bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(res_round2$res, res_round2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

final.calls_l2_multi <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls_l2_multi) <- c(names(round2.calls),neg.cells)

## Lane 3
cells <- unlist(strsplit(rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$LaneID == 3)], split="-3"))
bar.table <- barTable_l3_multi[cells,1:8]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round1 <- findThresh(call.list=bar.table_sweep.list)
round1.calls <- classifyCells(bar.table, q=findQ(res_round1$res, res_round1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]

bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(res_round2$res, res_round2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

final.calls_l3_multi <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls_l3_multi) <- c(names(round2.calls),neg.cells)

## Record results
seu_pbmc_clean@meta.data[,"MULTI"] <- rep("none")
seu_pbmc_clean@meta.data[names(final.calls_multi_l2),"MULTI"] <- final.calls_multi_l2
seu_pbmc_clean@meta.data[names(final.calls_multi_l3),"MULTI"] <- final.calls_multi_l3

## Convert to donor IDs
seu_pbmc_clean@meta.data[,"Donor"] <- rep("A")
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar1")] <- "A"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$LaneID %in% c(1,4))] <- "A"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar2")] <- "B"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar3")] <- "C"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar4")] <- "D"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar5")] <- "E"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar6")] <- "F"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar7")] <- "G"
seu_pbmc_clean@meta.data$Donor[which(seu_pbmc_clean@meta.data$MULTI == "Bar8")] <- "H"

## Step 2: HashTag classification for cells passing QC -------------------------------------------------------------------------------------------------------------------
## Lane 2
bar.table <- barTable_l2_hash[unlist(strsplit(rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$LaneID == 2)], split="-2")),9:12]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round1 <- findThresh(call.list=bar.table_sweep.list)
round1.calls <- classifyCells(bar.table, q=findQ(res_round1$res, res_round1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]

bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(res_round2$res, res_round2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round3 <- findThresh(call.list=bar.table_sweep.list)
round3.calls <- classifyCells(bar.table, q=findQ(res_round3$res, res_round3$extrema))

final.calls_l2_hash <- c(round3.calls, rep("Negative",length(neg.cells)))
names(final.calls_l2_hash) <- c(names(round3.calls),neg.cells)

## Lane 3
bar.table <- barTable_l3_hash[unlist(strsplit(rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc_clean@meta.data$LaneID == 3)], split="-3")),1:8]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round1 <- findThresh(call.list=bar.table_sweep.list)
round1.calls <- classifyCells(bar.table, q=findQ(res_round1$res, res_round1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]

bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round2 <- findThresh(call.list=bar.table_sweep.list)
round2.calls <- classifyCells(bar.table, q=findQ(res_round2$res, res_round2$extrema))
neg.cells <- c(neg.cells, names(round2.calls)[which(round2.calls == "Negative")])

bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]
bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}
res_round3 <- findThresh(call.list=bar.table_sweep.list)
round3.calls <- classifyCells(bar.table, q=findQ(res_round3$res, res_round3$extrema))

final.calls_l3_hash <- c(round3.calls, rep("Negative",length(neg.cells)))
names(final.calls_l3_hash) <- c(names(round3.calls),neg.cells)

## Add as meta data
seu_pbmc_clean@meta.data[,"HASH"] <- rep("none")
seu_pbmc_clean@meta.data[names(final.calls_l2_hash),"HASH"] <- final.calls_l2_hash
seu_pbmc_clean@meta.data[names(final.calls_l3_hash),"HASH"] <- final.calls_l3_hash


## Step 3: Compute barcode space for lane 3 cells ---------------------------------------------------------------------------------------------------
barTSNE_multi <- barTSNE(barTable_l3_multi[names(final.calls_l3_multi), 1:8])
barTSNE_hash <- barTSNE(barTable_l3_hash[names(final.calls_l3_multi), 1:8])

## Add classification results
barTSNE_multi[,"MULTI"] <- final.calls_l3_multi
barTSNE_hash[,"MULTI"] <- final.calls_l3_hash

## Step 4: souporcell classification ----------------------------------------------------------------------------------------------------------------
## Specify lane 3 cells
cellIDs_souporcell <- rownames(seu_pbmc_clean@meta.data)[which(seu_pbmc@meta.data$LaneID == 3)]
cellIDs_souporcell <- paste(unlist(strsplit(cellIDs_souporcell, split="-3")), "-1", sep="")
write.csv(cellIDs_souporcell, file="cellIDs_souporcell.tsv", quote=F, row.names = F)

## Run souporcell
singularity exec souporcell.sif souporcell_pipeline.py -i SLNR_L3_possorted_genome_bam.bam -b cellIDs_souporcell.tsv -f hg19_3.0.0_genome.fa -t 16 -o SLNR_PBMC_L3_souporcell -k 8

## Read results into R and record
souporcell <- read.table("./souporcell/clusters.tsv", sep="\t", header=T)
souporcell$barcode <- unlist(strsplit(as.character(souporcell$barcode), split="-1"))
rownames(souporcell) <- souporcell$barcode
barTSNE_multi[,"SOUP"] <- souporcell[rownames(barTSNE_multi),"assignment"]
barTSNE_hash[,"SOUP"] <- souporcell[rownames(barTSNE_hash),"assignment"]



## Step 4: demuxEM classification ----------------------------------------------------------------------------------------------------------------
barTable_multi_EM <- t(barTable_multi_EM) ## rows = BCs, cols = cells
barTable_hash_EM <- t(barTable_hash_EM)
write.csv(barTable_multi_EM, file="barTable_multi_EM.csv", quote=F, row.names=T)
write.csv(barTable_hash_EM, file="barTable_hash_EM.csv", quote=F, row.names=T)

...


barTSNE_multi[,"DEMUX"] <- ...
barTSNE_hash[,"DEMUX"] <- ...


####################
## Visualizations ##
####################
## Fig. 2A: Barcode tSNEs colored by (i) deMULTIplex, (ii) demuxEM, and (iii) souporcell
ggplot(barTSNE_multi, aes(x=TSNE1, y=TSNE2, color=MULTI)) + 
  geom_point(size=0.5) +
  scale_color_manual(values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","black","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","black")) +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

## Fig. 2B: Cell Hashing negative proportion aacross PBMC cell types
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

## Fig. 2C: Distribution of Cell Hashing negatives in CD4+ T-cell space
seu_CD4T_ficoll[,"Fig1D"] <- seu_CD4T_ficoll@meta.data$CD4TSubset
seu_CD4T_ficoll$Fig1D[which(seu_CD4T_ficoll@meta.data$HASH == "Negative")] <- "negative"
DimPlot(seu_CD4T_ficoll, group.by="Fig1D", cols=c("grey","beige","grey","black","red"), order="negative") + NoLegend()

