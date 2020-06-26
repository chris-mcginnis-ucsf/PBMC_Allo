############################################################
## In silico genotyping and MULTI-seq/SCMK Demultiplexing ##
## for 8-donor, 7-donor, and Zheng et al PBMC datasets #####
## Chris McGinnis, Gartner Lab, UCSF, 06/26/2020 ###########
############################################################

library(deMULTIplex)

##################
## 8-Donor PBMC ##
##################
## Step 1: MULTI-seq classification for cells passing QC --------------------------------------------------------------------------------------------------------------------------
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


## Step 2: Record results, convert to donor IDs  ----------------------------------------------------------------------------------------------------------------------------------
seu_pbmc_clean@meta.data[,"MULTI"] <- rep("none")
seu_pbmc_clean@meta.data[names(final.calls_multi_l2),"MULTI"] <- final.calls_multi_l2
seu_pbmc_clean@meta.data[names(final.calls_multi_l3),"MULTI"] <- final.calls_multi_l3

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


## Step 3: SCMK classification for cells passing QC -------------------------------------------------------------------------------------------------------------------------------
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


## Step 4: Compute barcode space for lane 3 cells ---------------------------------------------------------------------------------------------------------------------------------
barTSNE_multi <- barTSNE(barTable_l3_multi[names(final.calls_l3_multi), 1:8])
barTSNE_hash <- barTSNE(barTable_l3_hash[names(final.calls_l3_multi), 1:8])
barTSNE_multi[,"MULTI"] <- final.calls_l3_multi
barTSNE_hash[,"MULTI"] <- final.calls_l3_hash

## Step 5: souporcell classification ----------------------------------------------------------------------------------------------------------------------------------------------
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

##################
## 7-Donor PBMC ##
##################
## Step 1: Read in BD SCMK counts -------------------------------------------------------------------------------------------------------------------------------------------------
bd.raw <- Read10X('.')
bd.raw <- bd.raw[,rownames(seu_gh@meta.data)]
bd.raw <- t(as.data.frame(bd.raw))
bd.raw <- bd.raw[,1:7]
library(deMULTIplex)
bd.tsne <- barTSNE(bd.raw)
ggplot(bd.tsne, aes(x=TSNE1, y=TSNE2)) + geom_point() + theme_void() + scale_color_virids()

## Step 2: Perform Classification -------------------------------------------------------------------------------------------------------------------------------------------------
bar.table <- bd.raw
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

ggplot(data=res_round1$res, aes(x=q, y=Proportion, color=Subset)) + 
  geom_line() + 
  theme(legend.position = "none") +
  geom_vline(xintercept=res_round1$extrema, lty=2) + 
  scale_color_manual(values=c("red","black","blue"))

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

ggplot(data=res_round2$res, aes(x=q, y=Proportion, color=Subset)) + 
  geom_line() + 
  theme(legend.position = "none") +
  geom_vline(xintercept=res_round2$extrema, lty=2) + 
  scale_color_manual(values=c("red","black","blue"))

final.calls <- c(round2.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round2.calls),neg.cells)
final.calls <- final.calls[-which(duplicated(final.calls) == F)]

seu_gh@meta.data[,"class"] <- rep("Negative")
seu_gh@meta.data[names(final.calls),"class"] <- final.calls

######################
## Zheng et al PBMC ##
######################
cellIDs_XY <- rownames(seu_zheng_clean@meta.data)[which(seu_zheng_clean@meta.data$LaneID == 1)]
write.csv(cellIDs_XY, file="cellIDs_XY.tsv", quote=F, row.names = F)

singularity exec souporcell.sif souporcell_pipeline.py -i possorted_genome_bam.bam -b cellIDs_XY.tsv -f hg19-3.0.0-genome.fa -t 16 -o BC5050 -k 2

souporcell_XY <- read.table("./soup/clusters.tsv", sep="\t", header=T)
souporcell_XY$barcode <- as.character(souporcell_XY$barcode)
rownames(souporcell_XY) <- souporcell_XY$barcode
souporcell_XY$assignment <- as.character(souporcell_XY$assignment)
ind <- which(souporcell_XY$assignment %ni% c(0,1))
souporcell_XY$assignment[ind] <- "Doublet"
souporcell_XY[which(souporcell_XY$assignment == "1"), "assignment"] <- "X"
souporcell_XY[which(souporcell_XY$assignment == "0"), "assignment"] <- "Y"
souporcell_XY <- souporcell_XY[-1, ]

seu_zheng_clean@meta.data[,"Donor"] <- rep("X")
seu_zheng_clean@meta.data$Donor[which(seu_zheng_clean@meta.data$LaneID == 3)] <- "Y"
seu_zheng_clean@meta.data[rownames(souporcell_XY), "Donor"] <- souporcell_XY$assignment