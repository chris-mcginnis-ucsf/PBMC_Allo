#######################################################
## Jensen Shannon Divergence allorectivity analyses, ##
## 8-donor (Ficoll) and Zheng et al PBMC datasets #####
## Chris McGinnis, Gartner Lab, UCSF 06/26/2020 #######
#######################################################

library(Seurat)
library(MASS)
library(philentropy)

preProcess <- function(seu, cells, PCs) {
  seu <- SubsetData(seu, cells=cells)
  seu <- SCTransform(seu)
  seu <- RunPCA(seu, verbose=FALSE)
  seu <- RunUMAP(seu, dims=1:PCs, verbose=FALSE)
  return(seu)
}

##################
## 8-Donor PBMC ##
##################
## Step 1: Iteratively compute (n=100) JSD UMAPs for each PBMC cell type ------------------------------------------------------------------------------------------------------------------------------
## CD4T
cells_CD4T <- cbind(rownames(seu_CD4T_ficoll@meta.data), seu_CD4T_ficoll@meta.data$Donor_Lane_Type)
temp <- table(cells_CD4T[,2])
x_mem <- min(temp[grep("Memory",names(temp))])
x_act <- min(temp[grep("Activated",names(temp))])
x_naive <- min(temp[grep("Naive",names(temp))])
cells_CD4T_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_CD4T[,2])) {
    ind <- which(cells_CD4T[,2] == j)
    if (length(grep("Memory", j)) == 1) { cells_temp <- c(cells_temp, cells_CD4T[sample(ind, x_mem, replace=F),1]) }
    if (length(grep("Activated", j)) == 1) { cells_temp <- c(cells_temp, cells_CD4T[sample(ind, x_act, replace=F),1]) }
    if (length(grep("Naive", j)) == 1) { cells_temp <- c(cells_temp, cells_CD4T[sample(ind, x_naive, replace=F),1]) }
  }
  cells_CD4T_list[[i]] <- cells_temp
}

umap.list_CD4T <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_CD4T_ficoll, cells=cells_CD4T_list[[i]], PCs = 25)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Lane"] <- seu_temp@meta.data$Donor_Lane
  umap.list_CD4T[[i]] <- umap_temp
}

## CD8T
cells_CD8T <- which(seu_pbmc_ficoll@meta.data$CellType == "CD8T")
cells_CD8T <- cbind(rownames(seu_pbmc_ficoll@meta.data)[cells_CD8T], seu_pbmc_ficoll@meta.data$Donor_Lane[cells_CD8T])
x <- min(table(cells_CD8T[,2]))
cells_CD8T_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_CD8T[,2])) {
    ind <- which(cells_CD8T[,2] == j)
    cells_temp <- c(cells_temp, cells_CD8T[sample(ind, x, replace=F), 1])
  }
  cells_CD8T_list[[i]] <- cells_temp
}

umap.list_CD8T <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_pbmc_ficoll, cells=cells_CD8T_list[[i]], PCs = 15)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Lane"] <- seu_temp@meta.data$Donor_Lane
  umap.list_CD8T[[i]] <- umap_temp
}

## CD14Mono
cells_CD14Mono <- which(seu_pbmc_ficoll@meta.data$CellType == "CD14Mono")
cells_CD14Mono <- cbind(rownames(seu_pbmc_ficoll@meta.data)[cells_CD14Mono], seu_pbmc_ficoll@meta.data$Donor_Lane[cells_CD14Mono])
x <- min(table(cells_CD14Mono[,2]))
cells_CD14Mono_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_CD14Mono[,2])) {
    ind <- which(cells_CD14Mono[,2] == j)
    cells_temp <- c(cells_temp, cells_CD14Mono[sample(ind, x, replace=F), 1])
  }
  cells_CD14Mono_list[[i]] <- cells_temp
}

umap.list_CD14Mono <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_pbmc_ficoll, cells=cells_CD14Mono_list[[i]], PCs = 17)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Lane"] <- seu_temp@meta.data$Donor_Lane
  umap.list_CD14Mono[[i]] <- umap_temp
}

## CD16Mono
cells_CD16Mono <- which(seu_pbmc_ficoll@meta.data$CellType == "CD16Mono")
cells_CD16Mono <- cbind(rownames(seu_pbmc_ficoll@meta.data)[cells_CD16Mono], seu_pbmc_ficoll@meta.data$Donor_Lane[cells_CD16Mono])
x <- min(table(cells_CD16Mono[,2]))
cells_CD16Mono_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_CD16Mono[,2])) {
    ind <- which(cells_CD16Mono[,2] == j)
    cells_temp <- c(cells_temp, cells_CD16Mono[sample(ind, x, replace=F), 1])
  }
  cells_CD16Mono_list[[i]] <- cells_temp
}

umap.list_CD16Mono <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_pbmc_ficoll, cells=cells_CD16Mono_list[[i]], PCs = 15)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Lane"] <- seu_temp@meta.data$Donor_Lane
  umap.list_CD16Mono[[i]] <- umap_temp
}

## NK
cells_NK <- which(seu_pbmc_ficoll@meta.data$CellType == "NK")
cells_NK <- cbind(rownames(seu_pbmc_ficoll@meta.data)[cells_NK], seu_pbmc_ficoll@meta.data$Donor_Lane[cells_NK])
x <- min(table(cells_NK[,2]))
cells_NK_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_NK[,2])) {
    ind <- which(cells_NK[,2] == j)
    cells_temp <- c(cells_temp, cells_NK[sample(ind, x, replace=F), 1])
  }
  cells_NK_list[[i]] <- cells_temp
}

umap.list_NK <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_pbmc_ficoll, cells=cells_NK_list[[i]], PCs = 14)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Lane"] <- seu_temp@meta.data$Donor_Lane
  umap.list_NK[[i]] <- umap_temp
}

## Step 2: Compute JSD for each UMAP ------------------------------------------------------------------------------------------------------------------------------------------------------------------
JSD.div_list_CD4T <- list()
JSD.div_list_CD8T <- list()
JSD.div_list_CD14Mono <- list()
JSD.div_list_CD16Mono <- list()
JSD.div_list_NK <- list()

for (iter in 1:100) {
  print(iter)
  ## Etract embeddings
  umap_CD4T_temp <- umap.list_CD4T[[iter]]
  umap_CD8T_temp <- umap.list_CD8T[[iter]]
  umap_CD14Mono_temp <- umap.list_CD14Mono[[iter]]
  umap_CD16Mono_temp <- umap.list_CD16Mono[[iter]]
  umap_NK_temp <- umap.list_NK[[iter]]
  
  p <- ggplot_build(ggplot(umap_CD4T_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD4T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD4T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD8T_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD8T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD8T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD14Mono_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD14Mono <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD14Mono <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD16Mono_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD16Mono <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD16Mono <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_NK_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_NK <- p$layout$panel_scales_x[[1]]$range$range
  yrange_NK <- p$layout$panel_scales_y[[1]]$range$range
  
  ## Initialize JSD
  kde2d_CD4T_temp <- list()
  kde2d_CD8T_temp <- list()
  kde2d_CD14Mono_temp <- list()
  kde2d_CD16Mono_temp <- list()
  kde2d_NK_temp <- list()
  
  JSD.div_CD4T_temp <- data.frame()
  JSD.div_CD8T_temp <- data.frame()
  JSD.div_CD14Mono_temp <- data.frame()
  JSD.div_CD16Mono_temp <- data.frame()
  JSD.div_NK_temp<- data.frame()
  
  ## Compute GKDE and JSD
  for (i in 1:8){
    for (j in 1:8){
      kde2d_CD4T_temp[[i]] <- kde2d(umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Lane == groups[i])], 1], 
                                    umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Lane == groups[i])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      kde2d_CD4T_temp[[j]] <- kde2d(umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Lane == groups[j])], 1], 
                                    umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Lane == groups[j])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      JSD.div_CD4T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD4T_temp[[i]]$z), as.vector(kde2d_CD4T_temp[[j]]$z)))
      
      kde2d_CD8T_temp[[i]] <- kde2d(umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Lane == groups[i])], 1], 
                                    umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Lane == groups[i])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      kde2d_CD8T_temp[[j]] <- kde2d(umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Lane == groups[j])], 1], 
                                    umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Lane == groups[j])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      JSD.div_CD8T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD8T_temp[[i]]$z), as.vector(kde2d_CD8T_temp[[j]]$z)))
      
      kde2d_CD14Mono_temp[[i]] <- kde2d(umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Lane == groups[i])], 1], 
                                        umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Lane == groups[i])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      kde2d_CD14Mono_temp[[j]] <- kde2d(umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Lane == groups[j])], 1], 
                                        umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Lane == groups[j])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      JSD.div_CD14Mono_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD14Mono_temp[[i]]$z), as.vector(kde2d_CD14Mono_temp[[j]]$z)))
      
      kde2d_CD16Mono_temp[[i]] <- kde2d(umap_CD16Mono_temp[rownames(umap_CD16Mono_temp)[which(umap_CD16Mono_temp$Donor_Lane == groups[i])], 1], 
                                        umap_CD16Mono_temp[rownames(umap_CD16Mono_temp)[which(umap_CD16Mono_temp$Donor_Lane == groups[i])], 2],
                                        n=500, lims=c(xrange_CD16Mono,yrange_CD16Mono))
      kde2d_CD16Mono_temp[[j]] <- kde2d(umap_CD16Mono_temp[rownames(umap_CD16Mono_temp)[which(umap_CD16Mono_temp$Donor_Lane == groups[j])], 1], 
                                        umap_CD16Mono_temp[rownames(umap_CD16Mono_temp)[which(umap_CD16Mono_temp$Donor_Lane == groups[j])], 2],
                                        n=500, lims=c(xrange_CD16Mono,yrange_CD16Mono))
      JSD.div_CD16Mono_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD16Mono_temp[[i]]$z), as.vector(kde2d_CD16Mono_temp[[j]]$z)))
      
      kde2d_NK_temp[[i]] <- kde2d(umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Lane == groups[i])], 1], 
                                  umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Lane == groups[i])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      kde2d_NK_temp[[j]] <- kde2d(umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Lane == groups[j])], 1], 
                                  umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Lane == groups[j])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      JSD.div_NK_temp[i,j] <- JSD(rbind(as.vector(kde2d_NK_temp[[i]]$z), as.vector(kde2d_NK_temp[[j]]$z)))
    }
  }
  
  ## Store results
  colnames(JSD.div_CD4T_temp) <- groups
  rownames(JSD.div_CD4T_temp) <- groups
  colnames(JSD.div_CD8T_temp) <- groups
  rownames(JSD.div_CD8T_temp) <- groups
  colnames(JSD.div_CD14Mono_temp) <- groups
  rownames(JSD.div_CD14Mono_temp) <- groups
  colnames(JSD.div_CD16Mono_temp) <- groups
  rownames(JSD.div_CD16Mono_temp) <- groups
  colnames(JSD.div_NK_temp) <- groups
  rownames(JSD.div_NK_temp) <- groups
  
  JSD.div_CD4T_temp <- apply(JSD.div_CD4T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD4T_temp))))
  JSD.div_CD8T_temp <- apply(JSD.div_CD8T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD8T_temp))))
  JSD.div_CD14Mono_temp <- apply(JSD.div_CD14Mono_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD14Mono_temp))))
  JSD.div_CD16Mono_temp <- apply(JSD.div_CD16Mono_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD16Mono_temp))))
  JSD.div_NK_temp <- apply(JSD.div_NK_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_NK_temp))))
  
  JSD.div_list_CD4T[[iter]] <- JSD.div_CD4T_temp
  JSD.div_list_CD8T[[iter]] <- JSD.div_CD8T_temp
  JSD.div_list_CD14Mono[[iter]] <- JSD.div_CD14Mono_temp
  JSD.div_list_CD16Mono[[iter]] <- JSD.div_CD16Mono_temp
  JSD.div_list_NK[[iter]] <- JSD.div_NK_temp
  
}   

## Step 3: Repeat JSD computation after permuting donor A labels --------------------------------------------------------------------------------------------------------------------------------------
groups <- c("A_1","A_2","A_3","A_4","B_2","B_3","C_2","C_3")
JSD.div_list_CD4T_permute <- list()
JSD.div_list_CD8T_permute <- list()
JSD.div_list_CD14Mono_permute <- list()
JSD.div_list_CD16Mono_permute <- list()
JSD.div_list_NK_permute <- list()

for (iter in 1:100) {
  print(iter)
  # Randomly select UMAP embedding
  umap_CD4T_jsd <- umap.list_CD4T[[sample(1:100, 1)]]
  umap_CD8T_jsd <- umap.list_CD8T[[sample(1:100, 1)]]
  umap_CD14Mono_jsd <- umap.list_CD14Mono[[sample(1:100, 1)]]
  umap_CD16Mono_jsd <- umap.list_CD16Mono[[sample(1:100, 1)]]
  umap_NK_jsd <- umap.list_NK[[sample(1:100, 1)]]
  
  p <- ggplot_build(ggplot(umap_CD4T_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD4T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD4T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD8T_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD8T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD8T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD14Mono_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD14Mono <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD14Mono <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD16Mono_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD16Mono <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD16Mono <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_NK_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_NK <- p$layout$panel_scales_x[[1]]$range$range
  yrange_NK <- p$layout$panel_scales_y[[1]]$range$range
  
  # Permute labels
  umap_CD4T_jsd[,"permute"] <- umap_CD4T_jsd$Donor_Lane
  umap_CD8T_jsd[,"permute"] <- umap_CD8T_jsd$Donor_Lane
  umap_CD14Mono_jsd[,"permute"] <- umap_CD14Mono_jsd$Donor_Lane
  umap_CD16Mono_jsd[,"permute"] <- umap_CD16Mono_jsd$Donor_Lane
  umap_NK_jsd[,"permute"] <- umap_NK_jsd$Donor_Lane
  
  ind_cd4t <- which(umap_CD4T_jsd$permute %in% c("A_1","A_2","A_3","A_4"))
  ind_cd8t <- which(umap_CD8T_jsd$permute %in% c("A_1","A_2","A_3","A_4"))
  ind_cd14 <- which(umap_CD14Mono_jsd$permute %in% c("A_1","A_2","A_3","A_4"))
  ind_cd16 <- which(umap_CD16Mono_jsd$permute %in% c("A_1","A_2","A_3","A_4"))
  ind_nk <- which(umap_NK_jsd$permute %in% c("A_1","A_2","A_3","A_4"))
  
  umap_CD4T_jsd$permute[ind_cd4t] <- sample(umap_CD4T_jsd$permute[ind_cd4t],replace=F,length(ind_cd4t))
  umap_CD8T_jsd$permute[ind_cd8t] <- sample(umap_CD8T_jsd$permute[ind_cd8t],replace=F,length(ind_cd8t))
  umap_CD14Mono_jsd$permute[ind_cd14] <- sample(umap_CD14Mono_jsd$permute[ind_cd14],replace=F,length(ind_cd14))
  umap_CD16Mono_jsd$permute[ind_cd16] <- sample(umap_CD16Mono_jsd$permute[ind_cd16],replace=F,length(ind_cd16))
  umap_NK_jsd$permute[ind_nk] <- sample(umap_NK_jsd$permute[ind_nk],replace=F,length(ind_nk))
  
  # compute JSD
  kde2d_CD4T_temp <- list()
  kde2d_CD8T_temp <- list()
  kde2d_CD14Mono_temp <- list()
  kde2d_CD16Mono_temp <- list()
  kde2d_NK_temp <- list()
  
  JSD.div_CD4T_temp <- data.frame()
  JSD.div_CD8T_temp <- data.frame()
  JSD.div_CD14Mono_temp <- data.frame()
  JSD.div_CD16Mono_temp <- data.frame()
  JSD.div_NK_temp <- data.frame()
  
  for (i in 1:8){
    for (j in 1:8){
      kde2d_CD4T_temp[[i]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[i])], 1], 
                                    umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[i])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      kde2d_CD4T_temp[[j]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[j])], 1], 
                                    umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[j])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      JSD.div_CD4T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD4T_temp[[i]]$z), as.vector(kde2d_CD4T_temp[[j]]$z)))
      
      kde2d_CD8T_temp[[i]] <- kde2d(umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[i])], 1], 
                                    umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[i])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      kde2d_CD8T_temp[[j]] <- kde2d(umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[j])], 1], 
                                    umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[j])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      JSD.div_CD8T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD8T_temp[[i]]$z), as.vector(kde2d_CD8T_temp[[j]]$z)))
      
      kde2d_CD14Mono_temp[[i]] <- kde2d(umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[i])], 1], 
                                        umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[i])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      kde2d_CD14Mono_temp[[j]] <- kde2d(umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[j])], 1], 
                                        umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[j])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      JSD.div_CD14Mono_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD14Mono_temp[[i]]$z), as.vector(kde2d_CD14Mono_temp[[j]]$z)))
      
      kde2d_CD16Mono_temp[[i]] <- kde2d(umap_CD16Mono_jsd[rownames(umap_CD16Mono_jsd)[which(umap_CD16Mono_jsd$permute == groups[i])], 1], 
                                        umap_CD16Mono_jsd[rownames(umap_CD16Mono_jsd)[which(umap_CD16Mono_jsd$permute == groups[i])], 2],
                                        n=500, lims=c(xrange_CD16Mono,yrange_CD16Mono))
      kde2d_CD16Mono_temp[[j]] <- kde2d(umap_CD16Mono_jsd[rownames(umap_CD16Mono_jsd)[which(umap_CD16Mono_jsd$permute == groups[j])], 1], 
                                        umap_CD16Mono_jsd[rownames(umap_CD16Mono_jsd)[which(umap_CD16Mono_jsd$permute == groups[j])], 2],
                                        n=500, lims=c(xrange_CD16Mono,yrange_CD16Mono))
      JSD.div_CD16Mono_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD16Mono_temp[[i]]$z), as.vector(kde2d_CD16Mono_temp[[j]]$z)))
      
      kde2d_NK_temp[[i]] <- kde2d(umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[i])], 1], 
                                  umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[i])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      kde2d_NK_temp[[j]] <- kde2d(umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[j])], 1], 
                                  umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[j])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      JSD.div_NK_temp[i,j] <- JSD(rbind(as.vector(kde2d_NK_temp[[i]]$z), as.vector(kde2d_NK_temp[[j]]$z)))
    }
  }
  
  # store results
  colnames(JSD.div_CD4T_temp) <- groups
  rownames(JSD.div_CD4T_temp) <- groups
  colnames(JSD.div_CD8T_temp) <- groups
  rownames(JSD.div_CD8T_temp) <- groups
  colnames(JSD.div_CD14Mono_temp) <- groups
  rownames(JSD.div_CD14Mono_temp) <- groups
  colnames(JSD.div_CD16Mono_temp) <- groups
  rownames(JSD.div_CD16Mono_temp) <- groups
  colnames(JSD.div_NK_temp) <- groups
  rownames(JSD.div_NK_temp) <- groups
  
  JSD.div_CD4T_temp <- apply(JSD.div_CD4T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD4T_temp))))
  JSD.div_CD8T_temp <- apply(JSD.div_CD8T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD8T_temp))))
  JSD.div_CD14Mono_temp <- apply(JSD.div_CD14Mono_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD14Mono_temp))))
  JSD.div_CD16Mono_temp <- apply(JSD.div_CD16Mono_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD16Mono_temp))))
  JSD.div_NK_temp <- apply(JSD.div_NK_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_NK_temp))))
  
  JSD.div_list_CD4T_permute[[iter]] <- JSD.div_CD4T_temp
  JSD.div_list_CD8T_permute[[iter]] <- JSD.div_CD8T_temp
  JSD.div_list_CD14Mono_permute[[iter]] <- JSD.div_CD14Mono_temp
  JSD.div_list_CD16Mono_permute[[iter]] <- JSD.div_CD16Mono_temp
  JSD.div_list_NK_permute[[iter]] <- JSD.div_NK_temp
}

## Step 4: Record results -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
JSD_og_list <- list()
JSD_og_list[[1]] <- JSD.div_list_CD4T
JSD_og_list[[2]] <- JSD.div_list_CD4T_permute
JSD_og_list[[3]] <- JSD.div_list_CD8T
JSD_og_list[[4]] <- JSD.div_list_CD8T_permute
JSD_og_list[[5]] <- JSD.div_list_CD14Mono
JSD_og_list[[6]] <- JSD.div_list_CD14Mono_permute
JSD_og_list[[7]] <- JSD.div_list_CD16Mono
JSD_og_list[[8]] <- JSD.div_list_CD16Mono_permute
JSD_og_list[[9]] <- JSD.div_list_NK
JSD_og_list[[10]] <- JSD.div_list_NK_permute
names(JSD_og_list) <- c("CD4T","CD4T_permute","CD8T","CD8T_permute","CD14Mono","CD14Mono_permute","CD16Mono","CD16Mono_permute","NK","NK_permute")
save(JSD_og_list, file="JSD_og_list.Robj")

## Step 5: Compute JSD summary statistics -------------------------------------------------------------------------------------------------------------------------------------------------------------
jsd.summary <- as.data.frame(matrix(0L, nrow=25, ncol=4))
colnames(jsd.summary) <- c("CellType", "Group", "JSD", "SD")
jsd.summary$CellType <- rep(c("CD4T","CD8T","CD14Mono","CD16Mono","NK"), each=5)
jsd.summary$Group <- rep(c("Donor","Tech","Mix","Delta","Permute"), 5)

for (celltype in c("CD4T","CD8T","CD14Mono","CD16Mono","NK")) {
  jsd.temp <- get(paste("JSD.div_list_",celltype,sep=""), envir = .GlobalEnv)
  permute.temp <- get(paste("JSD.div_list_",celltype,"_permute",sep=""), envir = .GlobalEnv)
  
  donor.temp <- mean(unlist(lapply(jsd.temp, function(x) mean(x[c("B_2","B_3","C_2","C_3"),"A_1"]))))
  tech.temp <- mean(unlist(lapply(jsd.temp, function(x) mean(c(x["A_4","A_1"], x["A_2","A_3"], x["B_2","B_3"], x["C_2","C_3"])))))
  mix.temp <- mean(unlist(lapply(jsd.temp, function(x) mean(x[c("A_2","A_3"),c("A_1","A_4")]))))
  delta.temp <- abs(tech.temp-mix.temp)
  perm.temp <- mean(unlist(lapply(permute.temp, function(x) mean(x[c("A_2","A_3","A_4"),"A_1"]))))
  
  jsd.summary[which(jsd.summary$CellType == celltype), "JSD"] <- c(donor.temp, tech.temp, mix.temp, delta.temp, perm.temp)
  
  donor.sd.temp <- sd(unlist(lapply(jsd.temp, function(x) mean(x[c("B_2","B_3","C_2","C_3"),"A_1"]))))      
  tech.sd.temp <- sd(unlist(lapply(jsd.temp, function(x) mean(c(x["A_4","A_1"], x["A_2","A_3"], x["B_2","B_3"], x["C_2","C_3"])))))
  mix.sd.temp <- sd(unlist(lapply(jsd.temp, function(x) mean(x[c("A_2","A_3"),c("A_1","A_4")]))))
  delta.sd.temp <- 0
  perm.sd.temp <- sd(unlist(lapply(permute.temp, function(x) mean(x[c("A_2","A_3","A_4"),"A_1"]))))
  
  jsd.summary[which(jsd.summary$CellType == celltype), "SD"] <- c(donor.sd.temp, tech.sd.temp, mix.sd.temp, delta.sd.temp, perm.sd.temp)
  
}

jsd.summary$Group <- factor(jsd.summary$Group, levels=c("Donor","Mix","Tech","Delta","Permute"))
jsd.summary$CellType <- factor(jsd.summary$CellType, levels=c("CD4T","CD8T","CD14Mono","CD16Mono","NK"))


######################
## Zheng et al PBMC ##
######################
## Step 1: Iteratively compute (n=100) JSD UMAPs for each PBMC cell type ------------------------------------------------------------------------------------------------------------------------------
## CD4T
cells_CD4T <- cbind(rownames(seu_zheng_CD4T@meta.data), seu_zheng_CD4T@meta.data$Donor_Mix_Subtype)
temp <- table(cells_CD4T[,2])
x_mem <- min(temp[grep("memory",names(temp))])
x_act <- min(temp[grep("activated",names(temp))])
x_naive <- min(temp[grep("naive",names(temp))])
cells_CD4T_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_CD4T[,2])) {
    ind <- which(cells_CD4T[,2] == j)
    if (length(grep("memory", j)) == 1) { cells_temp <- c(cells_temp, cells_CD4T[sample(ind, x_mem, replace=F),1]) }
    if (length(grep("activated", j)) == 1) { cells_temp <- c(cells_temp, cells_CD4T[sample(ind, x_act, replace=F),1]) }
    if (length(grep("naive", j)) == 1) { cells_temp <- c(cells_temp, cells_CD4T[sample(ind, x_naive, replace=F),1]) }
  }
  cells_CD4T_list[[i]] <- cells_temp
}

umap.list_CD4T <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_zheng_clean, cells=cells_CD4T_list[[i]], PCs = 15)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Mix"] <- seu_temp@meta.data$Donor_Mix
  umap.list_CD4T[[i]] <- umap_temp
}

## CD8T
cells_CD8T <- which(seu_zheng_clean@meta.data$CellType == "CD8T")
cells_CD8T <- cbind(rownames(seu_zheng_clean@meta.data)[cells_CD8T], seu_zheng_clean@meta.data$Donor_Mix[cells_CD8T])
x <- min(table(cells_CD8T[,2]))
cells_CD8T_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_CD8T[,2])) {
    ind <- which(cells_CD8T[,2] == j)
    cells_temp <- c(cells_temp, cells_CD8T[sample(ind, x, replace=F), 1])
  }
  cells_CD8T_list[[i]] <- cells_temp
}

umap.list_CD8T <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_zheng_clean, cells=cells_CD8T_list[[i]], PCs = 13)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Mix"] <- seu_temp@meta.data$Donor_Mix
  umap.list_CD8T[[i]] <- umap_temp
}

## CD14Mono
cells_CD14Mono <- which(seu_zheng_clean@meta.data$CellType == "CD14Mono")
cells_CD14Mono <- cbind(rownames(seu_zheng_clean@meta.data)[cells_CD14Mono], seu_zheng_clean@meta.data$Donor_Mix[cells_CD14Mono])
x <- min(table(cells_CD14Mono[,2]))
cells_CD14Mono_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_CD14Mono[,2])) {
    ind <- which(cells_CD14Mono[,2] == j)
    cells_temp <- c(cells_temp, cells_CD14Mono[sample(ind, x, replace=F), 1])
  }
  cells_CD14Mono_list[[i]] <- cells_temp
}

umap.list_CD14Mono <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_zheng_clean, cells=cells_CD14Mono_list[[i]], PCs = 17)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Mix"] <- seu_temp@meta.data$Donor_Mix
  umap.list_CD14Mono[[i]] <- umap_temp
}

## B
cells_B <- which(seu_zheng_clean@meta.data$CellType == "B")
cells_B <- cbind(rownames(seu_zheng_clean@meta.data)[cells_B], seu_zheng_clean@meta.data$Donor_Mix[cells_B])
x <- min(table(cells_B[,2]))
cells_B_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_B[,2])) {
    ind <- which(cells_B[,2] == j)
    cells_temp <- c(cells_temp, cells_B[sample(ind, x, replace=F), 1])
  }
  cells_B_list[[i]] <- cells_temp
}

umap.list_B <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_zheng_clean, cells=cells_B_list[[i]], PCs = 16)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Mix"] <- seu_temp@meta.data$Donor_Mix
  umap.list_B[[i]] <- umap_temp
}

## NK
cells_NK <- which(seu_zheng_clean@meta.data$CellType == "NK")
cells_NK <- cbind(rownames(seu_zheng_clean@meta.data)[cells_NK], seu_zheng_clean@meta.data$Donor_Mix[cells_NK])
x <- min(table(cells_NK[,2]))
cells_NK_list <- list()
for (i in 1:100) {
  cells_temp <- NULL
  for (j in unique(cells_NK[,2])) {
    ind <- which(cells_NK[,2] == j)
    cells_temp <- c(cells_temp, cells_NK[sample(ind, x, replace=F), 1])
  }
  cells_NK_list[[i]] <- cells_temp
}

umap.list_NK <- list()
for (i in 1:100) {
  print(paste("Iteration ",i,"...",sep=""))
  seu_temp <- preProcess(seu_zheng_clean, cells=cells_NK_list[[i]], PCs = 14)
  umap_temp <- as.data.frame(seu_temp@reductions$umap@cell.embeddings[,1:2])
  umap_temp[,"Donor_Mix"] <- seu_temp@meta.data$Donor_Mix
  umap.list_NK[[i]] <- umap_temp
}


## Step 2: Compute JSD for each UMAP ------------------------------------------------------------------------------------------------------------------------------------------------------------------
groups <- c("X_mix","X_unmix","Y_mix","Y_unmix")
JSD.div_list_zheng_CD4T <- list()
JSD.div_list_zheng_CD8T <- list()
JSD.div_list_zheng_CD14Mono <- list()
JSD.div_list_zheng_B <- list()
JSD.div_list_zheng_NK <- list()

for (iter in 1:100) {
  print(iter)
  ## Etract embeddings
  umap_CD4T_temp <- umap.list_CD4T[[iter]]
  umap_CD8T_temp <- umap.list_CD8T[[iter]]
  umap_CD14Mono_temp <- umap.list_CD14Mono[[iter]]
  umap_B_temp <- umap.list_B[[iter]]
  umap_NK_temp <- umap.list_NK[[iter]]
  
  p <- ggplot_build(ggplot(umap_CD4T_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD4T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD4T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD8T_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD8T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD8T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD14Mono_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD14Mono <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD14Mono <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_NK_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_NK <- p$layout$panel_scales_x[[1]]$range$range
  yrange_NK <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_B_temp, aes(x=UMAP_1, y=UMAP_2)))
  xrange_B <- p$layout$panel_scales_x[[1]]$range$range
  yrange_B <- p$layout$panel_scales_y[[1]]$range$range
  
  ## Initialize JSD
  kde2d_CD4T_temp <- list()
  kde2d_CD8T_temp <- list()
  kde2d_CD14Mono_temp <- list()
  kde2d_B_temp <- list()
  kde2d_NK_temp <- list()
  
  JSD.div_CD4T_temp <- data.frame()
  JSD.div_CD8T_temp <- data.frame()
  JSD.div_CD14Mono_temp <- data.frame()
  JSD.div_B_temp <- data.frame()
  JSD.div_NK_temp<- data.frame()
  
  ## Compute GKDE and JSD
  for (i in 1:4){
    for (j in 1:4){
      kde2d_CD4T_temp[[i]] <- kde2d(umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Mix == groups[i])], 1], 
                                    umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Mix == groups[i])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      kde2d_CD4T_temp[[j]] <- kde2d(umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Mix == groups[j])], 1], 
                                    umap_CD4T_temp[rownames(umap_CD4T_temp)[which(umap_CD4T_temp$Donor_Mix == groups[j])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      JSD.div_CD4T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD4T_temp[[i]]$z), as.vector(kde2d_CD4T_temp[[j]]$z)))
      
      kde2d_CD8T_temp[[i]] <- kde2d(umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Mix == groups[i])], 1], 
                                    umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Mix == groups[i])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      kde2d_CD8T_temp[[j]] <- kde2d(umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Mix == groups[j])], 1], 
                                    umap_CD8T_temp[rownames(umap_CD8T_temp)[which(umap_CD8T_temp$Donor_Mix == groups[j])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      JSD.div_CD8T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD8T_temp[[i]]$z), as.vector(kde2d_CD8T_temp[[j]]$z)))
      
      kde2d_CD14Mono_temp[[i]] <- kde2d(umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Mix == groups[i])], 1], 
                                        umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Mix == groups[i])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      kde2d_CD14Mono_temp[[j]] <- kde2d(umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Mix == groups[j])], 1], 
                                        umap_CD14Mono_temp[rownames(umap_CD14Mono_temp)[which(umap_CD14Mono_temp$Donor_Mix == groups[j])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      JSD.div_CD14Mono_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD14Mono_temp[[i]]$z), as.vector(kde2d_CD14Mono_temp[[j]]$z)))
      
      kde2d_B_temp[[i]] <- kde2d(umap_B_temp[rownames(umap_B_temp)[which(umap_B_temp$Donor_Mix == groups[i])], 1], 
                                 umap_B_temp[rownames(umap_B_temp)[which(umap_B_temp$Donor_Mix == groups[i])], 2],
                                 n=500, lims=c(xrange_B,yrange_B))
      kde2d_B_temp[[j]] <- kde2d(umap_B_temp[rownames(umap_B_temp)[which(umap_B_temp$Donor_Mix == groups[j])], 1], 
                                 umap_B_temp[rownames(umap_B_temp)[which(umap_B_temp$Donor_Mix == groups[j])], 2],
                                 n=500, lims=c(xrange_B,yrange_B))
      JSD.div_B_temp[i,j] <- JSD(rbind(as.vector(kde2d_B_temp[[i]]$z), as.vector(kde2d_B_temp[[j]]$z)))
      
      kde2d_NK_temp[[i]] <- kde2d(umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Mix == groups[i])], 1], 
                                  umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Mix == groups[i])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      kde2d_NK_temp[[j]] <- kde2d(umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Mix == groups[j])], 1], 
                                  umap_NK_temp[rownames(umap_NK_temp)[which(umap_NK_temp$Donor_Mix == groups[j])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      JSD.div_NK_temp[i,j] <- JSD(rbind(as.vector(kde2d_NK_temp[[i]]$z), as.vector(kde2d_NK_temp[[j]]$z)))
    }
  }
  
  ## Store results
  colnames(JSD.div_CD4T_temp) <- groups
  rownames(JSD.div_CD4T_temp) <- groups
  colnames(JSD.div_CD8T_temp) <- groups
  rownames(JSD.div_CD8T_temp) <- groups
  colnames(JSD.div_CD14Mono_temp) <- groups
  rownames(JSD.div_CD14Mono_temp) <- groups
  colnames(JSD.div_B_temp) <- groups
  rownames(JSD.div_B_temp) <- groups
  colnames(JSD.div_NK_temp) <- groups
  rownames(JSD.div_NK_temp) <- groups
  
  JSD.div_CD4T_temp <- apply(JSD.div_CD4T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD4T_temp))))
  JSD.div_CD8T_temp <- apply(JSD.div_CD8T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD8T_temp))))
  JSD.div_CD14Mono_temp <- apply(JSD.div_CD14Mono_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD14Mono_temp))))
  JSD.div_B_temp <- apply(JSD.div_B_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_B_temp))))
  JSD.div_NK_temp <- apply(JSD.div_NK_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_NK_temp))))
  
  JSD.div_list_zheng_CD4T[[iter]] <- JSD.div_CD4T_temp
  JSD.div_list_zheng_CD8T[[iter]] <- JSD.div_CD8T_temp
  JSD.div_list_zheng_CD14Mono[[iter]] <- JSD.div_CD14Mono_temp
  JSD.div_list_zheng_B[[iter]] <- JSD.div_B_temp
  JSD.div_list_zheng_NK[[iter]] <- JSD.div_NK_temp
  
}    

## Step 3: Repeat JSD computation after permuting donor A labels --------------------------------------------------------------------------------------------------------------------------------------
groups <- c("X_unmix", "X_mix", "Y_unmix", "Y_mix")
JSD.div_list_zheng_CD4T_permute <- list()
JSD.div_list_zheng_CD8T_permute <- list()
JSD.div_list_zheng_CD14Mono_permute <- list()
JSD.div_list_zheng_B_permute <- list()
JSD.div_list_zheng_NK_permute <- list()

for (iter in 1:100) {
  print(iter)
  # Randomly select UMAP embedding
  umap_CD4T_jsd <- umap.list_zheng_CD4T[[sample(1:100, 1)]]
  umap_CD8T_jsd <- umap.list_zheng_CD8T[[sample(1:100, 1)]]
  umap_CD14Mono_jsd <- umap.list_zheng_CD14Mono[[sample(1:100, 1)]]
  umap_B_jsd <- umap.list_zheng_B[[sample(1:100, 1)]]
  umap_NK_jsd <- umap.list_zheng_NK[[sample(1:100, 1)]]
  
  p <- ggplot_build(ggplot(umap_CD4T_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD4T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD4T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD8T_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD8T <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD8T <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_CD14Mono_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_CD14Mono <- p$layout$panel_scales_x[[1]]$range$range
  yrange_CD14Mono <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_B_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_B <- p$layout$panel_scales_x[[1]]$range$range
  yrange_B <- p$layout$panel_scales_y[[1]]$range$range
  
  p <- ggplot_build(ggplot(umap_NK_jsd, aes(x=UMAP_1, y=UMAP_2)))
  xrange_NK <- p$layout$panel_scales_x[[1]]$range$range
  yrange_NK <- p$layout$panel_scales_y[[1]]$range$range
  
  # Permute labels
  umap_CD4T_jsd[,"permute"] <- umap_CD4T_jsd$Donor_Mix
  umap_CD8T_jsd[,"permute"] <- umap_CD8T_jsd$Donor_Mix
  umap_CD14Mono_jsd[,"permute"] <- umap_CD14Mono_jsd$Donor_Mix
  umap_B_jsd[,"permute"] <- umap_B_jsd$Donor_Mix
  umap_NK_jsd[,"permute"] <- umap_NK_jsd$Donor_Mix
  
  ind_cd4t_X <- grep("^X", umap_CD4T_jsd$permute)
  ind_cd4t_Y <- grep("^Y", umap_CD4T_jsd$permute)
  ind_cd8t_X <- grep("^X", umap_CD8T_jsd$Donor_Mix)
  ind_cd8t_Y <- grep("^Y", umap_CD8T_jsd$Donor_Mix)
  ind_cd14_X <- grep("^X", umap_CD14Mono_jsd$Donor_Mix)
  ind_cd14_Y <- grep("^Y", umap_CD14Mono_jsd$Donor_Mix)
  ind_b_X <- grep("^X", umap_B_jsd$Donor_Mix)
  ind_b_Y <- grep("^Y", umap_B_jsd$Donor_Mix)
  ind_nk_X <- grep("^X", umap_NK_jsd$Donor_Mix)
  ind_nk_Y <- grep("^Y", umap_NK_jsd$Donor_Mix)
  
  umap_CD4T_jsd$permute[ind_cd4t_X] <- sample(umap_CD4T_jsd$permute[ind_cd4t_X],replace=F,length(ind_cd4t_X))
  umap_CD4T_jsd$permute[ind_cd4t_Y] <- sample(umap_CD4T_jsd$permute[ind_cd4t_Y],replace=F,length(ind_cd4t_Y))
  umap_CD8T_jsd$permute[ind_cd8t_X] <- sample(umap_CD8T_jsd$permute[ind_cd8t_X],replace=F,length(ind_cd8t_X))
  umap_CD8T_jsd$permute[ind_cd8t_Y] <- sample(umap_CD8T_jsd$permute[ind_cd8t_Y],replace=F,length(ind_cd8t_Y))
  umap_CD14Mono_jsd$permute[ind_cd14_X] <- sample(umap_CD14Mono_jsd$permute[ind_cd14_X],replace=F,length(ind_cd14_X))
  umap_CD14Mono_jsd$permute[ind_cd14_Y] <- sample(umap_CD14Mono_jsd$permute[ind_cd14_Y],replace=F,length(ind_cd14_Y))
  umap_B_jsd$permute[ind_b_X] <- sample(umap_B_jsd$permute[ind_b_X],replace=F,length(ind_b_X))
  umap_B_jsd$permute[ind_b_Y] <- sample(umap_B_jsd$permute[ind_b_Y],replace=F,length(ind_b_Y))
  umap_NK_jsd$permute[ind_nk_X] <- sample(umap_NK_jsd$permute[ind_nk_X],replace=F,length(ind_nk_X))
  umap_NK_jsd$permute[ind_nk_Y] <- sample(umap_NK_jsd$permute[ind_nk_Y],replace=F,length(ind_nk_Y))
  
  # compute JSD
  kde2d_CD4T_temp <- list()
  kde2d_CD8T_temp <- list()
  kde2d_CD14Mono_temp <- list()
  kde2d_B_temp <- list()
  kde2d_NK_temp <- list()
  
  JSD.div_CD4T_temp <- data.frame()
  JSD.div_CD8T_temp <- data.frame()
  JSD.div_CD14Mono_temp <- data.frame()
  JSD.div_B_temp <- data.frame()
  JSD.div_NK_temp<- data.frame()
  
  for (i in 1:4){
    for (j in 1:4){
      kde2d_CD4T_temp[[i]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[i])], 1], 
                                    umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[i])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      kde2d_CD4T_temp[[j]] <- kde2d(umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[j])], 1], 
                                    umap_CD4T_jsd[rownames(umap_CD4T_jsd)[which(umap_CD4T_jsd$permute == groups[j])], 2],
                                    n=500, lims=c(xrange_CD4T,yrange_CD4T))
      JSD.div_CD4T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD4T_temp[[i]]$z), as.vector(kde2d_CD4T_temp[[j]]$z)))
      
      kde2d_CD8T_temp[[i]] <- kde2d(umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[i])], 1], 
                                    umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[i])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      kde2d_CD8T_temp[[j]] <- kde2d(umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[j])], 1], 
                                    umap_CD8T_jsd[rownames(umap_CD8T_jsd)[which(umap_CD8T_jsd$permute == groups[j])], 2],
                                    n=500, lims=c(xrange_CD8T,yrange_CD8T))
      JSD.div_CD8T_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD8T_temp[[i]]$z), as.vector(kde2d_CD8T_temp[[j]]$z)))
      
      kde2d_CD14Mono_temp[[i]] <- kde2d(umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[i])], 1], 
                                        umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[i])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      kde2d_CD14Mono_temp[[j]] <- kde2d(umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[j])], 1], 
                                        umap_CD14Mono_jsd[rownames(umap_CD14Mono_jsd)[which(umap_CD14Mono_jsd$permute == groups[j])], 2],
                                        n=500, lims=c(xrange_CD14Mono,yrange_CD14Mono))
      JSD.div_CD14Mono_temp[i,j] <- JSD(rbind(as.vector(kde2d_CD14Mono_temp[[i]]$z), as.vector(kde2d_CD14Mono_temp[[j]]$z)))
      
      kde2d_B_temp[[i]] <- kde2d(umap_B_jsd[rownames(umap_B_jsd)[which(umap_B_jsd$permute == groups[i])], 1], 
                                 umap_B_jsd[rownames(umap_B_jsd)[which(umap_B_jsd$permute == groups[i])], 2],
                                 n=500, lims=c(xrange_B,yrange_B))
      kde2d_B_temp[[j]] <- kde2d(umap_B_jsd[rownames(umap_B_jsd)[which(umap_B_jsd$permute == groups[j])], 1], 
                                 umap_B_jsd[rownames(umap_B_jsd)[which(umap_B_jsd$permute == groups[j])], 2],
                                 n=500, lims=c(xrange_B,yrange_B))
      JSD.div_B_temp[i,j] <- JSD(rbind(as.vector(kde2d_B_temp[[i]]$z), as.vector(kde2d_B_temp[[j]]$z)))
      
      kde2d_NK_temp[[i]] <- kde2d(umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[i])], 1], 
                                  umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[i])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      kde2d_NK_temp[[j]] <- kde2d(umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[j])], 1], 
                                  umap_NK_jsd[rownames(umap_NK_jsd)[which(umap_NK_jsd$permute == groups[j])], 2],
                                  n=500, lims=c(xrange_NK,yrange_NK))
      JSD.div_NK_temp[i,j] <- JSD(rbind(as.vector(kde2d_NK_temp[[i]]$z), as.vector(kde2d_NK_temp[[j]]$z)))
    }
  }
  
  # store results
  colnames(JSD.div_CD4T_temp) <- groups
  rownames(JSD.div_CD4T_temp) <- groups
  colnames(JSD.div_CD8T_temp) <- groups
  rownames(JSD.div_CD8T_temp) <- groups
  colnames(JSD.div_CD14Mono_temp) <- groups
  rownames(JSD.div_CD14Mono_temp) <- groups
  colnames(JSD.div_B_temp) <- groups
  rownames(JSD.div_B_temp) <- groups
  colnames(JSD.div_NK_temp) <- groups
  rownames(JSD.div_NK_temp) <- groups
  
  JSD.div_CD4T_temp <- apply(JSD.div_CD4T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD4T_temp))))
  JSD.div_CD8T_temp <- apply(JSD.div_CD8T_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD8T_temp))))
  JSD.div_CD14Mono_temp <- apply(JSD.div_CD14Mono_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_CD14Mono_temp))))
  JSD.div_B_temp <- apply(JSD.div_B_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_B_temp))))
  JSD.div_NK_temp <- apply(JSD.div_NK_temp, 1, function(x) rescale(x, to=c(0,1),from=c(0,max(JSD.div_NK_temp))))
  
  JSD.div_list_zheng_CD4T_permute[[iter]] <- JSD.div_CD4T_temp
  JSD.div_list_zheng_CD8T_permute[[iter]] <- JSD.div_CD8T_temp
  JSD.div_list_zheng_CD14Mono_permute[[iter]] <- JSD.div_CD14Mono_temp
  JSD.div_list_zheng_B_permute[[iter]] <- JSD.div_B_temp
  JSD.div_list_zheng_NK_permute[[iter]] <- JSD.div_NK_temp
}


## Step 4: Record results -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
JSD_zheng_list <- list()
JSD_zheng_list[[1]] <- JSD.div_list_zheng_CD4T
JSD_zheng_list[[2]] <- JSD.div_list_zheng_CD4T_permute
JSD_zheng_list[[3]] <- JSD.div_list_zheng_CD8T
JSD_zheng_list[[4]] <- JSD.div_list_zheng_CD8T_permute
JSD_zheng_list[[5]] <- JSD.div_list_zheng_CD14Mono
JSD_zheng_list[[6]] <- JSD.div_list_zheng_CD14Mono_permute
JSD_zheng_list[[7]] <- JSD.div_list_zheng_B
JSD_zheng_list[[8]] <- JSD.div_list_zheng_B_permute
JSD_zheng_list[[9]] <- JSD.div_list_zheng_NK
JSD_zheng_list[[10]] <- JSD.div_list_zheng_NK_permute
names(JSD_zheng_list) <- c("CD4T","CD4T_permute","CD8T","CD8T_permute","CD14Mono","CD14Mono_permute","B","B_permute","NK","NK_permute")
save(JSD_zheng_list, file="JSD_zheng_list.Robj")

## Step 5: Compute JSD summary statistics -------------------------------------------------------------------------------------------------------------------------------------------------------------
jsd.summary_zheng <- as.data.frame(matrix(0L, nrow=15, ncol=4))
colnames(jsd.summary_zheng) <- c("CellType", "Group", "JSD", "SD")
jsd.summary_zheng$CellType <- rep(c("CD4T","CD8T","CD14Mono","B","NK"), each=3)
jsd.summary_zheng$Group <- rep(c("Donor","Mix","Permute"), 5)

for (celltype in c("CD4T","CD8T","CD14Mono","B","NK")) {
  jsd.temp <- get(paste("JSD.div_list_zheng_",celltype,sep=""), envir = .GlobalEnv)
  permute.temp <- get(paste("JSD.div_list_zheng_",celltype,"_permute",sep=""), envir = .GlobalEnv)
  temp <- c(mean(unlist(lapply(jsd.temp, function(x) mean(x["X_mix","Y_mix"],x["X_unmix","Y_unmix"])))),
            mean(unlist(lapply(jsd.temp, function(x) mean(x["X_mix","X_unmix"],x["Y_mix","Y_unmix"])))),
            mean(unlist(lapply(permute.temp, function(x) mean(x["X_mix","X_unmix"],x["Y_mix","Y_unmix"])))))
  
  jsd.summary_zheng[which(jsd.summary_zheng$CellType == celltype), "JSD"] <- temp
  temp <- c(sd(unlist(lapply(jsd.temp, function(x) mean(x["X_mix","Y_mix"],x["X_unmix","Y_unmix"])))),
            sd(unlist(lapply(jsd.temp, function(x) mean(x["X_mix","X_unmix"],x["Y_mix","Y_unmix"])))),
            sd(unlist(lapply(permute.temp, function(x) mean(x["X_mix","X_unmix"],x["Y_mix","Y_unmix"])))))
  jsd.summary_zheng[which(jsd.summary_zheng$CellType == celltype), "SD"] <- temp
}

jsd.summary_zheng$Group <- factor(jsd.summary_zheng$Group, levels=c("Donor","Mix","Permute"))
jsd.summary_zheng$CellType <- factor(jsd.summary_zheng$CellType, levels=c("CD4T","CD8T","CD14Mono","NK","B"))