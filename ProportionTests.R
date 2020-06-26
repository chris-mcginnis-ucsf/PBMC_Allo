#################################################################
## PBMC cell type and CD4+ T-cell subtype proportion analyses, ##
## (un)mixed cells in 8-donor and Zheng et al PBMC datasets #####
## Chris McGinnis, Gartner Lab, UCSF, 06/26/2020 ################
#################################################################

library(ggplot2)

## Step 1: Compute cell type frequencies, 8-donor PBMC ----------------------------------------------------------------------------------------------------------------------------
anno_freq <- as.data.frame(matrix(0L, nrow=56, ncol=3))
colnames(anno_freq) <- c("DonorLane","CellType","Frequency")
anno_freq$DonorLane <- rep(unique(seu_pbmc_ficoll@meta.data$Donor_Lane),7)
anno_freq$CellType <- rep(unique(seu_pbmc_ficoll@meta.data$CellType), each=8)
for (i in 1:nrow(anno_freq)) {
  dl.temp <- anno_freq$DonorLane[i]
  ct.temp <- anno_freq$CellType[i]
  x <- length(which(seu_pbmc_ficoll@meta.data$Donor_Lane == dl.temp & seu_pbmc_ficoll@meta.data$CellType == ct.temp))
  y <- length(which(seu_pbmc_ficoll@meta.data$Donor_Lane == dl.temp))
  anno_freq$Frequency[i] <- x/y
}


## Step 2: Compute cell type frequencies, Zheng et al PBMC ------------------------------------------------------------------------------------------------------------------------
anno_freq_zheng <- as.data.frame(matrix(0L, nrow=32, ncol=3))
colnames(anno_freq_zheng) <- c("Donor_Mix","CellType","Frequency")
temp.meta <- seu_zheng_clean@meta.data[which(seu_zheng_clean@meta.data$CellType != "LowFrequency"), ]
anno_freq_zheng$DonorLane <- rep(unique(temp.meta$Donor_Mix),8)
anno_freq_zheng$CellType <- rep(unique(temp.meta$CellType), each=4)
for (i in 1:nrow(anno_freq_zheng)) {
  dm.temp <- anno_freq_zheng$Donor_Mix[i]
  ct.temp <- anno_freq_zheng$CellType[i]
  x <- length(which(temp.meta$Donor_Mix == dm.temp & temp.meta$CellType == ct.temp))
  y <- length(which(stemp.meta$Donor_Mix == dm.temp))
  anno_freq_zheng$Frequency[i] <- x/y
}


## Step 3: Compute CD4+ T-cell subtype frequencies, 8-donor PBMC ------------------------------------------------------------------------------------------------------------------
anno_freq_CD4T <- as.data.frame(matrix(0L, nrow=24, ncol=3))
colnames(anno_freq_CD4T) <- c("DonorLane","Subtype","Frequency")
temp.meta <- seu_CD4T_ficoll@meta.data[which(seu_CD4T_ficoll@meta.data$CD4TSubset != 0), ]
anno_freq_CD4T$DonorLane <- rep(unique(temp.meta$Donor_Lane),3)
anno_freq_CD4T$Subtype <- rep(unique(temp.meta$CD4TSubset), each=8)
for (i in 1:nrow(anno_freq_CD4T)) {
  dl.temp <- anno_freq_CD4T$DonorLane[i]
  ct.temp <- anno_freq_CD4T$Subtype[i]
  x <- length(which(temp.meta$Donor_Lane == dl.temp & temp.meta$CD4TSubset == ct.temp))
  y <- length(which(temp.meta$Donor_Lane == dl.temp))
  anno_freq_CD4T$Frequency[i] <- x/y
}


## Step 4: Compute CD4+ T-cell subtype frequencies, Zheng et al PBMC --------------------------------------------------------------------------------------------------------------
anno_freq_zheng_CD4T <- as.data.frame(matrix(0L, nrow=12, ncol=3))
colnames(anno_freq_zheng_CD4T) <- c("Donor_Mix","Subtype","Frequency")
temp.meta <- seu_zheng_CD4T@meta.data
anno_freq_zheng_CD4T$Donor_Mix <- rep(unique(temp.meta$Donor_Mix),3)
anno_freq_zheng_CD4T$Subtype <- rep(unique(temp.meta$Subtype), each=4)
for (i in 1:nrow(anno_freq_zheng_CD4T)) {
  dm.temp <- anno_freq_zheng_CD4T$Donor_Mix[i]
  ct.temp <- anno_freq_zheng_CD4T$Subtype[i]
  x <- length(which(temp.meta$Donor_Mix == dm.temp & temp.meta$Subtype == ct.temp))
  y <- length(which(temp.meta$Donor_Mix == dm.temp))
  anno_freq_zheng_CD4T$Frequency[i] <- x/y
}


## Step 5: PBMC cell type pairwise proportion test, 8-donor PBMC ------------------------------------------------------------------------------------------------------------------
xtab_list <- list()
temp <- seu_pbmc_ficoll@meta.data
x <- 0
for (i in unique(temp$CellType)) {
  print(i)
  xtab <- as.table(rbind(c(length(which(temp$Donor_Lane == "A_1" & temp$CellType != i)), length(which(temp$Donor_Lane == "A_2" & temp$CellType != i)),
                           length(which(temp$Donor_Lane == "A_3" & temp$CellType != i)), length(which(temp$Donor_Lane == "A_4" & temp$CellType != i)),
                           length(which(temp$Donor_Lane == "B_2" & temp$CellType != i)), length(which(temp$Donor_Lane == "B_3" & temp$CellType != i)),
                           length(which(temp$Donor_Lane == "C_2" & temp$CellType != i)), length(which(temp$Donor_Lane == "C_3" & temp$CellType != i))),
                         c(length(which(temp$Donor_Lane == "A_1" & temp$CellType == i)), length(which(temp$Donor_Lane == "A_2" & temp$CellType == i)),
                           length(which(temp$Donor_Lane == "A_3" & temp$CellType == i)), length(which(temp$Donor_Lane == "A_4" & temp$CellType == i)),
                           length(which(temp$Donor_Lane == "B_2" & temp$CellType == i)), length(which(temp$Donor_Lane == "B_3" & temp$CellType == i)),
                           length(which(temp$Donor_Lane == "C_2" & temp$CellType == i)), length(which(temp$Donor_Lane == "C_3" & temp$CellType == i)))))
  
  dimnames(xtab) <- list(CellType = c(paste("Not ",i,sep=""), i),
                         Class = c("A_1","A_2","A_3","A_4","B_2","B_3","C_2","C_3"))
  x <- x + 1
  xtab_list[[x]] <- xtab
  names(xtab_list)[x] <- i
}

proptest_data <- as.data.frame(matrix(0L, nrow=28*7, ncol=4))
colnames(proptest_data) <- c("group1","group2","p.adj","CellType")
for (i in 1:length(xtab_list)) { 
  temp <- as.data.frame(pairwise_prop_test(xtab_list[[i]]))
  proptest_data[(1+(i-1)*28):(i*28), c("group1","group2","p.adj")] <- temp[,c("group1","group2","p.adj")]
  proptest_data[(1+(i-1)*28):(i*28), "CellType"] <- rep(names(xtab_list)[i], 28)
}
proptest_data$p.adj <- -log10(proptest_data$p.adj)
proptest_data[,"Bin"] <- rep(0L)
proptest_data$Bin[which(proptest_data$p.adj > 2)] <- 1
proptest_data$Bin[which(proptest_data$p.adj >= 5)] <- 2
proptest_data$Bin[which(proptest_data$p.adj >= 10)] <- 3
proptest_data$Bin[which(proptest_data$p.adj >= 15)] <- 4
proptest_data$Bin[which(proptest_data$p.adj >= 20)] <- 5


## Step 6: CD4+ T-cell subtype pairwise proportion test, 8-donor PBMC -------------------------------------------------------------------------------------------------------------
xtab_list <- list()
temp.meta <- seu_CD4T_ficoll@meta.data[which(seu_CD4T_ficoll@meta.data$CD4TSubset != 0), ]
x <- 0
for (i in unique(temp.meta$CD4TSubset)) {
  print(i)
  xtab <- as.table(rbind(c(length(which(temp.meta$Donor_Lane == "A_1" & temp.meta$CD4TSubset != i)), length(which(temp.meta$Donor_Lane == "A_2" & temp.meta$CD4TSubset != i)),
                           length(which(temp.meta$Donor_Lane == "A_3" & temp.meta$CD4TSubset != i)), length(which(temp.meta$Donor_Lane == "A_4" & temp.meta$CD4TSubset != i)),
                           length(which(temp.meta$Donor_Lane == "B_2" & temp.meta$CD4TSubset != i)), length(which(temp.meta$Donor_Lane == "B_3" & temp.meta$CD4TSubset != i)),
                           length(which(temp.meta$Donor_Lane == "C_2" & temp.meta$CD4TSubset != i)), length(which(temp.meta$Donor_Lane == "C_3" & temp.meta$CD4TSubset != i))),
                         c(length(which(temp.meta$Donor_Lane == "A_1" & temp.meta$CD4TSubset == i)), length(which(temp.meta$Donor_Lane == "A_2" & temp.meta$CD4TSubset == i)),
                           length(which(temp.meta$Donor_Lane == "A_3" & temp.meta$CD4TSubset == i)), length(which(temp.meta$Donor_Lane == "A_4" & temp.meta$CD4TSubset == i)),
                           length(which(temp.meta$Donor_Lane == "B_2" & temp.meta$CD4TSubset == i)), length(which(temp.meta$Donor_Lane == "B_3" & temp.meta$CD4TSubset == i)),
                           length(which(temp.meta$Donor_Lane == "C_2" & temp.meta$CD4TSubset == i)), length(which(temp.meta$Donor_Lane == "C_3" & temp.meta$CD4TSubset == i)))))
  
  dimnames(xtab) <- list(CD4TSubset = c(paste("Not ",i,sep=""), i),
                         Class = c("A_1","A_2","A_3","A_4","B_2","B_3","C_2","C_3"))
  x <- x + 1
  xtab_list[[x]] <- xtab
  names(xtab_list)[x] <- i
}

proptest_CD4T_data <- as.data.frame(matrix(0L, nrow=28*3, ncol=4))
colnames(proptest_CD4T_data) <- c("group1","group2","p.adj","Subset")
for (i in 1:length(xtab_list)) { 
  temp <- as.data.frame(pairwise_prop_test(xtab_list[[i]]))
  proptest_CD4T_data[(1+(i-1)*28):(i*28), c("group1","group2","p.adj")] <- temp[,c("group1","group2","p.adj")]
  proptest_CD4T_data[(1+(i-1)*28):(i*28), "Subset"] <- rep(names(xtab_list)[i], 28)
}
proptest_CD4T_data$p.adj <- -log10(proptest_CD4T_data$p.adj)
proptest_CD4T_data[,"Bin"] <- rep(0L)
proptest_CD4T_data$Bin[which(proptest_CD4T_data$p.adj > 2)] <- 1
proptest_CD4T_data$Bin[which(proptest_CD4T_data$p.adj >= 4)] <- 2
proptest_CD4T_data$Bin[which(proptest_CD4T_data$p.adj >= 6)] <- 3
proptest_CD4T_data$Bin[which(proptest_CD4T_data$p.adj >= 8)] <- 4
proptest_CD4T_data$Bin[which(proptest_CD4T_data$p.adj >= 10)] <- 5


## Step 7: PBMC cell type pairwise proportion test, Zheng et al PBMC --------------------------------------------------------------------------------------------------------------
xtab_zheng_list <- list()
temp <- seu_zheng_clean@meta.data
x <- 0
for (i in unique(temp$CellType)) {
  print(i)
  xtab <- as.table(rbind(c(length(which(temp$Donor_Mix == "X_mix" & temp$CellType != i)), length(which(temp$Donor_Mix == "X_unmix" & temp$CellType != i)),
                           length(which(temp$Donor_Mix == "Y_mix" & temp$CellType != i)), length(which(temp$Donor_Mix == "Y_unmix" & temp$CellType != i))),
                         c(length(which(temp$Donor_Mix == "X_mix" & temp$CellType == i)), length(which(temp$Donor_Mix == "X_unmix" & temp$CellType == i)),
                           length(which(temp$Donor_Mix == "Y_mix" & temp$CellType == i)), length(which(temp$Donor_Mix == "Y_unmix" & temp$CellType == i)))))
  
  dimnames(xtab) <- list(CellType = c(paste("Not ",i,sep=""), i),
                         Class = c("X_mix","X_unmix","Y_mix","Y_unmix"))
  x <- x + 1
  xtab_zheng_list[[x]] <- xtab
  names(xtab_zheng_list)[x] <- i
}

proptest_zheng_data <- as.data.frame(matrix(0L, nrow=6*8, ncol=4))
colnames(proptest_zheng_data) <- c("group1","group2","p.adj","CellType")
for (i in 1:length(xtab_zheng_list)) { 
  temp <- as.data.frame(pairwise_prop_test(xtab_zheng_list[[i]]))
  proptest_zheng_data[(1+(i-1)*6):(i*6), c("group1","group2","p.adj")] <- temp[,c("group1","group2","p.adj")]
  proptest_zheng_data[(1+(i-1)*6):(i*6), "CellType"] <- rep(names(xtab_zheng_list)[i], 6)
}
proptest_zheng_data$p.adj <- -log10(proptest_zheng_data$p.adj)
proptest_zheng_data[,"Bin"] <- rep(0L)
proptest_zheng_data$Bin[which(proptest_zheng_data$p.adj > 2)] <- 1
proptest_zheng_data$Bin[which(proptest_zheng_data$p.adj >= 10)] <- 2
proptest_zheng_data$Bin[which(proptest_zheng_data$p.adj >= 100)] <- 3
proptest_zheng_data$Bin[which(proptest_zheng_data$p.adj >= 200)] <- 4
proptest_zheng_data$Bin[which(proptest_zheng_data$p.adj >= 300)] <- 5


## Step 8: CD4+ T-cell subtype pairwise proportion test, Zheng et al PBMC ---------------------------------------------------------------------------------------------------------
xtab_zheng_list <- list()
temp <- seu_zheng_CD4T@meta.data
x <- 0
for (i in unique(temp$Subtype)) {
  print(i)
  xtab <- as.table(rbind(c(length(which(temp$Donor_Mix == "X_mix" & temp$Subtype != i)), length(which(temp$Donor_Mix == "X_unmix" & temp$Subtype != i)),
                           length(which(temp$Donor_Mix == "Y_mix" & temp$Subtype != i)), length(which(temp$Donor_Mix == "Y_unmix" & temp$Subtype != i))),
                         c(length(which(temp$Donor_Mix == "X_mix" & temp$Subtype == i)), length(which(temp$Donor_Mix == "X_unmix" & temp$Subtype == i)),
                           length(which(temp$Donor_Mix == "Y_mix" & temp$Subtype == i)), length(which(temp$Donor_Mix == "Y_unmix" & temp$Subtype == i)))))
  
  dimnames(xtab) <- list(Subtype = c(paste("Not ",i,sep=""), i),
                         Class = c("X_mix","X_unmix","Y_mix","Y_unmix"))
  x <- x + 1
  xtab_zheng_list[[x]] <- xtab
  names(xtab_zheng_list)[x] <- i
}

proptest_zheng_CD4T_data <- as.data.frame(matrix(0L, nrow=6*3, ncol=4))
colnames(proptest_zheng_CD4T_data) <- c("group1","group2","p.adj","Subtype")
for (i in 1:length(xtab_zheng_list)) { 
  temp <- as.data.frame(pairwise_prop_test(xtab_zheng_list[[i]]))
  proptest_zheng_CD4T_data[(1+(i-1)*6):(i*6), c("group1","group2","p.adj")] <- temp[,c("group1","group2","p.adj")]
  proptest_zheng_CD4T_data[(1+(i-1)*6):(i*6), "Subtype"] <- rep(names(xtab_zheng_list)[i], 6)
}
proptest_zheng_CD4T_data$p.adj <- -log10(proptest_zheng_CD4T_data$p.adj)
proptest_zheng_CD4T_data[,"Bin"] <- rep(0L)
proptest_zheng_CD4T_data$Bin[which(proptest_zheng_CD4T_data$p.adj > 2)] <- 1
proptest_zheng_CD4T_data$Bin[which(proptest_zheng_CD4T_data$p.adj >= 5)] <- 2
proptest_zheng_CD4T_data$Bin[which(proptest_zheng_CD4T_data$p.adj >= 10)] <- 3
proptest_zheng_CD4T_data$Bin[which(proptest_zheng_CD4T_data$p.adj >= 15)] <- 4
proptest_zheng_CD4T_data$Bin[which(proptest_zheng_CD4T_data$p.adj >= 20)] <- 5