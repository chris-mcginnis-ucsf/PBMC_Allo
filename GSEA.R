##################################################################
## GSEA for (un)mixed PBMCs, 8-donor and Zheng et al PBMC data ###
## Chris McGinnis, Gartner Lab, UCSF, 06/26/2020 #################
##################################################################

library(Seurat)
library(XML)

## Step 1: Compute signed p-values for all expressed genes via (un)mixed DEG analysis -----------------------------------------------------------------------------------------------------------------
# 8-donor PBMC
seu_pbmc_ficoll@meta.data[,"CellType_Donor_Mix"] <- paste(seu_pbmc_ficoll@meta.data$CellType, seu_pbmc_ficoll@meta.data$Donor_Mix, sep="_")
seu_pbmc_ficoll <- SetIdent(seu_pbmc_ficoll, value=seu_pbmc_ficoll@meta.data$CellType_Donor_Mix)

allo.markers_CD8T <- FindMarkers(seu_pbmc_ficoll, ident.1="CD8T_A_unmix", ident.2="CD8T_A_mix", test.use = "MAST", logfc.threshold = 0)
allo.markers_CD14Mono <- FindMarkers(seu_pbmc_ficoll, ident.1="CD14Mono_A_unmix", ident.2="CD14Mono_A_mix", test.use = "MAST", logfc.threshold = 0)
allo.markers_CD16Mono <- FindMarkers(seu_pbmc_ficoll, ident.1="CD16Mono_A_unmix", ident.2="CD16Mono_A_mix", test.use = "MAST", logfc.threshold = 0)
allo.markers_NK <- FindMarkers(seu_pbmc_ficoll, ident.1="NK_A_unmix", ident.2="NK_A_mix", test.use = "MAST", logfc.threshold = 0)
allo.markers_B <- FindMarkers(seu_pbmc_ficoll, ident.1="B_A_unmix", ident.2="B_A_mix", test.use = "MAST", logfc.threshold = 0)
allo.markers_DC <- FindMarkers(seu_pbmc_ficoll, ident.1="DC_A_unmix", ident.2="DC_A_mix", test.use = "MAST", logfc.threshold = 0)

seu_CD4T_ficoll@meta.data[,"Subtype_Donor_Mix"] <- paste(seu_CD4T_ficoll@meta.data$CD4TSubset, seu_CD4T_ficoll@meta.data$Donor_Mix, sep="_")
seu_CD4T_ficoll <- SetIdent(seu_CD4T_ficoll, value=seu_CD4T_ficoll@meta.data$Subtype_Donor_Mix)

allo.markers_CD4T_act <- FindMarkers(seu_CD4T_ficoll, ident.1="Activated_A_unmix", ident.2="Activated_A_mix", test.use = "MAST", logfc.threshold = 0)
allo.markers_CD4T_mem <- FindMarkers(seu_CD4T_ficoll, ident.1="Memory_A_unmix", ident.2="Memory_A_mix", test.use = "MAST", logfc.threshold = 0)
allo.markers_CD4T_naive <- FindMarkers(seu_CD4T_ficoll, ident.1="Naive_A_unmix", ident.2="Naive_A_mix", test.use = "MAST", logfc.threshold = 0)

write.csv(allo.markers_CD8T, file="allo.markers_CD8T.csv")
write.csv(allo.markers_CD14Mono, file="allo.markers_CD14Mono.csv")
write.csv(allo.markers_CD16Mono, file="allo.markers_CD16Mono.csv")
write.csv(allo.markers_NK, file="allo.markers_NK.csv")
write.csv(allo.markers_B, file="allo.markers_B.csv")
write.csv(allo.markers_DC, file="allo.markers_DC.csv")
write.csv(allo.markers_CD4T_act, file="allo.markers_CD4T_act.csv")
write.csv(allo.markers_CD4T_mem, file="allo.markers_CD4T_mem.csv")
write.csv(allo.markers_CD4T_naive, file="allo.markers_CD4T_naive.csv")

# Zheng et al PBMC
seu_zheng_clean <- SetIdent(seu_zheng_clean, value=seu_zheng_clean@meta.data$Donor_Mix_CellType)

allo.markers_zheng_X_CD8T <- FindMarkers(seu_zheng_clean, ident.1="X_CD8T_ISO", ident.2="X_CD8T_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_CD8T <- FindMarkers(seu_zheng_clean, ident.1="Y_CD8T_ISO", ident.2="Y_CD8T_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_CD14Mono <- FindMarkers(seu_zheng_clean, ident.1="X_CD14Mono_ISO", ident.2="X_CD14Mono_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_CD14Mono <- FindMarkers(seu_zheng_clean, ident.1="Y_CD14Mono_ISO", ident.2="Y_CD14Mono_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_CD16Mono <- FindMarkers(seu_zheng_clean, ident.1="X_CD16Mono_ISO", ident.2="X_CD16Mono_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_CD16Mono <- FindMarkers(seu_zheng_clean, ident.1="Y_CD16Mono_ISO", ident.2="Y_CD16Mono_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_NK <- FindMarkers(seu_zheng_clean, ident.1="X_NK_ISO", ident.2="X_NK_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_NK <- FindMarkers(seu_zheng_clean, ident.1="Y_NK_ISO", ident.2="Y_NK_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_DC <- FindMarkers(seu_zheng_clean, ident.1="X_DC_ISO", ident.2="X_DC_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_DC <- FindMarkers(seu_zheng_clean, ident.1="Y_DC_ISO", ident.2="Y_DC_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_B <- FindMarkers(seu_zheng_clean, ident.1="X_B_ISO", ident.2="X_B_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_B <- FindMarkers(seu_zheng_clean, ident.1="Y_B_ISO", ident.2="Y_B_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_Platelet <- FindMarkers(seu_zheng_clean, ident.1="X_Platelet_ISO", ident.2="X_Platelet_MIX", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_Platelet <- FindMarkers(seu_zheng_clean, ident.1="Y_Platelet_ISO", ident.2="Y_Platelet_MIX", test.use = "MAST", logfc.threshold = 0)

seu_zheng_CD4T <- SetIdent(seu_zheng_CD4T, value=seu_zheng_CD4T@meta.data$Donor_Mix_CellType)

allo.markers_zheng_X_CD4T_act <- FindMarkers(seu_zheng_CD4T, ident.1="X_unmix_activated", ident.2="X_mix_activated", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_CD4T_act <- FindMarkers(seu_zheng_CD4T, ident.1="Y_unmix_activated", ident.2="Y_mix_activated", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_CD4T_mem <- FindMarkers(seu_zheng_CD4T, ident.1="X_unmix_memory", ident.2="X_mix_memory", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_CD4T_mem <- FindMarkers(seu_zheng_CD4T, ident.1="Y_unmix_memory", ident.2="Y_mix_memory", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_X_CD4T_naive <- FindMarkers(seu_zheng_CD4T, ident.1="X_unmix_naive", ident.2="X_mix_naive", test.use = "MAST", logfc.threshold = 0)
allo.markers_zheng_Y_CD4T_naive <- FindMarkers(seu_zheng_CD4T, ident.1="Y_unmix_naive", ident.2="Y_mix_naive", test.use = "MAST", logfc.threshold = 0)

write.csv(allo.markers_zheng_X_CD8T, file="allo.markers_zheng_X_CD8T.csv")
write.csv(allo.markers_zheng_Y_CD8T, file="allo.markers_zheng_Y_CD8T.csv")
write.csv(allo.markers_zheng_X_CD14Mono, file="allo.markers_zheng_X_CD14Mono.csv")
write.csv(allo.markers_zheng_Y_CD14Mono, file="allo.markers_zheng_Y_CD14Mono.csv")
write.csv(allo.markers_zheng_X_CD16Mono, file="allo.markers_zheng_X_CD16Mono.csv")
write.csv(allo.markers_zheng_Y_CD16Mono, file="allo.markers_zheng_Y_CD16Mono.csv")
write.csv(allo.markers_zheng_X_NK, file="allo.markers_zheng_X_NK.csv")
write.csv(allo.markers_zheng_Y_NK, file="allo.markers_zheng_Y_NK.csv")
write.csv(allo.markers_zheng_X_DC, file="allo.markers_zheng_X_DC.csv")
write.csv(allo.markers_zheng_Y_DC, file="allo.markers_zheng_Y_DC.csv")
write.csv(allo.markers_zheng_X_B, file="allo.markers_zheng_X_B.csv")
write.csv(allo.markers_zheng_Y_B, file="allo.markers_zheng_Y_B.csv")
write.csv(allo.markers_zheng_X_Platelet, file="allo.markers_zheng_X_Platelet.csv")
write.csv(allo.markers_zheng_Y_Platelet, file="allo.markers_zheng_Y_Platelet.csv")
write.csv(allo.markers_zheng_X_CD4T_act, file="allo.markers_zheng_X_CD4T_act.csv")
write.csv(allo.markers_zheng_Y_CD4T_act, file="allo.markers_zheng_Y_CD4T_act.csv")
write.csv(allo.markers_zheng_X_CD4T_mem, file="allo.markers_zheng_X_CD4T_mem.csv")
write.csv(allo.markers_zheng_Y_CD4T_mem, file="allo.markers_zheng_Y_CD4T_mem.csv")
write.csv(allo.markers_zheng_X_CD4T_naive, file="allo.markers_zheng_X_CD4T_naive.csv")
write.csv(allo.markers_zheng_Y_CD4T_naive, file="allo.markers_zheng_Y_CD4T_naive.csv")


## Step 2: Read GSEA outputs into R, selecting signficiant gene sets only -----------------------------------------------------------------------------------------------------------------------------
for (i in list.files()[-1]) {
  temp <- readHTMLTable(i, header=F)[[1]][,c(2,7,8)]
  colnames(temp) <- c("GS","nomP","fdrQ")
  temp$nomP <- as.numeric(as.character(temp$nomP))
  temp$fdrQ <- as.numeric(as.character(temp$fdrQ))
  temp <- temp[which(temp$nomP <= 0.05 & temp$fdrQ <= 0.05), ]
  if (nrow(temp) == 0) { next }
  i <- gsub(x=i, pattern="Up", replacement = "MIX")
  i <- gsub(x=i, pattern="Down", replacement = "ISO")
  i <- gsub(x=i, pattern=".html", replacement = "")
  print(i)
  assign(x=i, value=temp, envir = .GlobalEnv)
}

## Step 3: Compile gene sets increased in iso conditions, 8-donor -------------------------------------------------------------------------------------------------------------------------------------
CD4T_Activated_Original_ISO <- cbind(rep("CD4T_act",nrow(CD4T_Activated_Original_ISO)), CD4T_Activated_Original_ISO)
colnames(CD4T_Activated_Original_ISO)[1] <- "CellType"
DC_Original_ISO <- cbind(rep("DC",nrow(DC_Original_ISO)), DC_Original_ISO)
colnames(DC_Original_ISO)[1] <- "CellType"
gsea_og_iso <- rbind(CD4T_Activated_Original_ISO, DC_Original_ISO)
## Result -- Activated CD4T-cells enriched in "HUMORAL_IMMUNE_RESPONSE", DCs enriched in cell killing, epigenetic regulation (e.g., lysine methylation)

## Step 3: Compile gene sets increased in mixed conditions, 8-donor -----------------------------------------------------------------------------------------------------------------------------------
## Result -- No enriched gene sets in MIX conditions
  
## Step 4: Compile donor-constant gene sets increased in iso conditions, Zheng ------------------------------------------------------------------------------------------------------------------------
## Result -- No enriched gene sets in ISO conditions shared between donors

## Step 5: Compile donor-constant gene sets increased in mixed conditions, Zheng ----------------------------------------------------------------------------------------------------------------------
B_X_Zheng_MIX <- cbind(rep("B",nrow(B_X_Zheng_MIX)), rep("X",nrow(B_X_Zheng_MIX)), B_X_Zheng_MIX)
colnames(B_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
B_Y_Zheng_MIX <- cbind(rep("B",nrow(B_Y_Zheng_MIX)), rep("Y",nrow(B_Y_Zheng_MIX)), B_Y_Zheng_MIX)
colnames(B_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD16Mono_X_Zheng_MIX <- cbind(rep("CD16Mono",nrow(CD16Mono_X_Zheng_MIX)), rep("X",nrow(CD16Mono_X_Zheng_MIX)), CD16Mono_X_Zheng_MIX)
colnames(CD16Mono_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD16Mono_Y_Zheng_MIX <- cbind(rep("CD16Mono",nrow(CD16Mono_Y_Zheng_MIX)), rep("Y",nrow(CD16Mono_Y_Zheng_MIX)), CD16Mono_Y_Zheng_MIX)
colnames(CD16Mono_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD4T_Activated_X_Zheng_MIX <- cbind(rep("CD4T_act",nrow(CD4T_Activated_X_Zheng_MIX)), rep("X",nrow(CD4T_Activated_X_Zheng_MIX)), CD4T_Activated_X_Zheng_MIX)
colnames(CD4T_Activated_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD4T_Activated_Y_Zheng_MIX <- cbind(rep("CD4T_act",nrow(CD4T_Activated_Y_Zheng_MIX)), rep("Y",nrow(CD4T_Activated_Y_Zheng_MIX)), CD4T_Activated_Y_Zheng_MIX)
colnames(CD4T_Activated_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD4T_Memory_X_Zheng_MIX <- cbind(rep("CD4T_mem",nrow(CD4T_Memory_X_Zheng_MIX)), rep("X",nrow(CD4T_Memory_X_Zheng_MIX)), CD4T_Memory_X_Zheng_MIX)
colnames(CD4T_Memory_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD4T_Memory_Y_Zheng_MIX <- cbind(rep("CD4T_mem",nrow(CD4T_Memory_Y_Zheng_MIX)), rep("Y",nrow(CD4T_Memory_Y_Zheng_MIX)), CD4T_Memory_Y_Zheng_MIX)
colnames(CD4T_Memory_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD4T_Naive_X_Zheng_MIX <- cbind(rep("CD4T_naive",nrow(CD4T_Naive_X_Zheng_MIX)), rep("X",nrow(CD4T_Naive_X_Zheng_MIX)), CD4T_Naive_X_Zheng_MIX)
colnames(CD4T_Naive_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD4T_Naive_Y_Zheng_MIX <- cbind(rep("CD4T_naive",nrow(CD4T_Naive_Y_Zheng_MIX)), rep("Y",nrow(CD4T_Naive_Y_Zheng_MIX)), CD4T_Naive_Y_Zheng_MIX)
colnames(CD4T_Naive_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD8T_X_Zheng_MIX <- cbind(rep("CD8T",nrow(CD8T_X_Zheng_MIX)), rep("X",nrow(CD8T_X_Zheng_MIX)), CD8T_X_Zheng_MIX)
colnames(CD8T_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
CD8T_Y_Zheng_MIX <- cbind(rep("CD8T",nrow(CD8T_Y_Zheng_MIX)), rep("Y",nrow(CD8T_Y_Zheng_MIX)), CD8T_Y_Zheng_MIX)
colnames(CD8T_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
DC_X_Zheng_MIX <- cbind(rep("DC",nrow(DC_X_Zheng_MIX)), rep("X",nrow(DC_X_Zheng_MIX)), DC_X_Zheng_MIX)
colnames(DC_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
DC_Y_Zheng_MIX <- cbind(rep("DC",nrow(DC_Y_Zheng_MIX)), rep("Y",nrow(DC_Y_Zheng_MIX)), DC_Y_Zheng_MIX)
colnames(DC_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
NK_X_Zheng_MIX <- cbind(rep("NK",nrow(NK_X_Zheng_MIX)), rep("X",nrow(NK_X_Zheng_MIX)), NK_X_Zheng_MIX)
colnames(NK_X_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")
NK_Y_Zheng_MIX <- cbind(rep("NK",nrow(NK_Y_Zheng_MIX)), rep("Y",nrow(NK_Y_Zheng_MIX)), NK_Y_Zheng_MIX)
colnames(NK_Y_Zheng_MIX)[c(1,2)] <- c("CellType","Donor")

# Step 6: Select only genes consistent between donors -------------------------------------------------------------------------------------------------------------------------------------------------
B_X_Zheng_MIX <- B_X_Zheng_MIX[which(B_X_Zheng_MIX$GS %in% B_Y_Zheng_MIX$GS), ]
B_Y_Zheng_MIX <- B_Y_Zheng_MIX[which(B_Y_Zheng_MIX$GS %in% B_X_Zheng_MIX$GS), ]
CD16Mono_X_Zheng_MIX <- CD16Mono_X_Zheng_MIX[which(CD16Mono_X_Zheng_MIX$GS %in% CD16Mono_Y_Zheng_MIX$GS), ]
CD16Mono_Y_Zheng_MIX <- CD16Mono_Y_Zheng_MIX[which(CD16Mono_Y_Zheng_MIX$GS %in% CD16Mono_X_Zheng_MIX$GS), ]
CD4T_Activated_X_Zheng_MIX <- CD4T_Activated_X_Zheng_MIX[which(CD4T_Activated_X_Zheng_MIX$GS %in% CD4T_Activated_Y_Zheng_MIX$GS), ]
CD4T_Activated_Y_Zheng_MIX <- CD4T_Activated_Y_Zheng_MIX[which(CD4T_Activated_Y_Zheng_MIX$GS %in% CD4T_Activated_X_Zheng_MIX$GS), ]
CD4T_Memory_X_Zheng_MIX <- CD4T_Memory_X_Zheng_MIX[which(CD4T_Memory_X_Zheng_MIX$GS %in% CD4T_Memory_Y_Zheng_MIX$GS), ]
CD4T_Memory_Y_Zheng_MIX <- CD4T_Memory_Y_Zheng_MIX[which(CD4T_Memory_Y_Zheng_MIX$GS %in% CD4T_Memory_X_Zheng_MIX$GS), ]
CD4T_Naive_X_Zheng_MIX <- CD4T_Naive_X_Zheng_MIX[which(CD4T_Naive_X_Zheng_MIX$GS %in% CD4T_Naive_Y_Zheng_MIX$GS), ]
CD4T_Naive_Y_Zheng_MIX <- CD4T_Naive_Y_Zheng_MIX[which(CD4T_Naive_Y_Zheng_MIX$GS %in% CD4T_Naive_X_Zheng_MIX$GS), ]
CD8T_X_Zheng_MIX <- CD8T_X_Zheng_MIX[which(CD8T_X_Zheng_MIX$GS %in% CD8T_Y_Zheng_MIX$GS), ]
CD8T_Y_Zheng_MIX <- CD8T_Y_Zheng_MIX[which(CD8T_Y_Zheng_MIX$GS %in% CD8T_X_Zheng_MIX$GS), ]
DC_X_Zheng_MIX <- DC_X_Zheng_MIX[which(DC_X_Zheng_MIX$GS %in% DC_Y_Zheng_MIX$GS), ]
DC_Y_Zheng_MIX <- DC_Y_Zheng_MIX[which(DC_Y_Zheng_MIX$GS %in% DC_X_Zheng_MIX$GS), ]
NK_X_Zheng_MIX <- NK_X_Zheng_MIX[which(NK_X_Zheng_MIX$GS %in% NK_Y_Zheng_MIX$GS), ]
NK_Y_Zheng_MIX <- NK_Y_Zheng_MIX[which(NK_Y_Zheng_MIX$GS %in% NK_X_Zheng_MIX$GS), ]

gsea_zheng_mix <- rbind(B_X_Zheng_MIX, B_Y_Zheng_MIX,
                        CD16Mono_X_Zheng_MIX, CD16Mono_Y_Zheng_MIX,
                        CD4T_Activated_X_Zheng_MIX, CD4T_Activated_Y_Zheng_MIX,
                        CD4T_Memory_X_Zheng_MIX, CD4T_Memory_Y_Zheng_MIX,
                        CD4T_Naive_X_Zheng_MIX, CD4T_Naive_Y_Zheng_MIX,
                        CD8T_X_Zheng_MIX, CD8T_Y_Zheng_MIX,
                        DC_X_Zheng_MIX, DC_Y_Zheng_MIX,
                        NK_X_Zheng_MIX, NK_Y_Zheng_MIX)

gsea_zheng_mix <- gsea_zheng_mix[which(gsea_zheng_mix$Donor == "X"), ]
gsea_zheng_mix <- gsea_zheng_mix[,-2]
## Result -- B-cells enriched in protein trafficking, translation, nonsense mediated decay, and viral gene expression
## CD16Mono enriched in protein trafficking, translation, nonsense mediated decay, and viral gene expression
## Activated CD4Ts enriched in protein trafficking, translation, nonsense mediated decay, viral gene expression, amino acid metabolism and "INTERSPECIES_INTERACTION_BETWEEN_ORGANISMS"
## Memory CD4Ts enriched in protein trafficking, translation, nonsense mediated decay, viral gene expression, and amino acid metabolism 
## Naive CD4Ts enriched in protein trafficking, translation, nonsense mediated decay, viral gene expression, amino acid metabolism and "INTERSPECIES_INTERACTION_BETWEEN_ORGANISMS"
## CD8Ts enriched in protein trafficking
## DCs enriched in protein trafficking and nonsense mediated decay
## NKs enriched in protein trafficking, translation, nonsense mediated decay, and viral gene expression