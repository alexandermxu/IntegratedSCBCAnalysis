library(Seurat)
library(dplyr)
library(viridis)
library(reshape2)
library(extrafont)
setwd("C:/Users/alexm/Documents/git/Protein Analysis/")
mingeneappearancethreshold <- 5
lowUMIpercellthreshold <- 500
lowgenepercellthreshold <- 100

# load("AX206genes")
# AX206 <- Genes
# load("AX207genes")
# AX207 <- Genes
# load("AX208genes")
# AX208 <- Genes
# load("AX206Redogenes")
# AX206Redo <- Genes
# load("AX208Redogenes")
# AX208Redo <- Genes
# load("AX218genes")
# AX218 <- Genes
# load("AX219genes")
# AX219 <- Genes

# colnames(IntegratedData) <- gsub("X", "AX219X", colnames(IntegratedData))
# colnames(NoCellIntegratedData) <- gsub("X", "AX219X", colnames(NoCellIntegratedData))
# save(list = c("IntegratedData", "NoCellIntegratedData"), file = "AX219alldata")

print("Loading data")
load("AX206alldata")
AX206all <- IntegratedData
AX206NoCell <- NoCellIntegratedData
load("AX207alldata")
AX207all <- IntegratedData
AX207NoCell <- NoCellIntegratedData
load("AX208alldata")
AX208all <- IntegratedData
AX208NoCell <- NoCellIntegratedData
load("AX206Redoalldata")
AX206Redoall <- IntegratedData
AX206RedoNoCell <- NoCellIntegratedData
load("AX208Redoalldata")
AX208Redoall <- IntegratedData
AX208RedoNoCell <- NoCellIntegratedData
load("AX218alldata")
AX218all <- IntegratedData
AX218NoCell <- NoCellIntegratedData
load("AX219alldata")
AX219all <- IntegratedData
AX219NoCell <- NoCellIntegratedData
# save(list=c("AX206Vals","AX207Vals","AX208Vals","AX218Vals","AX219Vals","AX206Zeros","AX207Zeros","AX208Zeros","AX218Zeros","AX219Zeros","AX206RedoVals","AX208RedoVals","AX206RedoZeros","AX208RedoZeros"), file = "AllProteinValues")
load("AllProteinValues")
# AX206Vals <- data.frame(t(ProteinsPerBeads))
# rownames(AX219Vals) <- gsub("X","AX219X",rownames(AX219Vals))
# AX219Zeros <- AX219Vals[which(AX219Vals[,4]==0),]
print("Applying normalization and background subtraction")
AX206Background <- apply(AX206Zeros,2,mean)[1:3]
AX207Background <- apply(AX207Zeros,2,mean)[1:3]
AX208Background <- apply(AX208Zeros,2,mean)[1:3]
AX218Background <- apply(AX218Zeros,2,mean)[1:3]
AX219Background <- apply(AX219Zeros,2,mean)[1:3]
AX206RedoBackground <- apply(AX206RedoZeros,2,mean)[1:3]
AX208RedoBackground <- apply(AX208RedoZeros,2,mean)[1:3]
AX206ConversionFactors <- AX206Background/100
AX207ConversionFactors <- AX207Background/100
AX208ConversionFactors <- AX208Background/100
AX218ConversionFactors <- AX218Background/100
AX219ConversionFactors <- AX219Background/100
AX206RedoConversionFactors <- AX206RedoBackground/100
AX208RedoConversionFactors <- AX208RedoBackground/100

AX206NormalizedProteins <- AX206Vals
AX206NormalizedProteins[,1:3] <- t(apply(AX206Vals[,1:3],1,function(x) (x-AX206Background)/AX206ConversionFactors))
AX207NormalizedProteins <- AX207Vals
AX207NormalizedProteins[,1:3] <- t(apply(AX207Vals[,1:3],1,function(x) (x-AX207Background)/AX207ConversionFactors))
AX208NormalizedProteins <- AX208Vals
AX208NormalizedProteins[,1:3] <- t(apply(AX208Vals[,1:3],1,function(x) (x-AX208Background)/AX208ConversionFactors))
AX218NormalizedProteins <- AX218Vals
AX218NormalizedProteins[,1:3] <- t(apply(AX218Vals[,1:3],1,function(x) (x-AX218Background)/AX218ConversionFactors))
AX219NormalizedProteins <- AX219Vals
AX219NormalizedProteins[,1:3] <- t(apply(AX219Vals[,1:3],1,function(x) (x-AX219Background)/AX219ConversionFactors))
AX206RedoNormalizedProteins <- AX206RedoVals
AX206RedoNormalizedProteins[,1:3] <- t(apply(AX206RedoVals[,1:3],1,function(x) (x-AX206RedoBackground)/AX206RedoConversionFactors))
AX208RedoNormalizedProteins <- AX208RedoVals
AX208RedoNormalizedProteins[,1:3] <- t(apply(AX208RedoVals[,1:3],1,function(x) (x-AX208RedoBackground)/AX208RedoConversionFactors))
# AX206NormalizedProteins[,"Chip"] <- "AX206"
# AX207NormalizedProteins[,"Chip"] <- "AX207"
# AX208NormalizedProteins[,"Chip"] <- "AX208"
# AX218NormalizedProteins[,"Chip"] <- "AX218"
# AX219NormalizedProteins[,"Chip"] <- "AX219"

AX206 <- AX206all[-((nrow(AX206all)-11):nrow(AX206all)),]
AX207 <- AX207all[-((nrow(AX207all)-11):nrow(AX207all)),]
AX208 <- AX208all[-((nrow(AX208all)-11):nrow(AX208all)),]
AX206Redo <- AX206Redoall[-((nrow(AX206Redoall)-11):nrow(AX206Redoall)),]
AX208Redo <- AX208Redoall[-((nrow(AX208Redoall)-11):nrow(AX208Redoall)),]
AX218 <- AX218all[-((nrow(AX218all)-11):nrow(AX218all)),]
AX219 <- AX219all[-((nrow(AX219all)-11):nrow(AX219all)),]

print("Creating Seurat objects")
AX206S <- CreateSeuratObject(raw.data=AX206, project="AX206", min.cells=mingeneappearancethreshold)
AX206S@meta.data$celltype <- "U87"
AX206S <- FilterCells(AX206S, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
AX206S <- NormalizeData(AX206S, display.progress=F)
AX206S <- ScaleData(AX206S, display.progress=F)
AX206S <- FindVariableGenes(AX206S, do.plot = F, display.progress=F)
# AX206S <- SetAssayData(AX206S, assay.type = "SCBC", slot = "raw.data", new.data = AX206all[((nrow(AX206all)-3):(nrow(AX206all)-1)),])
# AX206S <- NormalizeData(AX206S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# AX206S <- ScaleData(AX206S, assay.type = "SCBC", display.progress = F)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AX206S@data), value = TRUE)
percent.mito <- Matrix::colSums(AX206S@raw.data[mito.genes, ])/Matrix::colSums(AX206S@raw.data)
AX206S <- AddMetaData(object = AX206S, metadata = percent.mito, col.name = "percent.mito")

AX207S <- CreateSeuratObject(raw.data=AX207, project="AX207", min.cells=mingeneappearancethreshold)
AX207S@meta.data$celltype <- "HEK"
AX207S <- FilterCells(AX207S, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
AX207S <- NormalizeData(AX207S, display.progress=F)
AX207S <- ScaleData(AX207S, display.progress=F)
AX207S <- FindVariableGenes(AX207S, do.plot = F, display.progress=F)
# AX207S <- SetAssayData(AX207S, assay.type = "SCBC", slot = "raw.data", new.data = AX207all[((nrow(AX207all)-3):(nrow(AX207all)-1)),])
# AX207S <- NormalizeData(AX207S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# AX207S <- ScaleData(AX207S, assay.type = "SCBC", display.progress = F)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AX207S@data), value = TRUE)
percent.mito <- Matrix::colSums(AX207S@raw.data[mito.genes, ])/Matrix::colSums(AX207S@raw.data)
AX207S <- AddMetaData(object = AX207S, metadata = percent.mito, col.name = "percent.mito")

AX208S <- CreateSeuratObject(raw.data=AX208, project="AX208", min.cells=mingeneappearancethreshold)
AX208S@meta.data$celltype <- "HEK"
AX208S <- FilterCells(AX208S, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
AX208S <- NormalizeData(AX208S, display.progress=F)
AX208S <- ScaleData(AX208S, display.progress=F)
AX208S <- FindVariableGenes(AX208S, do.plot = F, display.progress=F)
# AX208S <- SetAssayData(AX208S, assay.type = "SCBC", slot = "raw.data", new.data = AX208all[((nrow(AX208all)-3):(nrow(AX208all)-1)),])
# AX208S <- NormalizeData(AX208S, assay.type = "SCBC", normalization.method = "genesCLR")
# AX208S <- ScaleData(AX208S, assay.type = "SCBC", display.progress = F)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AX208S@data), value = TRUE)
percent.mito <- Matrix::colSums(AX208S@raw.data[mito.genes, ])/Matrix::colSums(AX208S@raw.data)
AX208S <- AddMetaData(object = AX208S, metadata = percent.mito, col.name = "percent.mito")

AX218S <- CreateSeuratObject(raw.data=AX218, project="AX218", min.cells=mingeneappearancethreshold)
AX218S@meta.data$celltype <- "U87"
AX218S <- FilterCells(AX218S, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
AX218S <- NormalizeData(AX218S, display.progress=F)
AX218S <- ScaleData(AX218S, display.progress=F)
AX218S <- FindVariableGenes(AX218S, do.plot = F, display.progress=F)
# AX218S <- SetAssayData(AX218S, assay.type = "SCBC", slot = "raw.data", new.data = AX218all[((nrow(AX218all)-3):(nrow(AX218all)-1)),])
# AX218S <- NormalizeData(AX218S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# AX218S <- ScaleData(AX218S, assay.type = "SCBC", display.progress = F)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AX218S@data), value = TRUE)
percent.mito <- Matrix::colSums(AX218S@raw.data[mito.genes, ])/Matrix::colSums(AX218S@raw.data)
AX218S <- AddMetaData(object = AX218S, metadata = percent.mito, col.name = "percent.mito")

AX219S <- CreateSeuratObject(raw.data=AX219, project="AX219", min.cells=mingeneappearancethreshold)
AX219S@meta.data$celltype <- "U87"
AX219S <- FilterCells(AX219S, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
AX219S <- NormalizeData(AX219S, display.progress=F)
AX219S <- ScaleData(AX219S, display.progress=F)
AX219S <- FindVariableGenes(AX219S, do.plot = F, display.progress=F)
# AX219S <- SetAssayData(AX219S, assay.type = "SCBC", slot = "raw.data", new.data = AX219all[((nrow(AX219all)-3):(nrow(AX219all)-1)),])
# AX219S <- NormalizeData(AX219S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# AX219S <- ScaleData(AX219S, assay.type = "SCBC", display.progress = F)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AX219S@data), value = TRUE)
percent.mito <- Matrix::colSums(AX219S@raw.data[mito.genes, ])/Matrix::colSums(AX219S@raw.data)
AX219S <- AddMetaData(object = AX219S, metadata = percent.mito, col.name = "percent.mito")

AX206RedoS <- CreateSeuratObject(raw.data=AX206Redo, project="AX206Redo", min.cells=mingeneappearancethreshold)
AX206RedoS@meta.data$celltype <- "U87"
AX206RedoS <- FilterCells(AX206RedoS, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
AX206RedoS <- NormalizeData(AX206RedoS, display.progress=F)
AX206RedoS <- ScaleData(AX206RedoS, display.progress=F)
AX206RedoS <- FindVariableGenes(AX206RedoS, do.plot = F, display.progress=F)
# AX206RedoS <- SetAssayData(AX206RedoS, assay.type = "SCBC", slot = "raw.data", new.data = AX206Redoall[((nrow(AX206Redoall)-3):(nrow(AX206Redoall)-1)),])
# AX206RedoS <- NormalizeData(AX206RedoS, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# AX206RedoS <- ScaleData(AX206RedoS, assay.type = "SCBC", display.progress = F)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AX206RedoS@data), value = TRUE)
percent.mito <- Matrix::colSums(AX206RedoS@raw.data[mito.genes, ])/Matrix::colSums(AX206RedoS@raw.data)
AX206RedoS <- AddMetaData(object = AX206RedoS, metadata = percent.mito, col.name = "percent.mito")

AX208RedoS <- CreateSeuratObject(raw.data=AX208Redo, project="AX208Redo", min.cells=mingeneappearancethreshold)
AX208RedoS@meta.data$celltype <- "HEK"
AX208RedoS <- FilterCells(AX208RedoS, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
AX208RedoS <- NormalizeData(AX208RedoS, display.progress=F)
AX208RedoS <- ScaleData(AX208RedoS, display.progress=F)
AX208RedoS <- FindVariableGenes(AX208RedoS, do.plot = F, display.progress=F)
# AX208RedoS <- SetAssayData(AX208RedoS, assay.type = "SCBC", slot = "raw.data", new.data = AX208Redoall[((nrow(AX208Redoall)-3):(nrow(AX208Redoall)-1)),])
# AX208RedoS <- NormalizeData(AX208RedoS, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# AX208RedoS <- ScaleData(AX208RedoS, assay.type = "SCBC", display.progress = F)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = AX208RedoS@data), value = TRUE)
percent.mito <- Matrix::colSums(AX208RedoS@raw.data[mito.genes, ])/Matrix::colSums(AX208RedoS@raw.data)
AX208RedoS <- AddMetaData(object = AX208RedoS, metadata = percent.mito, col.name = "percent.mito")

# U871 <- read.csv("GSM2794663_U87_con_1_Genes_ReadCount.txt", sep = "\t", row.names = 1)
# colnames(U871) <- "U87Control1"
# U872 <- read.csv("GSM2794664_U87_con_2_Genes_ReadCount.txt", sep = "\t", row.names = 1)
# colnames(U872) <- "U87Control2"
# 
# HEKCombinedSingleCell <- CombinedGenesbyMerge@raw.data[,CombinedGenesbyMerge@meta.data$celltype=="HEK"]
# U87CombinedSingleCell <- CombinedGenesbyMerge@raw.data[,CombinedGenesbyMerge@meta.data$celltype=="U87"]
# U87CombinedSingleCell <- apply(U87CombinedSingleCell,1,mean)
# HEKCombinedSingleCell <- apply(HEKCombinedSingleCell,1,mean)
# 
# BulkComp <- data.frame(cbind(U87CombinedSingleCell, HEKCombinedSingleCell))

# CombinedS <- CreateSeuratObject(raw.data=BulkComp, project="CombinedCells")
# CombinedS@meta.data$celltype <- "U87"
# CombinedS <- FilterCells(CombinedS, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
# CombinedS <- NormalizeData(CombinedS, display.progress=F)
# CombinedS <- ScaleData(CombinedS, display.progress=F)
# CombinedS <- FindVariableGenes(CombinedS, do.plot = F, display.progress=F)
# # AX208RedoS <- SetAssayData(AX208RedoS, assay.type = "SCBC", slot = "raw.data", new.data = AX208Redoall[((nrow(AX208Redoall)-3):(nrow(AX208Redoall)-1)),])
# # AX208RedoS <- NormalizeData(AX208RedoS, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# # AX208RedoS <- ScaleData(AX208RedoS, assay.type = "SCBC", display.progress = F)
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = CombinedS@data), value = TRUE)
# percent.mito <- Matrix::colSums(CombinedS@raw.data[mito.genes, ])/Matrix::colSums(CombinedS@raw.data)
# CombinedS <- AddMetaData(object = CombinedS, metadata = percent.mito, col.name = "percent.mito")

# # GSM2794664
# U87S <- CreateSeuratObject(raw.data=U87BulkControls, project="U87Control1")
# U87S@meta.data$celltype <- "U87"
# U87S <- FilterCells(U87S, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
# U87S <- NormalizeData(U87S, display.progress=F)
# U87S <- ScaleData(U87S, display.progress=F)
# U87S <- FindVariableGenes(U87S, do.plot = F, display.progress=F)
# # AX208RedoS <- SetAssayData(AX208RedoS, assay.type = "SCBC", slot = "raw.data", new.data = AX208Redoall[((nrow(AX208Redoall)-3):(nrow(AX208Redoall)-1)),])
# # AX208RedoS <- NormalizeData(AX208RedoS, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# # AX208RedoS <- ScaleData(AX208RedoS, assay.type = "SCBC", display.progress = F)
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = U87S@data), value = TRUE)
# percent.mito <- Matrix::colSums(U87S@raw.data[mito.genes, ])/Matrix::colSums(U87S@raw.data)
# U87S <- AddMetaData(object = U87S, metadata = percent.mito, col.name = "percent.mito")
# 
# # GSM2599702
# HEKS <- CreateSeuratObject(raw.data=UMI_count, project="HEK")
# HEKS@meta.data$celltype <- "HEK"
# HEKS <- FilterCells(HEKS, subset.names=c("nUMI","nGene"), low.thresholds=c(lowUMIpercellthreshold,lowgenepercellthreshold), high.thresholds=c(Inf,Inf))
# HEKS <- NormalizeData(HEKS, display.progress=F)
# HEKS <- ScaleData(HEKS, display.progress=F)
# HEKS <- FindVariableGenes(HEKS, do.plot = F, display.progress=F)
# # AX208RedoS <- SetAssayData(AX208RedoS, assay.type = "SCBC", slot = "raw.data", new.data = AX208Redoall[((nrow(AX208Redoall)-3):(nrow(AX208Redoall)-1)),])
# # AX208RedoS <- NormalizeData(AX208RedoS, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
# # AX208RedoS <- ScaleData(AX208RedoS, assay.type = "SCBC", display.progress = F)
# mito.genes <- grep(pattern = "^MT-", x = rownames(x = HEKS@data), value = TRUE)
# percent.mito <- Matrix::colSums(HEKS@raw.data[mito.genes, ])/Matrix::colSums(HEKS@raw.data)
# HEKS <- AddMetaData(object = HEKS, metadata = percent.mito, col.name = "percent.mito")

print("Adding protein values to Seurat")
AX206NormalizedProteins <- AX206NormalizedProteins[rownames(AX206NormalizedProteins) %in% AX206S@cell.names,]
AX206AllProts <- AX206NormalizedProteins
AX206NormalizedProteins[,1:3] <- AX206NormalizedProteins[,1:3]/AX206NormalizedProteins[,4]

AX207NormalizedProteins <- AX207NormalizedProteins[rownames(AX207NormalizedProteins) %in% AX207S@cell.names,]
AX207AllProts <- AX207NormalizedProteins
AX207NormalizedProteins[,1:3] <- AX207NormalizedProteins[,1:3]/AX207NormalizedProteins[,4]

AX208NormalizedProteins <- AX208NormalizedProteins[rownames(AX208NormalizedProteins) %in% AX208S@cell.names,]
AX208AllProts <- AX208NormalizedProteins
AX208NormalizedProteins[,1:3] <- AX208NormalizedProteins[,1:3]/AX208NormalizedProteins[,4]

AX218NormalizedProteins <- AX218NormalizedProteins[rownames(AX218NormalizedProteins) %in% AX218S@cell.names,]
AX218AllProts <- AX218NormalizedProteins
AX218NormalizedProteins[,1:3] <- AX218NormalizedProteins[,1:3]/AX218NormalizedProteins[,4]

AX219NormalizedProteins <- AX219NormalizedProteins[rownames(AX219NormalizedProteins) %in% AX219S@cell.names,]
AX219AllProts <- AX219NormalizedProteins
AX219NormalizedProteins[,1:3] <- AX219NormalizedProteins[,1:3]/AX219NormalizedProteins[,4]

AX206RedoNormalizedProteins <- AX206RedoNormalizedProteins[rownames(AX206RedoNormalizedProteins) %in% AX206RedoS@cell.names,]
AX206RedoAllProts <- AX206RedoNormalizedProteins
AX206RedoNormalizedProteins[,1:3] <- AX206RedoNormalizedProteins[,1:3]/AX206RedoNormalizedProteins[,4]

AX208RedoNormalizedProteins <- AX208RedoNormalizedProteins[rownames(AX208RedoNormalizedProteins) %in% AX208RedoS@cell.names,]
AX208RedoAllProts <- AX208RedoNormalizedProteins
AX208RedoNormalizedProteins[,1:3] <- AX208RedoNormalizedProteins[,1:3]/AX208RedoNormalizedProteins[,4]

AX206S <- SetAssayData(AX206S, assay.type = "SCBC", slot = "raw.data", new.data = t(AX206NormalizedProteins[,1:3]))
AX206S <- AddMetaData(object = AX206S, metadata = AX206NormalizedProteins[,4:5], col.name = c("cells","beads"))
AX206S <- NormalizeData(AX206S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
AX206S <- ScaleData(AX206S, assay.type = "SCBC", display.progress = F)
AX207S <- SetAssayData(AX207S, assay.type = "SCBC", slot = "raw.data", new.data = t(AX207NormalizedProteins[,1:3]))
AX207S <- AddMetaData(object = AX207S, metadata = AX207NormalizedProteins[,4:5], col.name = c("cells","beads"))
AX207S <- NormalizeData(AX207S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
AX207S <- ScaleData(AX207S, assay.type = "SCBC", display.progress = F)
AX208S <- SetAssayData(AX208S, assay.type = "SCBC", slot = "raw.data", new.data = t(AX208NormalizedProteins[,1:3]))
AX208S <- AddMetaData(object = AX208S, metadata = AX208NormalizedProteins[,4:5], col.name = c("cells","beads"))
AX208S <- NormalizeData(AX208S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
AX208S <- ScaleData(AX208S, assay.type = "SCBC", display.progress = F)
AX218S <- SetAssayData(AX218S, assay.type = "SCBC", slot = "raw.data", new.data = t(AX218NormalizedProteins[,1:3]))
AX218S <- AddMetaData(object = AX218S, metadata = AX218NormalizedProteins[,4:5], col.name = c("cells","beads"))
AX218S <- NormalizeData(AX218S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
AX218S <- ScaleData(AX218S, assay.type = "SCBC", display.progress = F)
AX219S <- SetAssayData(AX219S, assay.type = "SCBC", slot = "raw.data", new.data = t(AX219NormalizedProteins[,1:3]))
AX219S <- AddMetaData(object = AX219S, metadata = AX219NormalizedProteins[,4:5], col.name = c("cells","beads"))
AX219S <- NormalizeData(AX219S, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
AX219S <- ScaleData(AX219S, assay.type = "SCBC", display.progress = F)
AX206RedoS <- SetAssayData(AX206RedoS, assay.type = "SCBC", slot = "raw.data", new.data = t(AX206RedoNormalizedProteins[,1:3]))
AX206RedoS <- AddMetaData(object = AX206RedoS, metadata = AX206RedoNormalizedProteins[,4:5], col.name = c("cells","beads"))
AX206RedoS <- NormalizeData(AX206RedoS, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
AX206RedoS <- ScaleData(AX206RedoS, assay.type = "SCBC", display.progress = F)
AX208RedoS <- SetAssayData(AX208RedoS, assay.type = "SCBC", slot = "raw.data", new.data = t(AX208RedoNormalizedProteins[,1:3]))
AX208RedoS <- AddMetaData(object = AX208RedoS, metadata = AX208RedoNormalizedProteins[,4:5], col.name = c("cells","beads"))
AX208RedoS <- NormalizeData(AX208RedoS, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
AX208RedoS <- ScaleData(AX208RedoS, assay.type = "SCBC", display.progress = F)

# AX206SGeneNames <- head(rownames(AX206S@hvg.info), 1000)
# AX207SGeneNames <- head(rownames(AX207S@hvg.info), 1000)
# AX208SGeneNames <- head(rownames(AX208S@hvg.info), 1000)
# AX218SGeneNames <- head(rownames(AX218S@hvg.info), 1000)
# AX219SGeneNames <- head(rownames(AX219S@hvg.info), 1000)
# AX206RedoSGeneNames <- head(rownames(AX206RedoS@hvg.info), 1000)
# AX208RedoSGeneNames <- head(rownames(AX208RedoS@hvg.info), 1000)
print("Integrating multiple chips")
AX206NormalizedProteins[,"Chip"] <- "AX206"
AX207NormalizedProteins[,"Chip"] <- "AX207"
AX208NormalizedProteins[,"Chip"] <- "AX208"
AX218NormalizedProteins[,"Chip"] <- "AX218"
AX219NormalizedProteins[,"Chip"] <- "AX219"
AX206RedoNormalizedProteins[,"Chip"] <- "AX206Redo"
AX208RedoNormalizedProteins[,"Chip"] <- "AX208Redo"
Allprotsnormalized <- rbind(AX206NormalizedProteins,AX207NormalizedProteins,AX208NormalizedProteins,AX218NormalizedProteins,AX219NormalizedProteins, AX206RedoNormalizedProteins, AX208RedoNormalizedProteins)
Allprotsall <- rbind(AX206AllProts,AX207AllProts,AX208AllProts,AX218AllProts,AX219AllProts, AX206RedoAllProts, AX208RedoAllProts)
Allprotsall["Chip"] <- gsub(pattern = "*X(.*)", replacement="", x=gsub(pattern = "AX",replacement="A", x=rownames(Allprotsall)))
colnames(Allprotsall)[1:3] <- c("PKM2","c-MYC","PDHK1")
AllprotsallPlot <- melt(Allprotsall, id=c("Cells","Beads", "Chip"))
AllprotsallPlot[,1] <- as.factor(AllprotsallPlot[,1])
Allprotsnormalizedplot <- melt(Allprotsnormalized, id=c("Cells", "Beads", "Chip"))
Allprotsnormalizedplot$Chip <- as.factor(Allprotsnormalizedplot$Chip)
Allprotsnormalizedplot$Cells <- as.factor(Allprotsnormalizedplot$Cells)
Allprotsnormalizedplot["Celltype"] <- NA
Allprotsnormalizedplot$Celltype[Allprotsnormalizedplot$Chip %in% c("AX206", "AX206Redo", "AX218", "AX219")] <- "U87"
Allprotsnormalizedplot$Celltype[Allprotsnormalizedplot$Chip %in% c("AX208", "AX208Redo", "AX207")] <- "HEK"

ggplot(Allprotsnormalizedplot) + geom_boxplot(aes(x=variable, y=value, fill=Chip), outlier.shape = 3)+geom_point(aes(x=variable, y=value, fill=Chip, size=Cells), position=position_dodge(width = 0.75), alpha=0.5)+scale_size_discrete(range = c(1,5))

print("Choosing variable genes")
AX206SGeneNames <- AX206S@var.genes
AX207SGeneNames <- AX207S@var.genes
AX208SGeneNames <- AX208S@var.genes
AX218SGeneNames <- AX218S@var.genes
AX219SGeneNames <- AX219S@var.genes
AX206RedoSGeneNames <- AX206RedoS@var.genes
AX208RedoSGeneNames <- AX208RedoS@var.genes

GenestoUse <- unique(c(AX206SGeneNames, AX207SGeneNames, AX208SGeneNames, AX206RedoSGeneNames, AX208RedoSGeneNames, AX218SGeneNames, AX219SGeneNames))
GenestoUse <- intersect(GenestoUse, rownames(AX206S@raw.data))
# GenestoUse <- intersect(GenestoUse, rownames(AX207S@raw.data))
GenestoUse <- intersect(GenestoUse, rownames(AX208S@raw.data))
GenestoUse <- intersect(GenestoUse, rownames(AX208RedoS@raw.data))
GenestoUse <- intersect(GenestoUse, rownames(AX206RedoS@raw.data))
GenestoUse <- intersect(GenestoUse, rownames(AX218S@raw.data))
GenestoUse <- intersect(GenestoUse, rownames(AX219S@raw.data))

HEKOnly <- MergeSeurat(AX207S, AX208S)
HEKOnly <- MergeSeurat(HEKOnly, AX208RedoS)

U87Only <- MergeSeurat(AX206S, AX206RedoS)
U87Only <- MergeSeurat(U87Only, AX218S)
U87Only <- MergeSeurat(U87Only, AX219S)

CombinedGenesbyMerge <- MergeSeurat(AX206S, AX207S)
CombinedGenesbyMerge <- MergeSeurat(CombinedGenesbyMerge, AX208S)
CombinedGenesbyMerge <- MergeSeurat(CombinedGenesbyMerge, AX218S)
CombinedGenesbyMerge <- MergeSeurat(CombinedGenesbyMerge, AX219S)
CombinedGenesbyMerge <- MergeSeurat(CombinedGenesbyMerge, AX206RedoS)
CombinedGenesbyMerge <- MergeSeurat(CombinedGenesbyMerge, AX208RedoS)

# Allprotein <- cbind(AX206all[((nrow(AX206all)-3):(nrow(AX206all)-1)),], AX207all[((nrow(AX207all)-3):(nrow(AX207all)-1)),], AX208all[((nrow(AX208all)-3):(nrow(AX208all)-1)),], 
#                     AX218all[((nrow(AX218all)-3):(nrow(AX218all)-1)),], AX219all[((nrow(AX219all)-3):(nrow(AX219all)-1)),], AX206Redoall[((nrow(AX206Redoall)-3):(nrow(AX206Redoall)-1)),], 
#                     AX208Redoall[((nrow(AX208Redoall)-3):(nrow(AX208Redoall)-1)),])

CombinedGenesbyMerge <- SetAssayData(CombinedGenesbyMerge, assay.type = "SCBC", slot = "raw.data", new.data = t(Allprotsnormalized[,1:3]))
CombinedGenesbyMerge <- NormalizeData(CombinedGenesbyMerge, assay.type = "SCBC", normalization.method = "genesCLR", display.progress = F)
CombinedGenesbyMerge <- ScaleData(CombinedGenesbyMerge, assay.type = "SCBC", display.progress = F)
print("Analyzing combined data")
source("BulkComp.R")
# CombinedGenesbyMerge@var.genes <- GenestoUse
# CombinedGenesbyMerge@var.genes <- rownames(CombinedGenesbyMerge@raw.data)[rownames(CombinedGenesbyMerge@raw.data) %in% rownames(resOrdered)]
CombinedGenesbyMerge@var.genes <- TestBulkvar
CombinedGenesbyMerge <- NormalizeData(CombinedGenesbyMerge, display.progress = F)
CombinedGenesbyMerge <- ScaleData(CombinedGenesbyMerge, vars.to.regress = c("nUMI"), display.progress = F)
CombinedGenesbyMerge <- RunPCA(object = CombinedGenesbyMerge, pc.genes = CombinedGenesbyMerge@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)



# CombinedGenesbyMergePlusBulks <- MergeSeurat(CombinedS, U87S)
# # CombinedGenesbyMergePlusBulks <- MergeSeurat(CombinedGenesbyMergePlusBulks, HEKS)
# CombinedGenesbyMergePlusBulks@var.genes <- GenestoUse
# CombinedGenesbyMergePlusBulks <- ScaleData(CombinedGenesbyMergePlusBulks, vars.to.regress = c("nUMI", "orig.ident"))
# CombinedGenesbyMergePlusBulks <- RunPCA(object = CombinedGenesbyMergePlusBulks, pc.genes = CombinedGenesbyMergePlusBulks@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
# CombinedGenesbyMergePlusBulks <- ProjectPCA(object = CombinedGenesbyMergePlusBulks)
# CombinedGenesbyMergePlusBulks <- JackStraw(object = CombinedGenesbyMergePlusBulks, num.replicate = 50, display.progress = FALSE)
# CombinedGenesbyMergePlusBulks <- FindClusters(object = CombinedGenesbyMergePlusBulks, reduction.type = "pca", dims.use = 1:20, resolution = 1.1, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
# CombinedGenesbyMergePlusBulks <- RunTSNE(object = CombinedGenesbyMergePlusBulks, dims.use = 1:20, do.fast = TRUE)
# cluster1.markers <- FindMarkers(object = CombinedGenesbyMergePlusBulks, ident.1 = 1, min.pct = 0.25)
# CombinedGenesbyMergePlusBulks.markers <- FindAllMarkers(object = CombinedGenesbyMergePlusBulks, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

VizPCA(object = CombinedGenesbyMerge, pcs.use = 1:2)
PCAPlot(object = CombinedGenesbyMerge, dim.1 = 1, dim.2 = 2, group.by = "celltype")
CombinedGenesbyMerge <- ProjectPCA(object = CombinedGenesbyMerge)
PCHeatmap(object = CombinedGenesbyMerge, pc.use = 1, do.balanced = TRUE, label.columns = FALSE)
CombinedGenesbyMerge <- JackStraw(object = CombinedGenesbyMerge, num.replicate = 50, display.progress = FALSE)
# JackStrawPlot(object = CombinedGenesbyMerge, PCs = 1:20)
# PCElbowPlot(object = CombinedGenesbyMerge)
CombinedGenesbyMerge <- FindClusters(object = CombinedGenesbyMerge, reduction.type = "pca", dims.use = 1:20, resolution = 1.1, print.output = 0, save.SNN = TRUE, force.recalc=TRUE)
PrintFindClustersParams(object = CombinedGenesbyMerge)
CombinedGenesbyMerge <- RunTSNE(object = CombinedGenesbyMerge, dims.use = 1:20, do.fast = TRUE)
cluster1.markers <- FindMarkers(object = CombinedGenesbyMerge, ident.1 = 1, min.pct = 0.25)
print(x = head(x = cluster1.markers, n = 5))
CombinedGenesbyMerge.markers <- FindAllMarkers(object = CombinedGenesbyMerge, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
CombinedGenesbyMerge.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

Metadata <- CombinedGenesbyMerge@meta.data
Metadata[,"GeneCellRatio"] <- Metadata[,1]/Metadata[,6]
Metadata[,"GeneBeadRatio"] <- Metadata[,1]/Metadata[,7]
Metadata[,"CellBeadRatio"] <- Metadata[,7]/Metadata[,6]

NoCellIncludedMetadata <- data.frame(t(cbind(tail(AX206NoCell,6),tail(AX206RedoNoCell,6),tail(AX207NoCell,6),tail(AX208NoCell,6),tail(AX208RedoNoCell,6),tail(AX218NoCell,6),tail(AX219NoCell,6))))

TSNEPlot(object = CombinedGenesbyMerge, group.by = "orig.ident", pt.size = 3)

# Figure 3 B 
TSNEPlot(object = CombinedGenesbyMerge, group.by = "celltype", pt.size = 4, colors.use = c(NineColScheme[1], NineColScheme[6]), no.legend = TRUE)


RidgePlot(CombinedGenesbyMerge, features.plot = c("B","C","D"), nCol = 2, group.by = "celltype")
ggplot(Metadata, aes(x=Beads, y=nGene))+geom_point()+geom_smooth(method='lm',formula=y~x) + scale_x_continuous(breaks=seq(0,11,1)) + coord_fixed(ratio = 11/4000) + theme(text=element_text(family="Calibri"))
ggplot(Metadata, aes(x=Cells, y=nGene))+geom_point()+geom_smooth(method='lm',formula=y~x) + scale_x_continuous(breaks=seq(0,11,1)) + coord_fixed(ratio = 9/4000) + theme(text=element_text(family="Calibri"))

# cbmc_cite <- RunPCA(CombinedGenesbyMerge, pc.genes = c("B","C","D"), assay.type = "SCBC", pcs.print = 0, pcs.compute = 1:5)
# PCAPlot(cbmc_cite, pt.size = 3, group.by="celltype")
FileName <- "AllCells"
GenesofInterest <- list()
ProteinNames <- c()
U87cells <- CombinedGenesbyMerge@meta.data[,"celltype"]=="U87"
HEKcells <- CombinedGenesbyMerge@meta.data[,"celltype"]=="HEK"
# IntegratedSeuratDataset <- data.frame(as.matrix(t(rbind(CombinedGenesbyMerge@scale.data[CombinedGenesbyMerge@var.genes,U87cells], CombinedGenesbyMerge@assay$SCBC@raw.data[,U87cells]))))
IntegratedSeuratDataset <- data.frame(as.matrix(t(rbind(CombinedGenesbyMerge@scale.data, CombinedGenesbyMerge@assay$SCBC@scale.data))))
for (n in 1:3)
{
  Target <- colnames(Allprotsnormalized[,1:3])[n]
  print(Target)
  # ProteinNames <- c(ProteinNames, Target)
  
  PairwiseMatrixLinearRegression <- apply(IntegratedSeuratDataset[ , 1:(ncol(IntegratedSeuratDataset)-3)], 2, 
                                          function(x) lm(x ~ IntegratedSeuratDataset[ , ncol(IntegratedSeuratDataset)-3+n], 
                                                         data = IntegratedSeuratDataset))
  assign(paste0(Target,"PairwiseLinearRegression"), PairwiseMatrixLinearRegression)
  
  Coefficients <- sapply(PairwiseMatrixLinearRegression,coef)
  assign(paste0(Target,"Coefficients"), Coefficients)
  
  Rsquared <- sapply(PairwiseMatrixLinearRegression,summary)[8,,drop=FALSE]
  assign(paste0(Target,"Rsquared"), Rsquared)
  
  assign(paste0(FileName,Target,"LinearModel"), t(rbind(Coefficients,unlist(Rsquared))))
  
  SpearmanMatrix <- apply(IntegratedSeuratDataset[ , 1:(ncol(IntegratedSeuratDataset)-3)], 2, 
                          function(x) cor.test(x,IntegratedSeuratDataset[ , ncol(IntegratedSeuratDataset)-3+n], method="spearman"))
  assign(paste0(FileName,Target,"Spearman"), SpearmanMatrix)
  
  SpearmanPValues <- sapply(SpearmanMatrix, function(x) x$p.value)
  
  PearsonMatrix <- apply(IntegratedSeuratDataset[ , 1:(ncol(IntegratedSeuratDataset)-3)], 2, 
                         function(x) cor.test(x,IntegratedSeuratDataset[ , ncol(IntegratedSeuratDataset)-3+n], method="pearson"))
  assign(paste0(FileName,Target,"Pearson"), PearsonMatrix)                          
  
  PearsonPValues <- sapply(PearsonMatrix, function(x) x$p.value)   
  
  SignificanceTable <- data.frame(cbind(Rsquared=unlist(Rsquared), SpearmanPValues, PearsonPValues))
  SignificanceTable <- cbind(SignificanceTable, RsquaredThres=SignificanceTable[,"Rsquared"]>0.4,
                             SpearmanPValuesThres=SignificanceTable[,"SpearmanPValues"]<0.05,
                             PearsonPValuesThres=SignificanceTable[,"PearsonPValues"]<0.05)
  SignificanceTable <- cbind(SignificanceTable, SoftHit=SignificanceTable[,"RsquaredThres"]|SignificanceTable[,"SpearmanPValuesThres"]|SignificanceTable[,"PearsonPValuesThres"],
                             HardHit=SignificanceTable[,"RsquaredThres"]&SignificanceTable[,"SpearmanPValuesThres"]&SignificanceTable[,"PearsonPValuesThres"])
  assign(paste0(FileName,Target,"SignificanceTable"), SignificanceTable)
  
  SoftHits <- rownames(SignificanceTable[which(SignificanceTable["SoftHit"]==1),])
  names(SoftHits) <- SoftHits
  SoftHits <- list(data.frame(t(SoftHits)))
  GenesofInterest <- c(GenesofInterest, SoftHits)
}
library(plyr)
GenesofInterest <- t(do.call(rbind.fill, GenesofInterest))
colnames(GenesofInterest) <- ProteinNames
GenesofInterest[is.na(GenesofInterest)] <- ""
library(xlsx)
write.xlsx(GenesofInterest, paste0(FileName, "GenesofInterest.xlsx"), row.names = FALSE)
ggplot(Metadata) + 
  geom_violin(aes(x="nUMI", y=nUMI), width=0.7, fill="red") + 
  geom_jitter(aes(x="nUMI", y=nUMI), width=0.2, size=4, alpha=0.6) +
  geom_violin(aes(x="nGene", y=nGene), width=0.7) + 
  geom_jitter(aes(x="nGene", y=nGene), width=0.2, size=4, alpha=0.6) +
  theme(text=element_text(family="Calibri")) +
  labs(x = "Counts", y = "Metric")
ggplot(IntegratedSeuratDataset, aes(x=B, y=ITGA10))+geom_point()+geom_smooth(method='lm',formula=y~x)
ggplot(AllprotsallPlot, aes(x=variable, y=value, color=Cells)) + 
  geom_jitter(width=0.3, size=4, alpha=0.6) +
  scale_color_manual(values=rev(viridis(9))) +
  ggtitle("Proteins") +
  ylab("Fluorescence (arbitrary units)") +
  xlab("Protein") +
  theme(legend.position = c(0.8,0.8), text=element_text(family="Calibri"))
ggplot(Allprotsnormalizedplot) + 
  geom_boxplot(aes(x=variable, y=value, fill=Chip), outlier.shape = 3) +
  # geom_point(aes(x=variable, y=value, fill=Chip, size=Cells), position=position_dodge(width = 0.75), alpha=0.5) +
  scale_size_discrete(range = c(1,5)) +
  theme(text=element_text(family="Calibri"))
# viridis(9)

AllprotsnormalizedNoRep <- rbind(AX206NormalizedProteins,AX207NormalizedProteins,AX208NormalizedProteins,AX218NormalizedProteins,AX219NormalizedProteins)
# Allprotsall <- rbind(AX206AllProts,AX207AllProts,AX208AllProts,AX218AllProts,AX219AllProts, AX206RedoAllProts, AX208RedoAllProts)
# Allprotsall["Chip"] <- gsub(pattern = "*X(.*)", replacement="", x=gsub(pattern = "AX",replacement="A", x=rownames(Allprotsall)))
# colnames(Allprotsall)[1:3] <- c("PKM2","c-MYC","PDHK1")
# AllprotsallPlot <- melt(Allprotsall, id=c("Cells","Beads", "Chip"))
# AllprotsallPlot[,1] <- as.factor(AllprotsallPlot[,1])
colnames(AllprotsnormalizedNoRep)[1:3] <- c("PKM2", "c-MYC", "PDHK1")
AllprotsnormalizedplotNoRep <- melt(AllprotsnormalizedNoRep, id=c("Cells", "Beads", "Chip"))
AllprotsnormalizedplotNoRep$Chip <- as.factor(AllprotsnormalizedplotNoRep$Chip)
AllprotsnormalizedplotNoRep$Cells <- as.factor(AllprotsnormalizedplotNoRep$Cells)
AllprotsnormalizedplotNoRep["Celltype"] <- NA
AllprotsnormalizedplotNoRep$Celltype[AllprotsnormalizedplotNoRep$Chip %in% c("AX206", "AX218", "AX219")] <- "U87"
AllprotsnormalizedplotNoRep$Celltype[AllprotsnormalizedplotNoRep$Chip %in% c("AX208", "AX207")] <- "HEK"
BProts <- AllprotsnormalizedplotNoRep[AllprotsnormalizedplotNoRep$variable=="PKM2",]
CProts <- AllprotsnormalizedplotNoRep[AllprotsnormalizedplotNoRep$variable=="c-MYC",]
DProts <- AllprotsnormalizedplotNoRep[AllprotsnormalizedplotNoRep$variable=="PDHK1",]
t.test(BProts$value ~ CProts$Celltype)
t.test(CProts$value ~ CProts$Celltype)
t.test(DProts$value ~ DProts$Celltype)



# Figure 3 A Set width to 500

ggplot(AllprotsnormalizedplotNoRep) + 
  geom_boxplot(aes(x=variable, y=value, fill=Celltype), outlier.shape = 3, width = 0.5) +
  scale_size_discrete(range = c(1,5)) +
  scale_fill_manual(values=c(NineColScheme[1], NineColScheme[6])) + 
  theme(text = element_text(family = "Arial"), legend.position = c(0, 0.9)) + 
  coord_fixed(ratio = 1/80) +
  labs(x="Protein", y="Fluorescence (a.u.)") +
  theme(text=element_text(family="Arial", size = 15))



# Supfig 3 A 

ggplot(AllprotsnormalizedplotNoRep[AllprotsnormalizedplotNoRep$Celltype=="U87",]) +
geom_boxplot(aes(x=variable, y=value, fill=Chip), outlier.shape = 3, width = 0.5) +
scale_fill_manual(values=c(NineColScheme[1],NineColScheme[5],NineColScheme[6])) +
theme(text = element_text(family = "Arial"), legend.position = "none") +
scale_y_continuous(limits = c(-8,240)) +
coord_fixed(ratio = 1/100) +
labs(x="Protein", y="Fluorescence (a.u.)") +
theme(text=element_text(family="Arial", size = 15)) 

t.test(CProts$value ~ CProts$Celltype)

  
ggplot(AllprotsnormalizedplotNoRep[AllprotsnormalizedplotNoRep$Celltype=="HEK",]) +
geom_boxplot(aes(x=variable, y=value, fill=Chip), outlier.shape = 3, width = 0.5) +
scale_fill_manual(values=c(NineColScheme[1],NineColScheme[6])) +
theme(text = element_text(family = "Arial"), legend.position = "none") +
scale_y_continuous(limits = c(-5,240)) +
coord_fixed(ratio = 1/80) +
labs(x="Protein", y="Fluorescence (a.u.)") +
theme(text=element_text(family="Arial", size = 15))


Allcellrawmeta <- rbind(AX206Vals, AX207Vals, AX208Vals, AX218Vals, AX219Vals, AX206RedoVals, AX208RedoVals)
ggplot(Allcellrawmeta) + geom_jitter(aes(x=Cells, y=Beads), width = 0.2, height = 0.1, alpha = 0.3)

# Figure 2 b. Reduce width by 33%

NoCellIncludedMetadata <- NoCellIncludedMetadata[NoCellIncludedMetadata$Cells<7,]
NoCellIncludedMetadata <- NoCellIncludedMetadata[NoCellIncludedMetadata$Cells<7,]
NoCellIncludedMetadata$Cells <- factor(NoCellIncludedMetadata$Cells)
ggplot(NoCellIncludedMetadata, aes(x=Cells, fill=Cells, y=TotalReads/Beads)) + 
  geom_boxplot(aes(group = Cells), alpha = 0.4, outlier.color = NA) + 
  geom_jitter(width=0.1) +
  scale_fill_manual(values = NineColScheme) +
  labs(x="Cells", y="Reads per bead") +
  theme(text=element_text(family="Arial", size = 15), legend.position = "none")

CountHeatmap <- data.frame(as.matrix(NoCellIncludedMetadata[,c("Cells", "Beads")] %>% table)) #%>% group_by(Digital, Physical)
CountHeatmap$Cells <- as.numeric(as.character(CountHeatmap$Cells))
CountHeatmap$Beads <- as.numeric(as.character(CountHeatmap$Beads))

ggplot(CountHeatmap, aes(y=Cells, x=Beads, color=Freq))+geom_point(size = 9)+scale_color_gradientn(colors = c("#FFFFFF",NineColScheme[1:5]))+scale_x_continuous(breaks=0:max(CountHeatmap$Beads),limits = c(0,max(CountHeatmap$Beads)))+scale_y_continuous(breaks = 0:max(CountHeatmap$Cells), limits = c(0,max(CountHeatmap$Cells)))+theme(text=element_text(family="Arial", size = 15))+coord_fixed(ratio=1)
ggplot(NoCellIncludedMetadata) + geom_point(aes(x=Beads, y=TotalReads))

# SupFig 2 B
TestS <- AX206RedoS
colnames(TestS@scale.data) <- gsub("AX206Redo", "", colnames(TestS@scale.data))
colnames(TestS@scale.data) <- gsub("-.", " ", colnames(TestS@scale.data))
colnames(TestS@scale.data) <- gsub("Y", "Y-", colnames(TestS@scale.data))
colnames(TestS@scale.data) <- gsub("X", "X-", colnames(TestS@scale.data))
heatmap.2(as.matrix(TestS@scale.data), trace = "none", margins = c(5,2), labRow = FALSE)
heatmap.2(as.matrix(TestS@scale.data[TestS@var.genes,]), trace = "none", margins = c(5,2), labRow = FALSE)

heatmap.2(as.matrix(AX206RedoS@scale.data), trace="none", margins = c(8,5), labRow = FALSE)
VlnPlot(object = CombinedGenesbyMerge, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "orig.ident", y.lab.rot = TRUE)

FeaturePlot(CombinedGenesbyMerge, features.plot = c("B","C","D"), cols.use = c("lightgrey", "blue"), pt.size = 2, nCol = 1)

# SupFig 2

ggplot(Metadata, aes(fill=celltype)) + 
  geom_violin(aes(x="Genes", y=nGene), scale = "count") +
  geom_violin(aes(x="Transcripts", y=nUMI), scale = "count") +
  coord_fixed(ratio = 1/10000) +
  scale_fill_manual(values = c(NineColScheme[1],NineColScheme[6])) +
  labs(y="Counts") +
  theme(text=element_text(family="Arial", size = 15), legend.position = "none")

VlnPlot(object = CombinedGenesbyMerge, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3, group.by = "orig.ident", y.lab.rot = TRUE)

#SupFig3 

FeaturePlot(CombinedGenesbyMerge, features.plot = c("D"), cols.use = c("lightgrey", NineColScheme[6]), pt.size = 4)