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
print("Loading data")
load("AX206alldata")
AX206all <- IntegratedData
load("AX207alldata")
AX207all <- IntegratedData
load("AX208alldata")
AX208all <- IntegratedData
load("AX206Redoalldata")
AX206Redoall <- IntegratedData
load("AX208Redoalldata")
AX208Redoall <- IntegratedData
load("AX218alldata")
AX218all <- IntegratedData
load("AX219alldata")
AX219all <- IntegratedData

load("AllProteinValues")
# AX206Vals <- data.frame(t(ProteinsPerBeads)), rownames(AX219Vals) <- gsub("X","AX219X",rownames(AX219Vals)), AX219Zeros <- AX219Vals[which(AX219Vals[,4]==0),]
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
CombinedGenesbyMerge@var.genes <- GenestoUse
CombinedGenesbyMerge <- ScaleData(CombinedGenesbyMerge, vars.to.regress = c("nUMI", "orig.ident"))
CombinedGenesbyMerge <- RunPCA(object = CombinedGenesbyMerge, pc.genes = CombinedGenesbyMerge@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)


VizPCA(object = CombinedGenesbyMerge, pcs.use = 1:2)
PCAPlot(object = CombinedGenesbyMerge, dim.1 = 1, dim.2 = 2)
CombinedGenesbyMerge <- ProjectPCA(object = CombinedGenesbyMerge)
PCHeatmap(object = CombinedGenesbyMerge, pc.use = 1, do.balanced = TRUE, label.columns = FALSE)
CombinedGenesbyMerge <- JackStraw(object = CombinedGenesbyMerge, num.replicate = 50, display.progress = FALSE)
JackStrawPlot(object = CombinedGenesbyMerge, PCs = 1:20)
PCElbowPlot(object = CombinedGenesbyMerge)
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

ggplot(Metadata, aes(x=Beads, y=nGene))+geom_point()+geom_smooth(method='lm',formula=y~x) + scale_x_continuous(breaks=seq(0,11,1))

TSNEPlot(object = CombinedGenesbyMerge, group.by = "orig.ident", pt.size = 3)

RidgePlot(CombinedGenesbyMerge, features.plot = c("B","C","D"), 
          nCol = 2, group.by = "celltype")

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
                             SpearmanPValuesThres=SignificanceTable[,"SpearmanPValues"]<0.15,
                             PearsonPValuesThres=SignificanceTable[,"PearsonPValues"]<0.15)
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
ggplot(Metadata) +geom_violin(aes(x="nUMI", y=nUMI), width=0.7, fill="red")+ geom_jitter(aes(x="nUMI", y=nUMI), width=0.2, size=4, alpha=0.6)+geom_violin(aes(x="nGene", y=nGene), width=0.7)+ geom_jitter(aes(x="nGene", y=nGene), width=0.2, size=4, alpha=0.6)+theme(text=element_text(family="Calibri"))+ylab(label = "Counts")+xlab(label = "Metric")
ggplot(IntegratedSeuratDataset, aes(x=B, y=ITGA10))+geom_point()+geom_smooth(method='lm',formula=y~x)
ggplot(AllprotsallPlot, aes(x=variable, y=value, color=Cells)) + geom_jitter(width=0.3, size=4, alpha=0.6)+scale_color_manual(values=rev(viridis(9)))+ggtitle("Proteins")+ylab("Fluorescence (arbitrary units)")+xlab("Protein")+theme(legend.position = c(0.8,0.8), text=element_text(family="Calibri"))
ggplot(Allprotsnormalizedplot) + geom_boxplot(aes(x=variable, y=value, fill=Chip), outlier.shape = 3)+geom_point(aes(x=variable, y=value, fill=Chip, size=Cells), position=position_dodge(width = 0.75), alpha=0.5)+scale_size_discrete(range = c(1,5))+theme(text=element_text(family="Calibri"))
# viridis(9)