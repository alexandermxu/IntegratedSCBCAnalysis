library(xlsx)
library(plyr)
library(tidyverse)
library(Matrix)
library(gplots)

GeneAnalysis <- TRUE
GeneSource <- "C:/Users/alexm/Dropbox (Personal)/AX206-8/AX218/"

setwd("C:/Users/alexm/Documents/ProteinChips/")

MYCUp <- read.table("MYCUP.txt", as.is = TRUE)
MYCDown <- read.table("MYCDOWN.txt", as.is = TRUE)
PyruvateTargets <- read.table("PYRUVATE.txt", as.is = TRUE)
# Pfirst
# FileName <- "AX204"
# Cells <- c(1,0,0,0,0,1,1,3,1,1,0,0,
#            0,5,6,2,2,2,2,1,0,1,1,3,
#            2,2,1,0,1,2,0,0,0,1,0,0)
# Beads <- c(1,1,4,0,2,2,3,8,6,4,2,5,
#            5,8,0,6,7,4,5,2,2,6,9,3,
#            5,3,1,3,1,7,3,2,1,3,2,0)
# FileName <- "AX206"
# Cells <- c(2,1,0,0,0,1,4,0,3,7,0,3,
#            7,8,6,9,0,4,1,1,1,4,1,3,
#            0,5,1,1,0,0,0,1,0,2,0,0)
# Beads <- c(5,4,4,1,3,6,4,1,3,4,5,3,
#            9,2,0,1,1,2,2,4,4,6,2,3,
#            2,6,6,2,3,2,1,4,1,1,0,3)
# FileName <- "AX207"
# Cells <- c(0,0,0,2,0,0,1,0,1,0,0,3,
#            0,0,0,0,1,2,0,0,0,0,0,0,
#            0,1,1,0,0,0,0,0,0,0,0,0)
# Beads <- c(6,9,7,2,4,7,7,5,4,2,12,3,
#            5,1,2,3,2,1,2,16,1,2,2,2,
#            2,2,2,0,2,3,0,0,2,1,6,1)
# FileName <- "AX208"
# Cells <- c(0,0,1,1,1,1,6,1,1,5,2,1,
#            5,2,0,4,1,1,0,0,1,0,0,0,
#            0,0,0,0,4,2,0,0,1,0,0,0)
# Beads <- c(0,2,1,1,1,2,0,4,4,3,7,2,
#            1,4,3,5,2,3,0,0,0,0,0,0,
#            1,7,10,1,1,2,0,0,0,0,0,0)


# Psecond
# FileName <- "AX212"
# Cells <- c(0,0,1,0,0,0,3,4,0,1,0,0,
#            0,1,0,1,1,2,0,0,1,0,0,0,
#            0,0,0,0,0,0,0,0,0,0,0,0)
# Beads <- c(0,0,1,2,0,0,3,3,1,0,0,0,
#            1,0,1,4,2,1,0,5,2,0,1,2,
#            2,3,3,1,0,2,2,0,1,0,2,2)
# FileName <- "AX218"
# Cells <- c(0,0,1,1,0,1,0,2,1,3,5,0,
#            0,8,3,4,0,2,1,2,2,4,4,3,
#            0,0,0,0,1,1,0,1,1,0,0,2)
# Beads <- c(2,1,2,0,2,1,3,6,6,11,10,0,
#            7,9,6,7,3,2,5,2,6,3,4,10,
#            5,5,3,4,4,2,1,2,1,3,3,5)
FileName <- "AX219"
Cells <- c(1,3,1,0,1,2,5,0,2,2,5,4,
           2,4,2,6,2,6,0,6,0,2,4,0,
           1,6,2,1,1,0,2,0,0,2,1,0) #maybe 1 in first of this block
Beads <- c(1,6,8,2,4,7,9,8,4,6,2,5,
           6,7,2,5,4,7,4,6,7,3,3,5,
           13,5,4,4,3,4,9,8,11,9,7,3)

Barcodes <- c("B", "C", "D")
File <- read.xlsx(paste0(FileName,".xlsx"), sheetIndex = 2)
Blocks <- unique(File$Block)
Group1 <- c(1:6, 13:18, 25:30)
Group2 <- c(7:12, 19:24, 31:36)
# #PFirst
# Registry1 <- c(NA,NA,NA,NA,"D","C","B",NA,NA,NA,NA)
# Registry2 <- c(NA,NA,NA,"B","C","D",NA,NA,NA,NA,NA)
# PSecond
# Registry1 <- c(NA,NA,NA,"B","C","D",NA,NA,NA,NA,NA)
# Registry2 <- c(NA,NA,NA,NA,"D","C","B",NA,NA,NA,NA)

# I at position 9 first row
Registry1 <- c(NA,NA,NA,NA,"B","C","D",NA,NA,NA,NA)
Registry2 <- c(NA,NA,NA,"D","C","B",NA,NA,NA,NA,NA)
A1TopRight <- c(36,25,24,13,12,1,35,26,23,14,11,2,34,27,22,15,10,3,33,28,21,16,9,4,32,29,20,17,8,5,31,30,19,18,7,6)
# F1TopRight <- c(6,7,18,19,30,31,5,8,17,20,29,32,4,9,16,21,28,33,3,10,15,22,27,34,2,11,14,23,26,35,1,12,13,24,25,36)
KeepData <- c("Block", "Row", "F635.Mean", "F532.Mean")
Entries <- length(KeepData)

Output <- array(dim = c(length(Blocks)*length(Barcodes), Entries))
Output <- data.frame(Output)
colnames(Output) <- KeepData

FileExtract <- File[,KeepData]
FileExtract[FileExtract$Block %in% Group1,Entries+1] <- Registry1[FileExtract[FileExtract$Block %in% Group1,]$Row]
FileExtract[FileExtract$Block %in% Group2,Entries+1] <- Registry2[FileExtract[FileExtract$Block %in% Group2,]$Row]
ReducedFileExtract <- na.omit(FileExtract)
colnames(ReducedFileExtract) <- c("Block", "Row", "Mean635", "Mean532", "Barcode")
Values <- ReducedFileExtract[,-c(1:2,Entries+1)]
i=1
Output[c("Chamber", "Cells", "Beads")] <- NA

for (n in 1:length(Blocks))
{
  for (m in 1:length(Barcodes))
  {
    Output[i,1] <- Blocks[n]
    Output[i,2] <- Barcodes[m]
    Output[i,3:Entries] <- colMeans(Values[ReducedFileExtract$Block==Blocks[n] & ReducedFileExtract$Barcode==Barcodes[m],])
    i=i+1
  }
}

Output[,Entries+1] <- A1TopRight[as.numeric(as.character(Output[,1]))]

Output[,Entries+2] <- Cells[as.numeric(as.character(Output[,Entries+1]))]

Output[,Entries+3] <- Beads[as.numeric(as.character(Output[,Entries+1]))]

colnames(Output) <- c("Block", "Barcode", "Mean635", "Mean532", "Chamber", "Cells", "Beads")
# colnames(Output) <- c("Block", "Barcode", "Mean635", "Mean594", "Mean532", "Mean488", "Chamber", "Cells")

Output["Normalized635"] <- Output[,"Mean635"]/Output[,"Cells"]

ZeroCellChambers <- Output[Output["Cells"]==0,][c("Barcode", "Mean635")]

BackgroundValues <- ddply(ZeroCellChambers, .(Barcode), numcolwise(mean))


RawPlot <- ggplot(Output, aes(x=Barcode,y=Mean635,label=Block)) + 
  geom_jitter(aes(color=factor(Cells)), width=0.3, size=4, alpha=0.4) + 
  labs(title=FileName)
assign(paste0(FileName, "RawPlot"), RawPlot)

NormOutput <- Output


for(i in 1:length(Barcodes))
{
  Label <- Barcodes[i]
  Suboutput <- Output[which(Output$Barcode == Label),]
  Max635 <- max(Suboutput$Mean635)
  Min635 <- min(Suboutput$Mean635)
  Norm635 <- (Suboutput$Mean635-Min635)/(Max635-Min635)
  NormOutput[which(NormOutput$Barcode == Label),"Mean635"] <- Norm635
}
NormPlot <- ggplot(NormOutput, aes(x=Barcode,y=Mean635,label=Block)) + geom_jitter(aes(color=factor(Cells)), width=0.3, size=4, alpha=0.4) + labs(title=FileName)
assign(paste0(FileName, "NormPlot"), NormPlot)

write.csv(Output, file = paste0(FileName,"Output.csv"))

if(GeneAnalysis==TRUE)
{
  # A1TopRightBarcodeConversion <- c("X6-FY1-1","X5-EY1-1","X4-DY1-1","X3-CY1-1","X2-BY1-1","X1-AY1-1",
  #                                "X6-FY2-2","X5-EY2-2","X4-DY2-2","X3-CY2-2","X2-BY2-2","X1-AY2-2",
  #                                "X6-FY3-3","X5-EY3-3","X4-DY3-3","X3-CY3-3","X2-BY3-3","X1-AY3-3",
  #                                "X6-FY4-4","X5-EY4-4","X4-DY4-4","X3-CY4-4","X2-BY4-4","X1-AY4-4",
  #                                "X6-FY5-5","X5-EY5-5","X4-DY5-5","X3-CY5-5","X2-BY5-5","X1-AY5-5",
  #                                "X6-FY6-6","X5-EY6-6","X4-DY6-6","X3-CY6-6","X2-BY6-6","X1-AY6-6")

A1TopRightBarcodeConversion <- c("X6-FY1-1","X5-EY1-1","X4-DY1-1","X3-CY1-1","X2-BY1-1","X1-AY1-1",
                                 "X1-AY2-2","X2-BY2-2","X3-CY2-2","X4-DY2-2","X5-EY2-2","X6-FY2-2",
                                 "X6-FY3-3","X5-EY3-3","X4-DY3-3","X3-CY3-3","X2-BY3-3","X1-AY3-3",
                                 "X1-AY4-4","X2-BY4-4","X3-CY4-4","X4-DY4-4","X5-EY4-4","X6-FY4-4",
                                 "X6-FY5-5","X5-EY5-5","X4-DY5-5","X3-CY5-5","X2-BY5-5","X1-AY5-5",
                                 "X1-AY6-6","X2-BY6-6","X3-CY6-6","X4-DY6-6","X5-EY6-6","X6-FY6-6")
setwd(GeneSource)

Genes <- read.table("output.dge.txt.gz", header=TRUE, row.names=1)

RefGenes <- Genes

Location_BeadRegistry <- read.csv(paste0(FileName, "LocationCompendium.csv"))

colnames(Genes) <- Location_BeadRegistry$GroupNumberLabel[sapply(colnames(Genes), function(x) which(Location_BeadRegistry$Reference==x))]

MaxGenes <- apply(Genes,1,max)>5

TotalGenes <- apply(Genes,1,sum)>15

VarianceGenes <- apply(Genes,1,var)>2

PresentGenes <- apply(Genes,1,median)>0

GenesSelected <- MaxGenes & VarianceGenes & TotalGenes #& PresentGenes

GenesNormalizedByExpression <- Genes / apply(Genes, 2, sum)[col(Genes)]*10000

LogGenesNormalizedByExpression <- log(GenesNormalizedByExpression+1)

SDNormalizedLogGenes <- LogGenesNormalizedByExpression/apply(LogGenesNormalizedByExpression,1,sd)

GenesForLinearModel <- SDNormalizedLogGenes[GenesSelected,]

TotalReads <- colSums(Genes)

Location_BeadRegistry["Chamber"] <- NA

XYChamberConversion <- vapply(Location_BeadRegistry$GroupNumberLabel,function(x) which(A1TopRightBarcodeConversion==x), FUN.VALUE = 1)

Location_BeadRegistry$Chamber <- as.factor(XYChamberConversion)

ProteinsPerBeads <- data.frame(row.names=c(Barcodes, "Cells", "Beads"))

ProteinsPerBeads[names(Genes)] <- NA



for(n in 1:length(unique(Location_BeadRegistry$Replacements)))
{
  CurrentBead <- unique(Location_BeadRegistry$Replacements)[n]
  CurrentChamber <- Location_BeadRegistry$Chamber[match(CurrentBead,Location_BeadRegistry$Replacements)]
  ProteinLevels <- Output$Mean635[Output$Chamber==CurrentChamber]
  CellCount <- Output$Cells[Output$Chamber==CurrentChamber][1]
  BeadCount <- Output$Beads[Output$Chamber==CurrentChamber][1]
  ColumnNumber <- which(names(RefGenes)==CurrentBead)
  ProteinsPerBeads[ColumnNumber] <- c(ProteinLevels, CellCount, BeadCount)
}

CellBeadPositiveBarcodes <- ProteinsPerBeads["Cells",]>0
HitThresholdChambers <- apply(Genes,2,sum)>500

KeptChambers <- CellBeadPositiveBarcodes & HitThresholdChambers
GenesWithCells <- Genes[,KeptChambers]


BGSubtractedProteins <- sweep(ProteinsPerBeads[Barcodes,KeptChambers],1,BackgroundValues[,2])
rownames(BGSubtractedProteins) <- paste0(Barcodes, "minusBG")
BGSubtractedProteins[BGSubtractedProteins<0] <- 1

CellNormalizedBGSProteins <- BGSubtractedProteins / ProteinsPerBeads["Cells",KeptChambers][col(BGSubtractedProteins)]
rownames(CellNormalizedBGSProteins) <- paste0(Barcodes, "minusBGpercell")

LogProteins <- log(CellNormalizedBGSProteins)

SDNormalizedNoLogProteins <- CellNormalizedBGSProteins/apply(CellNormalizedBGSProteins,1,sd)
SDNormalizedProteins <- LogProteins/apply(LogProteins,1,sd)

rownames(SDNormalizedProteins) <- c("PKMprotein", "CMYCprotein", "PDHK1protein")
rownames(SDNormalizedNoLogProteins) <- c("PKMprotein", "CMYCprotein", "PDHK1protein")

IntegratedData <- rbind(GenesWithCells, ProteinsPerBeads[,KeptChambers], BGSubtractedProteins, CellNormalizedBGSProteins, "TotalReads"=TotalReads[KeptChambers])
NoCellIntegratedData <- rbind(Genes, ProteinsPerBeads, "TotalReads"=TotalReads)

SingleCellChambers <- IntegratedData["Cells",]==1
# IntegratedDataPositiveCells <- IntegratedData[,CellBeadPositiveBarcodes]
write.csv(x=IntegratedData, file=paste0(FileName,"FullData.csv"))




# ProteinsandGenes <- rbind(apply(LogGenesNormalizedByExpression[rownames(IntegratedDataPositiveCells) %in% MYCUp[,1],CellBeadPositiveBarcodes],2,mean),
#                           apply(LogGenesNormalizedByExpression[rownames(IntegratedDataPositiveCells) %in% MYCDown[,1],CellBeadPositiveBarcodes],2,mean),
#                           apply(LogGenesNormalizedByExpression[rownames(IntegratedDataPositiveCells) %in% PyruvateTargets[,1],CellBeadPositiveBarcodes],2,mean),
#                           IntegratedDataPositiveCells[c("B -BG per Cell", "C -BG per Cell", "D -BG per Cell"),])
# 
# rownames(ProteinsandGenes)[1:3] <- c("MYCUP", "MYCDOWN", "PYRUVATE")
# 
# PKM_PyruvateGeneRatio <- ProteinsandGenes["C -BG per Cell",]/ProteinsandGenes["PYRUVATE",]
# PDHK1_PyruvateGeneRatio <- ProteinsandGenes["D -BG per Cell",]/ProteinsandGenes["PYRUVATE",]
# CMYC_CMYYCUpRatio <- ProteinsandGenes["C -BG per Cell",]/ProteinsandGenes["MYCUP",]
# CMYC_CMYYCDownRatio <- ProteinsandGenes["C -BG per Cell",]/ProteinsandGenes["MYCDOWN",]
# ProteinsandGenes <- rbind(ProteinsandGenes, PKM_PyruvateGeneRatio,PDHK1_PyruvateGeneRatio,CMYC_CMYYCUpRatio, CMYC_CMYYCDownRatio)
# 
# rownames(ProteinsandGenes)[4:10] <- c("PKM", "CMYC", "PDHK1", "PKM/PYRUVATE", "PDHK1/PYRUVATE","CMYC/CMYCUP", "CMYC/CMYCDOWN")
# 
# 

LinearModelInitial <- rbind(GenesForLinearModel[,KeptChambers], SDNormalizedProteins)

LinearModelTranspose <- data.frame(t(LinearModelInitial))

LinearModelData <- (LinearModelTranspose / apply(LinearModelTranspose,2,max)[col(LinearModelTranspose)])

# LinearModelData <- LinearModelDatawNA[,colnames(LinearModelDatawNA)[!is.na(apply(LinearModelDatawNA,2,max))]]
GenesofInterest <- list()
ProteinNames <- c()
for (n in 1:3)
{
  Target <- rownames(SDNormalizedProteins)[n]
  print(Target)
  ProteinNames <- c(ProteinNames, Target)
  
  PairwiseMatrixLinearRegression <- apply(LinearModelData[ , 1:(ncol(LinearModelData)-3)], 2, 
                          function(x) lm(x ~ LinearModelData[ , ncol(LinearModelData)-3+n], 
                                         data = LinearModelData))
  assign(paste0(Target,"PairwiseLinearRegression"), PairwiseMatrixLinearRegression)
  
  Coefficients <- sapply(PairwiseMatrixLinearRegression,coef)
  assign(paste0(Target,"Coefficients"), Coefficients)
  
  Rsquared <- sapply(PairwiseMatrixLinearRegression,summary)[8,,drop=FALSE]
  assign(paste0(Target,"Rsquared"), Rsquared)
  
  assign(paste0(FileName,Target,"LinearModel"), t(rbind(Coefficients,unlist(Rsquared))))
  
  SpearmanMatrix <- apply(LinearModelData[ , 1:(ncol(LinearModelData)-3)], 2, 
                          function(x) cor.test(x,LinearModelData[ , ncol(LinearModelData)-3+n], method="spearman"))
  assign(paste0(FileName,Target,"Spearman"), SpearmanMatrix)
  
  SpearmanPValues <- sapply(SpearmanMatrix, function(x) x$p.value)
                          
  PearsonMatrix <- apply(LinearModelData[ , 1:(ncol(LinearModelData)-3)], 2, 
                          function(x) cor.test(x,LinearModelData[ , ncol(LinearModelData)-3+n], method="pearson"))
  assign(paste0(FileName,Target,"Pearson"), PearsonMatrix)                          
           
  PearsonPValues <- sapply(PearsonMatrix, function(x) x$p.value)   
  
  SignificanceTable <- data.frame(cbind(Rsquared=unlist(Rsquared), SpearmanPValues, PearsonPValues))
  SignificanceTable <- cbind(SignificanceTable, RsquaredThres=SignificanceTable[,"Rsquared"]>0.4,
                             SpearmanPValuesThres=SignificanceTable[,"SpearmanPValues"]<0.10,
                             PearsonPValuesThres=SignificanceTable[,"PearsonPValues"]<0.10)
  SignificanceTable <- cbind(SignificanceTable, SoftHit=SignificanceTable[,"RsquaredThres"]|SignificanceTable[,"SpearmanPValuesThres"]|SignificanceTable[,"PearsonPValuesThres"],
                             HardHit=SignificanceTable[,"RsquaredThres"]&SignificanceTable[,"SpearmanPValuesThres"]&SignificanceTable[,"PearsonPValuesThres"])
  assign(paste0(FileName,Target,"SignificanceTable"), SignificanceTable)
  
  SoftHits <- rownames(SignificanceTable[which(SignificanceTable["SoftHit"]==1),])
  names(SoftHits) <- SoftHits
  SoftHits <- list(data.frame(t(SoftHits)))
  GenesofInterest <- c(GenesofInterest, SoftHits)
}

GenesofInterest <- t(do.call(rbind.fill, GenesofInterest))
colnames(GenesofInterest) <- ProteinNames
GenesofInterest[is.na(GenesofInterest)] <- ""

write.xlsx(GenesofInterest, paste0(FileName, "GenesofInterest.xlsx"), row.names = FALSE)

# ggplot(LinearModelData, aes(x=CMYCprotein, y=PDHK1protein))+geom_point()+geom_smooth(method='lm',formula=y~x)
# heatmap.2(as.matrix(t(LinearModelData)), trace="none", colRow = )
}
