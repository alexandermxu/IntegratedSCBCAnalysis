# Paper Figs
# Fig 2

# run ProteinChipProcessing on AX206
Output$Cells <- as.factor(Output$Cells)
RawPlot <- ggplot(Output, aes(x=Barcode,y=Mean635,label=Block)) +
geom_jitter(aes(color=Cells), width=0.3, size=4, alpha=0.4) +
coord_fixed(ratio=1/400) +
scale_y_continuous(breaks = seq(0, 1500, by=300), limits=c(0,1500)) +
scale_x_discrete(labels = c("PKM2", "c-MYC", "PDHK1")) +
labs(x="Protein", y="Fluorescence (a.u.)") +
theme(text=element_text(family="Arial", size = 15), legend.position = c(0.75, 0.7)) +
scale_color_manual(values = NineColScheme)



# Run ProteinChipProcessing on Chip 1 
OldGeneData <- Genes[,KeptChambers]

# Run Protein Chip Processing on Chip 2
NewGeneData <- Genes[,KeptChambers]

# Run ReplicateTest

# read.csv("SpeciesMixingOnChip")
# speciesmix <- read.xlsx("SpeciesMixingOnChip.xlsx", header = TRUE, sheetIndex = 1)
# speciesmix <- speciesmix[-73,]
# speciesmix <- speciesmix[-1,]
# speciesmix <- speciesmix[,-4]
# colnames(speciesmix) <- c("Bead","Human Genes","Human Transcripts", "Mouse Genes", "Mouse Transcripts")
# speciesmix[,2] <- as.numeric(levels(speciesmix[,2])[speciesmix[,2]])
# speciesmix[,3] <- as.numeric(levels(speciesmix[,3])[speciesmix[,3]])
# speciesmix[,4] <- as.numeric(levels(speciesmix[,4])[speciesmix[,4]])
# speciesmix[,5] <- as.numeric(levels(speciesmix[,5])[speciesmix[,5]])
# ggplot(speciesmix) + geom_point(aes(x=`Human Genes`, y=`Mouse Genes`), size=2) + ylim(0,1400) + xlim(0,2100)+theme(text=element_text(family="Calibri"))

# SupFig 2 
speciesmix <- read.xlsx("SpeciesMixingOnChip.xlsx", header = TRUE, sheetIndex = 2)
colnames(speciesmix) <- c("Bead","Human Genes","Human Transcripts", "Mouse Genes", "Mouse Transcripts")
ggplot(speciesmix) + geom_point(aes(x=`Human Genes`, y=`Mouse Genes`), size=2) + ylim(0,1400) + xlim(0,2100)+theme(text=element_text(family="Calibri"))

# Figure 2 C Change height to 500

locationcomp <- read.xlsx(file = "SequencingLocationCompression.xlsx", sheetIndex = 2)
colnames(locationcomp) <- c("Bead","Genes", "Transcripts", "Clustering Method")
ggplot(locationcomp, aes(x=`Clustering Method`, y=Transcripts, fill=`Clustering Method`)) + 
  geom_violin(alpha = 0.4) + 
  geom_jitter(width = 0.1, show.legend= FALSE) + 
  scale_fill_manual(values=c(NineColScheme[1], NineColScheme[6])) + 
  scale_color_manual(values = c(NineColScheme[1], NineColScheme[6])) +
  theme(text = element_text(family = "Arial"), legend.position = "none") + 
  coord_fixed(ratio = 1/5000)


PCA1 <- c("LGALS1", "VIM", "B2M", "S100A6", "TMSB4X", "RPS27L", "MT2A", "AKR1B1", "ANXA2","GNG11", "FTL", "FTH1", "ZFP36L1", "ANXA1", "CD63", "SCD", "HSP90B1",  "SNHG3","SH3BGRL3", "UCHL1", "RPL22L1", "CALR", "CSTB","NEAT1","PRDX5","RAB13","SKP1","TMSB10","TFPI2","S100A10","XIST","CKB","RPS3AP47","RP11-234A1.1","RPS2","RPL5","PTMA", "RP11-478C6.4","RPS10","RPL38","HMGB1","RPL5P1","RPL39", "RPL21","RPS6","RPL22","CTD-3035D6.1","PARP1","RPL24","RPS29","COX7C","RPL31","RPS8","RPS13P2","RPL34","HNRNPA1","NCL","RPS17","RPS4X","RPL10A")
  # c("PSMA3",  "PSMA7" , "TARS" ,  "ARPC3" , "SERF2" , "RPLP2" , "ALDOA" , "SRP14",  "EEF1D" , "ACTG1" , "PSMB1" ,"OAZ1"  , "RPLP1",  "RPL9"   ,"ATP5G3" ,"EEF1A1", "YWHAB" , "ATPIF1" ,"MATR3" , "RPS19"  ,"ENY2"  , "RPL37" , "GTF2H5" ,"ATP5E",  "EEF2" ,  "RPL13"  ,"SNRPD2" ,"ATF4" ,  "XRCC5"  ,"RPS5"  , "PARP1"   ,   "RPL39"   ,   "RBMX"  ,     "AC004453.8" ,"NOLC1"  ,    "RPL10P6" ,   "RPS10"  ,    "HMGB1" , "MT-ATP8"  ,"SERBP1"  ,   "MT-ND3"   ,  "SET"     ,   "MTATP6P1"   ,"MT-RNR1" ,   "MT-ND2"   ,  "PNN"      ,  "MT-ATP6" ,   "NASP"    ,   "HNRNPU"  ,   "SMC1A"  ,    "PSMA4"   ,   "HMGB2"   ,   "MT-CYB"   ,  "SNRPD1"   , "NHP2" ,      "HNRNPA1"   , "NCL"  ,     "NDUFS5"  ,   "MT-ND6"    , "MRPS21" )
Reorderbycelltypevector <- c(which(CombinedGenesbyMerge@meta.data$celltype=="HEK"), which(CombinedGenesbyMerge@meta.data$celltype=="U87"))
PCA3toPlotbyCellType <- CombinedGenesbyMerge@scale.data[PCA3,Reorderbycelltypevector]
PCA3toPlotbyPCALoad <- CombinedGenesbyMerge@scale.data[PCA3,]
HEKPCA3 <- PCA3toPlot[,1:19]
U87PCA3 <- PCA3toPlot[,20:ncol(PCA3toPlot)]
PCA3CellLoadingOrder <- order(CombinedGenesbyMerge@dr$pca@cell.embeddings[,3], decreasing=TRUE)
PCA3ReorderByPCA <- match(CombinedGenesbyMerge@cell.names,PCA3CellLoadingNames)
PCA3CelltypeOrderbyPCALoad[PCA3CelltypeOrderbyPCALoad=="HEK"] <- "blue"
PCA3CelltypeOrderbyPCALoad[PCA3CelltypeOrderbyPCALoad=="U87"] <- "red"

ReorderedCellVector <- c(hclust(dist(t(HEKPCA3)))$order, hclust(dist(t(U87PCA3)))$order)
PCA3CellLoading <- sort(CombinedGenesbyMerge@dr$pca@cell.embeddings[,3], decreasing=TRUE)
PCA3CellLoadingNames <- names(PCA3CellLoading)
PCA3GenesNewOrder <- c(31:43, seq(13,1,-1))
PCA3toPlotbyPCALoad[PCA3toPlotbyPCALoad<=-2.5] <- -2.5
PCA3toPlotbyPCALoad[PCA3toPlotbyPCALoad>=2.5] <- 2.5
heatmap.2(PCA3toPlotbyPCALoad[PCA3GenesNewOrder, PCA3CellLoadingOrder], trace="none", dendrogram="none", col = PurpleAndYellow(), Rowv = NULL, Colv=NULL, ColSideColors = PCA3CelltypeOrderbyPCALoad, key.title = NA, key.ylab = NA, labCol = NA)

# Focus on AX206 
Keeps <- tail(AX206RedoNoCell)[,c(1,4,19,28)]
Keeps <- data.frame(t(Keeps))
Keeps[c("X","Y")] <- cbind(c(6,5,6,5), c(3,3,2,2))
KeepsB <- Keeps[,c(1,4:8)]
KeepsC <- Keeps[,c(2,4:8)]
KeepsD <- Keeps[,c(3,4:8)]
KeepsB["Protein"] <- "B"
KeepsC["Protein"] <- "C"
KeepsD["Protein"] <- "D"
colnames(KeepsB)[1] <- "Value"
colnames(KeepsC)[1] <- "Value"
colnames(KeepsD)[1] <- "Value"
NewKeeps <- rbind(KeepsB, KeepsC, KeepsD)
g <- ggplot(NewKeeps)
g+geom_bar(aes(y=Value, fill=Protein, x=Protein), stat="identity")+facet_grid(facets = Y ~ X, labeller = label_bquote(cols = X - .(X), rows = Y - .(Y)), switch = "y")+scale_y_continuous(position = "right")
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
ggplot(LocCompress)+geom_col(aes(x=1, y=Beads, fill=Chamber), width=0.3) + geom_col(aes(x=2, y=Ones, fill=Chamber), width=0.3) +scale_fill_manual(values = col_vector)
brewercycle <- c(brewer.pal(n = 9, name = "Set1"),brewer.pal(n = 8, name = "Dark2"),brewer.pal(n = 8, name = "Accent"),brewer.pal(n = 8, name = "Pastel2"))
ggplot(LocCompress)+geom_col(aes(x=1, y=Beads, fill=Chamber), width=0.3) + geom_col(aes(x=2, y=Ones, fill=Chamber), width=0.3) +scale_fill_manual(values = brewercycle)

LOCSCom <- read.csv("AX206LocationCompendium.csv")
BeadsPutBack <- table(LOCSCom$GroupNumberLabel)
BeadCaptureRate <- data.frame("Digital" = cbind(as.vector(BeadsPutBack), "Physical" = Beads[match(names(BeadsPutBack), A1TopRightBarcodeConversion)], "Chip" = "AX206"))
rownames(BeadCaptureRate) <- names(BeadsPutBack)
TotalBeadCapture <- BeadCaptureRate
TotalBeadCapture <- cbind(TotalBeadCapture, BeadCaptureRate)


LOCSCom <- read.csv("AX206LocationCompendium.csv")
BeadsPutBack <- table(LOCSCom$GroupNumberLabel)
BeadCaptureRate <- data.frame(cbind(as.vector(BeadsPutBack), Beads[match(names(BeadsPutBack), A1TopRightBarcodeConversion)]))
rownames(BeadCaptureRate) <- names(BeadsPutBack)
colnames(BeadCaptureRate) <- c("Digital", "Physical")
#colors: red #bf, orange #f7, yellow #f9ed32, green 079247, light blue #c3e7ea, dark blue 1b75bc, purple 662d, pink da1c5c, light red d3757f

NineColScheme <- c("#BF1E2D", "#f79421", "#f9ed32", "#8dc63f", "#079247", "#27aae1", "#00a79d", "#1B75BC", "#662D91", "#da1c5c")


# SupFig Location tracking load("DigitalPhysicalBeads")
AX219BeadCapture <- TotalBeadCapture[TotalBeadCapture$Chip=="AX219",1:2]

AllBeadsCapture <- TotalBeadCapture[,1:2]


CountHeatmap <- data.frame(as.matrix(AX219BeadCapture %>% table)) %>% group_by(Digital, Physical)
CountHeatmap$Digital <- as.numeric(as.character(CountHeatmap$Digital))
CountHeatmap$Physical <- as.numeric(as.character(CountHeatmap$Physical))

ggplot(CountHeatmap, aes(x=Physical, y=Digital, color=Freq))+
  geom_point(size = 15)+
  scale_color_manual(values = c("#FFFFFF",NineColScheme))+
  scale_y_continuous(breaks=0:max(CountHeatmap$Digital),limits = c(0,max(CountHeatmap$Digital)))+
  scale_x_continuous(breaks = 0:max(CountHeatmap$Physical), limits = c(0,max(CountHeatmap$Physical)))+
  theme(text=element_text(family="Arial", size = 15))+
  coord_fixed(ratio=1)
ggplot(CountHeatmap, aes(x=Physical, y=Digital, color=Freq))+
  geom_point(size = 12)+
  scale_color_gradientn(colors = c("#FFFFFF",NineColScheme[1:5]))+
  scale_y_continuous(breaks=0:max(CountHeatmap$Digital),limits = c(0,max(CountHeatmap$Digital)))+
  scale_x_continuous(breaks = 0:max(CountHeatmap$Physical), limits = c(0,max(CountHeatmap$Physical)))+
  theme(text=element_text(family="Arial", size = 15))+
  coord_fixed(ratio=1)


# SupFig ReadDependence on beads vs cells
NoCellIncludedAverages <- group_by(NoCellIncludedMetadata, Cells, Beads)
Readaveragesbychambercontents <- summarize(NoCellIncludedAverages, averages = mean(TotalReads), beadaverages = mean(TotalReads/Beads))

ggplot(Readaveragesbychambercontents, aes(x=Beads, y=Cells, size=beadaverages, color = beadaverages))+
  geom_point()+
  scale_color_gradientn(colors = c("#FFFFFF",NineColScheme[1:5]))+
  theme(text=element_text(family="Arial", size = 15))+
  coord_fixed(ratio=1)+
  guides(color = guide_legend(), size = guide_legend())+
  scale_size_continuous(breaks = seq(0,3500,500))+
  scale_color_continuous(breaks = seq(0,3500,500))

#SupFig3 Bulk
ggplot(BulkProteinRatios[BulkProteinRatios$Chip=="Chip 2",]) +
  geom_col(aes(x=Protein, y=Value, group=Celltype, fill = Celltype), width = 0.5, position="dodge") +
  scale_fill_manual(values=c(NineColScheme[1],NineColScheme[6])) +
  theme(text = element_text(family = "Arial"), legend.position = c(0.1,0.9)) +
  coord_fixed(ratio = 1/2000) +
  labs(x="Protein", y="Fluorescence (a.u.)") +
  theme(text=element_text(family="Arial", size = 15))