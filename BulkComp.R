HEKBulk <- read.csv("GSE79133_NS20160217_dge_ed1.txt", sep = "", row.names = 1)
colnames(HEKBulk) <- c("HEK1", "HEK2", "HEK3")
U87Bulk <- read.csv("GSM2794663_U87_con_1_Genes_ReadCount.txt", sep = "\t", row.names=1)
colnames(U87Bulk) <- "U871"
U87Bulk1 <- read.csv("GSM2794664_U87_con_2_Genes_ReadCount.txt", sep = "\t", row.names=1)
colnames(U87Bulk1) <- "U872"
HekBulk1 <- read.csv("GSM2486332_HEK_Population.txt", sep = "", row.names=1)
HekBulk1 <- subset(HekBulk1,select=-genes)
# GSE89164
colnames(HekBulk1) <- "HEK4"
HEKBulk2 <- read.csv("gene_counts_hek.txt", sep = "", row.names = 1)
colnames(HEKBulk2) <- "HEK5"
U87Bulk2 <- read.csv("GSM2333485_U87Control1_raw_count.txt", sep = "\t")
U87Bulk2 <- U87Bulk2[-which(duplicated(U87Bulk2$external_gene_id)==TRUE),1:7]
rownames(U87Bulk2) <- U87Bulk2$external_gene_id
U87Bulk2 <- subset(U87Bulk2, select = U87Control1_raw_count)
colnames(U87Bulk2) <- "U873"

GenesforBulk <- rownames(CombinedGenesbyMerge@raw.data)

BulkDatasets <- data.frame(row.names=GenesforBulk)
BulkDatasets[colnames(HEKBulk)] <- HEKBulk[GenesforBulk,]
BulkDatasets[colnames(HekBulk1)] <- HekBulk1[GenesforBulk,]
BulkDatasets[colnames(HEKBulk2)] <- HEKBulk2[GenesforBulk,]
BulkDatasets[colnames(U87Bulk)] <- U87Bulk[GenesforBulk,]
BulkDatasets[colnames(U87Bulk1)] <- U87Bulk1[GenesforBulk,]
BulkDatasets[colnames(U87Bulk2)] <- U87Bulk2[GenesforBulk,]
BulkDatasets[is.na(BulkDatasets)] <- 0

# RepMetadata <- cbind(c("HEK", "HEK", "HEK", "HEK", "HEK", "U87", "U87", "U87"), c(1:8))
# colnames(RepMetadata) <- c("celltype", "rep")
# rownames(RepMetadata) <- c("HEK1", "HEK2", "HEK3", "HEK4", "HEK5", "U871", "U872", "U873")

# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
# library("DESeq2")
# dds <- DESeqDataSetFromMatrix(countData = BulkDatasets,
#                               colData = RepMetadata,
#                               design = ~ celltype)
# dds <- DESeq(dds)
# res <- results(dds)
# resLFC <- lfcShrink(dds, coef="celltype_U87_vs_HEK", type="apeglm")
# resOrdered <- res[order(res$pvalue),]
# res05 <- results(dds, alpha=0.05)

BulkS <- CreateSeuratObject(raw.data=BulkDatasets, project = "Bulks")
BulkS <- NormalizeData(BulkS, display.progress = F)
BulkS <- ScaleData(BulkS, display.progress = F)

BulkS <- FindVariableGenes(BulkS, x.low.cutoff = 0.3, x.high.cutoff = 8, y.cutoff = 1.5, do.plot = F, display.progress = F)


TestBulkvar <- BulkS@var.genes
CombinedGenesbyMerge@var.genes <- TestBulkvar

# [1] "PC1"
# [1] "RPL39"         "RP11-466H18.1" "RPS10"         "RPS29"         "RPS17"         "RPL26"        
# [7] "RP6-24A23.7"   "RPL27"         "NDUFA13"       "MT-ATP8"       "NDUFC2"        "RPS12"        
# [13] "LAMB1"         "TRIM33"        "PRKAG2"        "GJA1"          "SH3YL1"        "NDUFA11"      
# [19] "TRAM2"         "NAV1"          "COL5A1"        "SEZ6L2"        "KMT2B"         "SLC9A3R1"     
# [25] "ITGA2"         "FOSL2"         "AFF3"          "NCS1"          "HTRA1"         "GFPT2"        
# [1] ""
# [1] "S100A6"      "TMSB4X"      "LGALS1"      "VIM"         "GNG11"       "ANXA1"       "TFPI2"      
# [8] "ANXA2"       "SH3BGRL3"    "CALU"        "DCBLD2"      "RPL22L1"     "CAV1"        "MAP1B"      
# [15] "S100A16"     "CALD1"       "CAPN2"       "PTRF"        "DST"         "EGFR"        "TM4SF1"     
# [22] "CD44"        "MINOS1"      "HSPA5"       "LRRC75A-AS1" "TLN1"        "STC1"        "NUPR1"      
# [29] "MLPH"        "LMNA"   

BulkProteinRatios <- data.frame("Protein" = c("PKM2", "c-MYC", "PDHK1","PKM2", "c-MYC", "PDHK1", "PKM2", "c-MYC", "PDHK1","PKM2", "c-MYC", "PDHK1"),
                                "Chip" = c("Chip 1","Chip 1","Chip 1","Chip 2","Chip 2","Chip 2","Chip 1","Chip 1","Chip 1","Chip 2","Chip 2","Chip 2"),
                                "Value" = c(178.9583333, 477.5, 883.5416667, 1104.270833, 3523.1875, 1242, 192.3333333, 1741.770833, 358.75, 1114.9375, 5498.8125, 439.375),
                                "Celltype" = c("U87","U87","U87","U87","U87","U87","HEK","HEK","HEK","HEK","HEK","HEK"))
BulkValues <- read.csv("BulkChipValues.txt", sep = "", header = FALSE)
Bulk1 <- BulkValues[8:155, 2:4]
Bulk1["celltype"]="U87"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="1"
TotalBulkDataFrame <- Bulk1
BulkValues <- BulkValues[-(1:155),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="U87"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="2"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="U87"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="3"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[9:156, 2:4]
Bulk1["celltype"]="HEK"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="1"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:156),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="HEK"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="2"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[8:155, 2:4]
Bulk1["celltype"]="Blank"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="1"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:155),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="Blank"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="2"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="Blank"
Bulk1["chip"]="AX210"
Bulk1["replicate"]="3"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[10:157, 2:4]
Bulk1["celltype"]="U87"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="1"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:157),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="U87"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="2"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="U87"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="3"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[8:155, 2:4]
Bulk1["celltype"]="HEK"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="1"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:155),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="HEK"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="2"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[6:153, 2:4]
Bulk1["celltype"]="HEK"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="3"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:153),]
Bulk1 <- BulkValues[8:155, 2:4]
Bulk1["celltype"]="Blank"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="1"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
BulkValues <- BulkValues[-(1:155),]
Bulk1 <- BulkValues[7:154, 2:4]
Bulk1["celltype"]="Blank"
Bulk1["chip"]="AX209"
Bulk1["replicate"]="2"
TotalBulkDataFrame <- cbind(TotalBulkDataFrame, Bulk1)
# BulkValues <- BulkValues[-(1:154),]
# Bulk1 <- BulkValues[10:157, 2:4]
# Bulk1["celltype"]="U87"
# Bulk1["chip"]="AX201"
# Bulk1["replicate"]="1"
# TotalBulkDataFrame <- rbind(TotalBulkDataFrame, Bulk1)
# BulkValues <- BulkValues[-(1:157),]
# Bulk1 <- BulkValues[6:153, 2:4]
# Bulk1["celltype"]="U87"
# Bulk1["chip"]="AX201"
# Bulk1["replicate"]="2"
# TotalBulkDataFrame <- rbind(TotalBulkDataFrame, Bulk1)
# BulkValues <- BulkValues[-(1:153),]
# Bulk1 <- BulkValues[6:153, 2:4]
# Bulk1["celltype"]="U87"
# Bulk1["chip"]="AX201"
# Bulk1["replicate"]="3"
# TotalBulkDataFrame <- rbind(TotalBulkDataFrame, Bulk1)
# BulkValues <- BulkValues[-(1:153),]