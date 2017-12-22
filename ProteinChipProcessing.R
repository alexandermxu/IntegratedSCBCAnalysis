library(xlsx)
library(tidyverse)
File <- read.xlsx("AX194_Adjusted.xlsx", sheetIndex = 2)
Cells <- c(1,0,0,0,0,0,2,0,1,0,0,0,
           0,1,3,1,2,2,2,1,2,0,4,3,
           0,0,0,1,3,1,0,0,0,0,0,0)
Barcodes <- c("B", "C", "D")
Blocks <- unique(File$Block)
Group1 <- c(1:6, 13:18, 25:30)
Group2 <- c(7:12, 19:24, 31:36)
Registry1 <- c(NA,NA,NA,NA,"B","C","D",NA,NA,NA,NA)
Registry2 <- c(NA,NA,NA,NA,NA,"D","C","B",NA,NA,NA)
A1TopRight <- c(36,25,24,13,12,1,35,26,23,14,11,2,34,27,22,15,10,3,33,28,21,16,9,4,32,29,20,17,8,5,31,30,19,18,7,6)
KeepData <- c("Block", "Row", "F635.Mean", "F532.Mean")
Entries <- length(KeepData)
Output <- array(dim = c(length(Blocks)*length(Barcodes), Entries+2))
FileExtract <- File[,KeepData]
FileExtract[FileExtract$Block %in% Group1,Entries+1] <- Registry1[FileExtract[FileExtract$Block %in% Group1,]$Row]
FileExtract[FileExtract$Block %in% Group2,Entries+1] <- Registry2[FileExtract[FileExtract$Block %in% Group2,]$Row]
ReducedFileExtract <- na.omit(FileExtract)
colnames(ReducedFileExtract) <- c("Block", "Row", "Mean635", "Mean532", "Barcode")
Values <- ReducedFileExtract[,-c(1:2,Entries+1)]
i=1
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
colnames(Output) <- c("Block", "Barcode", "Mean635", "Mean532", "Chamber", "Cells")
# colnames(Output) <- c("Block", "Barcode", "Mean635", "Mean594", "Mean532", "Mean488", "Chamber", "Cells")

Output <- data.frame(Output)
for (l in c(1,3:(Entries+2)))
{Output[,l] <- as.numeric(as.character(Output[,l]))}

Plot1 <- ggplot(Output, aes(x=Barcode,y=Mean635,label=Block)) + geom_jitter(aes(color=factor(Cells)))
Plot1
