barcodes <- c("B", "C", "D")
blocks <- unique(AX179T$Block)
output <- array(dim = c(length(blocks)*length(barcodes), ncol(AX179T)-2))
A1TopRight <- c(36,25,24,13,12,1,35,26,23,14,11,2,34,27,22,15,10,3,33,28,21,16,9,4,32,29,20,17,8,5,31,30,19,18,7,6)

numvalues <- length(colnames(AX179T)[-c(1:3,length(colnames(AX179T)))])

values <- AX179T[colnames(AX179T)[-c(1:3,length(colnames(AX179T)))]]
i=1
for (n in 1:length(blocks))
{
  for (m in 1:length(barcodes))
  {
    output[i,1] <- blocks[n]
    output[i,2] <- barcodes[m]
    output[i,3:18] <- colMeans(values[AX179T$Block==blocks[n] & AX179T$Barcode==barcodes[m],])
    i=i+1
  }
}

colnames(output) <- colnames(AX179T)[c(-2, -ncol(AX179T))]
TrueBlock <- A1TopRight[as.double(as.character(output1$Block))]

output <- data.frame(output)
