library(gplots)
library(ggdendro)
OldGenes <- rownames(OldGeneData)
NewGenes <- rownames(NewGeneData)

ReplicateGenes <- OldGeneData
colnames(ReplicateGenes) <- paste0(colnames(ReplicateGenes), ": 1")

ReplicateGenes[paste0(colnames(NewGeneData),": 2")] <- as.integer(0)
colnames(ReplicateGenes) <- gsub("-."," ",colnames(ReplicateGenes))
OldNewMatch <- sapply(NewGenes, function(x) match(x, OldGenes))

ReplicateGenes[which(rownames(ReplicateGenes) %in% names(OldNewMatch[!is.na(OldNewMatch)==TRUE])),(ncol(OldGeneData)+1):(ncol(OldGeneData)+ncol(NewGeneData))] <- NewGeneData[!is.na(OldNewMatch),]

ReplicateGenes[(nrow(ReplicateGenes)+1):(nrow(ReplicateGenes)+nrow(NewGeneData[is.na(OldNewMatch),])),] <- 0

ReplicateGenes[(nrow(ReplicateGenes)-nrow(NewGeneData[is.na(OldNewMatch),])+1):nrow(ReplicateGenes),(ncol(OldGeneData)+1):(ncol(OldGeneData)+ncol(NewGeneData))] <- NewGeneData[is.na(OldNewMatch),]

rownames(ReplicateGenes)[(nrow(ReplicateGenes)-nrow(NewGeneData[is.na(OldNewMatch),])+1):nrow(ReplicateGenes)] <- rownames(NewGeneData[is.na(OldNewMatch),])

VarRep <- apply(ReplicateGenes,1,var)

MaxRepGenes <- apply(ReplicateGenes,1,max)>5 

RepGenesNormalized <- ReplicateGenes / apply(ReplicateGenes,2,sum)[col(ReplicateGenes)]*10000

LogRepGenes <- log(ReplicateGenes+1)

SDNormLogRepGenes <- LogRepGenes/apply(LogRepGenes,1,sd)

RepSelectGenes <- SDNormLogRepGenes[MaxRepGenes,]

hc <- hclust(dist(t(RepSelectGenes)))
dend <- as.dendrogram(hc)
dend_data <- dendro_data(dend, type = "rectangle")

p <- ggplot(dend_data$segments) +
geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size=1, lineend = "square") +
geom_text(family="Calibri", data = dend_data$labels, aes(x, y, label = label), hjust = 1.2, angle = 45, size = 4) + 
ylim(-8, 85) + xlim(-1,41) + labs(title="Clustering after resampling") +
theme(text=element_text(family="Calibri"), axis.title=element_blank(), axis.text=element_blank(), axis.line = element_blank(), axis.ticks = element_blank())

pinset <- ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size=2, lineend = "square") +
  geom_text(family="Calibri", fontface = "bold", data = dend_data$labels, aes(x=x+0.2, y, label = label), hjust = 1.2, angle = 45, size = 7) + 
  ylim(-8, 33) + xlim(25.3,39.2) + labs(title="Clustering after resampling") +
  theme(text=element_text(family="Calibri"), axis.title=element_blank(), axis.text=element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
# heatmap.2(as.matrix(RepSelectGenes),trace="none")

# plot(hclust(dist(t(RepSelectGenes))), hang=-1)