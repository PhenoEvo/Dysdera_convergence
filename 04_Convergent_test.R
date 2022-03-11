
library(phytools)
library(convevol)
library(geomorph)

tr <- read.nexus("Dysdtree_1001.tre")
tr100<-sample(tr, 100)

load("Dysdera_data.bin")
clustD <- read.table("spp_clusters.txt", header = T, sep = "\t")

### Mean shape for species and extract the first two PC and join the results with the cheliceral type of each species
gr.means <- lapply(gr.data2, mshape)
gr.means <- simplify2array(gr.means)
pca.gr <- gm.prcomp(gr.means)
clustD$pc1 <- pca.gr$x[,1] ; clustD$pc2 <- pca.gr$x[,2]
rownames(clustD) <- clustD$sp


### Convergent test for each cheliceral type
gA <-clustD[clustD$type %in% "A",]
convtipsA<- as.character(gA$sp)

gB <-clustD[clustD$type %in% "B",]
convtipsB<- as.character(gB$sp)

gG <-clustD[clustD$type %in% "G",]
convtipsG<- as.character(gG$sp)

#Matrix with the different C statistic values (CsA/B/G), cutoffs (EsA/B/G) and the p-values (PvA/B/G)

CsA <-matrix(NA, ncol = 4, nrow = length(tr100))
colnames(CsA) <- c("C1", "C2", "C3", "C4")

PvA <-matrix(NA, ncol = 4, nrow = length(tr100))
colnames(PvA) <- c("C1", "C2", "C3", "C4")

EsA <- matrix(NA, ncol=4, nrow= length(tr100))
colnames(EsA) <- c("C1", "C2", "C3", "C4")

for (i in 1:length(tr100)) {
  is_tip <- tr100[[i]]$edge[,2] <= length(tr100[[i]]$tip.label)
  ordered_tips <- tr100[[i]]$edge[is_tip,2]
  lablesN<-tr100[[i]]$tip.label[ordered_tips]
  data_new1 <- clustD[match(lablesN, clustD$sp), ]
  data_new1 <- data_new1[,-c(1:2)]
  phD <- as.matrix(data_new1)
  
  answerD <- convratsig(tr100[[i]],phD,convtipsA,100)
  CsA[i,] <- answerD$ObservedCs
  PvA[i,] <- answerD$Pvals
  EsA[i,] <- answerD$Cutoffs
  
}


CsB <-matrix(NA, ncol = 4, nrow = length(tr100))
colnames(CsD) <- c("C1", "C2", "C3", "C4")

PvB <-matrix(NA, ncol = 4, nrow = length(tr100))
colnames(PvD) <- c("C1", "C2", "C3", "C4")

EsB <- matrix(NA, ncol=4, nrow= length(tr100))
colnames(EsB) <- c("C1", "C2", "C3", "C4")

for (i in 1:length(tr100)) {
  is_tip <- tr100[[i]]$edge[,2] <= length(tr100[[i]]$tip.label)
  ordered_tips <- tr100[[i]]$edge[is_tip,2]
  lablesN<-tr100[[i]]$tip.label[ordered_tips]
  data_new1 <- clustD[match(lablesN, clustD$sp), ]
  data_new1 <- data_new1[,-c(1:2)]
  phD <- as.matrix(data_new1)
  
  answerD <- convratsig(tr100[[i]],phD,convtipsB,100)
  CsB[i,] <- answerD$ObservedCs
  PvB[i,] <- answerD$Pvals
  EsB[i,] <- answerD$Cutoffs
  
}


CsG <-matrix(NA, ncol = 4, nrow = length(tr100))
colnames(CsG) <- c("C1", "C2", "C3", "C4")

PvG <-matrix(NA, ncol = 4, nrow = length(tr100))
colnames(PvG) <- c("C1", "C2", "C3", "C4")

EsG <- matrix(NA, ncol=4, nrow= length(tr100))
colnames(EsG) <- c("C1", "C2", "C3", "C4")

for (i in 1:length(tr100)) {
  is_tip <- tr100[[i]]$edge[,2] <= length(tr100[[i]]$tip.label)
  ordered_tips <- tr100[[i]]$edge[is_tip,2]
  lablesN<-tr100[[i]]$tip.label[ordered_tips]
  data_new1 <- clustD[match(lablesN, clustD$sp), ]
  data_new1 <- data_new1[,-c(1:2)]
  phD <- as.matrix(data_new1)
  
  answerD <- convratsig(tr100[[i]],phD,convtipsG,100)
  CsG[i,] <- answerD$ObservedCs
  PvG[i,] <- answerD$Pvals
  EsG[i,] <- answerD$Cutoffs
  
}
