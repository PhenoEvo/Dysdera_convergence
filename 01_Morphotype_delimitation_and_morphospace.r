
library(geomorph)
library(plyr)
library(ggplot2)

load("Dysdera_data.bin")
spec.selq1 <- read.table("Dysdera_data_names.txt", header = T, sep = " ")

gr.means <- lapply(gr.data2, mshape)
gr.means <- simplify2array(gr.means)


pca.gr <- gm.prcomp(gr.means)

#### Cluster methods test ####
D.gr <- dist(pca.gr$x)
methods <- c("ward.D", "ward.D2", "single", "complete", "average" , "mcquitty" , "median", "centroid")
correlations <- vector(length = 8)
for (i in 1:length(methods)) {
  hc <- hclust(D.gr, method = methods[i])
  correlations[i] <- cor(x = D.gr, cophenetic(hc))
  
}

names(correlations) <- methods
correlations
avr<-as.dendrogram(hclust(D.gr, method = "average"))


#### Anova Pairwise comparation between species & Headmap ####
sh <- array(unlist(gr.data2),dim=c(51,2,398))
spp <- vector()
for (i in 1:length(gr.data2)) {
  x <- dim(gr.data2[[i]])
  y <- rep.int(names(gr.data2[i]), x[[3]])
  spp<- append (spp, y)

}


dfc03 <- geomorph.data.frame(sh = sh, sp = spp)
prsh<-procD.lm(sh ~ sp, data = dfc03)
anova(prsh)

gp <- interaction(spp)

PW <- pairwise(prsh, groups = gp, covariate = NULL)
sum <- summary(PW, test.type = "dist", confidence = 0.95, stat.table = FALSE)

pvalues <- sum$pairwise.tables$P

pvalues[pvalues == 1.000] <- NA
pvalues[pvalues > 0.05] <- 1
pvalues[pvalues <= 0.05] <- 0


avr2<-rev(avr)
pvalues <- as.matrix(pvalues)
pvaluesLOG<-log(pvalues)
heatmap(pvaluesLOG, Rowv = avr, Colv=avr2, scale='none', col = "black")



#### Morphospace ####

#Load the file with each cheliceral morphotype defined from the Heatmap
clustD <- read.table("spp_clusters.txt", header = T, sep = "\t")
clustD[,c(3:4)] <- pca.gr$x[,1:2]

clustD <- clustD[,c(1,3,4,2)]
colnames(clustD) <- c( "sp", "pc1", "pc2", "type")



clustDM <- aggregate(clustD[, 2:3], list(clustD$type), mean)
clustD$mean1 <- clustD$type
clustD$mean2 <- clustD$type
clustD$mean1 <- revalue(clustD$mean1, 
                        c("A"=-0.040246048, "B"=0.022587225, "C"=-0.005136988,
                          "D"=-0.053226800, "E"=0.005683485	, "F"=-0.083433771,
                          "G"=0.083417905, "H"=0.081280612, "I"=0.112884318))
clustD$mean2 <- revalue(clustD$mean2, 
                        c("A"=-0.0189498363, "B"=-0.0087552859, "C"=-0.0316554577,
                          "D"=0.0489349712, "E"=0.0308343204	, "F"=0.0137312348,
                          "G"=0.0065201104, "H"=0.0362768092, "I"=0.0097695400))

find_hull <- function(clustD) clustD[chull(clustD$pc1, clustD$pc2), ]
hulls <- ddply(clustD, "type", find_hull)


cols <- c("red", "darkgoldenrod4", "green", "blue2", "darkcyan", "cyan3", "deeppink3", "purple", "pink")
plot <- ggplot(data = clustD, aes(x = pc1, y = pc2,  fill = type)) +
  geom_point()+
  geom_point(aes(x=-0.040246048, y=-0.0189498363), colour="red", size = 4)+ 
  geom_point(aes(x=0.022587225, y=-0.0087552859), colour="brown", size = 4 )+
  geom_point(aes(x=-0.005136988, y=-0.0316554577), colour="green", size = 4 )+
  geom_point(aes(x=-0.053226800, y=0.0489349712), colour="blue2", size = 4 )+
  geom_point(aes(x=0.005683485, y=0.0308343204), colour="deepskyblue3", size = 4 )+
  geom_point(aes(x=-0.083433771, y=0.0137312348), colour="purple", size = 4 )+
  geom_point(aes(x=0.083417905, y=0.0065201104), colour="palevioletred4", size = 4 )+
  geom_point(aes(x=0.081280612, y=0.0362768092), colour="yellow", size = 4 )+
  geom_point(aes(x=0.112884318, y=-0.0001073073), colour="lightsalmon", size = 4 )+
  geom_polygon(data = hulls, alpha = 0.6) +
  labs(x = "PC1", y = "PC2")

plot +  scale_fill_manual(breaks = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), 
                          values=c("red", "darkgoldenrod4", "green", "blue2", "darkcyan", "purple", "deeppink3", "yellow", "lightsalmon"))
