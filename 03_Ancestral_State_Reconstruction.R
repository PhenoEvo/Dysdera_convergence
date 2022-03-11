
library(phytools)
library(plyr)

tr <- read.nexus("Dysdtree_1001.tre")
clustD <- read.table("spp_clusters.txt", header = T, sep = "\t")

rownames(clustD)<-clustD$sp
clust_types <- clustD[,2]
clust_types<-as.matrix(clustD)[,2]

#### Model ER with the irreversibility constriction ####
ERC<-matrix(c(0,1,1,1,1,1,1,1,1, 
               0,0,1,0,1,1,1,1,1,
               0,1,0,0,1,1,1,1,1,
               1,1,1,0,1,1,1,1,1,
               0,1,1,0,0,1,1,1,1,
               0,1,1,0,1,0,1,1,1,
               0,1,1,0,1,1,0,1,1,
               0,1,1,0,1,1,1,0,1,
               0,1,1,0,1,1,1,1,0
), ncol = 9, byrow = T)

#### Simmap for each model (symmetric (SYM), all rates different (ARD), equal rates (ER) and the constrained ER (ERC)) ####
# for the 1001 tree from the posterior distribution #
mSYM <- make.simmap(tr, clust_types, model = "SYM", nsim=100)
mARD <- make.simmap(tr, clust_types, model = "ARD", nsim=100)
mER <- make.simmap(tr, clust_types, model = "ER", nsim=100)
mERC <- make.simmap(tr, clust_types, model = ERC, nsim=10)

objSYM <-describe.simmap(mSYM)
objARD <-describe.simmap(mARD)
objER <-describe.simmap(mER)
objERC<-describe.simmap(mERC)


colors<-setNames(c("red", "darkgoldenrod4", "green", "blue2", "darkcyan", "purple", "deeppink3", "yellow", "lightsalmon"),
                 c("A", "B", "C", "D", "E", "F", "G", "H", "I"))

plot(objSYM,colors=colors)
add.simmap.legend(colors = colors, vertical = T, prompt = FALSE, x = 82, y = 24)

plot(objARD,colors=colors)
add.simmap.legend(colors = colors, vertical = T, prompt = FALSE, x = 82, y = 24)

plot(objER,colors=colors)
add.simmap.legend(colors = colors, vertical = T, prompt = FALSE, x = 82, y = 24)

plot(objERC,colors=colors)
add.simmap.legend(colors = colors, vertical = T, prompt = FALSE, x = 82, y = 24)


#### AIC calculation ####
AIC_SYM <- matrix(NA, ncol = 1, nrow = length(mSYM))
for(i in 1:length(mSYM)){
  logT <- c(mSYM[[i]][["logL"]])
  AIC_SYM[i] <- logT
}
AIC_SYM <- as.data.frame(AIC_SYM)

AIC_ARD <- matrix(NA, ncol = 1, nrow = length(mARD))
for(i in 1:length(mARD)){
  logT <- c(mARD[[i]][["logL"]])
  AIC_ARD[i] <- logT
}
AIC_ARD <- as.data.frame(AIC_ARD)

AIC_ER <- matrix(NA, ncol = 1, nrow = length(mER))
for(i in 1:length(mER)){
  logT <- c(mER[[i]][["logL"]])
  AIC_ER[i] <- logT
}
AIC_ER <- as.data.frame(AIC_ER)

AIC_ERC <- matrix(NA, ncol = 1, nrow = length(mERC))
for(i in 1:length(mERC)){
  logT <- c(mERC[[i]][["logL"]])
  AIC_ERC[i] <- logT
}
AIC_ERC <- as.data.frame(AIC_ERC)


logL<-cbind(AIC_SYM, AIC_ARD, AIC_ER, AIC_ERC)
aic<-function(logL,k) 2*k-2*logL
aic.w<-function(aic){
  d.aic<-aic-min(aic)
  exp(-1/2*d.aic)/sum(exp(-1/2*d.aic))
}
AIC2<-mapply(aic,logL,c(36,72,1,1))
AIC.W<-aic.w(AIC2)

### Mean AIC and SD for each model
AIC2 <- as.data.frame(AIC2)
meansAICs <- matrix(NA, ncol = 3, nrow = 4)
meansAICs[,3] <- c("SYM","ARD","ER","ERC")
meansAICs[,1] <- c(mean(AIC2[,1]),mean(AIC2[,2]),mean(AIC2[,3]),mean(AIC2[,4]))
meansAICs[,2] <- c(sd(AIC2[,1]),sd(AIC2[,2]),sd(AIC2[,3]),sd(AIC2[,4]))
colnames(meansAICs) <- c("mean", "sd", "Model")
meansAICsM <- as.data.frame(meansAICs)
meansAICsM[,2] <- as.numeric(meansAICsM[,2])
meansAICsM[,1] <- as.numeric(meansAICsM[,1])

library(ggplot2)
ggplot(meansAICsM,aes(x=Model))+geom_boxplot(aes(lower=mean-sd,upper=mean+sd,middle=mean,ymin=mean-3*sd,ymax=mean+3*sd),stat="identity")
