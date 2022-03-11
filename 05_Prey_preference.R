
library(geomorph)
library(phytools)


# Feeding preferences
prey <- data.frame(readxl::read_xlsx("Prey Preference Acceptance.xlsx"))
prey[,3:7] <- log(prey[,3:7]/100 + 0.5)
rownames(prey) <- prey[,1]
prey <- prey[,3:7]

# Shape data
load("Dysdera_data.bin")
spec.selq1 <- read.table("Dysdera_data_names.txt", header = T, sep = " ")

gr.means <- lapply(gr.data2, mshape)
sh <- simplify2array(gr.means)
plotAllSpecimens(sh)

# Phylo
tr <- read.nexus("Dysdtree.tre")


# Match datasets
dim(sh)
dim(prey)
length(tr$tip.label)

all(dimnames(sh)[[3]]%in%tr$tip.label)
all(rownames(prey)%in%tr$tip.label)
rownames(prey)[14] <- "lea_t"
all(rownames(prey)%in%tr$tip.label)
all(rownames(prey)%in%dimnames(sh)[[3]])

drop.taxa <- dimnames(sh)[[3]][dimnames(sh)[[3]]%in%rownames(prey)==F]
tr.red <- drop.tip(tr, drop.taxa)
sh.red <- sh[,,tr.red$tip.label]
prey <- prey[tr.red$tip.label,]

# Phylogenetic PLS
pls.prey <- phylo.integration(sh.red, prey, tr.red)
summary(pls.prey)
plot(pls.prey)
text(pls.prey$XScores[,1], pls.prey$YScores[,1], labels = rownames(prey), pos = 2)


dt <- read.table("Prey_Pref_Detailed by exxADR.txt", header = T, sep ="\t")
dt <- dt[-which(rowSums(dt[,-1])==0),]
sp <- as.factor(dt$sp)
prey <- dt[,-1]

# Individual proportion of iso vs. non-iso prey consumed
iso.ind <- log(apply(prey, 1, function(x){sum(x[2:3])}) / rowSums(prey) + 0.5)
iso.ind <- apply(prey, 1, function(x){sum(x[2:3])}) / rowSums(prey) 
library(RRPP)
anova(lm.rrpp(iso.ind ~ sp))

isM <- as.data.frame(iso.ind)
isM$sp <- dt$sp

sp.means <- sort(tapply(iso.ind, sp, mean))
hist(sp.means)
median(sp.means)

spM <- as.data.frame(sp.means)
spM$isoG <- c("Generalist", "Generalist", "Generalist", "Generalist", "Generalist", "Generalist", 
              "Specialist", "Specialist", "Specialist", "Specialist", "Specialist", "Specialist", 
              "Specialist", "Specialist")


boxplot(sp.means ~ isoG, data = spM, 
        varwidth = TRUE,  las = 1, col = c("red", "blue"), 
        main = "Mean Prey Acceptance")

kmeans(sp.means, 2, iter.max = 10, nstart = 1,
       trace=FALSE)

# Proportion of iso consumed by species
isoN <- apply(prey[,2:3], 1, sum)
specialists <- as.factor(ifelse(sp%in%c("bre_t", "cri_t", "gom_g", "ins_t", "lev_t",
                                        "mac_t", "ram_g", "til_c"), 1, 0))
dtN <- data.frame(isoN = isoN, sp = sp, spec = specialists)
anova(glm(isoN ~ spec,
          data = dtN, family = poisson), test = "LRT")

ano <- anova(glm(isoN ~ spec,
                 data = dtN, family = poisson), test = "LRT")


