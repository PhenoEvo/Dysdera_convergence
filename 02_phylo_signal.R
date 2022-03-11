
library(phytools)
library(geomorph)

tr <- read.nexus("Dysdtree_1001.tre")
load("Dysdera_data.bin")

#### Mean shape for each species ####
gr.means <- lapply(gr.data2, mshape)
gr.means <- simplify2array(gr.means)


#### Phylogenetic signal for each tree from the posterior distribution ####
physig.geo <- matrix(NA, ncol = 2, nrow = length(tr))
colnames(physig.geo) <- c("K", "pval")
for(i in 1:length(tr)){
  ps <- physignal(gr.means, tr[[i]], iter = 999, seed = NULL, print.progress = FALSE)
  physig.geo[i,] <- c(ps$phy.signal, ps$pvalue)
}




#### Formula to test if the K value is significant lower than 1 ####
testKmult <- function(A, phy, iter = 999, seed = NULL, distr = TRUE){
  if (any(is.na(A))) 
    stop("Data matrix contains missing values. Estimate these first (see 'estimate.missing').\n", 
         call. = FALSE)
  if (length(dim(A)) == 3) {
    if (is.null(dimnames(A)[[3]])) 
      stop("Data array does not include taxa names as dimnames for 3rd dimension.\n", 
           call. = FALSE)
    x <- two.d.array(A)
  }
  if (length(dim(A)) == 2) {
    if (is.null(rownames(A))) 
      stop("Data matrix does not include taxa names as dimnames for rows.\n", 
           call. = FALSE)
    x <- A
  }
  if (is.vector(A)) {
    if (is.null(names(A))) 
      stop("Data vector does not include taxa names as names.\n", 
           call. = FALSE)
    x <- as.matrix(A)
  }
  if (!inherits(phy, "phylo")) 
    stop("tree must be of class 'phylo.'\n", call. = FALSE)
  N <- length(phy$tip.label)
  if (N != dim(x)[1]) 
    stop("Number of taxa in data matrix and tree are not not equal.\n", 
         call. = FALSE)
  if (length(match(rownames(x), phy$tip.label)) != N) 
    stop("Data matrix missing some taxa present on the tree.\n", 
         call. = FALSE)
  if (length(match(phy$tip.label, rownames(x))) != N) 
    stop("Tree missing some taxa in the data matrix.\n", 
         call. = FALSE)
  if (any(is.na(match(sort(phy$tip.label), sort(rownames(x)))))) 
    stop("Names do not match between tree and data matrix.\n", 
         call. = FALSE)
  if (is.null(dim(x))) 
    x <- matrix(x, dimnames = list(names(x)))
  
  x <- as.matrix(x[phy$tip.label, ])
  
  # Calculate observed rates matrix (for a single group)
  gp <- rep(1, N) 
  x <- scale(x, center = TRUE, scale = FALSE)
  phy.parts <- geomorph:::phylo.mat(x, phy)
  invC <- phy.parts$invC
  D.mat <- phy.parts$D.mat
  C <- phy.parts$C
  sigma.obs <- geomorph:::sigma.d(x, invC, D.mat, gp)
  
  # Simulate multivariate data with the empirical phylo, dims and rate matrix
  rate.mat <- sigma.obs$R
  diag(rate.mat) <- sigma.obs$sigma.d.all
  rate.mat <- geomorph:::makePD(rate.mat)
  x.sim <- geomorph:::sim.char.BM(phy = phy, par = rate.mat, 
                                  nsim = iter, seed = seed)
  
  nullK <- lapply(x.sim, function(x) physignal(x, phy = phy, iter = 2)$phy.signal)
  nullK <- unlist(nullK)
  
  # Compare to observed 
  physig.obs <- physignal(A, phy = phy)
  Kmult <- physig.obs$phy.signal
  P <- mean(abs(log(c(Kmult, nullK))) > abs(log(Kmult)))
  temp <- c(Kmult, P); names(temp) <- c("Kmult", "p-val")
  
  res <- list(testK1 = temp, nullK = nullK)
  if(distr){
    hist(nullK)
    abline(v = Kmult, lwd = 2, col = "red")
  }
  return(res)
}


#### Kmult test ####
Kmult <- matrix(NA, ncol = 2, nrow = length(tr))
colnames(Kmult) <- c("Kmult", "p-value")
for (i in 1:length(tr)) {
  fr <- testKmult(gr.means, tr[[1]])
  Kmult[i,] <- c(fr$testK1[1], fr$testK1[2])
  
}


