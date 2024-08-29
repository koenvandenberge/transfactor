library(matrixStats)
ksource <- function(x, ...) {
  library(knitr)
  source(purl(x, output = tempfile()), ...)
}
library(tidyverse)

simulateFromGRN_extended <- function(X,
                                     G=500,
                                     nTF=80,
                                     fracDETFs=0.25,
                                     nPerGroup=50,
                                     nCellTypes=5,
                                     shape_mu_t=1,
                                     scale_mu_t=2,
                                     shape_lambda_gc=1,
                                     scale_lambda_gc=scale_lambda_gc,
                                     magnitude_g=rep(100, G),
                                     seed=3,
                                     dirichletNoise = FALSE){

  set.seed(seed)

  ### Simulate mu_tc
  mu_tc <- matrix(rgamma(n=nTF, shape=shape_mu_t, scale=scale_mu_t),
                  nrow=nTF, ncol=nCellTypes, byrow=FALSE)
  ### add noise
  # mu_tc_noise <- mu_tc + matrix(rnorm(n=length(mu_tc), mean=mu_tc, sd=.5),
  #                               nrow=nTF, ncol=nCellTypes)
  # mu_tc_noise[mu_tc_noise <= 0] <- 1e-5
  mu_tc_noise <- mu_tc
  rownames(mu_tc) <- paste0("tf",1:nTF)
  rownames(mu_tc_noise) <- paste0("tf",1:nTF)

  ### simulate DE
  deTFs <- sample(nTF, size=nTF * fracDETFs)
  # simulate fold changes
  logFC <- log(rnorm(n=length(deTFs) * (nCellTypes-1), mean=2.5, sd=.3))
  # random sign
  x <- rbinom(n=length(deTFs) * (nCellTypes-1), size=1, prob=.5)
  # logFC <- log(rep(5, length(deTFs) * (nCellTypes-1)))
  # x <- rep(1, length(deTFs) * (nCellTypes-1))
  fcAll <- matrix(1, nrow = nTF, ncol = nCellTypes - 1)
  fcAll[deTFs, ] <- matrix(exp(ifelse(x == 1, 1, -1) * logFC),
                           nrow = length(deTFs), ncol = nCellTypes - 1)
  mu_tc_noise[, -1] <- mu_tc_noise[, -1] * fcAll

  ### mu_gc
  mu_gc_sim <- X %*% mu_tc_noise

  ### calculate mu_gtc, mu_gc
  alpha_gt <- list()
  for(gg in 1:G){
    id <- which(X[gg,] == 1)
    mu_tc_id <- mu_tc_noise[id,,drop=FALSE]
    curDimNames <- list(paste0(colnames(X)[id], ";gene",gg),
                        paste0("celltype",letters[1:nCellTypes]))

    curalpha_gt <- rowMeans(sweep(mu_tc_id, 2, colSums(mu_tc_id), "/")) * magnitude_g[gg]
    alpha_gt[[gg]] <- curalpha_gt
    curmu_gtc <- sweep(sweep(mu_tc_id, 2, colSums(mu_tc_id), "/"), 2, mu_gc_sim[gg,], "*")
    dimnames(curmu_gtc) <- curDimNames
    curmu_gc <- colSums(curmu_gtc)
    stopifnot(identical(round(unname(curmu_gc),2), round(mu_gc_sim[gg,],2)))
    if(!dirichletNoise){
      # noiseless:
      curpi_gtc <- sweep(curmu_gtc, 2, curmu_gc, "/")
    } else {
      # dirichlet noise:
      curpi_gtc <- rdirichlet(n=nCellTypes, alpha=curalpha_gt)
      if(length(id) > 1) curpi_gtc <- t(curpi_gtc)
      dimnames(curpi_gtc) <- curDimNames
    }
    if(gg == 1){
      mu_gtc <- curmu_gtc
      mu_gc <- curmu_gc
      pi_gtc <- curpi_gtc
    } else {
      mu_gtc <- rbind(mu_gtc, curmu_gtc)
      mu_gc <- rbind(mu_gc, curmu_gc)
      pi_gtc <- rbind(pi_gtc, curpi_gtc)
    }
  }
  rownames(mu_gc) <- paste0("gene",1:G)
  names(alpha_gt) <- paste0("gene",1:G)

  for(cc in 1:nCellTypes){
    curY <- matrix(rpois(n= nPerGroup * G, lambda=mu_gc_sim[,cc]),
                   nrow=nrow(mu_gc_sim), ncol= nPerGroup, byrow=FALSE)
    if(cc == 1){
      Y_igc <- curY
    } else {
      Y_igc <- cbind(Y_igc, curY)
    }
  }
  rownames(Y_igc) <- paste0("gene",1:G)

  return(list(Y=Y_igc,
              pi_gtc=pi_gtc,
              mu_tc=mu_tc,
              mu_tc_noise=mu_tc_noise,
              mu_gtc=mu_gtc,
              mu_gc=mu_gc,
              alpha_gt=alpha_gt,
              fcAll = fcAll,
              deTFs = deTFs
  ))
}

set.seed(123)
logit = function(x) log(x/(1-x))
expit = function(x) exp(x)/(1+exp(x))
set.seed(13)
## simulation parameters
G <- 500 #nr of genes
nTF <- 80 #nr of TFs
nPerGroup <- 50 #number of cells per cell type
shape_mu_t <- 1
scale_mu_t <- 2
shape_lambda_gc <- 1
scale_lambda_gc <- 1
nCellTypes <- 5
ct <- factor(rep(letters[1:nCellTypes], each=nPerGroup))
design <- model.matrix(~-1 + ct)
seed <- 3
magnitude_g <- rep(100,G)

## simulate GRN as TF regulation matrix
X <- matrix(0, nrow=G, ncol=nTF)
for(gg in 1:nrow(X)){
  nGenes <- rbinom(n=1, size=30, prob=.05)
  while(nGenes == 0) nGenes <- rbinom(n=1, size=30, prob=.05)
  id <- sample(ncol(X), nGenes)
  X[cbind(gg,id)] <- 1
}
rownames(X) <- paste0("gene",1:G)
colnames(X) <- paste0("tf",1:nTF)

simRes <- simulateFromGRN_extended(X,
                                   G=G,
                                   nTF=nTF,
                                   nPerGroup=nPerGroup,
                                   shape_mu_t=shape_mu_t,
                                   scale_mu_t=scale_mu_t,
                                   shape_lambda_gc=shape_lambda_gc,
                                   scale_lambda_gc=scale_lambda_gc,
                                   nCellTypes=nCellTypes,
                                   seed=3,
                                   magnitude_g=magnitude_g)


counts <- simRes$Y
ct <- factor(rep(letters[1:nCellTypes], each=nPerGroup))
U <- model.matrix(~ -1 + ct)

