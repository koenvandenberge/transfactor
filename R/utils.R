updateS <- function(sOld, mHat, p){
  # sOld is current precision estimate
  # mHat is current fraction estimate
  # p is the data: matrix of multinomial probabilities.
  # In p, each column is a class, each row is a data point.
  log_bar_p_k <- colMeans(log(p))
  K <- ncol(p)
  sNew <- (K-1)/((K-1)/sOld - digamma(sOld) + sum(mHat*digamma(sOld*mHat)) - sum(mHat*log_bar_p_k))
  return(sNew)
}


updateStepOLS <- function(mu_tc, mu_gc, X, n_c, Y_gc){
  mu_gc_scaled <- mu_gc/((X %*% mu_tc)+1e-10) # OK
  olsX <- t(n_c * (t(X) %*% diag(mu_gc_scaled[,1]))) # this is a TxG matrix with mu_gc_scaled as element whenever X=1.
  betaHat <- MASS::ginv(t(olsX) %*% olsX) %*% t(olsX) %*% Y_gc
  return(betaHat)
}


sparseInitialization_sufStats <- function(counts_suf, design,
                                          X, iterOLS=0,
                                          lassoFamily = "gaussian"){


  ## initialize mu_tc
  EZ_probOrig <- X / rowSums(X)
  ### mu_tc initialization: note that this is the mean across all genes a TF is regulating,
  ### assuming equal partitioning of reads among regulating TFs.
  mu_tc <- matrix(NA, nrow=ncol(X), ncol=ncol(counts_suf))
  for(tt in 1:ncol(X)){
    weights <- EZ_probOrig[,tt]
    hlp <- weights * counts_suf
    mu_tc[tt,] <- colMeans(hlp[weights > 0,,drop=FALSE])
  }
  rownames(mu_tc) <- colnames(X)

  betaMu <- matrix(NA, nrow=ncol(X), ncol=ncol(design))
  for(cc in 1:ncol(design)){
    cellId <- which(design[,cc] == 1)
    n_c <- sum(design[,cc])
    curMu_tc <- mu_tc[,cc]
    mu_gc <- counts_suf[,cc] / n_c
    Y_gc <- counts_suf[,cc]

    ## OLS update: really only useful if mu_tc initialization is bad (random).
    ## Defaults to zero iterations.
    ctrOLS <- 0
    while(ctrOLS < iterOLS){
      ctrOLS <- ctrOLS+1
      curMu_tc <- updateStepOLS(mu_tc = curMu_tc,
                                mu_gc = mu_gc,
                                X = X,
                                n_c = n_c,
                                Y_gc = Y_gc)
    }

    ## lasso
    mu_gc_scaled <- mu_gc/((X %*% curMu_tc)+1e-10)
    olsX <- t(n_c * (t(X) %*% diag(mu_gc_scaled[,1])))
    if(lassoFamily == "gaussian"){
      cvfit <- glmnet::cv.glmnet(olsX, Y_gc, lower.limits=0, intercept=FALSE, alpha=1,
                         family="gaussian", weights=1/(Y_gc+1), standardize=FALSE)
    } else if(lassoFamily == "poisson"){
      cvfit <- glmnet::cv.glmnet(olsX, Y_gc, lower.limits=0, intercept=FALSE, alpha=1,
                         family="poisson", standardize=FALSE)
    }
    betaHatLasso <- coef(cvfit, s = "lambda.min")[-1,]
    betaMu[,cc] <- betaHatLasso
  }
  dimnames(betaMu) <- list(colnames(X), colnames(design))
  return(betaMu)
}

thinningFactor <- function(counts,
                            X,
                            dfRepr,
                            U){

  ## Deleted checks on X rownames and counts rownames
  ## Deleted sorting of counts based on X
  repressedTFs <- as.character(unique(dfRepr$repressed))
  if(length(repressedTFs) == 0){
    message("No between-TF repressions present in the GRN.\n")
    rho_tc <- NULL
    return(rho_tc)
  }

  # design of experiment
  C <- ncol(U)
  lvl <- unlist(apply(U,1, function(row){
    which(row == 1)
  }))
  rho_tc <- matrix(NA, nrow=length(repressedTFs),
                   ncol=ncol(U))
  dimnames(rho_tc) <- list(repressedTFs, colnames(U))

  for(rr in 1:length(repressedTFs)){
    curTF <- as.character(repressedTFs[rr])
    repressingTFs <- as.character(dfRepr$repressor[dfRepr$repressed == curTF])
    repressingTFs <- repressingTFs[repressingTFs %in% rownames(counts)]
    if(length(repressingTFs) == 0) rho_tc[rr,] <- 1
    yRepressor <- colSums(counts[repressingTFs,, drop=FALSE])
    logX <- log1p(yRepressor)
    y <- counts[curTF,]+1
    ## power law model: Poisson
    mp <- glm(y ~ logX, family="poisson")
    # plot(x=logX, y=log(y) - logOffset)
    # summary(mp)
    # if X is all zero
    if(any(is.na(coef(mp)))){
      curRho <- 1
      rho_tc[rr,] <- curRho
      next
    }
    # if not negative or not significant
    if(coef(mp)[2] > 0 | summary(mp)$coefficients[2,4] > 0.05){
      curRho <- 1
      rho_tc[rr,] <- curRho
      next
    }
    grid <- seq(min(logX), max(logX), length=100)
    yhat <- predict(mp, newdata=data.frame(logX=grid))
    ## calculate thinning factor rho_tc
    ### X for each cell type
    colsumU <- colSums(U)
    meanLogX <- log1p((yRepressor %*% U) %*% diag(1/colsumU,
                                                  nrow=length(colsumU),
                                                  ncol=length(colsumU)))

    curRho <- predict(mp, newdata=data.frame(logX = c(0, meanLogX)))
    rhoBase <- curRho[1]
    curRho <- exp(curRho[-1]) / exp(rhoBase) #as compared to no expression of repressing TFs
    rho_tc[rr,] <- curRho
  }
  return(rho_tc)
}

pruneRepressingLinks <- function(counts, X){
  ## only keep repressions between TFs
  tfsInTargets <- colnames(X)[colnames(X) %in% rownames(X)]
  tfTargetRowId <- rownames(X) %in% tfsInTargets
  XNoTf <- X[!tfTargetRowId,]
  XNoTf[XNoTf == -1] <- 0
  if(any(rowSums(XNoTf) == 0)){
    message("Pruning ", sum(rowSums(XNoTf) == 0), " genes with only between-gene repressions.")
    pruneId <- rownames(XNoTf)[which(rowSums(XNoTf) == 0)]
    countsRep <- counts # for repressions (thinningFactor function)
    counts <- counts[!rownames(counts) %in% pruneId,]
    XNoTf <- XNoTf[!rownames(XNoTf) %in% pruneId,]
    X <- X[!rownames(X) %in% pruneId,]
  } else {
    countsRep <- counts
  }
  X[rownames(XNoTf),] <- XNoTf
  # now get all the repressions between TFs
  reprId <- arrayhelpers::vec2array(which(X==-1), dim(X))
  dfRepr <- data.frame(repressed=rownames(X)[reprId[,1]],
                       repressor=colnames(X)[reprId[,2]])
  # and finally make a GRN where repressions are set to zero.
  XPos <- X
  XPos[XPos == -1] <- 0
  keep <- rownames(XPos)[rowSums(XPos)>0]
  counts <- counts[keep,]
  XPos <- XPos[keep,]
  XPos <- XPos[,colSums(XPos)>0]
  return(list(XPos = XPos,
              countsRep = countsRep,
              dfRepr = dfRepr))
}

scaleAlpha <- function(Y_gc, XPos, alpha, alphaScale, design){
  message("Prior versus data weight is tuned to be ", alphaScale*100, "%.")
  alpha_gtcList <- list()
  for(cc in 1:ncol(design)){
    alpha_c <- alpha
    for(gg in 1:nrow(alpha)){
      curAlpha_gt <- alpha_c[gg,]
      # if all  positive alpha's are 1, don't change so it reduces back to no prior.
      if(all(curAlpha_gt[!curAlpha_gt==0] == 1)) next
      corFactor_g <- alphaScale * (Y_gc[gg,cc]+1) / sum(curAlpha_gt) #+1 to avoid alpha=0
      scaledAlpha_gt <- curAlpha_gt * corFactor_g
      alpha_c[gg,] <- scaledAlpha_gt
    }
    alpha_gtcList[[cc]] <- alpha_c
  }
  return(alpha_gtcList)
}

EStep <- function(XPos, mu_tc, Y_gc, tfNames){
  sumGene <- XPos %*% mu_tc
  countFracs <- Y_gc / (sumGene+1e-10)
  zList <- list()
  for(tt in 1:ncol(XPos)){
    id <- which(XPos[,tt] == 1)
    curZ_gtc <- sweep(countFracs[id,,drop=FALSE], 2, mu_tc[tt,], "*") # Z = Y * (mu / sum(mu))
    rownames(curZ_gtc) <- paste0(tfNames[tt],";",rownames(curZ_gtc)) #tf;gene
    zList[[tt]] <- curZ_gtc
  }
  Z_gtc <- do.call(rbind, zList)
  return(Z_gtc)
}




getPi_gtc_sufStats <- function(mu_gtc, counts, pt=NULL, qSteps=0.01, U=NULL){

  lvl <- unlist(apply(U,1, function(row){
    which(row == 1)
  }))
  design <- U

  tfRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 1))
  geneRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 2))
  tfUniq <- unique(tfRows)
  geneUniq <- unique(geneRows)
  rn <- rownames(mu_gtc)
  pi_gtc <- list()
  for(gg in 1:length(geneUniq)){
    curGene <- geneUniq[gg]
    id <- which(geneRows == curGene)
    curTFs <- tfRows[id]
    rowSel <- fastmatch::fmatch(paste0(curTFs,";",curGene), rn)
    curMu <- mu_gtc[rowSel,,drop=FALSE]
    curProbs <- sweep(curMu, 2, colSums(curMu)+1e-10, "/")
    pi_gtc[[gg]] <- curProbs
  }
  names(pi_gtc) <- geneUniq
  return(pi_gtc)
}


