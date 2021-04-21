#' @include utils.R

poissonEstimation <- function(counts,
                              X,
                              rho_t=NULL,
                              U=NULL,
                              maxIter=20,
                              plot=FALSE,
                              verbose=FALSE,
                              epsilon=1e-2,
                              iterOLS=0,
                              lassoFamily = "gaussian",
                              repressions=TRUE,
                              sparse=TRUE){

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


  counts <- counts[rownames(XPos),]

  # get design
  design <- U

  ## get sufficient statistics
  counts_suf <- counts %*% design

  tfNames <- colnames(XPos)
  ## initialize
  ## initialize
  if(sparse){
    mu_tc <- sparseInitialization_sufStats(counts_suf = counts_suf,
                                           design = design,
                                           X = XPos,
                                           iterOLS = iterOLS,
                                           lassoFamily = lassoFamily)
  } else {
    # initialize E(Z)
    EZ_probOrig <- XPos / rowSums(XPos)
    ### note that this is the mean across all genes a TF is regulating.
    ### we could consider other properties
    mu_tc <- matrix(NA, nrow=ncol(XPos), ncol=ncol(design))
    for(tt in 1:ncol(XPos)){
      weights <- EZ_probOrig[,tt]
      hlp <- weights * counts_suf
      mu_tc[tt,] <- colMeans(hlp[weights > 0,,drop=FALSE])
    }
    rownames(mu_tc) <- colnames(XPos)
    colnames(mu_tc) <- colnames(design)
  }

  ## incorporate repressions
  if(repressions){
    # if(!exists("countsRep")){
    #   rho_tc <- thinningFactor(counts=counts,
    #                          X=X,
    #                          U=U,
    #                          pt=pt,
    #                          qSteps=qSteps)
    # } else {
    rho_tc <- thinningFactor(counts=countsRep,
                             X=X,
                             U=U)
    #}
    rho_tc <- rho_tc[rownames(rho_tc) %in% rownames(mu_tc),]
    mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc
  }

  iter <- 0
  while(iter < maxIter){
    iter <- iter + 1
    if(verbose){
      if(iter == 1){
        message(paste0("iteration ", iter, "\n"))
      } else {
        message(paste0("iteration ", iter, ". Log-lik: ", round(tail(llAll,1), 3), "\n"))
      }
    }

    ## E-step: for a gene, select its regulating TFs, and normalize them
    sumGene <- XPos %*% mu_tc
    countFracs <- counts_suf / (sumGene+1e-10)
    zList <- list()
    for(tt in 1:ncol(XPos)){
      id <- which(XPos[,tt] == 1)
      curZ_gtc <- sweep(countFracs[id,,drop=FALSE], 2, mu_tc[tt,], "*") # Z = Y * (mu / sum(mu))
      rownames(curZ_gtc) <- paste0(tfNames[tt],";",rownames(curZ_gtc)) #tf;gene
      zList[[tt]] <- curZ_gtc
    }
    Z_gtc <- do.call(rbind, zList)

    ## M-step: estimate mean for each bin
    mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))

    ## log-likelihood (this actually takes quite some time...)
    llAllCur <- sum(log(dpois(x=round(Z_gtc), lambda=mu_gtc)+1e-16))

    ## update mu_tc
    tfRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 1))
    mu_tc <- t(sapply(unique(tfNames), function(curTF){
      tfid <- which(tfRows == curTF)
      colMeans(mu_gtc[tfid,,drop=FALSE])
    }))
    if(ncol(design) == 1) mu_tc <- t(mu_tc)

    ## incorporate thinning factor
    if(repressions){
      if(iter > 1){
        mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc
      }
    }

    if(iter == 1){
      llAll <- llAllCur
    } else {
      llAll <- c(llAll, llAllCur)
    }

    if(plot){
      par(mfrow=c(1,2))
      plot(x=1:iter, y=llAll, type="b", xlab="Iteration", ylab="Log likelihood")
      #hist(probVec, breaks=40)
    }
    ### check convergence
    if(iter > 1){
      if(abs(llAll[length(llAll)] - llAll[length(llAll)-1]) < epsilon){
        message("Converged.")

        # final Z_gtc
        # sumGene <- XPos %*% mu_tc
        # countFracs <- counts_suf / (sumGene+1e-10)
        # zList <- list()
        # for(tt in 1:ncol(XPos)){
        #   id <- which(XPos[,tt] == 1)
        #   curZ_gtc <- sweep(countFracs[id,,drop=FALSE], 2, mu_tc[tt,], "*") # Z = Y * (mu / sum(mu))
        #   rownames(curZ_gtc) <- paste0(tfNames[tt],";",rownames(curZ_gtc)) #tf;gene
        #   zList[[tt]] <- curZ_gtc
        # }
        # Z_gtc <- do.call(rbind, zList)
        # # final mu_gtc
        # mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
        # # final mu_tc
        # mu_tc <- t(sapply(unique(tfNames), function(curTF){
        #   tfid <- which(tfRows == curTF)
        #   colMeans(mu_gtc[tfid,,drop=FALSE])
        # }))
        # if(ncol(design) == 1) mu_tc <- t(mu_tc)


        # mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
        # mu_tc <- t(sapply(unique(tfNames), function(curTF){
        #   tfid <- which(tfRows == curTF)
        #   colMeans(mu_gtc[tfid,,drop=FALSE])
        # }))
        return(list(mu_tc=mu_tc,
                    mu_gtc=mu_gtc,
                    countsSufStats = counts_suf,
                    design = design))
      }
    }
  }

  # # final Z_gtc
  # sumGene <- XPos %*% mu_tc
  # countFracs <- counts_suf / (sumGene+1e-10)
  # zList <- list()
  # for(tt in 1:ncol(XPos)){
  #   id <- which(XPos[,tt] == 1)
  #   curZ_gtc <- sweep(countFracs[id,,drop=FALSE], 2, mu_tc[tt,], "*") # Z = Y * (mu / sum(mu))
  #   rownames(curZ_gtc) <- paste0(tfNames[tt],";",rownames(curZ_gtc)) #tf;gene
  #   zList[[tt]] <- curZ_gtc
  # }
  # Z_gtc <- do.call(rbind, zList)
  # # final mu_gtc
  # mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
  # # final mu_tc
  # mu_tc <- t(sapply(unique(tfNames), function(curTF){
  #   tfid <- which(tfRows == curTF)
  #   colMeans(mu_gtc[tfid,,drop=FALSE])
  # }))
  # if(ncol(design) == 1) mu_tc <- t(mu_tc)
  # mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
  # mu_tc <- t(sapply(unique(tfNames), function(curTF){
  #   tfid <- which(tfRows == curTF)
  #   colMeans(mu_gtc[tfid,,drop=FALSE])
  # }))
  return(list(mu_tc=mu_tc,
              mu_gtc=mu_gtc,
              countsSufStats = counts_suf,
              design = design))
}


#######################################################
############### DIRICHLET-MULTINOMIAL #################
#######################################################


dirMultEstimation <- function(counts,
                          X,
                          alpha = NULL,
                          rho_t=NULL,
                          U=NULL,
                          nIters=20,
                          plot=FALSE,
                          verbose=FALSE,
                          epsilon=1e-2,
                          iterOLS=0,
                          alphaScale=1,
                          repressions = TRUE,
                          sparse = TRUE,
                          lassoFamily = "gaussian"){
  # counts is nGenes x nCells count matrix of gene expression
  # X is TF regulation matrix with dims nGenes x nTranscriptionFactors
  # U is nCells x nVariables design matrix
  # nIters is number of iterations
  # qSteps are the quantile steps to make pseudotime bins
  # alphaScale is how you want to scale the prior vs the data. alphaScale=1 means we believe
  #       the prior as much as we believe the dat. alphaScale=1/2 means we believe the prior
  #       half of what we believe the data.
  # alphaScale = "none" means no effect of alpha and reduced to Poisson model.
  # Note that providing a very low alphaScale will result in many mu_tc=0, since due to the -1 in the
  #   calculation for mu_gtc, many will be negative and will be reset to zero.

  ## checks on X
  if(is.null(rownames(X))){
    stop("Please provide rownames for X.")
  }

  if(!all(rownames(X) %in% rownames(counts))){
    stop("Not all gene names in X are present in counts.")
  }

  # counts <- counts[rownames(X),]

  ## only keep repressions between TFs
  tfsInTargets <- colnames(X)[colnames(X) %in% rownames(X)]
  tfTargetRowId <- rownames(X) %in% tfsInTargets
  XNoTf <- X[!tfTargetRowId,]
  XNoTf[XNoTf == -1] <- 0
  if(any(rowSums(XNoTf) == 0)){
    message("Pruning ", sum(rowSums(XNoTf) == 0), " genes with only between-gene repressions.")
    pruneId <- rownames(XNoTf)[which(rowSums(XNoTf) == 0)]
    counts <- counts[!rownames(counts) %in% pruneId,]
    countsRep <- counts # for repressions (thinningFactor function)
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
  keep <- rowSums(XPos)>0
  counts <- counts[keep,]
  XPos <- XPos[keep,]
  XPos <- XPos[,colSums(XPos)>0]

  counts <- counts[rownames(XPos),]

  ## checks on alpha
  alpha <- alpha[rownames(XPos), colnames(XPos)]
  if(is.null(alpha)){
    message("alpha not provided, using X instead.")
    alpha <- XPos
  }
  if(!all(alpha[XPos==0] == 0)){
    id <- which(alpha[XPos==0] != 0)
    message("alpha_gt cannot be positive if there is no edge in the GRN,",
            "setting alpha for ", length(id),  " links to zero.")
    alpha[id] <- 0
  }

  if(any(alpha>0 & alpha<1)){
    smallAlpha <- which(alpha>0 & alpha<1)
    warning(paste0("alpha_gt must be 0 or greater than or equal to 1. Converting ",
                   length(smallAlpha)," values to 1."))
    alpha[smallAlpha] <- 1
  }


  if(is.null(colnames(XPos))){
    colnames(XPos) <- paste0("tf", 1:ncol(XPos))
  }

  # get design groups
  lvl <- unlist(apply(U,1, function(row){
    which(row == 1)
  }))
  design <- U

  ## get sufficient statistics
  Y_gc <- counts %*% design

  ## set alpha to corresponding scale
  if(alphaScale == "none"){
    alpha <- XPos
  } else {
    message("Prior versus data weight is tuned to be ", alphaScale*100, "%.")
    alpha_gtcList <- list()
    for(cc in 1:ncol(design)){
      alpha_c <- alpha
      for(gg in 1:nrow(counts)){
        curAlpha_gt <- alpha_c[gg,]
        # if all  positive alpha's are 1, don't change so it reduces back to no prior.
        if(all(curAlpha_gt[!curAlpha_gt==0] == 1)) next
        corFactor_g <- alphaScale * (Y_gc[gg,cc]+1) / sum(curAlpha_gt) #+1 to avoid alpha=0
        scaledAlpha_gt <- curAlpha_gt * corFactor_g
        alpha_c[gg,] <- scaledAlpha_gt
      }
      alpha_gtcList[[cc]] <- alpha_c
    }
  }

  tfNames <- colnames(XPos)

  ## initialize
  if(sparse){
    mu_tc <- sparseInitialization_sufStats(counts_suf = Y_gc,
                                           design = design,
                                           X = XPos,
                                           iterOLS = iterOLS,
                                           lassoFamily = lassoFamily)
  } else {
    # initialize E(Z)
    EZ_probOrig <- XPos / rowSums(XPos)
    ### note that this is the mean across all genes a TF is regulating.
    ### we could consider other properties
    mu_tc <- matrix(NA, nrow=ncol(XPos), ncol=ncol(design))
    for(tt in 1:ncol(XPos)){
      weights <- EZ_probOrig[,tt]
      hlp <- weights * Y_gc
      mu_tc[tt,] <- colMeans(hlp[weights > 0,,drop=FALSE])
    }
    rownames(mu_tc) <- colnames(XPos)
  }


  ## incorporate repressions
  if(repressions){
    # if(!exists("countsRep")){
    #   rho_tc <- thinningFactor(counts=counts,
    #                          X=X,
    #                          U=U,
    #                          pt=pt,
    #                          qSteps=qSteps)
    # } else {
    rho_tc <- thinningFactor(counts=countsRep,
                             X=X,
                             U=U)
    #}

    rho_tc <- rho_tc[rownames(rho_tc) %in% rownames(mu_tc),]
    mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc
  }

  iter <- 0
  while(iter < nIters){
    iter <- iter + 1
    if(verbose){
      if(iter == 1){
        message(paste0("iteration ", iter, "\n"))
      } else {
        message(paste0("iteration ", iter, ". Log-lik: ", round(tail(llAll,1), 3), "\n"))
      }
    }

    ## E-step: for a gene, select its regulating TFs, and normalize them
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

    if(iter == 1){
      alpha_gtc <- matrix(NA, nrow=nrow(Z_gtc), ncol=ncol(design))
      if(alphaScale == "none"){ # no effect of alpha => Poisson model
        alpha_gtc[] <- 1
      } else {
        for(cc in 1:ncol(design)){
          alpha_gt <- Z_gtc
          tfIk <- unlist(lapply(strsplit(rownames(Z_gtc), split=";"), "[[", 1))
          geneIk <- unlist(lapply(strsplit(rownames(Z_gtc), split=";"), "[[", 2))
          alpha_gt <- alpha_gtcList[[cc]][cbind(geneIk, tfIk)]
          alpha_gtc[,cc] <- alpha_gt
        }
      }
    }

    ## M-step: estimate mean for each bin
    # mu_gtc <- Z_gtc %*% diag(1/colSums(design)) ## Poisson
    mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design)) + alpha_gtc - 1
    if(any(mu_gtc<0)) mu_gtc[mu_gtc<0] <- 0

    ## log-likelihood (this actually takes quite some time...)
    llAllCur <- sum(log(dpois(x=round(Z_gtc), lambda=mu_gtc)+1e-16))

    ## update mu_tc
    tfRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 1))
    mu_tc <- t(sapply(unique(tfNames), function(curTF){
      tfid <- which(tfRows == curTF)
      colMeans(mu_gtc[tfid,,drop=FALSE])
    }))
    if(ncol(design) == 1) mu_tc <- t(mu_tc)

    ## incorporate thinning factor
    if(repressions) mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc

    if(iter == 1){
      llAll <- llAllCur
    } else {
      llAll <- c(llAll, llAllCur)
    }

    if(plot){
      par(mfrow=c(1,2))
      plot(x=1:iter, y=llAll, type="b", xlab="Iteration", ylab="Log likelihood")
      #hist(probVec, breaks=40)
    }
    ### check convergence
    if(iter > 1){
      if(abs(llAll[length(llAll)] - llAll[length(llAll)-1]) < epsilon){
        message("Converged.")
        ## get estimates without alpha as final update.
        # # final mu_tc
        # mu_tc <- t(sapply(unique(tfNames), function(curTF){
        #   tfid <- which(tfRows == curTF)
        #   colMeans(mu_gtc[tfid,,drop=FALSE])
        # }))
        # if(ncol(design) == 1) mu_tc <- t(mu_tc)
        # final Z_gtc
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
        # final mu_gtc
        mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
        # final mu_tc
        mu_tc <- t(sapply(unique(tfNames), function(curTF){
          tfid <- which(tfRows == curTF)
          colMeans(mu_gtc[tfid,,drop=FALSE])
        }))
        if(ncol(design) == 1) mu_tc <- t(mu_tc)
        return(list(mu_tc=mu_tc,
                    mu_gtc=mu_gtc,
                    countsSufStats = Y_gc,
                    design = design))
      }
    }
  }

  # # final mu_tc
  # mu_tc <- t(sapply(unique(tfNames), function(curTF){
  #   tfid <- which(tfRows == curTF)
  #   colMeans(mu_gtc[tfid,,drop=FALSE])
  # }))
  # if(ncol(design) == 1) mu_tc <- t(mu_tc)
  # final Z_gtc
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
  # final mu_gtc
  mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
  # final mu_tc
  mu_tc <- t(sapply(unique(tfNames), function(curTF){
    tfid <- which(tfRows == curTF)
    colMeans(mu_gtc[tfid,,drop=FALSE])
  }))
  if(ncol(design) == 1) mu_tc <- t(mu_tc)
  return(list(mu_tc=mu_tc,
              mu_gtc=mu_gtc,
              countsSufStats = Y_gc,
              design = design))
}

