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

  if(repressions){
    pruned <- pruneRepressingLinks(counts, X)
    XPos <- pruned$XPos
    countsRep <- pruned$countsRep
    dfRepr <- pruned$dfRepr
  } else {
    if(any(X == -1)){
      stop("Repressing links are present in GRN, but repressions is set to FALSE.")
    }
    XPos <- X
    countsRep <- counts
  }

  counts <- counts[rownames(XPos),]

  # get design
  design <- U

  ## get sufficient statistics
  counts_suf <- counts %*% design

  tfNames <- colnames(XPos)

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
    rho_tc <- thinningFactor(counts=countsRep,
                             dfRepr=dfRepr,
                             X=X,
                             U=U)
    rho_tc <- rho_tc[rownames(rho_tc) %in% rownames(mu_tc),]
    mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc
    if(is.null(rho_tc)) repressions <- FALSE
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
    Z_gtc <- EStep(XPos = XPos,
                   mu_tc = mu_tc,
                   Y_gc = counts_suf,
                   tfNames = tfNames)

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
  if(maxIter == 0){
    return(list(mu_tc=mu_tc,
                design = design))
  } else {
    return(list(mu_tc=mu_tc,
                mu_gtc=mu_gtc,
                countsSufStats = counts_suf,
                design = design))
  }
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

  if(repressions){
    pruned <- pruneRepressingLinks(counts, X)
    XPos <- pruned$XPos
    countsRep <- pruned$countsRep
    dfRepr <- pruned$dfRepr
  } else {
    if(any(X == -1)){
      stop("Repressing links are present in GRN, but repressions is set to FALSE.")
    }
    XPos <- X
    countsRep <- counts
  }

  counts <- counts[rownames(XPos),]
  alpha <- alpha[rownames(XPos),]

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
  design <- U

  ## get sufficient statistics
  Y_gc <- counts %*% design

  ## set alpha to corresponding scale
  if(alphaScale == "none"){
    alpha <- XPos
  } else {
    message("Prior versus data weight is tuned to be ", alphaScale*100, "%.")
    alpha_gtcList <- scaleAlpha(Y_gc, XPos, alpha, alphaScale, design)
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
    rho_tc <- thinningFactor(counts=countsRep,
                             X=X,
                             dfRepr=dfRepr,
                             U=U)
    rho_tc <- rho_tc[rownames(rho_tc) %in% rownames(mu_tc),]
    mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc
    if(is.null(rho_tc)) repressions <- FALSE
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
    Z_gtc <- EStep(XPos = XPos,
                   mu_tc = mu_tc,
                   Y_gc = Y_gc,
                   tfNames = tfNames)
    #### added on 24/01/23
    # not sure as this is (a) already piNew? (b) does not use prior.
    # pi_gtc <- Z_gtc / Y_gc[unlist(lapply(strsplit(rownames(Z_gtc),split=";"),"[[",2)),]


    if(iter == 1){
      alpha_gtc <- matrix(NA, nrow=nrow(Z_gtc), ncol=ncol(design))
      if(alphaScale == "none"){ # no effect of alpha => Poisson model
        alpha_gtc[] <- 1
      } else {
        for(cc in seq_len(ncol(design))){
          tfIk <- unlist(lapply(strsplit(rownames(Z_gtc), split=";"), "[[", 1))
          geneIk <- unlist(lapply(strsplit(rownames(Z_gtc), split=";"), "[[", 2))
          alpha_gt <- alpha_gtcList[[cc]][cbind(geneIk, tfIk)]
          alpha_gtc[,cc] <- alpha_gt
        }
      }
    }

    ## M-step: estimate mean for each bin
    mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design)) + alpha_gtc - 1
    ### adapted on 24/01/23
    # TODO: pi = (pi*Y + alpha - 1) / sum(pi*Y + alpha - 1)
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


#######################################################
############### DIRICHLET-MULTINOMIAL USING PI ########
#######################################################



dirMultEstimation2 <- function(counts,
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

  if(repressions){
    pruned <- pruneRepressingLinks(counts, X)
    XPos <- pruned$XPos
    countsRep <- pruned$countsRep
    dfRepr <- pruned$dfRepr
  } else {
    if(any(X == -1)){
      stop("Repressing links are present in GRN, but repressions is set to FALSE.")
    }
    XPos <- X
    countsRep <- counts
  }

  counts <- counts[rownames(XPos),]
  alpha <- alpha[rownames(XPos),]

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
  design <- U

  ## get sufficient statistics
  Y_gc <- counts %*% design

  ## set alpha to corresponding scale
  if(alphaScale == "none"){
    alpha <- XPos
  } else {
    # TODO: in gamma modelbased simulation this returnls alpha == 0
    # where they used to be 1.
    message("Prior versus data weight is tuned to be ", alphaScale*100, "%.")
    alpha_gtcList <- scaleAlpha(Y_gc, XPos, alpha, alphaScale, design)
  }

  tfNames <- colnames(XPos)

  ## initialize
  # initialize pi
  pi_gtc <- array(XPos / rowSums(XPos), dim=c(nrow(XPos), ncol(XPos), ncol(design)))
  if(sparse){
    mu_tc <- sparseInitialization_sufStats(counts_suf = Y_gc,
                                           design = design,
                                           X = XPos,
                                           iterOLS = iterOLS,
                                           lassoFamily = lassoFamily)
  } else {
    ### note that this is the mean across all genes a TF is regulating.
    ### we could consider other properties
    mu_tc <- matrix(NA, nrow=ncol(XPos), ncol=ncol(design))
    for(tt in 1:ncol(XPos)){
      weights <- pi_gtc[,tt,1] # same for all design columns
      hlp <- weights * Y_gc
      mu_tc[tt,] <- colMeans(hlp[weights > 0,,drop=FALSE])
    }
    rownames(mu_tc) <- colnames(XPos)
  }


  ## incorporate repressions
  if(repressions){
    rho_tc <- thinningFactor(counts=countsRep,
                             X=X,
                             dfRepr=dfRepr,
                             U=U)
    rho_tc <- rho_tc[rownames(rho_tc) %in% rownames(mu_tc),]
    mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc
    if(is.null(rho_tc)) repressions <- FALSE
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
    # Z_gtc <- EStep(XPos = XPos,
    #                mu_tc = mu_tc,
    #                Y_gc = Y_gc,
    #                tfNames = tfNames)
    #### added on 31/01/23
    ### for each design group: (make an array of nrow(X), ncol(X), ncol(design))
    # Z_gtc <- diag(Y_gc[,1]) %*% pi_gtc
    # Z_gtc <- array(Z_gtc, dim=c(nrow(Z_gtc), ncol(Z_gtc), ncol(design)))
    # TODO: check if works appropriately with multiple design columns
    Z_gtc <- EStep2(XPos = XPos,
                   mu_tc = mu_tc,
                   Y_gc = Y_gc,
                   design = design)


    if(iter == 1){
      alpha_gtc <- array(NA, dim=c(nrow(Z_gtc), ncol(Z_gtc), ncol(design)))
      if(alphaScale == "none"){ # no effect of alpha => Poisson model
        alpha_gtc[] <- 1
      } else {
        for(cc in seq_len(ncol(design))){
          alpha_gtc[,,cc] <- alpha_gtcList[[cc]]
          # alpha should minimum be 1.
          alpha_gtc[,,cc][XPos==1][which(alpha_gtc[,,cc][XPos==1]<1)] <- 1
        }
      }
    }

    ## M-step: estimate mean for each bin
    # TODO: this should be outside of the loop
    mu_gc <- Y_gc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
    # Original: mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design)) + alpha_gtc - 1
    mu_gtc <- array(0, dim=c(nrow(XPos), ncol(XPos), ncol(design)))
    ### adapted on 31/01/23
    ## TODO make efficient
    ### for each design group:
    for(cc in 1:ncol(design)){
      # denom <- (diag(Z_gtc[,,cc] %*% t(alpha_gtc[,,cc])) - rowSums(XPos))
      for(gg in 1:nrow(XPos)){
        posid <- which(XPos[gg,]>0)
        # pi_gtc[gg,posid] <- (Z_gtc[gg,posid,cc] + alpha_gtc[gg,posid,cc] - 1) / denom[gg]
        curPi <- ((Z_gtc[gg,posid,cc] + alpha_gtc[gg,posid,cc] - 1)+1e-10) / (sum(Z_gtc[gg,posid,cc] + alpha_gtc[gg,posid,cc] - 1)+1e-10)
        if(all(curPi==1)) curPi <- rep(1/length(posid), length(posid)) # when all Z=0
        stopifnot(all.equal.numeric(sum(curPi),1))
        pi_gtc[gg,posid,cc] <- curPi
      }
      mu_gtc[,,cc] <- diag(mu_gc[,cc]) %*% pi_gtc[,,cc]
    }
    # stopifnot(all(rowSums(pi_gtc)==1))
    if(any(mu_gtc<0)) mu_gtc[mu_gtc<0] <- 0

    ## update mu_tc
    # tfRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 1))
    # mu_tc <- t(sapply(unique(tfNames), function(curTF){
    #   tfid <- which(tfRows == curTF)
    #   colMeans(mu_gtc[tfid,,drop=FALSE])
    # }))
    # if(ncol(design) == 1) mu_tc <- t(mu_tc)
    ### adapted on 31/01/23
    #TODO: make efficient
    for(cc in 1:ncol(design)){
      for(tt in 1:ncol(XPos)){
        mu_tc[tt,cc] <- mean(mu_gtc[which(XPos[,tt]==1),tt,cc])
      }
    }


    ## log-likelihood (this actually takes quite some time...)
    llAllCur <- sum(log(dpois(x=round(Z_gtc), lambda=mu_gtc)+1e-16))

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
        # TODO: I guess we should now be returning estimates including alpha.
        # ## get estimates without alpha as final update.
        # # # final mu_tc
        # # mu_tc <- t(sapply(unique(tfNames), function(curTF){
        # #   tfid <- which(tfRows == curTF)
        # #   colMeans(mu_gtc[tfid,,drop=FALSE])
        # # }))
        # # if(ncol(design) == 1) mu_tc <- t(mu_tc)
        # # final mu_gtc
        # mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
        # # final mu_tc
        # mu_tc <- t(sapply(unique(tfNames), function(curTF){
        #   tfid <- which(tfRows == curTF)
        #   colMeans(mu_gtc[tfid,,drop=FALSE])
        # }))
        # if(ncol(design) == 1) mu_tc <- t(mu_tc)
        return(list(mu_tc=mu_tc,
                    mu_gtc=mu_gtc,
                    countsSufStats = Y_gc,
                    design = design))
      }
    }
  }

  # # # final mu_tc
  # # mu_tc <- t(sapply(unique(tfNames), function(curTF){
  # #   tfid <- which(tfRows == curTF)
  # #   colMeans(mu_gtc[tfid,,drop=FALSE])
  # # }))
  # # if(ncol(design) == 1) mu_tc <- t(mu_tc)
  # # final mu_gtc
  # mu_gtc <- Z_gtc %*% array(diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design)),
  #                           dim=c(ncol(design),ncol(design),ncol(design)))
  # # final mu_tc
  # mu_tc <- t(sapply(unique(tfNames), function(curTF){
  #   tfid <- which(tfRows == curTF)
  #   colMeans(mu_gtc[tfid,,drop=FALSE])
  # }))
  # if(ncol(design) == 1) mu_tc <- t(mu_tc)
  return(list(mu_tc=mu_tc,
              mu_gtc=mu_gtc,
              countsSufStats = Y_gc,
              design = design))
}


#######################################################
########### emp Bayes DIRICHLET-MULTINOMIAL ###########
#######################################################
dirMultEstimationAlpha <- function(counts,
                                   X,
                                   alpha=NULL,
                                   rho_t=NULL,
                                   pt=NULL,
                                   U=NULL,
                                   nIters=20,
                                   qSteps=0.01,
                                   plot=FALSE,
                                   verbose=FALSE,
                                   epsilon=1e-2,
                                   iterOLS=0,
                                   repressions = TRUE,
                                   sparse = TRUE,
                                   alphaScale = 1,
                                   lassoFamily = "gaussian"){
  # counts is nGenes x nCells count matrix of gene expression
  # X is TF regulation matrix with dims nGenes x nTranscriptionFactors
  # U is nCells x nVariables design matrix
  # nIters is number of iterations
  # qSteps are the quantile steps to make pseudotime bins
  # alphaScale is how you want to scale the prior vs the data. alphaScale=1 means we believe
  #       the prior as much as we believe the dat. alphaScale=1/2 means we believe the prior
  #       half of what we believe the data.


  ## checks on X
  if(is.null(rownames(X))){
    stop("Please provide rownames for X.")
  }

  if(!all(rownames(X) %in% rownames(counts))){
    stop("Not all gene names in X are present in counts.")
  }
  if(!is.null(alpha)){
    if(!all(alpha[X==0] == 0)){
      stop("alpha_gt cannot be positive if there is no edge in the GRN.")
    }
  }

  counts <- counts[rownames(X),]

  if(repressions){
    pruned <- pruneRepressingLinks(counts, X)
    XPos <- pruned$XPos
    countsRep <- pruned$countsRep
    dfRepr <- pruned$dfRepr
    if(is.null(alpha)) alpha <- XPos
  } else {
    if(any(X == -1)){
      stop("Repressing links are present in GRN, but repressions is set to FALSE.")
    }
    XPos <- X
    countsRep <- counts
    if(is.null(alpha)) alpha <- XPos
  }

  counts <- counts[rownames(XPos),]
  alpha <- alpha[rownames(XPos),]

  if(is.null(colnames(XPos))){
    colnames(XPos) <- paste0("tf", 1:ncol(XPos))
  }

  design <- U

  ## get sufficient statistics
  Y_gc <- counts %*% design

  ## set alpha to corresponding scale
  if(alphaScale == "none"){
    alpha <- XPos
  } else {
    message("Prior versus data weight is tuned to be ", alphaScale*100, "%.")
    alpha_gtcList <- scaleAlpha(Y_gc, XPos, alpha, alphaScale, design)
    alpharange <- range(alpha_gtcList)
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
    rho_tc <- thinningFactor(counts=countsRep,
                             X=X,
                             U=U,
                             dfRepr=dfRepr)
    rho_tc <- rho_tc[rownames(rho_tc) %in% rownames(mu_tc),]
    mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc
    if(is.null(rho_tc)) repressions <- FALSE
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
    Z_gtc <- EStep(XPos = XPos,
                   mu_tc = mu_tc,
                   Y_gc = Y_gc,
                   tfNames = tfNames)

    # if(iter == 1){
    # get right structure
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
    # } else if(iter > 1){
    #   # update alpha
    #   alpha_gtc <- matrix(NA, nrow=nrow(Z_gtc), ncol=ncol(design))
    #   if(alphaScale == "none"){ # no effect of alpha => Poisson model
    #     alpha_gtc[] <- 1
    #   } else {
    #     #### TODO: does not take into account sparsity of mu_tc in some c...
    #     #### 2nd TODO: should alpha estimation be focussed on non-zero pi?
    #     for(cc in 1:ncol(design)){
    #       alpha_gt <- alpha[cbind(geneIk, tfIk)]
    #       alpha_gtc[,cc] <- alpha_gt
    #     }
    #   }
    # }

    ## M-step I: estimate mean for each bin
    # mu_gtc <- Z_gtc %*% diag(1/colSums(design)) ## Poisson
    mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design)) + alpha_gtc - 1
    if(any(mu_gtc<0)) mu_gtc[mu_gtc<0] <- 0

    ## update mu_tc
    tfRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 1))
    mu_tc <- t(sapply(unique(tfNames), function(curTF){
      tfid <- which(tfRows == curTF)
      colMeans(mu_gtc[tfid,,drop=FALSE])
    }))
    if(ncol(design) == 1) mu_tc <- t(mu_tc)


    ### M-step II: estimate alpha_gt
    ## TODO: faster way?
    if(!alphaScale == "none"){
      for(gg in 1:nrow(XPos)){
        idTF <- which(XPos[gg,] == 1)
        if(length(idTF) == 1) next
        # if all mu's are zero, set alpha to zero.
        if(all(mu_tc[idTF,] == 0)){
          # pi_gtc <- mu_tc[idTF,,drop=FALSE]
          # pi_gtc[] <- 1/nrow(mu_tc[idTF,,drop=FALSE])
          alpha[cbind(gg, idTF)] <- 0
          next
        } else {
          pi_gtc <- sweep(mu_tc[idTF,,drop=FALSE]+1e-10, 2, colSums(mu_tc[idTF,,drop=FALSE])+1e-10, "/")
          ## make sure sums to 1
          pi_gtc <- sweep(pi_gtc, 2, colSums(pi_gtc),"/")
        }
        alpha_relative <- rowMeans(pi_gtc)
        # alpha_relative <- sapply(seq_len(nrow(pi_gtc)), function(rr){
        #   if(any(mu_tc[idTF[rr],] > 0)){
        #     return(mean(pi_gtc[rr, mu_tc[idTF[rr],] > 0]))
        #   } else {
        #     return(1e-10)
        #   }
        # })
        alpha_relative <- alpha_relative / sum(alpha_relative)
        # curAlpha <- sirt::dirichlet.mle(t(pi_gtc))$alpha
        ## if a gene is regulated by only one TF, no need to estimate magnitude
        if(!any(alpha_relative > .999)){
          magnitudeOld <- sum(alpha[gg,])
          if(magnitudeOld < 1e-10) magnitudeOld <- 1
          magnitudeNew <- updateS(magnitudeOld, alpha_relative, t(pi_gtc))
        } else {
          magnitudeNew <- nrow(pi_gtc)
        }
        alpha[cbind(gg, idTF)] <- alpha_relative*magnitudeNew
      }
      # reset where necessary
      alpha[alpha<0] <- 0
      # alpha[alpha > max(alpharange)] <- max(alpharange)
      alpha[is.infinite(alpha)] <- max(alpharange)
      # alpha[alpha < 1 & alpha > 0] <- 1
      # alpha[alpha > 1e4] <- 100
    }


    ## scale alpha
    ## set alpha to corresponding scale
    if(alphaScale == "none"){
      alpha <- XPos
    } else {
      message("Prior versus data weight is tuned to be ", alphaScale*100, "%.")
      alpha_gtcList <- scaleAlpha(Y_gc, XPos, alpha, alphaScale, design)
    }

    ## incorporate thinning factor
    if(repressions) mu_tc[rownames(rho_tc),] <- mu_tc[rownames(rho_tc),]*rho_tc

    ## log-likelihood (this actually takes quite some time...)
    llAllCur <- sum(log(dpois(x=round(Z_gtc), lambda=mu_gtc)+1e-16))

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
        mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
        mu_tc <- t(sapply(unique(tfNames), function(curTF){
          tfid <- which(tfRows == curTF)
          colMeans(mu_gtc[tfid,,drop=FALSE])
        }))
        return(list(mu_tc=mu_tc,
                    mu_gtc=mu_gtc,
                    countsSufStats = Y_gc,
                    design = design))
      }
    }
  }

  # final mu_gtc
  mu_gtc <- Z_gtc %*% diag(1/colSums(design), nrow=ncol(design), ncol=ncol(design))
  # final mu_tc
  mu_tc <- t(sapply(unique(tfNames), function(curTF){
    tfid <- which(tfRows == curTF)
    colMeans(mu_gtc[tfid,,drop=FALSE])
  }))
  return(list(mu_tc=mu_tc,
              mu_gtc=mu_gtc,
              countsSufStats = Y_gc,
              design = design))
}
