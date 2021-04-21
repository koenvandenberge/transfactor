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
