
#' @export
tfDistanceRanking <- function(emRes,
                              X,
                              counts,
                              U,
                              cellGroups=NULL,
                              distance="Euclidean",
                              scaleDistance=FALSE,
                              contrast="consecutive",
                              referenceGroup=NULL){
  # emRes is output from EM
  # X is GRN
  # counts is gene expression
  # U is design matrix
  # cellGroups is the pairwise comparison of interest
  # distance is the type of distance used (Euclidean vs L1)
  # contrast specifies the contrast to use for a global test.
  #   if consecutive, it tests each bin vs the next
  #   if reference, it tests referenceBin vs all others
  # for reference contrast, referenceGroup is the group being compared to all others.

  .distanceCalculation <- function(mu_gc,
                                   pi_gtc,
                                   cellID,
                                   distance="Euclidean",
                                   scaleDistance=FALSE){
    ## avoid division by zero
    if(scaleDistance){
      denomOrig <- (mu_gc[,cellID[2]] * pi_gtc[,cellID[2]] + mu_gc[,cellID[1]] * pi_gtc[,cellID[1]])
      denom <- denomOrig
      denom[denomOrig == 0] <- 1
    }

    if(distance == "Euclidean"){
      if(scaleDistance){
        ## avoid division by zero and set corresponding distances to 0
        curDist_g <- ((mu_gc[,cellID[2]] * pi_gtc[,cellID[2]] - mu_gc[,cellID[1]] * pi_gtc[,cellID[1]])^2) / denom
        curDist_g[denomOrig == 0] <- 0
        curDist <- sum(curDist_g)
      } else {
        curDist <- sum((mu_gc[,cellID[2]] * pi_gtc[,cellID[2]] - mu_gc[,cellID[1]] * pi_gtc[,cellID[1]])^2)
      }
    } else if(distance == "L1"){
      if(scaleDistance){
        curDist_g <- (abs(mu_gc[,cellID[2]] * pi_gtc[,cellID[2]] - mu_gc[,cellID[1]] * pi_gtc[,cellID[1]])) / denom
        curDist_g[denomOrig == 0] <- 0
        curDist <- sum(curDist_g)
      } else {
        curDist <- sum(abs(mu_gc[,cellID[2]] * pi_gtc[,cellID[2]] - mu_gc[,cellID[1]] * pi_gtc[,cellID[1]]))
      }
    }
    return(curDist)
  }

  # ignore repressions
  if(any(X == -1)){
    X[X==-1] <- 0
  }

  if(!is.null(cellGroups)){
    if(is.character(cellGroups)){
      cellGroups <- which(colnames(U) %in% cellGroups)
    }
  }

  # tf;gene notation in rownames
  pi_gtc <- getPi_gtc_sufStats(emRes$mu_gtc,
                               counts = as.matrix(counts),
                               U = U)
  mu_gc <- ((counts %*% U) %*% diag(1/colSums(U)))
  pi_gtc_mat <- do.call(rbind, pi_gtc)
  pi_tfNames <- unlist(lapply(strsplit(rownames(pi_gtc_mat), split=";"), "[[", 1))
  pi_geneNames <- unlist(lapply(strsplit(rownames(pi_gtc_mat), split=";"), "[[", 2))

  allTF <- colnames(X)
  tfDist <- vector(length=length(allTF))

  ## loop over all TFs
  for(tt in seq_len(length(allTF))){
    curTF <- allTF[tt]
    curTargets <- rownames(X)[which(X[,curTF] == 1)]
    rowSel <- fastmatch::fmatch(paste0(curTF,";",curTargets), rownames(pi_gtc_mat))
    curPi <- pi_gtc_mat[rowSel,,drop=FALSE]

    ## if two groups are compared
    if(!is.null(cellGroups)){
      curDist <- .distanceCalculation(mu_gc = mu_gc[curTargets,],
                                      pi_gtc = curPi,
                                      cellID = cellGroups,
                                      distance = distance,
                                      scaleDistance = scaleDistance)
    }

    ## global contrasts
    ### consecutive comparisons
    if(contrast == "consecutive"){
      conDist <- c()
      for(kk in seq_len(ncol(U)-1)){
        conDist[kk] <- .distanceCalculation(mu_gc = mu_gc[curTargets,,drop=FALSE],
                                            pi_gtc = curPi,
                                            cellID = c(kk, kk+1),
                                            distance = distance,
                                            scaleDistance = scaleDistance)
      }
      curDist <- sum(conDist)
      ### comparisons wrt reference group
    } else if(contrast == "reference"){
      varsToCompare <- seq_len(ncol(U))[-referenceGroup]
      for(kk in seq_len(length(varsToCompare))){
        conDist[kk] <- .distanceCalculation(mu_gc = mu_gc[curTargets,],
                                            pi_gtc = curPi,
                                            cellID = c(referenceGroup, varsToCompare[kk]),
                                            distance = distance,
                                            scaleDistance = scaleDistance)
      }
      curDist <- sum(conDist)
    }
    tfDist[tt] <- curDist
  }
  names(tfDist) <- allTF
  return(tfDist)
}
