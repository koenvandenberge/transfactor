#' @include utils.R


.distanceCalculation_original <- function(mu_gc,
                                 pi_gtc,
                                 cellID,
                                 distance="Euclidean",
                                 scaleDistance=FALSE){
  ## avoid division by zero
  if(scaleDistance){
    if(distance %in% c("Euclidean", "L1")){
      denomOrig <- (mu_gc[,cellID[2]] * pi_gtc[,cellID[2]] + mu_gc[,cellID[1]] * pi_gtc[,cellID[1]])
      denom <- denomOrig
      denom[denomOrig == 0] <- 1
    }
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

.distanceCalculation_newMethods <- function(mu_tc,
                                            cellID,
                                            distance="rank",
                                            scaleDistance=FALSE){

  if(distance == "rank") scaleDistance <- FALSE
  if(scaleDistance){
    denomOrig <- (mu_tc[,cellID[2]] + mu_tc[,cellID[1]])
    denom <- denomOrig
    denom[denomOrig == 0] <- 1
  }

  if(distance == "rank"){
    mu_tc <- mu_tc[,cellID]
    tfRanks <- apply(mu_tc, 2, rank)
    curDist <- abs(tfRanks[,2] - tfRanks[,1])
  } else if(distance == "EuclideanTF"){
    if(scaleDistance){
      ## avoid division by zero and set corresponding distances to 0
      curDist <- ((mu_tc[,cellID[2]] - mu_tc[,cellID[1]])^2) / denom
      curDist[denomOrig == 0] <- 0
    } else {
      curDist <- (mu_tc[,cellID[2]] - mu_tc[,cellID[1]])^2
    }
  } else if(distance == "L1TF"){
    if(scaleDistance){
      curDist <- (abs(mu_tc[,cellID[2]] - mu_tc[,cellID[1]])) / denom
      curDist[denomOrig == 0] <- 0
    } else {
      curDist <- abs(mu_tc[,cellID[2]] - mu_tc[,cellID[1]])
    }
  }
  return(curDist)
}

#' @title Distance-based ranking of transcription factors.
#' @description Rank transcription factors based on a number of distance metrics
#' and contrasts.
#'
#' @param activity Estimated TF activity. The output from \code{\link{estimateActivity}}.
#' @param X Gene regulatory network of dimensions G x T.
#' @param counts Gene expression counts, of dimenstions G x n.
#' @param U Design matrix, of dimensions n x p. The design matrix should not
#' contain an intercept. This is the same design matrix provided to \code{\link{estimateActivity}}.
#' @param cellGroups The columns of \code{U} to be compared, if a pairwise
#' comparison of sets of cells is of interest.
#' @param distance The distance function to use. Defaults to \code{"Euclidean"}.
#' Options are
#'  - \code{"Euclidean"}: Euclidean distance on mu_gtc.
#'  - \code{"L1"}: L1 distance on mu_gtc.
#'  - \code{"EuclideanTF"}: Euclidean distance on mu_tc.
#'  - \code{"L1TF"}: L1 distance on mu_tc.
#'  - \code{"rank"}: Difference of ranks on mu_tc.
#' @param scaleDistance Logical. Should the distances be scaled?
#' @param contrast The contrast of interest. If \code{"consecutive"},
#' it tests each group versus the next (relevant for datasets with trajectories),
#' if \code{"reference"}, it tests one set of cells (defined via \code{referenceGroup})
#' versus all others.
#' @param referenceGroup If \code{contrast} is \code{"reference"}, \code{referenceGroup}
#' is the set of cells being compared to all others.
#' @return A vector with distances for each TF.
#' @examples
#' counts <- matrix(rpois(n= 1e4, lambda=4), nrow=100, ncol=100)
#' X <- matrix(0, nrow=100, ncol=12)
#' X[cbind(sample(100, size=250, replace=TRUE), sample(12, size=250, replace=TRUE))] <- 1
#' rownames(X) <- rownames(counts) <- paste0("gene",1:100)
#' act <- estimateActivity(counts, X)
#' ict <- rep(1, ncol(counts))
#' cellType <- factor(rep(c("a", "b"), each=50))
#' U <- model.matrix(~ -1 + cellType)
#' act <- estimateActivity(counts, X, U=U)
#' tfDist <- tfDistance(act, X, counts, U)
#' @export
#' @rdname tfDistance
setMethod(f = "tfDistance",
          signature = c(activity = "list",
                        counts = "matrix",
                        X = "matrix",
                        U = "matrix"),
          definition = function(activity,
                       X,
                       counts,
                       U,
                       cellGroups=NULL,
                       distance="Euclidean",
                       scaleDistance=FALSE,
                       contrast="consecutive",
                       referenceGroup=NULL){

  if(is.null(colnames(X))){
    colnames(X) <- paste0("tf",1:ncol(X))
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

  if(distance %in% c("Euclidean", "L1")){
    # tf;gene notation in rownames
    pi_gtc <- getPi_gtc_sufStats(activity$mu_gtc,
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
        curDist <- .distanceCalculation_original(mu_gc = mu_gc[curTargets,],
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
          conDist[kk] <- .distanceCalculation_original(mu_gc = mu_gc[curTargets,,drop=FALSE],
                                                       pi_gtc = curPi,
                                                       cellID = c(kk, kk+1),
                                                       distance = distance,
                                                       scaleDistance = scaleDistance)
        }
        curDist <- sum(conDist)
        ### comparisons wrt reference group
      } else if(contrast == "reference"){
        refDist <- c()
        varsToCompare <- seq_len(ncol(U))[-referenceGroup]
        for(kk in seq_len(length(varsToCompare))){
          refDist[kk] <- .distanceCalculation_original(mu_gc = mu_gc[curTargets,],
                                                       pi_gtc = curPi,
                                                       cellID = c(referenceGroup, varsToCompare[kk]),
                                                       distance = distance,
                                                       scaleDistance = scaleDistance)
        }
        curDist <- sum(refDist)
      }
      tfDist[tt] <- curDist
    }
    names(tfDist) <- allTF

  } else if(distance %in% c("rank", "EuclideanTF", "L1TF")){
    ## if two groups are compared
    if(!is.null(cellGroups)){
      curDist <- .distanceCalculation_newMethods(mu_tc = activity$mu_tc,
                                               cellID = cellGroups,
                                               distance = distance,
                                               scaleDistance = scaleDistance)
      return(curDist)
    }

    ## global contrasts
    ### consecutive comparisons
    if(contrast == "consecutive"){
      conDist <- matrix(NA, nrow=ncol(X), ncol=ncol(U)-1)
      for(kk in seq_len(ncol(U)-1)){
        conDist[,kk] <- .distanceCalculation_newMethods(mu_tc = activity$mu_tc,
                                                     cellID = c(kk, kk+1),
                                                     distance = distance,
                                                     scaleDistance = scaleDistance)
      }
      tfDist <- rowSums(conDist)
      names(tfDist) <- rownames(activity$mu_tc)
      ### comparisons wrt reference group
    } else if(contrast == "reference"){
      varsToCompare <- seq_len(ncol(U))[-referenceGroup]
      refDist <- matrix(NA, nrow=ncol(X), ncol=ncol(U)-1)
      for(kk in seq_len(length(varsToCompare))){
        refDist[,kk] <- .distanceCalculation_newMethods(mu_tc = activity$mu_tc,
                                                     cellID = c(referenceGroup, varsToCompare[kk]),
                                                     distance = distance,
                                                     scaleDistance = scaleDistance)
      }
      tfDist <- rowSums(refDist)
      names(tfDist) <- rownames(activity$mu_tc)
    }
  }

  return(tfDist)
}
)
