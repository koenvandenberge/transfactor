# TODO: allow providing only pigtc or mugtc as input.

#' @title Get counts produced by each transcription factor.
#' @description Get counts produced by each transcription factor.
#'
#' @param mu_gtc Output from estimateActivity.
#' @param counts Gene expression counts.
#' @param design Design matrix.
#' @return A matrix of estimated counts, with rows corresponding to transcription
#' factors and columns to cells.
#' @examples
#' counts <- matrix(rpois(n= 1e4, lambda=4), nrow=100, ncol=100)
#' X <- matrix(0, nrow=100, ncol=12)
#' X[cbind(sample(100, size=250, replace=TRUE), sample(12, size=250, replace=TRUE))] <- 1
#' rownames(X) <- rownames(counts) <- paste0("gene",1:100)
#' act <- estimateActivity(counts, X)
#' Y_ti <- tfCounts(act$mu_gtc, counts)
#' @export
#' @rdname tfCounts
setMethod(f = "tfCounts",
          signature = c(counts = "matrix",
                        pi_gtc = NULL,
                        mu_gtc = "matrix"),
          definition = function(counts,
                                pi_gtc=NULL,
                                mu_gtc,
                                design=NULL){

   if(!(is.null(pi_gtc) | is.null(mu_gtc))){
     stop("Either provide pi_gtc or mu_gtc, not both.")
   }

  if(is.null(design)){
    message("No design matrix provided. Working with intercept only.")
    ict <- rep(1, length = ncol(counts))
    design <- stats::model.matrix(~ -1 + ict)
  }

  if(nrow(design) != ncol(counts)){
    stop("Dimensions of design matrix and count matrix don't match.")
  }

  tfRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 1))
  geneRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 2))
  tfUniq <- unique(tfRows)
  geneUniq <- unique(geneRows)
  lvl <- unlist(apply(design,1, function(row){
    which(row == 1)
  }))
  colnames(mu_gtc) <- paste0("ct",1:ncol(mu_gtc))
  mu_gtc_tibble <- suppressWarnings(tibble::as_tibble(mu_gtc))

  Y_ti <- matrix(0, nrow=length(tfUniq), ncol=ncol(counts),
                 dimnames=list(tfUniq, colnames(counts)))
  for(gg in 1:length(geneUniq)){
    curGene <- geneUniq[gg]
    id <- which(geneRows == curGene)
    curTFs <- tfRows[id]
    curMu <- as.matrix(mu_gtc_tibble[id,1:ncol(design)])
    #curMu <- mu_gtc[paste0(curTFs,";",curGene),,drop=FALSE]
    curProbs <- sweep(curMu, 2, colSums(curMu)+1e-10, "/")
    curProbs <- curProbs[,lvl,drop=FALSE]
    Y_ti[curTFs,] <- Y_ti[curTFs,] + sweep(curProbs, 2, counts[curGene,], "*")
  }
  return(Y_ti)
}
)


#' @name tfCounts
setMethod(f = "tfCounts",
          signature = c(counts = "matrix",
                        pi_gtc = "array",
                        mu_gtc = NULL),
          definition = function(counts,
                                pi_gtc,
                                mu_gtc=NULL,
                                design=NULL){

            if(!(is.null(pi_gtc) | is.null(mu_gtc))){
              stop("Either provide pi_gtc or mu_gtc, not both.")
            }

            if(is.null(design)){
              message("No design matrix provided. Working with intercept only.")
              ict <- rep(1, length = ncol(counts))
              design <- stats::model.matrix(~ -1 + ict)
            }

            if(nrow(design) != ncol(counts)){
              stop("Dimensions of design matrix and count matrix don't match.")
            }

            lvl <- unlist(apply(design,1, function(row){
              which(row == 1)
            }))

            Y_ti <- matrix(0, nrow=ncol(pi_gtc), ncol=ncol(counts))
            for(cc in 1:dim(pi_gtc)[3]){
              Y_ti[,lvl==cc] <- t(pi_gtc[,,cc]) %*% counts[,lvl==cc]
            }

            return(Y_ti)
          }
)
