
mugtcOldToNew <- function(mu_gtc, design, X){

  tfRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 1))
  geneRows <- unlist(lapply(strsplit(rownames(mu_gtc), split=";"), "[[", 2))
  nGenes <- length(unique(geneRows))
  nTFs <- length(unique(tfRows))

  mu_gtcNew <- array(0, dim=c(nGenes, nTFs, ncol(mu_gtc)),
                     dimnames = list(sort(unique(geneRows)),
                                     sort(unique(tfRows)),
                                     NULL))
  for(cc in 1:ncol(mu_gtc)){
    mu_gtcNew[,,cc][cbind(geneRows,tfRows)] <- mu_gtc[,cc]
  }
  return(mu_gtcNew)
}


elementToRowCol <- function(element, nrows, ncols){
  rowColMat <- matrix(NA, nrow = length(element), ncol = 2)
  colnames(rowColMat) <- c("row", "column")
  for (ee in 1:length(element)) {
    col <- ceiling(element[ee]/nrows)
    row <- element[ee] - (nrows * (col - 1))
    rowColMat[ee, ] <- c(row, col)
  }
  return(rowColMat)
}

mugtcNewToOld <- function(mu_gtc, design, X){

  if(is.null(colnames(X))){
    colnames(X) <- paste0("tf", 1:ncol(X))
  }
  id <- which(X==1)
  rc <- elementToRowCol(id, nrow(X), ncol(X))
  rcNames <- cbind(rownames(X)[rc[,1]], colnames(X)[rc[,2]])

  mu_gtcOld <- matrix(0, nrow=length(id), ncol=dim(mu_gtc)[3],
                      dimnames=list(paste0(rcNames[,2],";",rcNames[,1]),
                                    NULL))
  for(cc in 1:dim(mu_gtc)[3]){
    mu_gtcOld[,cc] <- mu_gtc[,,cc][cbind(rcNames[,1],rcNames[,2])]
  }
  return(mu_gtcOld)
}


