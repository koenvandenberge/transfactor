#' @include utils.R
#' @include activityEstimationFunctions.R


#' @title Estimate TF activity
#' @description Estimate TF activity given a gene regulatory network.
#'
#' @param counts Gene expression counts of dimensions G x n.
#' @param X Gene regulatory network of dimensions G x T.
#' @param model The model to use. Options are \code{"poisson"}, \code{"dirMult"},
#' and \code{"dirMultEmpBayes"} (not recommended). Defaults to \code{"poisson"}.
#' @param U Design matrix, of dimensions n x p. The design matrix should not
#' contain an intercept.
#' @param alpha Prior belief in GRN edges.
#' @param alphaScale Scaling factor to tune the belief of the prior versus the
#' observed dataset. Setting \code{alphaScale=1} (the default) means assuming equal belief to
#' the prior as to the observed dataset. Setting \code{alphaScale=1/2}
#' (\code{alphaScale=2}) means you belief the prior half (twice) as much as the
#' observed data. Note that providing a very low number will result in many
#' estimate TF activities to equal zero. Indeed, due to the -1 in the
#' calculation for mu_gtc, many will be negative and will hence be reset to zero.
#' @param maxIter Maximum number of iterations.
#' @param plot Logical: should log-likelihood be plotted across iterations?
#' @param verbose Logical: should log-likelihood and iteration be printed?
#' @param epsilon Log-likelihood convergence criterion.
#' @param iterOLS Number of OLS iterations to use to set up sparse initialization.
#' Defaults to 0.
#' @param lassoFamily The family for the sparse initialization. Defaults to \code{"gaussian"}.
#' @param repressions Logical: should repressions be taken into account?
#' @param rho_t Thinning factor to account for repressions. \code{NULL} by default.
#' @param sparse Logical: should sparse initialization be used?
#' @return A list of objects, including
#'  - mu_tc: TF activity for each TF t in each group c (as defined by columns in U).
#'  - mu_gtc: The expected number of molecules for gene g produced by TF t in group c.
#' @examples
#' counts <- matrix(rpois(n= 1e4, lambda=4), nrow=100, ncol=100)
#' X <- matrix(0, nrow=100, ncol=12)
#' X[cbind(sample(100, size=250, replace=TRUE), sample(12, size=250, replace=TRUE))] <- 1
#' rownames(X) <- rownames(counts) <- paste0("gene",1:100)
#' cellType <- factor(rep(c("a", "b"), each=50))
#' U <- model.matrix(~ -1 + cellType)
#' act <- estimateActivity(counts, X, U=U)
#' @export
#' @rdname estimateActivity
#' @import SingleCellExperiment
setMethod(f = "estimateActivity",
          signature = c(counts = "matrix",
                        X = "matrix"),
          definition = function(counts,
                                X,
                                U = NULL,
                                model = "poisson",
                                alpha = NULL,
                                alphaScale = 1,
                                maxIter = 500,
                                plot = FALSE,
                                verbose = FALSE,
                                epsilon = 1e-2,
                                iterOLS = 0,
                                sparse = TRUE,
                                lassoFamily = "gaussian",
                                repressions = TRUE,
                                rho_t = NULL){
## for dev:
# U = U
# model = "poisson"
# alpha = NULL
# alphaScale = 1
# maxIter = 500
# plot = FALSE
# verbose = FALSE
# epsilon = 1e-2
# iterOLS = 0
# lassoFamily = "gaussian"
# repressions = TRUE
# rho_t = NULL
# sparse = TRUE


            # TODO: restrict range on alpha (lower bound)?
            # TODO: check alpha scaling e.g. when all alpha=1 no scaling happens.

            ## checks on X
            if(is.null(rownames(X)) | is.null(rownames(counts))){
              stop("Please ensure rownames for both X and counts.")
            }
            if(!all(rownames(X) %in% rownames(counts))){
              stop("Not all genes in X are in the count matrix.")
            }
            if(is.null(colnames(X))){
              colnames(X) <- paste0("tf", 1:ncol(X))
            }
            if(any(rowSums(abs(X)) == 0)){
              noLinkGenes <- which(rowSums(abs(X)) == 0)
              message("Removing ", length(noLinkGenes), " genes from the GRN",
                      " without any links.")
              X <- X[-noLinkGenes,]
            }
            if(any(colSums(abs(X)) == 0)){
              noLinkTFs <- which(colSums(abs(X)) == 0)
              message("Removing ", length(noLinkTFs), "TFs from the GRN",
                      " without any links.")
              X <- X[,-noLinkTFs]
            }

            ## checks on U
            if(is.null(U)){
              message("No design matrix provided. Working with intercept only.")
              ict <- rep(1, length = ncol(counts))
              U <- stats::model.matrix(~ -1 + ict)
            }

            ## checks on arguments
            if(repressions){
              if(!any(X == -1)){
                message("No repressions in the GRN. Setting repressions to FALSE.")
                repressions <- FALSE
              }
            }

            if(model == "poisson"){
              res <- poissonEstimation(counts = counts,
                                       X = X,
                                       U = U,
                                       maxIter = maxIter,
                                       plot = plot,
                                       verbose = verbose,
                                       epsilon = epsilon,
                                       iterOLS = iterOLS,
                                       lassoFamily = lassoFamily,
                                       repressions = repressions,
                                       rho_t = rho_t,
                                       sparse = sparse)
            } else if(model == "dirMult"){
              res <- dirMultEstimation2(counts = counts,
                                       X = X,
                                       U = U,
                                       alpha = alpha,
                                       alphaScale = alphaScale,
                                       nIters = maxIter,
                                       plot = plot,
                                       verbose = verbose,
                                       epsilon = epsilon,
                                       iterOLS = iterOLS,
                                       lassoFamily = lassoFamily,
                                       repressions = repressions,
                                       rho_t = rho_t,
                                       sparse = sparse)
            } else if(model == "dirMultAlpha"){
              res <- dirMultEstimationAlpha2(counts = counts,
                                       X = X,
                                       U = U,
                                       alpha = alpha,
                                       alphaScale = alphaScale,
                                       nIters = maxIter,
                                       plot = plot,
                                       verbose = verbose,
                                       epsilon = epsilon,
                                       iterOLS = iterOLS,
                                       lassoFamily = lassoFamily,
                                       repressions = repressions,
                                       rho_t = rho_t,
                                       sparse = sparse)
            }

            return(res)

          }
)
