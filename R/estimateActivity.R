#' @include utils.R
#' @include activityEstimationFunctions.R


#' @title Estimate TF activity
#' @description Estimate TF activity given a gene regulatory network.
#'
#' @param counts Gene expression counts of dimensions G x n.
#' @param X Gene regulatory network of dimensions G x T.
#' @param model The model to use. Options are \code{"poisson"}, \code{"dirMult"},
#' and \code{"dirMultEmpBayes"}. Defaults to \code{"poisson"}.
#' @param U Design matrix, of dimensions n x p.
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
#' @examples
#' x=2
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
                                lassoFamily = "gaussian",
                                repressions = TRUE,
                                rho_t = NULL,
                                sparse = TRUE){

            # TODO: restrict range on alpha (lower bound)?

            ## checks on X
            if(!all(rownames(X) %in% rownames(counts))){
              stop("Not all genes in X are in the count matrix.")
            }
            if(is.null(colnames(X))){
              colnames(X) <- paste0("tf", 1:ncol(X))
            }

            ## checks on U
            if(is.null(U)){
              message("No design matrix provided. Working with intercept only.")
              ict <- rep(1, length = ncol(counts))
              U <- stats::model.matrix(~ -1 + ict)
            }

            if(model == "poisson"){
              res <- poissonEstimation(counts = counts,
                                       X = X,
                                       U = U,
                                       maxIter = maxIter,
                                       plot = plot,
                                       verbose = verbose,
                                       epsilon = epsilon-2,
                                       iterOLS = iterOLS,
                                       lassoFamily = lassoFamily,
                                       repressions = repressions,
                                       rho_t = rho_t,
                                       sparse = sparse)
            } else if(model == "dirMult"){
              res <- dirMultEstimation(counts = counts,
                                       X = X,
                                       U = U,
                                       alpha = alpha,
                                       alphaScale = alphaScale,
                                       maxIter = maxIter,
                                       plot = plot,
                                       verbose = verbose,
                                       epsilon = epsilon-2,
                                       iterOLS = iterOLS,
                                       lassoFamily = lassoFamily,
                                       repressions = repressions,
                                       rho_t = rho_t,
                                       sparse = sparse)
            }

            return(res)

          }
)
