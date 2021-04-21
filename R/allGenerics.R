#' @export
#' @name estimateActivity
#' @title estimateActivity
#' @param ... parameters including:
setGeneric(
  name = "estimateActivity",
  signature = c('counts', 'X'),
  def = function(counts, X, ...) {
    standardGeneric("estimateActivity")
  }
)
