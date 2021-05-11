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

#' @export
#' @name tfCounts
#' @title tfCounts
#' @param ... parameters including:
setGeneric(
  name = "tfCounts",
  def = function(mu_gtc = "matrix",
                 counts = "matrix",
                 ...) {
    standardGeneric("tfCounts")
  }
)

#' @export
#' @name tfDistance
#' @title tfDistance
#' @param ... parameters including:
setGeneric(
  name = "tfDistance",
  signature = c('activity',
                'X',
                'counts',
                'U'),
  def = function(activity = "list",
                 X = "matrix",
                 counts = "matrix",
                 U = "matrix",
                 ...) {
    standardGeneric("tfDistance")
  }
)
