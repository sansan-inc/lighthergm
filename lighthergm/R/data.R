#' A toy network to play `lighthergm` with.
#'
#' This network has a clear cluster structure.
#' The number of clusters is four, and which cluster each node belongs to is defined in the variable "block".
#'
#' @format A `statnet`'s network class object. It has three nodal features.
#' \describe{
#'   \item{block}{block membership of each node}
#'   \item{x}{a covariate. It has 10 labels.}
#'   \item{y}{a covariate. It has 10 labels.}
#'   ...
#' }
#' `x` and `y` are not variables with any particular meaning.
#'
"toyNet"
