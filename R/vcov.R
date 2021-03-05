#' Calculate Variance-Covariance Matrix for a Fitted Model Object
#'
#' Provides variance estimates based on the following three options
#' \itemize{
#'   \item The bootstrap, Algorithm 2 in Arkhangelsky et al.
#'   \item The jackknife, Algorithm 3 in Arkhangelsky et al.
#'   \item Placebo, Algorithm 4 in Arkhangelsky et al.
#' }
#'
#' The jackknife is not recommended for SC, see section 5 in Arkhangelsky et al.
#' "placebo" is the only option that works for only one treated unit.
#'
#' @param object A synthdid model
#' @param method, the CI method. The default is bootstrap (warning: this may be slow on large
#'  data sets, the jackknife option is the fastest, with the caveat that it is not recommended
#'  for SC).
#' @param weights, like attr(estimate, 'weights')
#' @param replications, the number of bootstrap replications
#' @param ... Additional arguments (currently ignored).
#'
#' @references Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
#'  "Synthetic Difference in Differences". arXiv preprint arXiv:1812.09970, 2019.
#'
#' @method vcov synthdid_estimate
#' @export
vcov.synthdid_estimate = function(
  object,
  method = c("bootstrap", "jackknife", "placebo"),
  weights = attr(object, 'weights'),
  replications = 200, ...) {

    matrix(synthdid_se(object, method, weights, replications)^2)
}
