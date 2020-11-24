#' Format a synthdid object
#' @param x The object to format
#' @param ... Additional arguments (currently ignored).
#' @method format synthdid_estimate
#' @export
format.synthdid_estimate = function(x, ...) {
  setup = attr(x, 'setup')
  weights = attr(x, 'weights')
  sprintf('synthdid: %1.3f. Effective N0/N0 = %1.1f/%d. Effective T0/T0 = %1.1f/%d. N1,T1 = %d,%d.',
    c(x), 1 / sum(weights$omega^2), setup$N0, 1 / sum(weights$lambda^2), setup$T0,
    nrow(setup$Y) - setup$N0, ncol(setup$Y) - setup$T0)
}
