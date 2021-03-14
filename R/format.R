#' Format a synthdid object
#' @param x The object to format
#' @param ... Additional arguments (currently ignored).
#' @method format synthdid_estimate
#' @export
format.synthdid_estimate = function(x, ...) {
  info = summary(x, fast=TRUE)
  d = as.list(info$dimensions)
  sprintf('synthdid: %1.3f +- %1.3f. Effective N0/N0 = %1.1f/%d~%1.1f. Effective T0/T0 = %1.1f/%d~%1.1f. N1,T1 = %d,%d.',
    c(x), 1.96*info$se, 
    d$N0.effective, d$N0, d$N0.effective/d$N0,
    d$T0.effective, d$T0, d$T0.effective/d$T0, 
    d$N1, d$T1)
}
