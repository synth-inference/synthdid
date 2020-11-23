#' Summarize a synthdid object
#' @param object The object to summarize
#' @param ... Additional arguments (currently ignored).
#' @method summary synthdid_estimate
#' @export
summary.synthdid_estimate = function(object, ...) {
    list(estimate = c(object),
         se = synthdid_se(object),
         controls = synthdid_controls(object))
}
