#' Print a synthdid object
#' @param x The object to print
#' @param ... Additional arguments (currently ignored).
#' @method print synthdid_estimate
#' @export
print.synthdid_estimate = function(x, ...) { cat(format(x, ...), "\n") }
