#' Outputs a table of important synthetic controls and their corresponding weights omega, sorted by weight.
#' The table is truncated to exclude synthetic controls that do not matter for any estimate ---
#' for each estimate, the truncated controls may have total weight no larger that 1-mass.
#' @param estimates, a list of estimates output by synthdid_estimate. Or a single estimate.
#' @param sort.by, the index of the estimate to sort by. Defaults to 1.
#' @param digits,  the number of digits of weight to display. Defaults to 3.
#' @param mass, which controls the length of the table. Defaults to 0.9.
#' @export synthdid_controls
synthdid_controls = function(estimates, sort.by=1, digits=3, mass=.9) {
    if(class(estimates) == 'synthdid_estimate') { estimates = list(estimates) }
    if(is.null(names(estimates))) { names(estimates) = sprintf('estimate %d', 1:length(estimates)) }

    omegas = do.call(cbind, lapply(estimates, function(est) { attr(est, 'weights')$omega }))
    if(is.null(dim(omegas))) { dim(omegas) = c(length(omegas), 1) }

    Y = attr(estimates[[1]], 'setup')$Y
    o = rev(order(omegas[,sort.by]))
    tab = round(omegas[o,,drop=FALSE], digits=digits)
    rownames(tab) = rownames(Y)[o]
    colnames(tab) = names(estimates)

    # truncate table to retain a weight sum of at least mass for each unit
    tab.len = max(apply(tab, 2, function(col) { Position(function(x){ x >= mass }, cumsum(col)) }))
    tab[1:tab.len, ]
}

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
