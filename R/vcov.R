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
#' @param replications, the number of bootstrap replications
#' @param ... Additional arguments (currently ignored).
#'
#' @references Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
#'  "Synthetic Difference in Differences". arXiv preprint arXiv:1812.09970, 2019.
#'
#' @method vcov synthdid_estimate
#' @export
vcov.synthdid_estimate = function(object,
  method = c("bootstrap", "jackknife", "placebo"),
  replications = 200, ...) {
    method = match.arg(method)
    if(method == 'bootstrap') { 
	se = bootstrap_se(object, replications)
    } else if(method == 'jackknife') {
	se = jackknife_se(object)
    } else if(method == 'placebo') {
	se = placebo_se(object, replications)
    }
    matrix(se^2)
}

#' Calculate the standard error of a synthetic diff in diff estimate. Deprecated. Use vcov.synthdid_estimate.
#' @export
synthdid_se = function(...) { sqrt(vcov(...)) }


# The bootstrap se: Algorithm 2 of Arkhangelsky et al.
bootstrap_se = function(estimate, replications) { sqrt((replications-1)/replications) * sd(bootstrap_sample(estimate, replications)) }
bootstrap_sample = function(estimate, replications) { 
    estimator = attr(estimate, 'estimator')
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
    theta = function(ind) {
	if(all(ind <= setup$N0) || all(ind > setup$N0)) { NA }
	else {do.call(estimator, c(list(Y=setup$Y[sort(ind),], N0=sum(ind <= setup$N0), T0=setup$T0, X=setup$X[sort(ind), ,]), opts))}
    }
    bootstrap.estimates = rep(NA, replications)
    count = 0
    while(count < replications) {
	bootstrap.estimates[count+1] = theta(sample(1:nrow(setup$Y), replace=TRUE))
	if(!is.na(bootstrap.estimates[count+1])) { count = count+1 }
    }
    bootstrap.estimates
}


# The fixed-weights jackknife estimate of variance: Algorithm 3 of Arkhangelsky et al.
# if weights = NULL is passed explicitly, calculates the usual jackknife estimate of variance.
jackknife_se = function(estimate, weights = attr(estimate, 'weights')) { 
    estimator = attr(estimate, 'estimator')
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
    sum_normalize = function(x) { x / sum(x) }
    theta = function(ind) {
	weights.jk = weights
	if (!is.null(weights)) { weights.jk$omega = sum_normalize(weights$omega[ind[ind <= setup$N0]]) }
	estimate.jk = do.call(estimator, 
	    c(list(Y=setup$Y[ind, ], N0=sum(ind <= setup$N0), T0=setup$T0, X = setup$X[ind, , ], weights = weights.jk), opts))
    }
    jackknife(1:nrow(setup$Y), theta)
}

# The placebo se: Algorithm 4 of Arkhangelsky et al.
placebo_se = function(estimate, replications) {
    estimator = attr(estimate, 'estimator')
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    N1 = nrow(setup$Y) - setup$N0
    if (setup$N0 <= N1) { stop('must have more controls than treated units to use the placebo se') }
    theta = function(ind) {
      do.call(estimator, c(list(Y=setup$Y[ind,], N0=length(ind)-N1,  T0=setup$T0,  X=setup$X[ind, ,]), opts))
    }
    sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:setup$N0))))
}


