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
#' @param ... Any valid arguments for vcov.synthdid_estimate
#' @export synthdid_se
synthdid_se = function(...) { sqrt(vcov(...)) }


# The bootstrap se: Algorithm 2 of Arkhangelsky et al.
bootstrap_se = function(estimate, replications) { sqrt((replications-1)/replications) * sd(bootstrap_sample(estimate, replications)) }
bootstrap_sample = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    if (setup$N0 == nrow(setup$Y) - 1) { return(NA) }
    theta = function(ind) {
	if(all(ind <= setup$N0) || all(ind > setup$N0)) { NA }
	else {
	    weights.boot = weights
	    weights.boot$omega = sum_normalize(weights$omega[sort(ind[ind <= setup$N0])])
	    do.call(synthdid_estimate, c(list(Y=setup$Y[sort(ind),], N0=sum(ind <= setup$N0), T0=setup$T0, X=setup$X[sort(ind), ,], weights=weights.boot), opts))
	}
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
# returns NA if there is one treated unit or, for the fixed-weights jackknife, one control with nonzero weight
jackknife_se = function(estimate, weights = attr(estimate, 'weights')) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    if (!is.null(weights)) {
      opts$update.omega = opts$update.lambda = FALSE
    }
    if (setup$N0 == nrow(setup$Y) - 1 || (!is.null(weights) && sum(weights$omega != 0) == 1)) { return(NA) }
    theta = function(ind) {
	weights.jk = weights
	if (!is.null(weights)) { weights.jk$omega = sum_normalize(weights$omega[ind[ind <= setup$N0]]) }
	estimate.jk = do.call(synthdid_estimate,
	    c(list(Y=setup$Y[ind, ], N0=sum(ind <= setup$N0), T0=setup$T0, X = setup$X[ind, , ], weights = weights.jk), opts))
    }
    jackknife(1:nrow(setup$Y), theta)
}

#' Jackknife standard error of function `theta` at samples `x`.
#' @param x vector of samples
#' @param theta a function which returns a scalar estimate
#' @importFrom stats var
#' @keywords internal
jackknife = function(x, theta) {
  n = length(x)
  u = rep(0, n)
  for (i in 1:n) {
    u[i] = theta(x[-i])
  }
  jack.se = sqrt(((n - 1) / n) * (n - 1) * var(u))

  jack.se
}



# The placebo se: Algorithm 4 of Arkhangelsky et al.
placebo_se = function(estimate, replications) {
    setup = attr(estimate, 'setup')
    opts = attr(estimate, 'opts')
    weights = attr(estimate, 'weights')
    N1 = nrow(setup$Y) - setup$N0
    if (setup$N0 <= N1) { stop('must have more controls than treated units to use the placebo se') }
    theta = function(ind) {
	N0 = length(ind)-N1
	weights.boot = weights
	weights.boot$omega = sum_normalize(weights$omega[ind[1:N0]])
        do.call(synthdid_estimate, c(list(Y=setup$Y[ind,], N0=N0,  T0=setup$T0,  X=setup$X[ind, ,], weights=weights.boot), opts))
    }
    sqrt((replications-1)/replications) * sd(replicate(replications, theta(sample(1:setup$N0))))
}

sum_normalize = function(x) {
    if(sum(x) != 0) { x / sum(x) }
    else { rep(1/length(x), length(x)) }
    # if given a vector of zeros, return uniform weights
    # this fine when used in bootstrap and placebo standard errors, where it is used only for initialization
    # for jackknife standard errors, where it isn't, we handle the case of a vector of zeros without calling this function.
}
