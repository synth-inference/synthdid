# collapse Y to an N0+1 x T0+1 vector by averaging the last N1=nrow(Y)-N0 rows and T1=ncol(Y)-T0 columns
collapsed.form = function(Y, N0, T0) {
  N = nrow(Y); T = ncol(Y)
  rbind(cbind(Y[1:N0, 1:T0, drop = FALSE], rowMeans(Y[1:N0, (T0 + 1):T, drop = FALSE])),
    cbind(t(colMeans(Y[(N0 + 1):N, 1:T0, drop = FALSE])), mean(Y[(N0 + 1):N, (T0 + 1):T, drop = FALSE])))
}

# return the component-wise sum of decreasing vectors in which NA is taken to mean that the vector has stopped decreasing
# and we can use the last non-na element. Where both are NA, leave as NA.
pairwise.sum.decreasing = function(x, y) {
  na.x = is.na(x)
  na.y = is.na(y)
  x[is.na(x)] = min(x[!na.x])
  y[is.na(y)] = min(y[!na.y])
  pairwise.sum = x + y
  pairwise.sum[na.x & na.y] = NA
  pairwise.sum
}

#' Convert a long (balanced) panel to a wide matrix
#'
#' Converts a data set in panel form to matrix format required by synthdid estimators.
#' A typical long panel date set looks like \[unit, time, outcome, treatment\]. Synthdid
#' requires a balanced panel with simultaneous adoption of treatment: each unit must be observed
#' at all times, and all treated units must begin treatment simultaneosly. This function
#' creates num.units x num.time.periods matrices Y and W of outcomes and treatment indicators.
#' In these matrices, columns are sorted by time, and by default (when treated.last=TRUE),
#' rows for control units appear before those of treated units. 
#'
#' @param panel A data.frame with columns consisting of units, time, outcome, and treatment indicator.
#' @param unit The column number/name corresponding to the unit identifier. Default is 1.
#' @param time The column number/name corresponding to the time identifier. Default is 2.
#' @param outcome The column number/name corresponding to the outcome identifier. Default is 3.
#' @param treatment The column number/name corresponding to the treatment status. Default is 4.
#' @param treated.last Should we sort the rows of Y and W so treated units are last. If FALSE, sort by unit number/name. Default is TRUE.
#' @return A list with entries `Y`: the data matrix, `N0`: the number of control units, `T0`:
#'  the number of time periods before treatment, `W`: the matrix of treatment indicators.
#'
#' @examples
#' \donttest{
#' # Load tobacco sales in long panel format.
#' data("california_prop99")
#' # Transform to N*T matrix format required for synthdid,
#' # where N is the number of units and T the time periods.
#' setup <- panel.matrices(california_prop99, unit = 1, time = 2, outcome = 3, treatment = 4)
#'
#' # Compute synthdid estimate
#' synthdid_estimate(setup$Y, setup$N0, setup$T0)
#' }
#'
#' @export
panel.matrices = function(panel, unit = 1, time = 2, outcome = 3, treatment = 4, treated.last = TRUE) {
  # TODO: add support for covariates X, i.e. could keep all other columns
  keep = c(unit, time, outcome, treatment)
  if (!all(keep %in% 1:ncol(panel) | keep %in% colnames(panel))) {
    stop("Column identifiers should be either integer or column names in `panel`.")
  }
  index.to.name = function(x) { if(x %in% 1:ncol(panel)) { colnames(panel)[x] } else { x } }
  unit = index.to.name(unit)
  time = index.to.name(time)
  outcome = index.to.name(outcome)
  treatment = index.to.name(treatment)
  keep = c(unit, time, outcome, treatment)

  panel = panel[keep]
  if (!is.data.frame(panel)){
    stop("Unsupported input type `panel.`")
  }
  if (anyNA(panel)) {
    stop("Missing values in `panel`.")
  }
  if (length(unique(panel[, treatment])) == 1) {
    stop("There is no variation in treatment status.")
  }
  if (!all(panel[, treatment] %in% c(0, 1))) {
    stop("The treatment status should be in 0 or 1.")
  }
  # Convert potential factor/date columns to character
  panel = data.frame(
    lapply(panel, function(col) {if (is.factor(col) || inherits(col, "Date")) as.character(col) else col}), stringsAsFactors = FALSE
  )
  val <- as.vector(table(panel[, unit], panel[, time]))
  if (!all(val == 1)) {
    stop("Input `panel` must be a balanced panel: it must have an observation for every unit at every time.")
  }

  panel = panel[order(panel[, unit], panel[, time]), ]
  num.years = length(unique(panel[, time]))
  num.units = length(unique(panel[, unit]))
  Y = matrix(panel[,outcome], num.units, num.years, byrow = TRUE,
             dimnames = list(unique(panel[,unit]), unique(panel[,time])))
  W = matrix(panel[,treatment], num.units, num.years, byrow = TRUE,
             dimnames = list(unique(panel[,unit]), unique(panel[,time])))
  w = apply(W, 1, any)                         # indicator for units that are treated at any time
  T0 = unname(which(apply(W, 2, any))[1]-1)    # last period nobody is treated
  N0 = sum(!w)

  if(! (all(W[!w,] == 0) && all(W[,1:T0] == 0) && all(W[w, (T0+1):ncol(Y)]==1))) {
    stop("The package cannot use this data. Treatment adoption is not simultaneous.")
  }

  unit.order = if(treated.last) { order(W[,T0+1], rownames(Y)) } else { 1:nrow(Y) }
  list(Y = Y[unit.order, ], N0 = N0, T0 = T0, W = W[unit.order, ])
}

#' Get timesteps from panel matrix Y
#'
#' timesteps are stored as colnames(Y), but column names cannot be Date objects.
#' Instead, we use strings. If they are strings convertible to dates, return that
#'
#' @param Y a matrix 
#' @return its column names interpreted as Dates if possible
#' @export
timesteps = function(Y) {
    tryCatch({
	as.Date(colnames(Y))
    }, error = function(e) { colnames(Y) })
}


## define some convenient accessors
setOldClass("synthdid_estimate")
get_slot = function(name) { function(object) { object[[name]] } }
setGeneric('weights')
setGeneric('Y',      get_slot('Y'))
setGeneric('lambda', get_slot('lambda'))
setGeneric('omega',  get_slot('omega'))
setMethod(weights, signature='synthdid_estimate',  definition=function(object) { attr(object, 'weights') })
setMethod(Y,       signature='synthdid_estimate',  definition=function(object) { attr(object, 'setup')$Y })
setMethod(lambda,  signature='synthdid_estimate',  definition=function(object) { lambda(weights(object)) })
setMethod(omega,   signature='synthdid_estimate',  definition=function(object) { omega(weights(object))  })


# A convenience function for generating data for unit tests.
random.low.rank = function() {
  n_0 <- 100
  n_1 <- 10
  T_0 <- 120
  T_1 <- 20
  n <- n_0 + n_1
  T <- T_0 + T_1
  tau <- 1
  sigma <- .5
  rank <- 2
  rho <- 0.7
  var <- outer(1:T, 1:T, FUN=function(x, y) rho^(abs(x-y)))

  W <- (1:n > n_0) %*% t(1:T > T_0)
  U <- matrix(rpois(rank * n, sqrt(sample(1:n)) / sqrt(n)), n, rank)
  V <- matrix(rpois(rank * T, sqrt(1:T) / sqrt(T)), T, rank)
  alpha <- outer(10*sample(1:n)/n, rep(1,T))
  beta <-  outer(rep(1,n), 10*(1:T)/T)
  mu <- U %*% t(V) + alpha + beta
  error <- mvtnorm::rmvnorm(n, sigma = var, method = "chol")
  Y <- mu + tau * W  + sigma * error
  rownames(Y) = 1:n
  colnames(Y) = 1:T
  list(Y=Y, L=mu, N0=n_0, T0=T_0)
}
