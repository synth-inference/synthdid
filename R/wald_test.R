#' Pre-Trends Wald Test for Synthetic Difference-in-Differences Estimates
#'
#' Tests the null hypothesis of no pre-treatment trends (H0: tau_\{g,t\} = 0 for all
#' t < g) using a joint Wald statistic. The statistic is asymptotically chi-squared
#' with degrees of freedom equal to the number of pre-treatment (g, t) pairs.
#'
#' Optionally accepts a contrast matrix \code{R} to test arbitrary linear
#' combinations of the pre-treatment ATTs (e.g., testing equality across cohorts).
#'
#' @section Covariance scaling convention:
#' This function requires the caller to be explicit about how \code{vcov_matrix_all_atts}
#' is scaled via \code{vcov_is_normalized}:
#' \itemize{
#'   \item \code{vcov_is_normalized = TRUE} (default): \code{vcov_matrix_all_atts} is the
#'     normalized asymptotic covariance, i.e. \code{Cov(att)}. The Wald statistic is
#'     \code{n * t(preatt) \%*\% solve(preV) \%*\% preatt}.
#'   \item \code{vcov_is_normalized = FALSE}: \code{vcov_matrix_all_atts} is already
#'     \code{n * Cov(att)} (the unnormalized form used in Callaway & Sant'Anna 2021
#'     internally). The Wald statistic is \code{t(preatt) \%*\% solve(preV) \%*\% preatt},
#'     and \code{n_scaling_factor} is ignored.
#' }
#' Passing the wrong convention will inflate or deflate the statistic by a factor of n.
#'
#' @param att_estimates Numeric vector of ATT(g, t) estimates for all (cohort, time) pairs.
#' @param cohort_treatment_time Numeric vector of cohort first-treatment times (g), same
#'   length as \code{att_estimates}.
#' @param calendar_time_period Numeric vector of calendar time periods (t), same length
#'   as \code{att_estimates}.
#' @param vcov_matrix_all_atts Numeric square matrix of dimension k x k (k =
#'   \code{length(att_estimates)}). See the Covariance scaling convention section.
#' @param std_errors_all_atts Numeric vector of standard errors, same length as
#'   \code{att_estimates}. Used only for pre-filtering degenerate cells; must be
#'   consistent with \code{diag(vcov_matrix_all_atts)}.
#' @param n_scaling_factor A single positive numeric value giving the sample size n
#'   used to scale the Wald statistic. Only used when \code{vcov_is_normalized = TRUE}.
#' @param vcov_is_normalized Logical. If \code{TRUE} (default), \code{vcov_matrix_all_atts}
#'   is the normalized \code{Cov(att)} and the statistic is multiplied by
#'   \code{n_scaling_factor}. If \code{FALSE}, \code{vcov_matrix_all_atts} is already
#'   n-scaled and no multiplication is applied.
#' @param R Optional numeric contrast matrix of dimension q x (number of pre-treatment
#'   indices). If \code{NULL} (default), tests H0: preatt = 0 directly (R = identity).
#'   Supply R to test arbitrary linear combinations, e.g. equality across cohorts.
#' @param r Optional numeric vector of length q giving null hypothesis values for
#'   \code{R \%*\% preatt}. Defaults to a zero vector.
#' @param rcond_threshold Threshold for the reciprocal condition number of the contrast
#'   covariance matrix below which it is considered numerically singular. Defaults to
#'   \code{.Machine$double.eps}.
#' @param verbose Logical. If \code{TRUE} (default), prints a status message with the
#'   test result. Set to \code{FALSE} to suppress output in production use.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{W}}{The Wald chi-squared statistic, or \code{NULL} if not computed.}
#'     \item{\code{Wpval}}{The p-value from \code{pchisq(W, df = q)}, or \code{NULL}.}
#'     \item{\code{q}}{Degrees of freedom (0 if test was not run due to an error).}
#'     \item{\code{message}}{A character string describing the outcome.}
#'   }
#'
#' @references
#' Callaway, B. and Sant'Anna, P.H.C. (2021). "Difference-in-Differences with Multiple
#' Time Periods." Journal of Econometrics, 225(2), 200-230.
#'
#' @importFrom stats pchisq
#' @export
calculate_wald_pre_test <- function(att_estimates,
                                    cohort_treatment_time,
                                    calendar_time_period,
                                    vcov_matrix_all_atts,
                                    std_errors_all_atts,
                                    n_scaling_factor,
                                    vcov_is_normalized = TRUE,
                                    R = NULL,
                                    r = NULL,
                                    rcond_threshold = .Machine$double.eps,
                                    verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # 1. Input validation
  # ---------------------------------------------------------------------------
  k <- length(att_estimates)

  if (length(cohort_treatment_time) != k ||
      length(calendar_time_period)  != k ||
      length(std_errors_all_atts)   != k ||
      nrow(vcov_matrix_all_atts)    != k ||
      ncol(vcov_matrix_all_atts)    != k) {
    stop("Input vectors/matrix have inconsistent lengths or dimensions.")
  }

  if (!is.logical(vcov_is_normalized) || length(vcov_is_normalized) != 1) {
    stop("'vcov_is_normalized' must be a single logical value (TRUE or FALSE).")
  }

  if (vcov_is_normalized) {
    if (missing(n_scaling_factor)) {
      stop("'n_scaling_factor' must be supplied when 'vcov_is_normalized = TRUE'.")
    }
    if (!is.numeric(n_scaling_factor) || length(n_scaling_factor) != 1 || n_scaling_factor <= 0) {
      stop("'n_scaling_factor' must be a single positive numeric value.")
    }
  }

  # Warn if std_errors_all_atts and sqrt(diag(vcov)) disagree — they should
  # be consistent since both come from the same estimation.
  # Use diag(vcov) as the authoritative source for all downstream filtering.
  vcov_se <- sqrt(pmax(diag(as.matrix(vcov_matrix_all_atts)), 0))
  se_discrepancy <- max(abs(std_errors_all_atts - vcov_se), na.rm = TRUE)
  if (se_discrepancy > sqrt(.Machine$double.eps) * 100) {
    warning(sprintf(
      paste0(
        "'std_errors_all_atts' and sqrt(diag(vcov_matrix_all_atts)) differ by up to %.2e. ",
        "Degenerate-cell filtering will use diag(vcov_matrix_all_atts) for consistency."
      ),
      se_discrepancy
    ))
  }

  # ---------------------------------------------------------------------------
  # 2. Identify pre-treatment indices: t < g  <=>  g > t
  # ---------------------------------------------------------------------------
  pre_indices <- which(cohort_treatment_time > calendar_time_period)

  # Initialise return values; q = 0 signals "test not run"
  W_stat <- NULL
  W_pval <- NULL
  q      <- 0
  status <- ""

  if (length(pre_indices) == 0) {
    status <- "No pre-treatment periods found to test."
    if (verbose) message(status)
    return(list(W = W_stat, Wpval = W_pval, q = q, message = status))
  }

  # ---------------------------------------------------------------------------
  # 3. Filter degenerate pre-treatment cells using diag(vcov) — not std_errors,
  #    which may have been computed separately and could disagree.
  # ---------------------------------------------------------------------------
  se_eps    <- sqrt(.Machine$double.eps) * 10
  bad_global <- which(is.na(vcov_se) | vcov_se <= se_eps)
  pre_indices <- pre_indices[!(pre_indices %in% bad_global)]

  if (length(pre_indices) == 0) {
    status <- "No valid pre-treatment periods after filtering for NA/zero SEs."
    if (verbose) message(status)
    return(list(W = W_stat, Wpval = W_pval, q = q, message = status))
  }

  # ---------------------------------------------------------------------------
  # 4. Extract pre-treatment ATTs and covariance sub-matrix
  # ---------------------------------------------------------------------------
  preatt <- as.matrix(att_estimates[pre_indices])                          # p x 1
  preV   <- as.matrix(vcov_matrix_all_atts[pre_indices, pre_indices,
                                            drop = FALSE])                 # p x p
  p      <- length(pre_indices)

  # ---------------------------------------------------------------------------
  # 5. Apply contrast matrix R  (default: identity — test preatt = 0 directly)
  # ---------------------------------------------------------------------------
  if (is.null(R)) {
    R <- diag(p)
  } else {
    R <- as.matrix(R)
    if (ncol(R) != p) {
      stop(sprintf(
        "'R' must have %d columns (one per pre-treatment period), but has %d.", p, ncol(R)
      ))
    }
  }

  if (is.null(r)) {
    r <- rep(0, nrow(R))
  } else {
    r <- as.numeric(r)
    if (length(r) != nrow(R)) {
      stop(sprintf(
        "'r' must have length %d (= nrow(R)), but has length %d.", nrow(R), length(r)
      ))
    }
  }

  q <- nrow(R)   # degrees of freedom

  # ---------------------------------------------------------------------------
  # 6. Feasibility checks on the contrast covariance  V_c = R preV R'
  # ---------------------------------------------------------------------------
  if (anyNA(preV)) {
    status <- "Not returning pre-test Wald statistic: NA values in pre-treatment covariance sub-matrix."
    warning(status)
    return(list(W = W_stat, Wpval = W_pval, q = 0, message = status))
  }

  V_contrast <- R %*% preV %*% t(R)   # q x q

  if (rcond(V_contrast) <= rcond_threshold) {
    status <- "Not returning pre-test Wald statistic: contrast covariance matrix (R preV R') is singular."
    warning(status)
    return(list(W = W_stat, Wpval = W_pval, q = 0, message = status))
  }

  V_inv <- try(solve(V_contrast), silent = TRUE)
  if (inherits(V_inv, "try-error")) {
    status <- "Not returning pre-test Wald statistic: failed to invert contrast covariance matrix."
    warning(status)
    return(list(W = W_stat, Wpval = W_pval, q = 0, message = status))
  }

  # ---------------------------------------------------------------------------
  # 7. Wald statistic
  #
  #   contrast = R %*% preatt - r                          (q x 1)
  #   W_base   = contrast' (R preV R')^{-1} contrast
  #
  #   vcov_is_normalized = TRUE  => preV = Cov(att),  W = n * W_base  ~ chi-sq(q)
  #   vcov_is_normalized = FALSE => preV = n*Cov(att), W = W_base     ~ chi-sq(q)
  # ---------------------------------------------------------------------------
  contrast <- as.numeric(R %*% preatt) - r
  W_stat   <- as.numeric(t(contrast) %*% V_inv %*% contrast)

  if (vcov_is_normalized) {
    W_stat <- n_scaling_factor * W_stat
  }

  W_pval <- 1 - pchisq(W_stat, df = q)

  status <- sprintf(
    "Wald pre-trends test: W = %.4f, df = %d, p-value = %.6g",
    W_stat, q, W_pval
  )
  if (verbose) message(status)

  return(list(W = W_stat, Wpval = W_pval, q = q, message = status))
}
