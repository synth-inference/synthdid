test_that("a simple workflow doesn't error", {
  setup = random.low.rank()
  tau.hat = synthdid_estimate(setup$Y,setup$N0,setup$T0)
  se = synthdid_se(tau.hat)

  print(paste0("point estimate: ", round(tau.hat, 2)))
  print(paste0("95% CI for tau: (", round(tau.hat - 1.96 * se, 2), ", ", round(tau.hat + 1.96 * se, 2), ")"))
  plot(tau.hat)

  expect_equal(1, 1)
})

test_that("adjustment for covariates works: random noise less influential if passed as covariate", {
  setup = random.low.rank()
  X = setup$Y - setup$L

  tau.hat = synthdid_estimate(setup$Y,setup$N0,setup$T0)
  tau.hat.noise = synthdid_estimate(setup$Y+X,setup$N0,setup$T0)
  tau.hat.cov = synthdid_estimate(setup$Y+X,setup$N0,setup$T0,X)
  se = synthdid_se(tau.hat)

  expect_lt(abs(tau.hat - tau.hat.cov), abs(tau.hat - tau.hat.noise))
})

test_that("default synthdid behavior has not changed", {
  setup = readRDS("setup.Rds")
  estimate = synthdid_estimate(setup$Y,setup$N0,setup$T0)
  weights = attr(estimate, 'weights')
  coef = c(estimate)

  weights.expected = readRDS("weights.expected.Rds")
  coef.expected = readRDS("coef.expected.Rds")

  expect_equal(coef, coef.expected)
  expect_equal(weights, weights.expected)

  # To update this test:
  # setup = random.low.rank()
  # saveRDS(setup, "setup.Rds")
  # saveRDS(coef, "coef.expected.Rds")
  # saveRDS(weights, "weights.expected.Rds")
})

test_that("invariances hold with default options", {
  # Test that three types of invariances hold, for details see
  # https://github.com/synth-inference/synthdid/issues/38
  estimators = list(sc_estimate, did_estimate, synthdid_estimate)
  setup = random.low.rank()

  # 1. Invariance to column fixed effect (all)
  # Re-mapping Yit <- Yit + bt for arbitrary bt doesn't change anything
  bt = 2 * matrix(1:ncol(setup$Y), nrow(setup$Y), ncol(setup$Y), byrow = TRUE)
  for (estimator in estimators) {
    estimate = estimator(setup$Y,
                         setup$N0,
                         setup$T0)
    estimate.se = synthdid_se(estimate)
    estimate.col.scaled = estimator(setup$Y + bt,
                                    setup$N0,
                                    setup$T0)
    estimate.se.col.scaled = synthdid_se(estimate.col.scaled)
    expect_equal(c(estimate), c(estimate.col.scaled))
    expect_equal(estimate.se, estimate.se.col.scaled)
  }

  # 2. Invariance to row fixed effect (exception: "SC")
  # Re-mapping Yit <- Yit + ai for arbitrary ai doesn't change anything
  ai = 2.5 * matrix(1:nrow(setup$Y), nrow(setup$Y), ncol(setup$Y), byrow = FALSE)
  for (estimator in estimators[-1]) {
    estimate = estimator(setup$Y,
                         setup$N0,
                         setup$T0)
    estimate.se = synthdid_se(estimate)
    estimate.row.scaled = estimator(setup$Y + ai,
                                    setup$N0,
                                    setup$T0)
    estimate.se.row.scaled = synthdid_se(estimate.row.scaled)
    expect_equal(c(estimate), c(estimate.row.scaled))
    expect_equal(estimate.se, estimate.se.row.scaled)
  }

  # 3. Invariance to scaling (all)
  # Re-mapping Yit <- c * Yit for c > 0 doesn't change weights and scales tau by c
  # 3.1: c is very large
  const = 10^6
  for (estimator in estimators) {
    estimate = estimator(setup$Y,
                         setup$N0,
                         setup$T0)
    estimate.se = synthdid_se(estimate)
    weights = attr(estimate, "weights")
    estimate.scaled = estimator(const * setup$Y,
                                setup$N0,
                                setup$T0)
    estimate.se.scaled = synthdid_se(estimate.scaled)
    weights.scaled = attr(estimate.scaled, "weights")
    expect_equal(c(const * estimate), c(estimate.scaled))
    expect_equal(const * estimate.se, estimate.se.scaled)
    expect_equal(weights[c("lambda", "omega")], weights.scaled[c("lambda", "omega")])
  }
  # 3.2: c is very small
  const = 10^-6
  for (estimator in estimators) {
    estimate = estimator(setup$Y,
                         setup$N0,
                         setup$T0)
    estimate.se = synthdid_se(estimate)
    weights = attr(estimate, "weights")
    estimate.scaled = estimator(const * setup$Y,
                                setup$N0,
                                setup$T0)
    estimate.se.scaled = synthdid_se(estimate.scaled)
    weights.scaled = attr(estimate.scaled, "weights")
    expect_equal(c(const * estimate), c(estimate.scaled))
    expect_equal(const * estimate.se, estimate.se.scaled)
    expect_equal(weights[c("lambda", "omega")], weights.scaled[c("lambda", "omega")])
  }
})
