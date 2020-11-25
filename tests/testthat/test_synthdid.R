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
