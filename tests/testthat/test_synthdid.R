test_that("a simple workflow doesn't error", {
  setup = random.low.rank()
  tau.hat = synthdid_estimate(setup$Y,setup$N0,setup$T0)
  se = sqrt(vcov(tau.hat, replications = 10))
  se.jackknife = sqrt(vcov(tau.hat, method='jackknife'))
  se.placebo   = sqrt(vcov(tau.hat, method='placebo', replications = 10))

  print(tau.hat)
  summary(tau.hat, fast = TRUE)
  plot(tau.hat)

  expect_equal(1, 1)
})

test_that("plotting doesn't error with (i) dates as colnames (ii) spaghetti units", {
  data(california_prop99)

  california_prop99$date = as.Date(sprintf('%04d/%02d/%02d', california_prop99$Year, 1, 1))
  setup = panel.matrices(california_prop99, time='date')
  estimate = synthdid_estimate(setup$Y, setup$N0, setup$T0)
  top.controls = synthdid_controls(estimate)[1:10, , drop=FALSE]
  plot(estimate, spaghetti.units=rownames(top.controls))

  expect_equal(1, 1)
})

test_that("adjustment for covariates works: random noise less influential if passed as covariate", {
  setup = random.low.rank()
  X = setup$Y - setup$L

  tau.hat = synthdid_estimate(setup$Y,setup$N0,setup$T0)
  tau.hat.noise = synthdid_estimate(setup$Y+X,setup$N0,setup$T0)
  tau.hat.cov = synthdid_estimate(setup$Y+X,setup$N0,setup$T0,X)
  se = synthdid_se(tau.hat, replications = 10)

  expect_lt(abs(tau.hat - tau.hat.cov), abs(tau.hat - tau.hat.noise))
})

test_that("column/row/scaling invariances hold with default options", {
  # Test that three types of invariances hold, for details see
  # https://github.com/synth-inference/synthdid/issues/38
  estimators = list(sc_estimate, did_estimate, synthdid_estimate)
  CI.methods = c("jackknife", "bootstrap", "placebo")
  setup = random.low.rank()
  seed = sample(1:1e6, 1)

  # 1. Invariance to column fixed effect (all)
  # Re-mapping Yit <- Yit + bt for arbitrary bt doesn't change anything
  bt = 2 * matrix(1:ncol(setup$Y), nrow(setup$Y), ncol(setup$Y), byrow = TRUE)
  for (estimator in estimators) {
    for (CI.method in CI.methods) {
      estimate = estimator(setup$Y, setup$N0, setup$T0)
      set.seed(seed); estimate.se = synthdid_se(estimate, method = CI.method, replications = 10)
      estimate.col.scaled = estimator(setup$Y + bt, setup$N0, setup$T0)
      set.seed(seed); estimate.se.col.scaled = synthdid_se(estimate.col.scaled, method = CI.method, replications = 10)
      expect_equal(c(estimate), c(estimate.col.scaled))
      expect_equal(estimate.se, estimate.se.col.scaled)
    }
  }

  # 2. Invariance to row fixed effect (exception: "SC")
  # Re-mapping Yit <- Yit + ai for arbitrary ai doesn't change anything
  ai = 2.5 * matrix(1:nrow(setup$Y), nrow(setup$Y), ncol(setup$Y), byrow = FALSE)
  for (estimator in estimators[-1]) {
    for (CI.method in CI.methods) {
      estimate = estimator(setup$Y, setup$N0, setup$T0)
      set.seed(seed); estimate.se = synthdid_se(estimate, method = CI.method, replications = 10)
      estimate.row.scaled = estimator(setup$Y + ai, setup$N0, setup$T0)
      set.seed(seed); estimate.se.row.scaled = synthdid_se(estimate.row.scaled, method = CI.method, replications = 10)
      expect_equal(c(estimate), c(estimate.row.scaled))
      expect_equal(estimate.se, estimate.se.row.scaled)
    }
  }

  # 3. Invariance to scaling (all)
  # Re-mapping Yit <- c * Yit for c > 0 doesn't change weights and scales tau by c
  # 3.1: c is very small
  c.small = 10^-6
  for (estimator in estimators) {
    for (CI.method in CI.methods) {
      estimate = estimator(setup$Y, setup$N0, setup$T0)
      set.seed(seed); estimate.se = synthdid_se(estimate, method = CI.method, replications = 10)
      weights = attr(estimate, "weights")
      estimate.scaled = estimator(c.small * setup$Y, setup$N0, setup$T0)
      set.seed(seed); estimate.se.scaled = synthdid_se(estimate.scaled, method = CI.method, replications = 10)
      weights.scaled = attr(estimate.scaled, "weights")
      expect_equal(c(c.small * estimate), c(estimate.scaled))
      expect_equal(c.small * estimate.se, estimate.se.scaled)
      expect_equal(weights[c("lambda", "omega")], weights.scaled[c("lambda", "omega")])
    }
  }
  # 3.2: c is very large
  c.large = 10^6
  for (estimator in estimators) {
    for (CI.method in CI.methods) {
      estimate = estimator(setup$Y, setup$N0, setup$T0)
      set.seed(seed); estimate.se = synthdid_se(estimate, method = CI.method, replications = 10)
      weights = attr(estimate, "weights")
      estimate.scaled = estimator(c.large * setup$Y, setup$N0, setup$T0)
      set.seed(seed); estimate.se.scaled = synthdid_se(estimate.scaled, method = CI.method, replications = 10)
      weights.scaled = attr(estimate.scaled, "weights")
      expect_equal(c(c.large * estimate), c(estimate.scaled))
      expect_equal(c.large * estimate.se, estimate.se.scaled)
      expect_equal(weights[c("lambda", "omega")], weights.scaled[c("lambda", "omega")])
    }
  }
})

test_that("treated effect shifts correctly with scalar shifts to the 4 blocks", {
  # Test that four types of invariances hold, for details see
  # https://github.com/synth-inference/synthdid/issues/43
  estimators = list(sc_estimate, did_estimate, synthdid_estimate)
  seed = sample(1:1e6, 1)
  setup = random.low.rank()
  T0 = setup$T0
  N0 = setup$N0
  T = ncol(setup$Y)
  exposed = 1:nrow(setup$Y) > N0
  Y.orig = setup$Y

  lambda = function(obj) { attr(obj, 'weights')$lambda }
  omega  = function(obj) { attr(obj, 'weights')$omega }

  for (c in c(1e-6, 0.25, 1e6)) {
    # 1.
    # Re-mapping Yit <- Yit + c for exposed t > T0 increases tau by c
    Y1 <- setup$Y
    Y1[exposed, (T0+1):T] <- c + Y1[exposed, (T0+1):T]
    for (estimator in estimators) {
      estimate = estimator(Y.orig, N0, T0)
      estimate.shift = estimator(Y1, N0, T0)
      expect_equal(lambda(estimate), lambda(estimate.shift))
      expect_equal(omega(estimate),  omega(estimate.shift))
      expect_equal(c(estimate.shift), c(estimate) + c)
    }

    # 2.
    # Re-mapping Yit <- Yit + c for exposed t < T0 decreases tau by c (exception: "SC")
    Y2 <- setup$Y
    Y2[exposed, 1:T0] <- c + Y2[exposed, 1:T0]
    for (estimator in estimators[-1]) {
      estimate = estimator(Y.orig, N0, T0)
      estimate.shift = estimator(Y2, N0, T0)
      expect_equal(lambda(estimate.shift), lambda(estimate))
      expect_equal(c(estimate.shift), c(estimate) - c, tol = 1e-10)
    }

    # 3.
    # Re-mapping Yit <- Yit + c for unexposed t < T0 increases tau by c (exception: "SC")
    Y3 <- setup$Y
    Y3[!exposed, 1:T0] <- c + Y3[!exposed, 1:T0]
    for (estimator in estimators[-1]) {
      estimate = estimator(Y.orig, N0, T0)
      estimate.shift = estimator(Y3, N0, T0)
      expect_equal(lambda(estimate.shift), lambda(estimate))
      expect_equal(c(estimate.shift), c(estimate) + c, tol = 1e-10)
    }

    # 4.
    # Re-mapping Yit <- Yit + c for unexposed t > T0 decreases tau by c
    Y4 <- setup$Y
    Y4[!exposed, (T0+1):T] <- c + Y4[!exposed, (T0+1):T]
    for (estimator in estimators[-1]) {     # shouldn't work for sc -- relies on unit fixed effects
      estimate = estimator(Y.orig, N0, T0)
      estimate.shift = estimator(Y4, N0, T0)
      expect_equal(omega(estimate.shift), omega(estimate))
      expect_equal(c(estimate.shift), c(estimate) - c, tol = 1e-10)
    }
  }
})
