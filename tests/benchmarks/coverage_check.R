# This script checks the coverage properties of all
# estimators for the random.low.rank DGP.
# R: 3.6
set.seed(42)
library(synthdid)

estimators = list("sc_estimate", "did_estimate", "synthdid_estimate")
CI.methods = c("jackknife", "bootstrap", "placebo")
res = 0
for (estimator in estimators) {
  for (CI.method in CI.methods) {
    cov = replicate(1000, {
      setup = synthdid:::random.low.rank()
      estimate = do.call(estimator, setup[-2])
      se = synthdid_se(estimate, method = CI.method)
      0 > estimate - 1.96*se && 0 < estimate + 1.96*se
    })
    res = c(res, c(paste0(estimator, CI.method, ": ", mean(cov))))
  }
}
matrix(res[-1], length(res[-1]), 1)

# Output at https://github.com/synth-inference/synthdid/pull/42 is:
