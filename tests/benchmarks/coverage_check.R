# This script checks the coverage properties of all
# estimators for the random.low.rank DGP.
# R: 3.6.0
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
      1 > estimate - 1.96*se && 1 < estimate + 1.96*se
    })
    res = c(res, c(paste0(estimator, CI.method, ": ", mean(cov))))
  }
}
matrix(res[-1], length(res[-1]), 1)

# Output at https://github.com/synth-inference/synthdid/pull/42 is:
# [,1]
# [1,] "sc_estimatejackknife: 1"
# [2,] "sc_estimatebootstrap: 0.954"
# [3,] "sc_estimateplacebo: 0.939"
# [4,] "did_estimatejackknife: 0.933"
# [5,] "did_estimatebootstrap: 0.92"
# [6,] "did_estimateplacebo: 0.951"
# [7,] "synthdid_estimatejackknife: 0.932"
# [8,] "synthdid_estimatebootstrap: 0.938"
# [9,] "synthdid_estimateplacebo: 0.937"
