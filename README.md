[![Build Status](https://dev.azure.com/grf-labs/synth-inference/_apis/build/status/synth-inference.synthdid?branchName=master)](https://dev.azure.com/grf-labs/synth-inference/_build/latest?definitionId=4&branchName=master)

# synthdid: Synthetic Difference in Differences Estimation

This package implements the synthetic difference in difference estimator (SDID) for the
average treatment effect in panel data, as proposed in Arkhangelsky et al (2019).
We consider a setting in which we observe a matrix `Y = L + tau W + noise` where W
is a matrix of indicators for treatment.  All treated units must begin treatment simultaneously,
so W indicates a treated block, i.e. `W[i,j] = 1` for `i > N_0, j > T_0` and is zero otherwise.
This applies, in particular, to the case of a single treated unit.

This package is currently in beta and the functionality and interface is subject to change.

To install this package in R, run the following commands:
```R
library(devtools)
install_github("synth-inference/synthdid")
```
Example usage:

```R
library(synthdid)

setup = synthdid:::random.low.rank()
tau.hat = synthdid_estimate(Y,n_0,T_0)
se = synthdid_se(tau.hat)

print(paste("true tau:", tau))
print(paste0("point estimate: ", round(tau.hat, 2)))
print(paste0("95% CI for tau: (", round(tau.hat - 1.96 * se, 2), ", ", round(tau.hat + 1.96 * se, 2), ")"))
plot(tau.hat)
```

#### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
<b>Synthetic Difference in Differences</b>
2019.
[<a href="https://arxiv.org/abs/1812.09970">arxiv</a>]
