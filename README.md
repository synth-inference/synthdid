# synthdid: Synthetic Difference in Differences Estimation

[![Build Status](https://dev.azure.com/grf-labs/synth-inference/_apis/build/status/synth-inference.synthdid?branchName=master)](https://dev.azure.com/grf-labs/synth-inference/_build/latest?definitionId=4&branchName=master)

This package implements the synthetic difference in difference estimator (SDID) for the
average treatment effect in panel data, as proposed in Arkhangelsky et al (2019).
We consider a setting in which we observe a matrix `Y = L + tau W + noise` where W
is a matrix of indicators for treatment.  All treated units must begin treatment simultaneously,
so W indicates a treated block, i.e. `W[i,j] = 1` for `i > N_0, j > T_0` and is zero otherwise.
This applies, in particular, to the case of a single treated unit.

This package is currently in beta and the functionality and interface is subject to change.

Some helpful links for getting started:

- The [R package documentation](https://synth-inference.github.io/synthdid/) contains usage examples and method reference.
- The [online vignettes](https://synth-inference.github.io/synthdid/articles/more-plotting.html) contains a gallery of plot examples.
- For community questions and answers around usage, see [Github issues page](https://github.com/synth-inference/synthdid/issues).

### Installation

The current development version can be installed from source using devtools.

```R
devtools::install_github("synth-inference/synthdid")
```

### Example

```R
library(synthdid)

# Estimate the effect of California Proposition 99 on cigarette consumption
data('california_prop99')
setup = panel.matrices(california_prop99)
tau.hat = synthdid_estimate(setup$Y, setup$N0, setup$T0)
se = sqrt(vcov(tau.hat, method='placebo'))
sprintf('point estimate: %1.2f', tau.hat)
sprintf('95%% CI (%1.2f, %1.2f)', tau.hat - 1.96 * se, tau.hat + 1.96 * se)
plot(tau.hat)
```

#### References
Dmitry Arkhangelsky, Susan Athey, David A. Hirshberg, Guido W. Imbens, and Stefan Wager.
<b>Synthetic Difference in Differences</b>, 2019.
[<a href="https://arxiv.org/abs/1812.09970">arxiv</a>]
