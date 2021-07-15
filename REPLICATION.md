# Description
This file provides information sufficient to replicate the tables and figures in the paper [Synthetic Difference in Differences](https://arxiv.org/abs/1812.09970). In short, if you have installed the version of the *synthdid* R package used at that time,
for example via the commands
```
install.packages('devtools')
devtools::install_github('synth-inference/synthdid', ref='sdid-paper')
```
then running the RMarkdown notebook [paper-results](vignettes/paper-results.Rmd) will do it, installing all other packages necessary
and writing figures and tables to the directory [vignettes/figures](vignettes/figures). We include [saved simulation results](vignettes/all-simulations.rds) 
in the package --- delete this file before running the notebook to rerun them. 
This will take a few days to run on an 8-core machine; we include [instructions](#using-a-cluster) for the use of a cluster with the slurm scheduler. 

With the exception of the package *MCPanel* (and its dependency *glmnet*), which implements the matrix completion estimator
of [Athey et al.](https://arxiv.org/abs/1710.10251), we do not fix the versions of other packages to the versions we used, as we anticipate that updates to them will be
backward compatible. Because *MCPanel* uses a platform-dependent random number generator in its cross-validation scheme, results involving the matrix completion estimator do vary slightly from platform to platform. The package versions and platforms that were used are described [here](#package-versions).

# Data and Code Availability

All code used is included in the *synthdid* package,
available on github as described above. All datasets used are provided in the package.
They can be loaded as follows.

```
data('california_prop99')
data('CPS')
data('PENN')
```

This data was sourced as follows.
1. The California cigarette tax (Proposition 99) dataset is available as a part of
[*Synth* matlab toolbox](https://web.stanford.edu/~jhain/Synth_Matlab/Synth_MATLAB.zip), in the file 'MLAB_data.txt'.
2. The CPS data is available on the [NBER website](https://data.nber.org/morg/annual). Each year's data is in a file 'morgYY.dta',
where YY is the year written as 2 digits. We use a subset of this data; see the appendix of [Synthetic Difference in Differences](https://arxiv.org/abs/1812.09970) for details.
Along with it, we use three state level policy indicators,
which are included with the CPS data in our package.
  - [Minimum wage](http://www.dol.gov/whd/state/stateMinWageHis.htm)
  - [Open-carry gun law](https://lawcenter.giffords.org/gun-laws/policy-areas/guns-in-public/open-carry)
  - [Partial birth abortion bans](https://www.guttmacher.org/state-policy/explore/overview-abortion-laws)
3. The Penn World Table data is available for download [here](www.ggdc.net/pwt). We use a subset of this data;
see the appendix of [Synthetic Difference in Differences](https://arxiv.org/abs/1812.09970) for details.

# Computational Requirements

The *synthdid* R package is implemented in pure R with few dependencies. It should run on R versions 3 and 4 on most platforms.
We've used the following.

1. A 2020 Macbook Air (M1) running R 4.0.5.
2. A cluster of x86_64 machines running R 4.0.2.

The first was used for the California Proposition 99 example and for aggregation and plotting of the results from the simulations;
the second was used to run the simulations themselves.

## Using a cluster

To run code on our cluster, which uses the slurm scheduler, we did the following.
1. The [notebook](vignettes/paper-results.Rmd) was converted to a script by calling
`knitr::purl('paper-results.Rmd')`.
2. This script, and our [slurm template file](vignettes/paper-results-details/batchtools.slurm.tmpl), were put in the same folder on the cluster and the script was run.
3. The simulation results, saved in the [simulations directory](vignettes/simulations) on the cluster,
were copied to the Macbook Air and the notebook was run there to generate plots and tables.


Simply running the notebook on a single computer, e.g. the Macbook Air, should yield essentially the same results.


While we did not run all the simulations on the Macbook Air, we ran some and checked them against what we got on the cluster using using [this script](vignettes/paper-results-details/test-cross-platform.R.).
Results were essentially identical, with differences on the order of 1e-13 or smaller, except for the MC estimator: the *MCPanel* library uses the C++ std::default_random_engine RNG to choose folds for cross-validation, which does not guarantee identical results on different platforms.
Among the 5500 simulations on which we compared results for the MC estimator on our two platforms, the absolute difference of point estimates had the following quantiles.

```
         0%          10%          20%          30%          40%          50%          60%          70%          80%          90%         100%


0.000000e+00 3.505241e-07 6.974247e-06 1.939810e-05 4.741519e-05 9.357476e-05 1.737794e-04 3.053189e-04 6.039954e-04 1.520949e-03 5.482976e-02
```


## Package versions
We give the output of sessionInfo() on these two platforms after running the code in paper-results.Rmd.

### Macbook Air

> sessionInfo()


R version 4.0.5 (2021-03-31)


Platform: aarch64-apple-darwin20.3.0 (64-bit)


Running under: macOS Big Sur 11.2.3







Matrix products: default


BLAS:   /opt/homebrew/Cellar/openblas/0.3.13/lib/libopenblasp-r0.3.13.dylib


LAPACK: /opt/homebrew/Cellar/r/4.0.5/lib/R/lib/libRlapack.dylib







Random number generation:


RNG:     L'Ecuyer-CMRG


Normal:  Inversion


Sample:  Rejection






locale:


[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8







attached base packages:


[1] stats     graphics  grDevices utils     datasets  methods   base    







other attached packages:


[1] future.batchtools_0.10.0 doFuture_0.12.0          foreach_1.5.1          


[4] future_1.21.0            rngtools_1.5             MCPanel_0.0            


[7] synthdid_0.0.9         







loaded via a namespace (and not attached):


[1] progress_1.2.2    shape_1.4.5       tidyselect_1.1.1  purrr_0.3.4     


[5] listenv_0.8.0     splines_4.0.5     lattice_0.20-41   latex2exp_0.5.0


[9] colorspace_2.0-0  vctrs_0.3.8       generics_0.1.0    utf8_1.2.1      


[13] survival_3.2-10   rlang_0.4.10      pillar_1.6.0      withr_2.4.2     


[17] glue_1.4.2        rappdirs_0.3.3    lifecycle_1.0.0   stringr_1.4.0   


[21] munsell_0.5.0     gtable_0.3.0      mvtnorm_1.1-2     codetools_0.2-18


[25] batchtools_0.9.15 parallel_4.0.5    fansi_0.4.2       Rcpp_1.0.6      


[29] scales_1.1.1      backports_1.2.1   checkmate_2.0.0   parallelly_1.25.0


[33] brew_1.0-6        ggplot2_3.3.3     hms_1.1.0         digest_0.6.27   


[37] stringi_1.5.3     dplyr_1.0.6       grid_4.0.5        tools_4.0.5     


[41] magrittr_2.0.1    base64url_1.4     glmnet_4.1-1      tibble_3.1.1    


[45] crayon_1.4.1      pkgconfig_2.0.3   ellipsis_0.3.2    Matrix_1.3-2    


[49] data.table_1.14.0 prettyunits_1.1.1 iterators_1.0.13  R6_2.5.0        


[53] globals_0.14.0    compiler_4.0.5  





### Cluster

sessionInfo()


R version 4.0.2 (2020-06-22)


Platform: x86_64-pc-linux-gnu (64-bit)


Running under: CentOS Linux 7 (Core)







Matrix products: default


BLAS/LAPACK: /share/software/user/open/openblas/0.2.19/lib/libopenblasp-r0.2.19.so







Random number generation:


RNG:     L'Ecuyer-CMRG


Normal:  Inversion


Sample:  Rejection






locale:


[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C             


[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8   


[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8  


[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                


[9] LC_ADDRESS=C               LC_TELEPHONE=C           


[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C      







attached base packages:


[1] stats     graphics  grDevices utils     datasets  methods   base    







other attached packages:


[1] future.batchtools_0.10.0 doFuture_0.12.0          foreach_1.5.1          


[4] future_1.21.0            rngtools_1.5             MCPanel_0.0            


[7] synthdid_0.0.9         







loaded via a namespace (and not attached):


[1] latex2exp_0.5.0   Rcpp_1.0.6        compiler_4.0.2    pillar_1.6.1    


[5] prettyunits_1.1.1 progress_1.2.2    iterators_1.0.13  tools_4.0.2     


[9] digest_0.6.27     checkmate_2.0.0   lifecycle_1.0.0   tibble_3.1.2    


[13] gtable_0.3.0      lattice_0.20-41   pkgconfig_2.0.3   rlang_0.4.11    


[17] Matrix_1.2-18     parallel_4.0.2    mvtnorm_1.1-2     withr_2.4.2     


[21] stringr_1.4.0     rappdirs_0.3.3    hms_1.1.0         globals_0.14.0  


[25] vctrs_0.3.8       glmnet_4.1-1      grid_4.0.2        data.table_1.14.0


[29] glue_1.4.2        listenv_0.8.0     R6_2.5.0          parallelly_1.25.0


[33] fansi_0.5.0       base64url_1.4     survival_3.1-12   ggplot2_3.3.3   


[37] magrittr_2.0.1    backports_1.2.1   batchtools_0.9.15 scales_1.1.1    


[41] codetools_0.2-16  ellipsis_0.3.2    splines_4.0.2     shape_1.4.6     


[45] colorspace_2.0-1  brew_1.0-6        utf8_1.2.1        stringi_1.6.2   


[49] munsell_0.5.0     crayon_1.4.1
