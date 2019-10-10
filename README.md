# losmix: Inference for Gaussian Location-Scale Mixed-Effects Models

*Martin Lysy*

---

### Description

Tools for fitting Gaussian linear mixed-effects models where both the mean and variance terms are random.  Marginal likelihoods are compiled in C++ with automatic differentiation using the [**TMB**](https://github.com/kaskr/adcomp/wiki) library, such that efficient gradient-based optimization algorithms can be readily employed.  The package provides C++ entry points for users to extend the basic model with minimal effort.

### Installation

Install the R package [**devtools**](https://CRAN.R-project.org/package=devtools) and run
```r
devtools::install_github("mlysy/losmix")
```

### Usage

Please see package [vignette](http://htmlpreview.github.io/?https://github.com/mlysy/losmix/master/doc/losmix.html).
