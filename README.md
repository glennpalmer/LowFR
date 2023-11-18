# LowFR

This repository contains code that generates the results from the paper: G. Palmer, A.H. Herring, D.B. Dunson (2023+) 'Low-rank longitudinal factor regression', as well as an implementation of LowFR for use on other data sets.

If you have any questions, find bugs, etc. please reach out to glenn.palmer@duke.edu.

# Fitting a LowFR model

Once you've cloned the repo, run the following code to import libraries and source the relevant R files.

```r
library(rstan)
library(miceadds)
source.all(path="R")
```
Then, given a vector of outcomes $y \in \mathbb{R}$ and a matrix of longitudinal exposures $X \in \mathbb{R}^{n \times pT}$, you can fit LowFR by running

```r
fit1 <- fit_LowFR(y_obs=y, X_obs=X, p, TT, k=NULL)
```
where `p` is the number of exposures, `TT` is the number of measurement times, and `k` is the number of latent factors for each measurement time. If `k` is left as `NULL`, it will be chosen automatically based on the singular value decomposition of the matrix of $x_{it}$s as described in the paper.

_Note: the rows of `X` must be organized so that_ $x_i = (x_{i11}, ..., x_{i1T}, ...., x_{ip1}, ... , x_{ipT})$ _where_ $x_{ijt}$ _is the_ $j\text{th}$ _exposure at measurement time_ $t$ _for patient_ $i$.
