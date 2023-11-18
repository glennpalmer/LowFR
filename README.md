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
where `p` is the number of exposures, `TT` is the number of measurement times, and `k` is the number of latent factors for each measurement time. If `k` is left as `NULL`, it will be chosen automatically based on the singular value decomposition of the matrix of $x_{it}\text{s}$ as described in the paper.

_Note: the rows of `X` must be organized so that_ $x_i = (x_{i11}, ..., x_{i1T}, ...., x_{ip1}, ... , x_{ipT})$ _where_ $x_{ijt}$ _is the_ $j\text{th}$ _exposure at measurement time_ $t$ _for patient_ $i$.

In the posterior samples saved to `fit1` in the above, the induced coefficients for $\mathbb{E}[y_i | x_i]$ are saved as `alpha_0` (intercept), `alpha` (linear terms), and `Gamma` (the symmetric interaction matrix), matching the notation from Corollary 3.2 in the paper.

For a concrete example, the functions in `R/simdata.R` generate simulated data under the three scenarios described in Section 4 of the paper. For example, the following code generates data under scenario 2 and then fits a LowFR model with a small number of posterior draws.

```r
data <- simulate_scenario2(random_seed=1234)
fit <- fit_LowFR(y_obs=data$y_obs, X_obs=data$X_obs, p=10, TT=3, burnin=100, samples=100)
```

Note: if `burnin` and `samples` are not specified, they each default to 1000 (which we used for generating the results in the paper, and recommend for general use). There is also a `chains` option in `fit_LowFR` that defaults to 4. Similarly, `p`, `TT`, `n`, and a number of other data-generation settings can be modified when calling the `simulate_scenarioX` functions.
