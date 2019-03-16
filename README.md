
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="man/figures/flexcm_hex.png" align="right" height=160/ display="block">

# `flexcm`: Flexible Models for Dispersed Count Data

[![Travis build
status](https://travis-ci.org/jreduardo/flexcm.svg?branch=master)](https://travis-ci.org/jreduardo/flexcm)

> [Eduardo E. R. Junior](http://leg.ufpr.br/~eduardojr) -
> <jreduardo@usp.br>, IME-USP

The `flexcm` package contains functions to fit flexible count models
that can handle equi-, over-, and underdispersion, namely we consider
COM-Poisson, Gamma-count, discrete Weibull, generalized Poisson, double
Poisson and Poisson-Tweedie models\[1\]. The normalizing constant for
double Poisson and COM-Poisson are written in C++.

Joint work with [Walmes M. Zeviani](www.leg.ufpr.br/~walmes/) and
[Clarice G.B.
Demétrio](http://ce.esalq.usp.br/equipe/clarice-garcia-borges-demetrio).

## Installation

You can install the development version of `flexcm` from
[GitHub](https://github.com/jreduardo/flexcm) with:

``` r

# install.packages("devtools")
devtools::install_github("jreduardo/flexcm")
```

## Usage and example

Basically, this package implements methods similar to those related to
glm objects. The main function is `flexcm(..., model)`.

``` r

library(flexcm)

# Fit models -----------------------------------------------------------

# Model families
families <- list("CMP" = "compoisson",
                 "GCT" = "gammacount",
                 "DWe" = "discreteweibull",
                 "GPo" = "generalizedpoisson",
                 "DPo" = "doublepoisson",
                 "PTw" = "poissontweedie")

# Fit models
models <- lapply(families, function(fam) {
  flexcm(ninsect ~ extract, model = fam, data = sitophilus)
})

# Methods --------------------------------------------------------------

vapply(models, coef, numeric(5))
#>                        CMP          GCT         DWe          GPo          DPo          PTw
#> log(nu)       -0.927155720 -0.927261394  1.03484864  0.019426056  0.865851300  0.347643101
#> (Intercept)    3.449732023  3.425483284 -9.88053024  3.449990747  3.450317139  3.449987546
#> extractLeaf   -0.006370969 -0.006605471 -0.04866513 -0.006369236 -0.006363878 -0.006369448
#> extractBranch -0.052157044 -0.053534485  0.07170600 -0.052129200 -0.052087403 -0.052129065
#> extractSeed   -3.255210977 -4.015742799  7.60692064 -3.354677303 -3.564366151 -3.354677366

vapply(models, logLik, numeric(1))
#>       CMP       GCT       DWe       GPo       DPo       PTw 
#> -121.6334 -121.6509 -128.7893 -122.2840 -121.7930 -121.8466

equitest(models[-3])
#> 
#> Likelihood ratio test for equidispersion 
#> 
#>                         Resid.df   Loglik LRT_stat LRT_df Pr(>LRT_stat)    
#> COM-Poisson         CMP       35 -121.633   17.899      1     2.330e-05 ***
#> Gamma-count         GCT       35 -121.651   17.864      1     2.373e-05 ***
#> Generalized Poisson GPo       35 -122.284   16.598      1     4.621e-05 ***
#> Double Poisson      DPo       35 -121.793   17.579      1     2.755e-05 ***
#> Poisson-Tweedie     PTw       34 -121.847   17.472      2     0.0001607 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

print(models[["PTw"]])
#> 
#> Poisson-Tweedie regression models
#> Call:  flexcm(formula = ninsect ~ extract, data = sitophilus, model = fam)
#> 
#> Mean coefficients:
#>   (Intercept)    extractLeaf  extractBranch    extractSeed  
#>      3.449988      -0.006369      -0.052129      -3.354677  
#> 
#> Dispersion coefficient: omega = 0.3476
#> Power coefficient (estimated): power = 1.403
#> 
#> Residual degrees of freedom: 34
#> Minus twice the log-likelihood: 243.6932

# Predict new data -----------------------------------------------------

newdf <- sitophilus[c(21, 31), -2, drop = FALSE]
purrr::map_dfr(models,
               .id = "model",
               .f = predict,
               newdata = newdf,
               type = "response",
               interval = "confidence",
               augment_data = TRUE)
#> # A tibble: 12 x 5
#>    model extract    fit    lwr   upr
#>    <chr> <fct>    <dbl>  <dbl> <dbl>
#>  1 CMP   Seed     1.21   0.765  1.93
#>  2 CMP   Control 31.5   26.5   37.5 
#>  3 GCT   Seed     1.10   0.563  2.43
#>  4 GCT   Control 31.5   26.5   37.5 
#>  5 DWe   Seed     1.50   1.08   2.03
#>  6 DWe   Control 29.3   23.4   36.6 
#>  7 GPo   Seed     1.10   0.602  2.01
#>  8 GPo   Control 31.5   26.4   37.6 
#>  9 DPo   Seed     0.892  0.310  2.57
#> 10 DPo   Control 31.5   26.6   37.4 
#> 11 PTw   Seed     1.1    0.552  2.19
#> 12 PTw   Control 31.5   26.6   37.4
```

Currently, the methods implemented for `"flexcm"` objects are

``` r

methods(class = "flexcm")
#>  [1] anova        coef         equitest     fitted       logLik       model.matrix predict     
#>  [8] print        summary      vcov        
#> see '?methods' for accessing help and source code
```

## Related projects

There are other R packages to deal with COM-Poisson models that have
somehow contributed to the writing of `flexcm`.

  - [`mcglm`](https://github.com/wbonat/mcglm): Routines for fitting the
    multivariate covariance genelralized linear models. Currently, this
    package is used to fit the Poisson-Tweedie models in `flexcm`.
  - [`DWreg`](https://cran.r-project.org/package=DWreg): Fit Discrete
    weibull models (allows to model \(q\) or \(\rho\) parameters).
  - [`gamlss.dist`](https://cran.r-project.org/package=gamlss.dist):
    Implements (among other) double Poisson and a restrictive
    generalized Poisson models.
  - [`glmmTMB`](https://github.com/glmmTMB/glmmTMB): Fit (among other)
    COM-Poisson, under a different mean-parametrization, and generalized
    Poisson models (includes zero-inflation, dispersion modeling and
    random effects).

## License

The `flexcm` package is licensed under the [GNU General Public License,
version 3](https://www.gnu.org/licenses/gpl-3.0.ht), see file
`LICENSE.md`, © 2019 E. E., Ribeiro Jr.

<!------------------------------------------- -->

<!-- Links -->

1.  Poisson-Tweedie models are fitted using
    [`mcglm`](https://github.com/wbonat/mcglm) package
