## Pump Model Sensitivity Example

This repo performs Markov Chain Monte Carlo (MCMC) and Variational Bayes
(VB) versions of the sensitivity analyses of [1] and [2].  Stan [3] is used
for the MCMC estimates as well as the derivative estimates needed for
VB estimation and sensitivity analysis.

You will need the R library `LRVBUtils`.  To install it, start R, and run

```R
> library(devtools)
> install_github("rgiordan/LRVBUtils")
```

Once installed, you should be able to run the scripts
* `fit_mcmc.R`
* `fit_vb.R`

...in that order.

### References

[1] Gustafson, Paul. "Local sensitivity of inferences to prior marginals." Journal of the American Statistical Association 91.434 (1996): 774-781.
APA

[2] Pérez, C. J., J. Martín, and M. J. Rufo. "Sensitivity estimations for Bayesian inference models solved by MCMC methods." Reliability Engineering & System Safety 91.10 (2006): 1310-1314.

[3] http://mc-stan.org/
