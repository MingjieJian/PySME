# Fitresults

The fitresults objects is a pretty straightforward
collection of fitresults.

Note that all results are purely based on the fitting and
do not include model uncertainties from e.g. the linelist.

The given uncertainties in punc are estimated from the residual
distribution and seem to be reasonable for some parameters,
but get problematic if the parameter only affects a small number
of points (e.g. individual abundances).

Here are the fields

- `maxiter`: Maximum number of iterations
- `chisq`: Final chi square
- `uncertainties`: Uncertainties of the fitparameters bases on SME statistics
- `fit_uncertainties`: Uncertainties of the fitparameters based solely on the least-squares fit
- `covar`: covariance matrix of the fitparameters
- `grad`: the gradient of the cost function for each parameter
- `pder`: The jacobian, i.e. the derivate at each point for each parameter
- `resid`: The residual between observation and model

## Uncertainties

As mentioned above we estimate the uncertainties using a special metric.
This metric is based on the destribution of the derivatives for each fitparameter.
We estimate the cumulative distribution function of the **generalized** normal distribution with:

```py
  x = residual / derivate
  y = abs(derivate) / uncs
  y = cumulative_sum(y)
  #normalize
  y /= y[-1]
```
![](https://upload.wikimedia.org/wikipedia/commons/thumb/c/ca/Normal_Distribution_CDF.svg/500px-Normal_Distribution_CDF.svg.png).

(Example of the cumulative distribution function for various values of sigma)


It is explained in detail in [Paper II (SME: Evolution)](https://ui.adsabs.harvard.edu/abs/2017A%26A...597A..16P/abstract).
