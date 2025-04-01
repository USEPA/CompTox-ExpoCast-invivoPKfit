# R Package "invivoPKfit"

invivoPKfit estimates the most likely parameter values for models describing the 
pharmacokinetics (PK, that is the absorption, distribution, metabolism, and 
elimination) of a compound given a set of in vivo concentration vs. time 
observations Cobs

This is achieved by creating a function that calculates probability of Cobs 
given a PK model M(x,t) and a statistical model S(σ)

This probability is described by a likelihood function L(Cobs, M, S)

We then use an optimizer function to optimize the values of x and σ such that 
L(Cobs, M(x), S(σ))

We then use the method of quadrature to estimate uncertainty in the estimates 
of x and σ as well as any functions of those parameters (for example, Cmax and 
AUC/area under the C(t) curve) 


## Background


## Censored Data

If the observation is that the concentration was below a limit of quantitation (LOQ) we call the observation “censored”
In this case we add to the likelihood all the probability from zero to the limit of quantitation
We use the cumulative distribution function of the log-normal distribution for this:

CDF(LOQ, ������, ������) =1/2+1/2 ������������������((ln⁡������������������−������)/(√2 ������))

We separate the observations into those above the LOQ and below the LOQ, above the LOQ we use the log-normal density, below the LOQ we use the CDF

## References

Anderson, D., and K. Burnham. "Model selection and multi-model inference." Second. NY: Springer-Verlag 63.2020 (2004): 10.

Akaike, H. "Maximum likelihood identification of Gaussian autoregressive moving average models." Biometrika 60.2 (1973): 255-265.

Boxenbaum, Harold G., Sidney Riegelman, and Robert M. Elashoff. "Statistical estimations in pharmacokinetics." Journal of pharmacokinetics and biopharmaceutics 2.2 (1974): 123-148.

Charles, Sandrine, Aude Ratier, and Christelle Lopes. "Generic solving of one-compartment toxicokinetic models." bioRxiv (2021).

Colburn, Wayne A. "Controversy III: To model or not to model." The Journal of Clinical Pharmacology 28.10 (1988): 879-888.

Dovì, V. G., O. Paladino, and A. P. Reverberi. "Some remarks on the use of the inverse hessian matrix of the likelihood function in the estimation of statistical properties of parameters." Applied Mathematics Letters 4.1 (1991): 87-90.

Sayre, Risa R., John F. Wambaugh, and Christopher M. Grulke. "Database of pharmacokinetic time-series data and parameters for 144 environmental chemicals." Scientific data 7.1 (2020): 1-10.

Sheiner, Lewis B. "Analysis of pharmacokinetic data using parametric models—1: Regression models." Journal of pharmacokinetics and biopharmaceutics 12.1 (1984): 93-117.

Sheiner, Lewis B. "Analysis of pharmacokinetic data using parametric models. II. Point estimates of an individual's parameters." Journal of pharmacokinetics and biopharmaceutics 13.5 (1985): 515-540.

Sheiner, Lewis B. "Analysis of pharmacokinetic data using parametric models. III. Hypothesis tests and confidence intervals." Journal of pharmacokinetics and biopharmaceutics 14.5 (1986): 539-555.

Wambaugh, John F., et al. "Evaluating in vitro-in vivo extrapolation of toxicokinetics." Toxicological Sciences 163.1 (2018): 152-169.

Yamaoka, Kiyoshi, Terumichi Nakagawa, and Toyozo Uno. "Application of Akaike's information criterion (AIC) in the evaluation of linear pharmacokinetic equations." Journal of pharmacokinetics and biopharmaceutics 6.2 (1978): 165-175.


## Getting Started

### Dependencies

* Users will need the freely available R statistical computing language: <https://www.r-project.org/>
* Users will likely want a development environment like RStudio: <https://posit.co/download/rstudio-desktop/>

### Installing

* Getting Started with R Package bayesmarker from the R command line
```
library(devtools)
install_github("USEPA/CompTox-ExpoCast-invivoPKfit")
```
* RStudio provides a menu ‘Install Packages’ under ‘Tools’ tab
* Load the invivoPKfit data and functions
```
library(invivoPKfit)
```
* Check what version you are using 
```
packageVersion(invivoPKfit)
```

## Authors

John Wambaugh
[@wambaugh.john@epa.gov]

Caroline Ring
[@ring.caroline@epa.gov]

Christopher Cook
[@Cook.Christopher@epa.gov]

Gilberto Padilla Mercado
[@padillamercado.gilberto@epa.gov]

## License

License: GPL-3 <https://www.gnu.org/licenses/gpl-3.0.en.html>
