
# flimo : Fixed Landscape Inference MethOd

[![CRAN
version](https://www.r-pkg.org/badges/version/flimo)](https://CRAN.R-project.org/package=flimo)

<!-- badges: start -->
<!-- badges: end -->

Fixed Landscape Inference MethOd allows likelihood free efficient
inference for stochastic models.

The framework is simply:

-   For each parameter value $\theta$, build a function that compare a
    certain number of (specific) simulations with data (e.g. for a
    scalar data, (mean(simulations($\theta$)-data))^2);

-   thanks to the specific management of the randomness, it’s possible
    to use determinsitic optimization algorithm to minimize that
    objective function.

This R-package works on its own and also with the Julia language for
more efficiency.

More details about the inference method and the package can be found in
the following paper : <https://doi.org/10.48550/arXiv.2210.06520>.
Please cite it if you use *flimo*.

## Requirements

The package uses the language Julia for some features. The link between
R and Julia is done with the R-package \[JuliaConnectoR\]
<https://github.com/stefan-m-lenz/JuliaConnectoR>.

Follow their recommendation:

“The package requires that [Julia (version ≥ 1.0) is
installed](https://julialang.org/downloads/) and that the Julia
executable is in the system search `PATH` or that the `JULIA_BINDIR`
environment variable is set to the `bin` directory of the Julia
installation.”

## Installation

You can install the released version of flimo from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("flimo")
```

## Overview

The *flimo* algorithm allows to infer parameters of continuous
stochastic models. It is based on simple simulations of the process,
with a specific randomness management.

The parameters of the model are grouped in a $\theta$ vector.

The user needs to define one number and to build two functions:

-   The number of random draws to be made for a single simulation,
    called $n_{draw}$. For example, a model with 10 runs, each based on
    5 binomial draws, will have $n_{draw} = 50$. An upper bound on this
    number is also appropriate.

-   a special simulator of the form *simulatorQ($\theta$, quantiles)*.
    The way to build it from a usual random simulator is detailed in
    section **Building simulatorQ with R** below;

-   a function of the form *dsumstats(data, simulations)* measuring the
    difference in summary statistics between the observed data and the
    simulations w.r.t. the tested parameters.

These two functions can also be combined directly by the user into an
objective function of the form: *objective($\theta$, quantiles)* to
minimize.

The use of *flimo* is then easy:

``` r
#flimoptim(ndraw, data, dsumstats, simulatorQ)
```

or

``` r
#flimoptim(ndraw, obj = objective)
```

The various options and features are documented in the package manual.
In particular, it is almost always necessary to adjust *nsim* (which
affects the speed and accuracy of the calculation), *lower* and *upper*
(the bounds of the parameter space).

## Example

Vignette provides two basic examples to learn how tu use *flimo*.

## Modes in the R-interface

The first mode (“R”) is basic R: the optimization functions are the ones
implemented in *optim* (package *stats*). Two methods are available:
“L-BFGS-B” and “Brent”.

The second mode (“Julia”) uses Julia functions wrapped in package
Jflimo.jl (see next section). It is designed to be up to 100 times
faster than R for non-trivial models. The interface is in R but the user
has to write two functions in Julia language: *simulatorQ* and
*dsumstats*. Both of these names are mandatory for the Julia functions.

See Vignette to see how to build them.

## Using Julia

Julia is relatively easy to learn. You can find a tutorial here:

<https://www.freecodecamp.org/news/learn-julia-programming-language/>

The Julia Jflimo package is available online
<https://metabarcoding.org/flimo/jflimo>. To install it and other
necessary packages once Julia is installed, the user can use the
following function (only the first time):

``` r
#julia_setup()
```

To write efficient code, you should follow these \[Performances
Tips\]<https://docs.julialang.org/en/v1/manual/performance-tips/>.
Package TimerOutputs is also a good tool.

The two examples of the vignette using Julia are presented at the end of
this page.

## Building simulatorQ with R

In concrete terms, each random draw must be replaced by a call to the
quantile function of the same distribution, taking into account the fact
that several simulations can be done at the same time.

The quantile functions are then applied to a matrix *quantiles* that the
user does not have to define himself, which is obtained in the package
*flimo* by:

``` r
ndraw <- 5 #random draws for one simulation, e.g. 5 cycles here
nsim <- 10 #number of simulations to average
quantiles <- matrix(runif(ndraw*nsim), nrow = nsim)
```

In the code of usual simulators, each line

*rdistrib(n,parameters)*

should be replaced by

*qdistrib(quantiles\[,i:(i+n-1)\], parameters)*

Example:

``` r
rnorm(5, mean = 0, sd = 1)
#> [1]  0.8672414 -1.5784402 -0.2918760  0.9841565 -1.0402933
```

becomes

``` r
qnorm(quantiles[,1:5], mean = 0, sd = 1)
#>              [,1]       [,2]        [,3]       [,4]        [,5]
#>  [1,] -0.34782224 -0.3967510 -1.30725379  0.5689414  0.16688808
#>  [2,]  0.12917020  0.2906169 -1.50181525  0.6019140  0.23814810
#>  [3,] -0.08018577  0.5334808  0.70101730 -0.3246339 -1.31112015
#>  [4,]  0.99915793 -0.4276736  0.80845545 -0.7295606 -0.83459642
#>  [5,]  2.18809576 -1.9612868  0.90043668  0.4904709 -0.15344521
#>  [6,]  0.30922129  1.0470538 -0.13561444  0.7846481 -0.56685544
#>  [7,]  1.52570792 -1.5697612 -0.07396329  0.7584592  0.41790599
#>  [8,]  0.08070942 -1.8478335  1.36840144  0.1754735  0.07035972
#>  [9,]  0.12500653 -0.7833601 -0.82447021 -0.3877886  0.14260788
#> [10,]  1.22803757  0.1919990 -1.08012009 -1.8801307  2.56684815
```

Each row is an independent simulation.

## Handling discrete models

The usual gradient-based optimization algorithm can’t covnerge in that
framework because the quantile functions of discrete distributions are
constant by pieces.

-   Your model has only one parameter: the Brent method is implemented
    for *flimo*;

-   replace every discrete distribution by a continuous one, e.g. :
    $Poisson(\theta) \leftarrow Normal(\mu = \theta, \sigma^2 = \theta)$;

-   (in project for multidimensional problems: adapted Nelder-Mead
    method.)

## Applications

Two applications are available in the \[applications
directory\]<https://metabarcoding.org/flimo/flimo/applications>:

-   Inferring parameters of a g-and-k distribution
-   Inferring a selective value for a Wright-Fisher model

## Examples of the Vignette with Julia

``` r
library(flimo)
library(JuliaConnectoR)
```

### Example 1: Poisson distribution

``` r

#Setup
set.seed(1)

#Create data

Theta_true1 <- 100 #data parameter
n1 <- 5 #data size

simulator1 <- function(Theta, n){
  #classical random simulator
  rpois(n, lambda = Theta)
}

data1 <- simulator1(Theta_true1, n1)

#Simulations with quantiles

ndraw1 <- n1 #number of random draws for one simulation

simulatorQ1 <- function(Theta, quantiles){
  qpois(quantiles, lambda = Theta)
}

#With Normal approximation
simulatorQ1N <- function(Theta, quantiles){
  qnorm(quantiles, mean = Theta, sd = sqrt(Theta))
}

nsim1 <- 10

#Sample Comparison: summary statistics

dsumstats1 <- function(simulations, data){
  #simulations : 2D array
  #data : 1D array
  mean_simu <- mean(rowMeans(simulations))
  mean_data <- mean(data)
  (mean_simu-mean_data)^2
}
```

There are optimization issues in R for normalized process and $\theta$
close to 0. Lower bound set to 1.

``` r
#Optimization with normal approximation of Poisson distribution

#Default mode: full R

system.time(optim1R <- flimoptim(ndraw1, data1, dsumstats1, simulatorQ1N,
                 ninfer = 10,
                 nsim = nsim1,
                 lower = 1, upper = 1000,
                 randomTheta0 = TRUE))
#> utilisateur     système      écoulé 
#>       0.044       0.003       0.052

summary(optim1R)
#> $Mode
#> [1] "R"
#> 
#> $method
#> [1] "L-BFGS-B"
#> 
#> $number_inferences
#> [1] 10
#> 
#> $number_converged
#> [1] 10
#> 
#> $minimizer
#>       par1       
#>  Min.   : 99.00  
#>  1st Qu.: 99.81  
#>  Median :100.71  
#>  Mean   :101.05  
#>  3rd Qu.:101.88  
#>  Max.   :105.06  
#> 
#> $minimum
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 1.140e-24 1.110e-22 2.989e-22 8.113e-21 6.841e-22 7.762e-20 
#> 
#> $median_time_inference
#> [1] 0.00235343

#Version with objective function provided

obj1 <- function(Theta, quantiles, data = data1){
  simulations <- simulatorQ1N(Theta, quantiles)
  dsumstats1(simulations, data)
}

system.time(optim1Rbis <- flimoptim(ndraw1,
                 obj = obj1,
                 ninfer = 10,
                 nsim = nsim1,
                 lower = 1, upper = 1000,
                 randomTheta0 = TRUE))
#> utilisateur     système      écoulé 
#>       0.008       0.000       0.009

#Second mode : full Julia

#Optimization with Automatic Differentiation
#Warning : you need to translate dsumstats and simulatorQ to Julia
#Both of these names for the Julia functions are mandatory !

#Most accurate config for complex problems : IPNewton with AD

if (juliaSetupOk()){
  julia_simulatorQ1N <-"
function simulatorQ(Theta, quantiles)
  quantile.(Normal.(Theta, sqrt.(Theta)), quantiles)
end
"
  
  julia_dsumstats1 <-"
function dsumstats(simulations, data)
  (mean(mean(simulations, dims = 2))-mean(data))^2
end
"
  
  system.time(optim1JAD <- flimoptim(ndraw1, data1, julia_dsumstats1, julia_simulatorQ1N,
                                     ninfer = 10,
                                     nsim = nsim1,
                                     lower = 0,
                                     upper = 1000,
                                     randomTheta0 = TRUE,
                                     mode = "Julia",
                                     load_julia = TRUE))

  summary(optim1JAD)
  
  #IPNewton without AD
  system.time(optim1J <- flimoptim(ndraw1, data1, julia_dsumstats1, julia_simulatorQ1N,
                                   ninfer = 10, nsim = nsim1, lower = 0, upper = 1000,
                                   randomTheta0 = TRUE, mode = "Julia", AD = FALSE))
  optim1J
  summary(optim1J)
  
  #Brent
  system.time(
    optim1JBrent <- flimoptim(ndraw1,
                              data1,
                              julia_dsumstats1,
                              julia_simulatorQ1N,
                              ninfer = 10,
                              nsim = nsim1,
                              lower = 0,
                              upper = 1000,
                              randomTheta0 = TRUE,
                              mode = "Julia",
                              method = "Brent")
  )

  summary(optim1JBrent)
}
#> Starting Julia ...
#> $Mode
#> [1] "Julia"
#> 
#> $method
#> [1] "Brent()"
#> 
#> $number_inferences
#> [1] 10
#> 
#> $number_converged
#> [1] 10
#> 
#> $minimizer
#>       par1      
#>  Min.   : 99.4  
#>  1st Qu.:100.5  
#>  Median :101.5  
#>  Mean   :101.3  
#>  3rd Qu.:102.1  
#>  Max.   :102.8  
#> 
#> $minimum
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 3.169e-16 1.429e-14 1.211e-13 2.065e-13 4.706e-13
```

### Example 2: Normal distribution

``` r
#Setup
set.seed(1)

#Create data

Theta_true2 <- c(3, 2) #data parameter
n2 <- 5 #data size

simulator2 <- function(Theta, n){
  #classical random simulator
  rnorm(n, mean = Theta[1], sd = Theta[2])
}

data2 <- simulator2(Theta_true2, n2)

#Simulations with quantiles

simulatorQ2 <- function(Theta, quantiles){
  qnorm(quantiles, mean = Theta[1], sd = Theta[2])
}

ndraw2 <- 5

dsumstats2 <-function(simulations, data){
  mean_simu <- mean(rowMeans(simulations))
  mean_data <- mean(data)
  sd_simu <-mean(apply(simulations, 1, sd))
  sd_data <- sd(data)
  (mean_simu-mean_data)^2+(sd_simu-sd_data)^2
}

nsim2 <- 10
```

``` r
#Optimization

#Default mode: full R

optim2R <- flimoptim(ndraw2, data2, dsumstats2, simulatorQ2,
                 ninfer = 10, nsim = nsim2,
                 lower = c(-5, 0), upper = c(10, 10),
                 randomTheta0 = TRUE)

summary(optim2R)
#> $Mode
#> [1] "R"
#> 
#> $method
#> [1] "L-BFGS-B"
#> 
#> $number_inferences
#> [1] 10
#> 
#> $number_converged
#> [1] 10
#> 
#> $minimizer
#>       par1            par2      
#>  Min.   :2.896   Min.   :1.824  
#>  1st Qu.:2.955   1st Qu.:1.979  
#>  Median :3.201   Median :2.168  
#>  Mean   :3.247   Mean   :2.128  
#>  3rd Qu.:3.422   3rd Qu.:2.229  
#>  Max.   :4.024   Max.   :2.393  
#> 
#> $minimum
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.000e+00 0.000e+00 1.400e-20 1.395e-15 4.208e-18 1.394e-14 
#> 
#> $median_time_inference
#> [1] 0.01571453

#Second mode : full Julia
#Optimization with Automatic Differentiation
#Warning : you need to translate dsumstats and simulatorQ to Julia
#Both of these names for the Julia functions are mandatory !

if (juliaSetupOk()){
  julia_simulatorQ2 <-"
function simulatorQ(Theta, quantiles)
  quantile.(Normal.(Theta[1], Theta[2]), quantiles)
end
"
  
  julia_dsumstats2 <-"
 function dsumstats(simulations, data)
  (mean(mean(simulations, dims = 2))-mean(data))^2+
  (mean(std(simulations, dims = 2))-std(data))^2
end
"
  
  #Most accurate config for complex problems : IPNewton with AD
  optim2JIPAD <- flimoptim(ndraw2, data2, julia_dsumstats2, julia_simulatorQ2,
                           ninfer = 10, nsim = nsim2,
                           lower = c(-5, 0), upper = c(10, 10),
                           randomTheta0 = TRUE,
                           mode = "Julia")

  summary(optim2JIPAD)
  
  #IPNewton without AD
  optim2JIP <- flimoptim(ndraw2, data2, julia_dsumstats2, julia_simulatorQ2,
                         ninfer = 10, nsim = nsim2,
                         lower = c(-5, 0), upper = c(10, 10),
                         randomTheta0 = TRUE,
                         mode = "Julia")

  summary(optim2JIP)
  
  #Nelder-Mead
  #No bounds allowed inside NM method: objective function is overridden
  
  julia_simulatorQ2NM <-"
function simulatorQ(Theta, quantiles)
  if Theta[2] < 0
    fill(Inf, size(quantiles))
  else
    quantile.(Normal.(Theta[1], Theta[2]), quantiles)
  end
end
"
  
  optim2JNM <- flimoptim(ndraw2,
                         data2, julia_dsumstats2, julia_simulatorQ2NM,
                         ninfer = 10,
                         nsim = nsim2,
                         lower = c(-5, 0),
                         upper = c(10, 10),
                         randomTheta0 = TRUE,
                         mode = "Julia",
                         method = "NelderMead")

  summary(optim2JNM)
}
#> $Mode
#> [1] "Julia"
#> 
#> $method
#> [1] "NelderMead{Optim.AffineSimplexer,Optim.AdaptiveParameters}"
#> 
#> $number_inferences
#> [1] 10
#> 
#> $number_converged
#> [1] 10
#> 
#> $minimizer
#>       par1            par2      
#>  Min.   :3.118   Min.   :1.917  
#>  1st Qu.:3.243   1st Qu.:2.124  
#>  Median :3.332   Median :2.174  
#>  Mean   :3.381   Mean   :2.189  
#>  3rd Qu.:3.547   3rd Qu.:2.259  
#>  Max.   :3.640   Max.   :2.580  
#> 
#> $minimum
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 2.032e-10 1.391e-09 3.262e-09 2.677e-09 3.843e-09 4.524e-09
```
