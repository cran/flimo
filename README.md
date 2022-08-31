
# flimo : Fixed Landscape Inference MethOd

[![CRAN
version](https://www.r-pkg.org/badges/version/flimo)](https://CRAN.R-project.org/package=flimo)

<!-- badges: start -->
<!-- badges: end -->

Fixed Landscape Inference MethOd allows likelihood free efficient
inference for stochastic models.

The framework is simply:

-   For each parameter value
    ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
    build a function that compare a certain number of (specific)
    simulations with data (e.g. for a scalar data,
    (mean(simulations(![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"))-data))^2);

-   thanks to the specific management of the randomness, it’s possible
    to use determinsitic optimization algorithm to minimize that
    objective function.

This R-package works on its own and also with the Julia language for
more efficiency.

More details about the inference method and the package can be found in
the following paper. Please cite it if you use *flimo*.

> Coming soon ! <https://doi.org/>…

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

The parameters of the model are grouped in a
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
vector.

The user needs to define one number and to build two functions:

-   The number of random draws to be made for a single simulation,
    called
    ![n\_{draw}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_%7Bdraw%7D "n_{draw}").
    For example, a model with 10 runs, each based on 5 binomial draws,
    will have
    ![n\_{draw} = 50](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n_%7Bdraw%7D%20%3D%2050 "n_{draw} = 50").
    An upper bound on this number is also appropriate.

-   a special simulator of the form
    *simulatorQ(![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
    quantiles)*. The way to build it from a usual random simulator is
    detailed in section **Building simulatorQ with R** below;

-   a function of the form *dsumstats(data, simulations)* measuring the
    difference in summary statistics between the observed data and the
    simulations w.r.t. the tested parameters.

These two functions can also be combined directly by the user into an
objective function of the form:
*objective(![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta"),
quantiles)* to minimize.

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
#> [1]  0.9580883 -0.5480208 -1.4569140 -0.2402742  1.8885880
```

becomes

``` r
qnorm(quantiles[,1:5], mean = 0, sd = 1)
#>             [,1]       [,2]        [,3]         [,4]       [,5]
#>  [1,] -0.1445568 -0.2805094 -0.80588279  2.039089952 -1.5477555
#>  [2,] -1.5422178 -1.7611784  0.05801355  0.822973238 -0.4517209
#>  [3,]  0.4708859 -1.9003471 -0.56537835  0.872456705  1.0743202
#>  [4,]  1.1132383  0.2197335  0.99611020  1.074424352 -0.5456380
#>  [5,] -1.0240074 -0.1747871 -1.35122350 -0.914113241  1.1102321
#>  [6,] -1.1946211  1.6320083  0.78535306  0.624859421  0.2370332
#>  [7,] -1.1195320 -2.7272738  0.14852302 -0.944825575  0.5690212
#>  [8,] -0.4388909  0.1362090 -0.26696724 -0.156150672 -1.5090528
#>  [9,]  0.6452006  0.0163198  0.69070983 -1.038505856  0.4817592
#> [10,]  0.6387688 -0.7962681  0.82967557  0.007502175 -0.3180369
```

Each row is an independent simulation.

## Handling discrete models

The usual gradient-based optimization algorithm can’t covnerge in that
framework because the quantile functions of discrete distributions are
constant by pieces.

-   Your model has only one parameter: the Brent method is implemented
    for *flimo*;

-   replace every discrete distribution by a continuous one, e.g. :
    ![Poisson(\theta) \leftarrow Normal(\mu = \theta, \sigma^2 = \theta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Poisson%28%5Ctheta%29%20%5Cleftarrow%20Normal%28%5Cmu%20%3D%20%5Ctheta%2C%20%5Csigma%5E2%20%3D%20%5Ctheta%29 "Poisson(\theta) \leftarrow Normal(\mu = \theta, \sigma^2 = \theta)");

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

There are optimization issues in R for normalized process and
![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")
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
#>       0.032       0.001       0.034

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
#> [1] 0.001893878

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
#>       0.009       0.000       0.011

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
#>  Min.   : 99.06  
#>  1st Qu.:100.69  
#>  Median :101.61  
#>  Mean   :101.37  
#>  3rd Qu.:102.38  
#>  Max.   :102.51  
#> 
#> $minimum
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 2.300e-19 8.400e-18 3.734e-16 2.619e-14 4.065e-14 9.542e-14
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
#> [1] 0.01606452

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
#>  Min.   :2.771   Min.   :1.545  
#>  1st Qu.:2.984   1st Qu.:1.794  
#>  Median :3.341   Median :1.878  
#>  Mean   :3.264   Mean   :1.858  
#>  3rd Qu.:3.495   3rd Qu.:1.902  
#>  Max.   :3.647   Max.   :2.206  
#> 
#> $minimum
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 1.388e-10 5.570e-10 1.803e-09 4.088e-09 6.211e-09 1.506e-08
```
