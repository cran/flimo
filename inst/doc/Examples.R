## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(flimo)
library(ggplot2)
theme_set(theme_bw())

## ----model_definition_1-------------------------------------------------------
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
#See README to know how to build this simulator

ndraw1 <- n1 #number of random draws for one simulation

check_simulator(simulator1, ndraw1, 0, 200)

simulatorQ1 <- function(Theta, quantiles){
  qpois(quantiles, lambda = Theta)
}
check_simulator(simulatorQ1, ndraw1, 0, 200)

#With Normal approximation
simulatorQ1N <- function(Theta, quantiles){
  qnorm(quantiles, mean = Theta, sd = sqrt(Theta))
}

check_simulator(simulatorQ1N, ndraw1, 0, 200)

#Quantile-simulator with Normal approximation

#First simulations

Theta11 <- 50
Theta21 <- 200

nsim1 <- 10

quantiles1 <- matrix(runif(ndraw1*nsim1), nrow = nsim1)

#Just one simulation
simulatorQ1(Theta11, quantiles1[1,])
#No random effect:
simulatorQ1(Theta11, quantiles1[1,]) 

#Independent values:
simulatorQ1(Theta11, quantiles1[2,])

#Matrix of nsim simulations
simu11 <- simulatorQ1(Theta11, quantiles1)
simu21 <- simulatorQ1(Theta21, quantiles1)

#Sample Comparison: summary statistics

dsumstats1 <- function(simulations, data){
  #simulations : 2D array
  #data : 1D array
  mean_simu <- mean(rowMeans(simulations))
  mean_data <- mean(data)
  (mean_simu-mean_data)^2
}

## ----objective_function_1-----------------------------------------------------
#Objective by parameter:
plot_objective(ndraw1, data1, dsumstats1, simulatorQ1,
               nsim = nsim1, lower = 0, upper = 200)

#Objective with Normal approximation :
plot_objective(ndraw1, data1, dsumstats1, simulatorQ1N,
               nsim = nsim1, lower = 0, upper = 200)

#both plots
#We use same quantiles for both

quantiles1 <- matrix(runif(ndraw1*nsim1), nrow = nsim1)

p <- plot_objective(data = data1,
                    dsumstats = dsumstats1,
                    simulatorQ = simulatorQ1,
                    lower = 0, upper = 200, quantiles = quantiles1)

plot_objective(data = data1,
               dsumstats = dsumstats1,
               simulatorQ = simulatorQ1N,
               lower = 0, upper = 200,
               visualize_min = FALSE,
               add_to_plot = p, quantiles = quantiles1)

## ----objective_function_1_zoom------------------------------------------------
p <- plot_objective(data = data1,
                    dsumstats = dsumstats1,
                    simulatorQ = simulatorQ1,
                    lower = 71, upper = 72,
                    visualize_min = FALSE, quantiles = quantiles1)

plot_objective(ndraw1, data1, dsumstats1, simulatorQ1N,
               lower = 71, upper = 72,
               nsim = nsim1,
               visualize_min = FALSE, add_to_plot = p,
               quantiles = quantiles1)

## ----optimization_1-----------------------------------------------------------
#Optimization with normal approximation of Poisson distribution

#Default mode: full R

system.time(optim1R <- flimoptim(ndraw1, data1, dsumstats1, simulatorQ1N,
                 ninfer = 10,
                 nsim = nsim1,
                 lower = 1, upper = 1000,
                 randomTheta0 = TRUE))
optim1R
summary(optim1R)
attributes(optim1R)

plot(optim1R)

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


## ----model_definition_2-------------------------------------------------------
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
#See README to know how to build this simulator

simulatorQ2 <- function(Theta, quantiles){
  qnorm(quantiles, mean = Theta[1], sd = Theta[2])
}

ndraw2 <- 5

check_simulator(simulatorQ2, ndraw2, c(0, 0), c(10, 10))
check_simulator(simulator2, ndraw2, c(0, 0), c(10, 10))


dsumstats2 <-function(simulations, data){
  mean_simu <- mean(rowMeans(simulations))
  mean_data <- mean(data)
  sd_simu <-mean(apply(simulations, 1, sd))
  sd_data <- sd(data)
  (mean_simu-mean_data)^2+(sd_simu-sd_data)^2
}

nsim2 <- 10

## ----objective_2--------------------------------------------------------------
plot_objective(ndraw2, data2, dsumstats2, simulatorQ2,
               index = 1, other_param = c(1, 2, 10),
               nsim = nsim2,
               lower = -5, upper = 10)

## ----optimization_2-----------------------------------------------------------
#Optimization

#Default mode: full R

optim2R <- flimoptim(ndraw2, data2, dsumstats2, simulatorQ2,
                 ninfer = 10, nsim = nsim2,
                 lower = c(-5, 0), upper = c(10, 10),
                 randomTheta0 = TRUE)

optim2R
summary(optim2R)
plot(optim2R, pairwise_par = TRUE, hist = TRUE, par_minimum = TRUE)

