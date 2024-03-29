% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{flimoptim_R}
\alias{flimoptim_R}
\title{Internal flimoptim function in R mode}
\usage{
flimoptim_R(
  ndraw,
  data = NULL,
  dsumstats = NULL,
  simulatorQ = NULL,
  obj = NULL,
  nsim = 10,
  ninfer = 1,
  lower = 0,
  upper = 1,
  Theta0 = (lower + upper)/2,
  randomTheta0 = FALSE,
  obj_threshold = Inf,
  method = "L-BFGS-B",
  number_tries = 1,
  maxit = 1000,
  factr = 1e+07,
  pgtol = 0,
  show_trace = FALSE,
  store_quantiles = FALSE,
  par_names = NULL
)
}
\arguments{
\item{ndraw}{Integer. Number of random variables to draw
for one simulation of the model.}

\item{data}{1D array containing the observations.}

\item{dsumstats}{Summary statistics to measure distance
between simulations and data.
R function of type dsumstats(simulations, data)
where simulations : 2D array and data : 1D array.
ncol(simulations) = length(data) mandatory.}

\item{simulatorQ}{Simulator of the stochastic process with fixed quantiles
(see README).
Theta is the parameter set for the simulations and
quantiles are drawn in U(0,1).}

\item{obj}{Objective function to minimize.
Default : is directly computed from dsumstats and simulatorQ.
R function of type objective(Theta, quantiles)}

\item{nsim}{Integer. Number of simulations to run for each step
of the optimization algorithm.
Computation time grows linearly with this number. Default to 10.}

\item{ninfer}{Integer. Number of independent inferences to run. Default to 1.}

\item{lower}{1D array. Lower bounds for parameters. Same length as upper.}

\item{upper}{1D array. Upper bounds for parameters. Same length as lower.}

\item{Theta0}{1D array. Initial values of the parameters.
Default : mean(lower, upper).}

\item{randomTheta0}{Boolean.
If True, Theta0 is randomly drawn between lower and upper bounds.}

\item{obj_threshold}{Float. Threshold score. If Final value of objective is
bigger, relaunch the inference if number_tries is not reached.
The purpose is to avoid local minima. Default to Inf (no threshold).}

\item{method}{String. Either "L-BFGS-B" (default) or any other method used by
the base function optim. Stochastic methods do not work with flimo.
If you want to provide bounds, you need to use L-BFGS-B or Brent.}

\item{number_tries}{Integer. Number of tries (inferences) for the objective
value to reach a point lower than obj_threshold. Default to 1.}

\item{maxit}{Integer. Max number of iterations during optimization.
Default to 1000.}

\item{factr}{Float. Control parameter for L-BFGS-B method in stats::optim.
Default to 1e7.}

\item{pgtol}{Float. Control parameter for L-BFGS-B method in stats::optim.
Default to 0.}

\item{show_trace}{Boolean. If true, shows standard trace. Default to false.}

\item{store_quantiles}{Boolean.
If true, stores every quantiles used for inference, to reproduce the results.}

\item{par_names}{vector of names for parameters.
Default is "par1", ..., "parn".}
}
\value{
Object of class flimo_result (list) containing every information
about convergence results.
}
\description{
Computes several parameter inferences with R optimizer.
}
