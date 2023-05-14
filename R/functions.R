
#' @import compiler
#' @import JuliaConnectoR
#' @import ggplot2
#' @importFrom stats median optim runif
#' @importFrom utils stack write.csv
#'
NULL

#_______________________________________________________________________________

#' @title Check Julia setup
#'
#' @description Checks installation of Julia and install the needed packages.
#' May take little time to run.Only run the first time you use Jflimo.
#'
#' @return Boolean. True if correct setup, False else.
#'

#' @export

julia_setup <- function(){
  requireNamespace("JuliaConnectoR")
  if (juliaSetupOk()){
    juliaEval('
       import Pkg
       Pkg.add(url = "https://git.metabarcoding.org/lecasofts/flimo/jflimo.git")
       Pkg.add("CSV")
       Pkg.add("DataFrames")
       Pkg.add("Distributions")
       Pkg.add("ForwardDiff")
       Pkg.add("Optim")
       Pkg.add("Random")
              ')
    return(TRUE)
  }
  else {
    warning("FAILURE :
            Uncomplete Julia setup - See JuliaConnectoR's documentation")
    return(FALSE)
  }
}

#_______________________________________________________________________________

#' @title Load Julia
#'
#' @description Load needed Julia packages. Run to use Jflimo.
#'
#' @return Boolean. True if load is done correctly
#'

#' @export

julia_load <- function(){
  tryCatch({
    juliaEval('
          using Jflimo
          using CSV
          using DataFrames
          using Distributions
          using ForwardDiff
          using Optim
          using Random
            ')
    return(TRUE)
  },
  error=function(cond) {
    warning("FAILURE :
            You should run flimo.julia_setup() once
            to complete installation and reload flimo.")
    return(FALSE)
  }
  )
}

#_______________________________________________________________________________

#' @title Objective function minimized by flimo
#'
#' @description Computes the summary statistics between simulations w.r.t.
#' Theta and data. This function is to be minimized by flimoptim.
#'
#' @param Theta 1D array. parameters for the simulations.
#' @param quantiles 2D array containing values drawn in U(0,1).
#' Row number = number of simulations.
#' Column number = number of random variables to draw in one simulation.
#' @param data 1D array containing the observations.
#' @param dsumstats Function computing the distance between simulations and data
#' of form dsumstats(simulations, data) where
#' simulations : 2D array and data : 1D array.
#' ncol(simulations) = length(data) mandatory.
#' @param simulatorQ Function of type simulatorQ(Theta, quantiles)
#' where Theta is the parameter set for the simulations
#' and quantiles are drawn in U(0,1).
#' See README for details.
#'
#' @return Numeric value. Distance between summary statistics of data
#' and simulations w.r.t. Theta.
#'
#' @examples
#'
#'quantiles <- matrix(runif(50), nrow = 10)
#'
#'data <- rep(100, 5)
#'
#'dsumstats <- function(simulations, data){
#' mean_simu <- mean(rowMeans(simulations))
#' mean_data <- mean(data)
#' (mean_simu-mean_data)^2
#'}
#'
#'simulatorQ <- function(Theta, quantiles){
#' qpois(quantiles, lambda = Theta)
#'}
#'
#'flimobjective(100, quantiles, data, dsumstats, simulatorQ)
#'
#' @export

flimobjective <- function(Theta, quantiles, data, dsumstats, simulatorQ){
  # simulations <- t(apply(quantiles, 1,
  #                        simulatorQ, Theta = Theta))
  first_sim <- simulatorQ(Theta, quantiles[1,])
  simulations <- matrix(0, nrow = nrow(quantiles), ncol = length(first_sim))
  simulations[1,] <- first_sim
  if (nrow(quantiles)>=2){
    for (i in 2:nrow(quantiles)){
      simulations[i,] <- simulatorQ(Theta, quantiles[i,])
    }
  }
  dsumstats(simulations, data)
}


#_______________________________________________________________________________

#' @title Main function to use flimo inference
#'
#' @description Computes several parameter inferences with R optimizer or
#' Julia optimizer in a full Julia mode.
#' In R mode (default) : L-BFGS-B optimization or other methods available for
#' the base::optim function.
#' In Julia mode : either IPNewton with or without Automatic Differentiation,
#' Nelder-Mead or Brent optimization.
#' Argument ndraw is mandatory.
#' You need either to provide data, dsumstats AND simulatorQ
#' OR obj.
#' @param ndraw Integer. Number of random variables to draw
#' for one simulation of the model.
#' @param data 1D array containing the observations.
#' @param dsumstats Summary statistics to measure distance
#' between simulations and data.
#' In R mode : R function of type dsumstats(simulations, data)
#' where simulations : 2D array and data : 1D array.
#' ncol(simulations) = length(data) mandatory.
#' In Julia mode : a string containing the script of the Julia function
#' dsumstats(simulations, data). The name "dsumstats" is mandatory.
#' @param simulatorQ Simulator of the stochastic process with fixed quantiles
#' (see README).
#  Either an R function of type simulatorQ(Theta, quantiles) (in mode "R")
#' or a string (in mode "Julia") containing the script of the Julia function
#' simulatorQ(Theta, quantiles).
#' In Julia mode, the name "simulatorQ" is mandatory.
#' Theta is the parameter set for the simulations and
#' quantiles are drawn in U(0,1).
#' @param obj Objective function to minimize.
#' Default : is directly computed from dsumstats and simulatorQ.
#' Either an R function of type objective(Theta, quantiles) (in mode "R")
#' or a string (in mode "Julia") containing the script of the Julia function
#' julia_obj(Theta, quantiles).
#' Warning : could be tricky if mode = "Julia" to call data.
#' In Julia mode, the name "julia_obj" is mandatory.
#' @param nsim Integer. Number of simulations to run for each step
#' of the optimization algorithm.
#' Computation time grows linearly with this number. Default to 10.
#' @param ninfer Integer. Number of independent inferences to run. Default to 1.
#' @param lower 1D array. Lower bounds for parameters. Same length as upper.
#' With Nelder-Mead in Julia mode: only used for starting point.
#' @param upper 1D array. Upper bounds for parameters. Same length as lower.
#' With Nelder-Mead in Julia mode: only used for starting point.
#' @param Theta0 1D array. Initial values of the parameters.
#' Default : mean(lower, upper).
#' @param randomTheta0 Boolean.
#' If True, Theta0 is randomly drawn between lower and upper bounds.
#' @param mode String. "R" (default) or "Julia". See README.
#' @param AD Boolean.
#' Only in Julia mode, uses Automatic Differentiation with IPNewton method.
#' Default to true.
#' @param method String.
#' In Julia mode, allows to choose the optimization method :
#' "IPNewton", "Brent" or "NelderMead". Default : IPNewton.
#' In R mode, allows to choose any of the optimization methods used by
#' base::optim. Default is L-BFGS-B. Random methods do not work with flimo.
#' Bounded methods are L-BFGS-B and Brent.
#' @param obj_threshold Float. Threshold score. If final value of objective is
#' bigger, relaunch the inference if number_tries is not reached.
#' The purpose is to avoid local minima. Default to Inf (no threshold).
#' @param number_tries Integer. Number of tries (inferences) for the objective
#' value to reach a point lower than obj_threshold. Default to 1.
#' @param maxit Integer. Max number of iterations during optimization.
#' Default to 1000.
#' @param time_limit Float. Time limit in second for each inference.
#' Default to no limit. Not available for R mode and Brent method in Julia mode.
#' @param factr Float.
#' In R-mode : control parameter for L-BFGS-B method in stats::optim.
#' Default to 1e7.
#' @param pgtol Float.
#' In R-mode : control parameter for L-BFGS-B method in stats::optim.
#' Default to 0.
#' @param xtol Float.
#' In Julia mode with IPNewton method : xtol option in Optim.Options.
#' Default to 0.
#' @param ftol Float.
#'In Julia mode with IPNewton method : ftol option in Optim.Options.
#Default to 0.
#' @param gtol Float.
#' In Julia mode with IPNewton method : gtol option in Optim.Options.
#' Default to 1e-8.
#' @param reltol Float.
#' In Julia mode with Brent method : reltol of Optim.optimize.
#' Default is sqrt(.Machine$double.eps), about 1e-8.
#' @param abstol Float.
#' In Julia mode with Brent method : abstol of Optim.optimize.
#' Default is .Machine$double.eps, about 1e-16.
#' @param show_trace Boolean. If true, shows standard trace. Default to false.
#' @param store_trace Boolean.
#' If true, stores standard trace as an array of strings.
#' Default to false. Not available for R mode.
#' @param store_quantiles Boolean.
#' If true, stores every quantiles used for inference, to reproduce the results.
#'Default to false.
#' @param par_names vector of names for parameters.
#' Default is "par1", ..., "parn".
#' @param load_julia Boolean. If true, run julia_load.
#' It can take few seconds. Default to False.
#'
#' @return Object of class flimo_result (list)
#' (converted from Julia object in Julia mode) containing every information
#' about convergence results.
#'
#' @examples
#'data <- rep(100, 5)
#'
#'simulatorQ <- function(Theta, quantiles){
#' qpois(quantiles, lambda = Theta)
#'}
#'dsumstats <- function(simulations, data){
#' mean_simu <- mean(rowMeans(simulations))
#' mean_data <- mean(data)
#' (mean_simu-mean_data)^2
#'}
#'
#' flimoptim(5, data, dsumstats, simulatorQ,
#' lower = 50,
#' upper = 150)
#'
#' @export

flimoptim <- function(ndraw,
                      data = NULL,
                      dsumstats = NULL,
                      simulatorQ = NULL,
                      obj = NULL,
                      nsim = 10,
                      ninfer = 1,
                      lower = 0,
                      upper = 1,
                      Theta0 = (lower+upper)/2,
                      randomTheta0 = FALSE,
                      mode = c("R", "Julia"),
                      AD = TRUE,
                      method = "",
                      obj_threshold = Inf,
                      number_tries = 1,
                      maxit = 1e3,
                      time_limit = NaN,
                      factr = 1e7,
                      pgtol = 0,
                      xtol = 0,
                      ftol = 0,
                      gtol = 1e-8,
                      reltol = sqrt(.Machine$double.eps),
                      abstol = .Machine$double.eps,
                      show_trace = FALSE,
                      store_trace = FALSE,
                      store_quantiles = FALSE,
                      par_names = NULL,
                      load_julia = FALSE){
  if (length(mode) > 1){
    mode <- mode[1] #default is R
  }
  if (mode == "Julia"){
    flimoptim_Julia(ndraw, data, dsumstats, simulatorQ,
                    julia_obj = obj,
                    nsim = nsim,
                    ninfer = ninfer,
                    lower = lower,
                    upper = upper,
                    Theta0 = Theta0,
                    randomTheta0 = randomTheta0,
                    AD = AD,
                    method = method,
                    obj_threshold = obj_threshold,
                    number_tries = number_tries,
                    maxit = maxit,
                    time_limit = time_limit,
                    xtol = xtol,
                    ftol = ftol,
                    gtol = gtol,
                    reltol = reltol,
                    abstol = abstol,
                    show_trace = show_trace,
                    store_trace = store_trace,
                    store_quantiles = store_quantiles,
                    par_names = par_names,
                    load_julia = load_julia)
  }
  else {
    if (mode != "R"){
      message("Default mode : full R optimization")
    }
    flimoptim_R(ndraw, data, dsumstats, simulatorQ,
                obj = obj,
                nsim = nsim,
                ninfer = ninfer,
                lower = lower,
                upper = upper,
                Theta0 = Theta0,
                randomTheta0 = randomTheta0,
                obj_threshold = obj_threshold,
                method = method,
                number_tries = number_tries,
                maxit = maxit,
                factr = factr,
                pgtol = pgtol,
                show_trace = show_trace,
                store_quantiles = store_quantiles,
                par_names = par_names)
  }
}

#_______________________________________________________________________________

#' @title Sample function with fixed quantiles
#'
#' @description Replace the sample function in the context of fixed quantiles.
#' Warning : first argument has less features.
#'
#' @param x a vector of one or more elements from which to choose.
#' @param quantiles 1D array containing values drawn in U(0,1).
#' Length has to be equal to size.
#' @param size a non-negative integer giving the number of items to choose.
#' @param replace should sampling be with replacement ?
#' @param prob a vector of probability weights for obtaining
#' the elements of the vector being sampled.
#'
#' @return a vector of length size with elements drawn from x.
#'
#' @examples
#'
#' set.seed(1)
#' quantiles <- runif(40)
#' sampleQ(1:10, 10, quantiles[1:10])
#' sampleQ(1:10, 10, quantiles[11:20])
#' sampleQ(11:20, 10, quantiles[1:10])
#' sampleQ(1:10, 20, quantiles[21:40], replace = TRUE)
#'
#' @export


sampleQ <- function(x, size, quantiles,
                    replace = FALSE, prob = NULL){
  res <- NULL
  if (is.null(prob)) prob <- rep(1, length(x))
  prob <- prob/sum(prob)
  if (replace){
    Sprob <- cumsum(prob)/sum(prob)
    fsample <- function(i){x[which(quantiles[i]<Sprob)[1]]}
    res <- rep(0, size)
    for (i in 1:size){
      res[i] <- fsample(i)
    }
    return(res)
  }
  else {
    for (n in 1:size){
      Sprob <- cumsum(prob)/sum(prob)
      index <- which(quantiles[n]<Sprob)[1]
      i <- x[index]
      x <- x[-index]
      prob <- prob[-index]
      res <- c(res, i)
    }
    return(res)
  }
}

#_______________________________________________________________________________

#' @title Internal flimoptim function in R mode
#'
#' @description Computes several parameter inferences with R optimizer.
#'
#' @param ndraw Integer. Number of random variables to draw
#' for one simulation of the model.
#' @param data 1D array containing the observations.
#' @param dsumstats Summary statistics to measure distance
#' between simulations and data.
#' R function of type dsumstats(simulations, data)
#' where simulations : 2D array and data : 1D array.
#' ncol(simulations) = length(data) mandatory.
#' @param simulatorQ Simulator of the stochastic process with fixed quantiles
#' (see README).
#  R function of type simulatorQ(Theta, quantiles)
#' Theta is the parameter set for the simulations and
#' quantiles are drawn in U(0,1).
#' @param obj Objective function to minimize.
#' Default : is directly computed from dsumstats and simulatorQ.
#' R function of type objective(Theta, quantiles)
#' @param nsim Integer. Number of simulations to run for each step
#' of the optimization algorithm.
#' Computation time grows linearly with this number. Default to 10.
#' @param ninfer Integer. Number of independent inferences to run. Default to 1.
#' @param lower 1D array. Lower bounds for parameters. Same length as upper.
#' @param upper 1D array. Upper bounds for parameters. Same length as lower.
#' @param Theta0 1D array. Initial values of the parameters.
#' Default : mean(lower, upper).
#' @param randomTheta0 Boolean.
#' If True, Theta0 is randomly drawn between lower and upper bounds.
#' @param obj_threshold Float. Threshold score. If Final value of objective is
#' bigger, relaunch the inference if number_tries is not reached.
#' The purpose is to avoid local minima. Default to Inf (no threshold).
#' @param method String. Either "L-BFGS-B" (default) or any other method used by
#' the base function optim. Stochastic methods do not work with flimo.
#' If you want to provide bounds, you need to use L-BFGS-B or Brent.
#' @param number_tries Integer. Number of tries (inferences) for the objective
#' value to reach a point lower than obj_threshold. Default to 1.
#' @param maxit Integer. Max number of iterations during optimization.
#' Default to 1000.
#' @param factr Float. Control parameter for L-BFGS-B method in stats::optim.
#'Default to 1e7.
#' @param pgtol Float. Control parameter for L-BFGS-B method in stats::optim.
#'Default to 0.
#' @param show_trace Boolean. If true, shows standard trace. Default to false.
#' @param store_quantiles Boolean.
#'If true, stores every quantiles used for inference, to reproduce the results.
# Default to False.
#' @param par_names vector of names for parameters.
#' Default is "par1", ..., "parn".
#'
#' @return Object of class flimo_result (list) containing every information
#' about convergence results.

flimoptim_R <- function(ndraw,
                        data = NULL,
                        dsumstats = NULL,
                        simulatorQ = NULL,
                        obj = NULL,
                        nsim = 10,
                        ninfer = 1,
                        lower = 0,
                        upper = 1,
                        Theta0 = (lower+upper)/2,
                        randomTheta0 = FALSE,
                        obj_threshold = Inf,
                        method = "L-BFGS-B",
                        number_tries = 1,
                        maxit = 1e3,
                        factr = 1e7,
                        pgtol = 0,
                        show_trace = FALSE,
                        store_quantiles = FALSE,
                        par_names = NULL){
  if (method == ""){method <- "L-BFGS-B"}
  minimizer <- matrix(rep(NA, ninfer*length(lower)), nrow = ninfer)
  if (is.null(par_names)) colnames(minimizer) <- paste0("par", 1:length(lower))
  else colnames(minimizer) <- par_names
  minimum <- rep(NA, ninfer)
  f_calls <- rep(NA, ninfer)
  g_calls <- rep(NA, ninfer)
  initial_x <- matrix(rep(NA, ninfer*length(lower)), nrow = ninfer)
  converged <- rep(NA, ninfer)
  message <- rep(NA, ninfer)
  if (store_quantiles){
    all_quantiles <- array(rep(NA, ninfer*nsim*ndraw),
                           dim=c(ninfer, nsim, ndraw))
  }
  time_run <- rep(NA, ninfer)

  for (infer in 1:ninfer){
    obj_value <- obj_threshold + 1
    if (is.finite(obj_threshold)) tries <- 0
    else tries <- number_tries - 1
    t <- 0 #time
    while ((obj_value >= obj_threshold) && (tries < number_tries)){
      quantiles <- matrix(runif(ndraw*nsim), nrow = nsim)
      if (is.null(obj)){
        intern_obj <- function(Theta){
          return(flimobjective(Theta, quantiles, data, dsumstats, simulatorQ))
        }
      }
      else {
        intern_obj <- function(Theta){
          obj(Theta, quantiles)
        }
      }
      if (randomTheta0){
        Theta0 <- runif(length(lower))*(upper-lower)+lower
      }
      start_time <- Sys.time()
      if (!any(method == c("L-BFGS-B", "Brent"))){
        opt <- stats::optim(par = Theta0, fn = intern_obj,
                            method = method,
                            control = list(trace = show_trace,
                                           maxit = maxit,
                                           factr = factr,
                                           pgtol = pgtol))
      }
      else {
        opt <- stats::optim(par = Theta0, fn = intern_obj,
                            method = method,
                            lower = lower, upper = upper,
                            control = list(trace = show_trace,
                                           maxit = maxit,
                                           factr = factr,
                                           pgtol = pgtol))
      }
      end_time <- Sys.time()
      t <- t + end_time - start_time
      obj_value <- opt$value
      tries <- tries + 1
    }
    time_run[infer] <- t
    initial_x[infer,] <- Theta0
    minimizer[infer,] <- opt$par
    minimum[infer] <- opt$value
    converged[infer] <- opt$convergence == 0
    f_calls[infer] <- opt$counts[1]
    g_calls[infer] <- opt$counts[2]
    if (!is.null(opt$message)){message[infer] <- opt$message}
    if (store_quantiles){
      all_quantiles[infer,,] <- quantiles
    }
  }
  optim_result <- NULL
  optim_result$mode <- "R"
  optim_result$method <- method
  optim_result$AD <- FALSE
  optim_result$minimizer <- minimizer
  optim_result$minimum <- minimum
  optim_result$converged <- converged
  optim_result$initial_x <- initial_x
  optim_result$f_calls <- f_calls
  optim_result$g_calls <- g_calls
  optim_result$message <- message
  if (store_quantiles){
    optim_result$quantiles <- all_quantiles
  }
  optim_result$time_run <- time_run

  class(optim_result) <- "flimo_result"

  optim_result
}

#_______________________________________________________________________________

#' @title Internal flimoptim function in Julia mode
#'
#' @description Computes several parameter inferences with Julia optimizer and
#' either IPNewton with or without Automatic Differentiation, Nelder-Mead
#' or Brent method.
#'
#' @param ndraw Integer. Number of random variables to draw
#' for one simulation of the model.
#' @param data 1D array containing the observations.
#' @param dsumstats Summary statistics to measure distance
#' between simulations and data.
#' String containing the script of the Julia function
#' dsumstats(simulations, data).
#' The name "dsumstats" is mandatory.
#' @param simulatorQ Simulator of the stochastic process
#' with fixed quantiles (see README).
#  String containing the script of the Julia function
#' simulatorQ(Theta, quantiles).
#' The name "simulatorQ" is mandatory.
#' Theta is the parameter set for the simulations and
#' quantiles are drawn in U(0,1).
#' @param julia_obj Objective function to minimize.
#' Default : is directly computed from dsumstats and simulatorQ.
#' String containing the script of the Julia function
#' julia_obj(Theta, quantiles).
#' Warning : can be tricky to call data.
#' The name "julia_obj" is mandatory.
#' @param nsim Integer. Number of simulations to run for each step
#' of the optimization algorithm.
#' Computation time grows linearly with this number. Default to 10.
#' @param ninfer Integer. Number of independent inferences to run. Default to 1.
#' @param lower 1D array. Lower bounds for parameters. Same length as upper.
#' @param upper 1D array. Upper bounds for parameters. Same length as lower.
#' @param Theta0 1D array. Initial values of the parameters.
#' Default : mean(lower, upper).
#' @param randomTheta0 Boolean.
#' If True, Theta0 is randomly drawn between lower and upper bounds.
#' @param AD Boolean.
#' Only in Julia mod, uses Automatic Differentiation with IPNewton method.
#' Default to true.
#' @param method String.
#' Allows to choose the optimization method :
#' "Brent", "IPNewton" or "NelderMead". Default : IPNewton.
#' @param obj_threshold Float. Threshold score. If Final value of objective is
#' bigger, relaunch the inference if number_tries is not reached.
#' The purpose is to avoid local minima. Default to Inf (no threshold).
#' @param number_tries Integer. Number of tries (inferences) for the objective
#' value to reach a point lower than obj_threshold. Default to 1.
#' @param maxit Integer. Max number of iterations during optimization.
#' Default to 1000.
#' @param time_limit Float. Time limit in second for each inference.
#' Default to no limit. Not available for Brent method.
#' @param xtol Float. With IPNewton method : xtol option in Optim.Options.
#' Default to 0.
#' @param ftol Float. With IPNewton method : ftol option in Optim.Options.
#' Default to 0.
#' @param gtol Float. With IPNewton method : gtol option in Optim.Options.
#' Default to 1e-8.
#' @param reltol Float. With Brent method : reltol of Optim.optimize.
#' Default is sqrt(.Machine$double.eps), about 1e-8.
#' @param abstol Float. With Brent method : abstol of Optim.optimize.
#' Default is .Machine$double.eps, about 1e-16.
#' @param show_trace Boolean. If true, shows standard trace. Default to false.
#' @param store_trace Boolean.
#' If true, stores standard trace as an array of strings.
#' Default to false. Not available for R mod.
#' @param store_quantiles Boolean.
#' If true, stores every quantiles used for inference, to reproduce the results.
#' Default to false.
#' @param par_names vector of names for parameters.
#' Default is "par1", ..., "parn".
#' @param load_julia Boolean. If true, run julia_load. It can take few seconds.
#' Default to False.
#'
#' @return Object of class flimo_result (list) converted from Julia object
#' containing every information about convergence results.
#'

flimoptim_Julia <- function(ndraw,
                            data = NULL,
                            dsumstats = NULL,
                            simulatorQ = NULL,
                            julia_obj = NULL,
                            nsim = 10,
                            ninfer = 1,
                            lower = 0,
                            upper = 1,
                            Theta0 = (lower+upper)/2,
                            randomTheta0 = FALSE,
                            AD = TRUE,
                            method = "",
                            obj_threshold = Inf,
                            number_tries = 1,
                            maxit = 1e3,
                            time_limit = NULL,
                            xtol = 0,
                            ftol = 0,
                            gtol = 1e-8,
                            reltol = sqrt(.Machine$double.eps),
                            abstol = .Machine$double.eps,
                            show_trace = FALSE,
                            store_trace = FALSE,
                            store_quantiles = FALSE,
                            par_names = NULL,
                            load_julia = FALSE){
  if (load_julia){
    julia_load()
  }

  if (!is.null(data)){
    file_data <- tempfile(pattern = "data_flimoptim", fileext = ".csv")
    write.csv(as.matrix(data), file_data,
              row.names = FALSE)
    cmd_julia <- ""
    juliaEval(paste0('julia_data = CSV.read("',file_data,'", DataFrame)[:,1]\n',
                     'julia_data = convert(Array{Float64,1}, julia_data)'))
    unlink(file_data)
  }
  else{
    cmd_julia <- 'julia_data = nothing'
  }
  cmd_julia <- paste0(cmd_julia,
                   paste0('julia_ndraw = ',ndraw,
                   '\njulia_nsim = ',nsim,
                   '\njulia_ninfer = ',ninfer,
                   '\njulia_method = "',method,'"',
                   '\njulia_obj_threshold = ', obj_threshold,
                   '\njulia_number_tries = ', number_tries,
                   '\njulia_maxit = ',maxit,
                   '\njulia_Theta0 = ',"[",paste(Theta0, collapse = ","),"]",
                   '\njulia_lower = ',"[",paste(lower, collapse = ","),"]",
                   '\njulia_upper = ',"[",paste(upper, collapse = ","),"]",
                   '\njulia_xtol = ',xtol,
                   '\njulia_ftol = ',ftol,
                   '\njulia_gtol = ',gtol,
                   '\njulia_reltol = ',reltol,
                   '\njulia_abstol = ',abstol),
                   sep = '\n')

  if (is.null(time_limit)){cmd_julia <- paste0(cmd_julia, 'julia_time_limit = NaN',
                                            sep = '\n')}
  else {cmd_julia <- paste0(cmd_julia, paste0('julia_time_limit = ',time_limit),
                         sep = '\n')}

  if (AD){cmd_julia <- paste0(cmd_julia, 'julia_AD = true', sep = '\n')}
  else {cmd_julia <- paste0(cmd_julia, 'julia_AD = false', sep = '\n')}

  if (randomTheta0){cmd_julia <- paste0(cmd_julia, 'julia_randomTheta0 = true',
                                     sep = '\n')}
  else {cmd_julia <- paste0(cmd_julia, 'julia_randomTheta0 = false', sep = '\n')}

  if (show_trace){cmd_julia <- paste0(cmd_julia, 'julia_show_trace = true',
                                   sep = '\n')}
  else {cmd_julia <- paste0(cmd_julia, 'julia_show_trace = false', sep = '\n')}

  if (store_trace){cmd_julia <- paste0(cmd_julia, 'julia_store_trace = true',
                                    sep = '\n')}
  else {cmd_julia <- paste0(cmd_julia, 'julia_store_trace = false', sep = '\n')}

  if (store_quantiles){cmd_julia <- paste0(cmd_julia, 'julia_store_quantiles = true',
                                        sep = '\n')}
  else {cmd_julia <- paste0(cmd_julia, 'julia_store_quantiles = false', sep = '\n')}

  if (is.null(julia_obj)){
    cmd_julia <- paste0(cmd_julia,
                     paste0('julia_obj = nothing','\n', dsumstats, '\n', simulatorQ),
                     sep = '\n')
  }
  else {
    cmd_julia <- paste0(cmd_julia,
                     paste0('julia_obj = nothing','\n', dsumstats, '\n', simulatorQ),
                     sep = '\n')
  }
  cmd_julia <- paste0(cmd_julia,
                   "
  julia_xtol = convert(Float64, julia_xtol)
  julia_ftol = convert(Float64, julia_ftol)
  julia_gtol = convert(Float64, julia_gtol)
  julia_reltol = convert(Float64, julia_reltol)
  julia_abstol = convert(Float64, julia_abstol)
  julia_Theta0 = convert(Array{Float64,1}, julia_Theta0)
  julia_lower = convert(Array{Float64,1}, julia_lower)
  julia_upper = convert(Array{Float64,1}, julia_upper)

  opt = Jflimo.flimoptim(julia_ndraw,
                   data = julia_data,
                   dsumstats = dsumstats,
                   simulatorQ = simulatorQ,
                   obj = julia_obj,
                   nsim = julia_nsim,
                   ninfer = julia_ninfer,
                   lower = julia_lower,
                   upper = julia_upper,
                   Theta0 = julia_Theta0,
                   randomTheta0 = julia_randomTheta0,
                   AD = julia_AD,
                   method = julia_method,
                   obj_threshold = julia_obj_threshold,
                   number_tries = julia_number_tries,
                   maxit = julia_maxit,
                   xtol = julia_xtol,
                   ftol = julia_ftol,
                   gtol = julia_gtol,
                   reltol = julia_reltol,
                   abstol = julia_abstol,
                   time_limit = julia_time_limit,
                   show_trace = julia_show_trace,
                   store_trace = julia_store_trace,
                   store_quantiles = julia_store_quantiles)", sep = '\n')

  optim_result <- juliaGet(juliaEval(cmd_julia))
  class(optim_result) <- "flimo_result"
  optim_result$mode <- "Julia"
  optim_result$AD <- AD
  if (is.null(par_names))
    colnames(optim_result$minimizer) <- paste0("par", 1:length(lower))
  else colnames(optim_result$minimizer) <- par_names
  names(optim_result)[names(optim_result) == "iteration_converged"] <-
    "iteration_limit_reached"

  optim_result
}

#_______________________________________________________________________________

#' @title Print flimo results
#'
#' @description Prints most important information about inference results.
#'
#' @param x Object of class flimo_result from any mod/method algorithm
#' of the flimo package.
#' @param ... optional args for generic method
#'
#' @return String containing most important information about argument
#' of class flimo_result.
#'
#' @export
#'

print.flimo_result <- function(x, ...){
  value <- 'optimization Result\n'
  value <- paste0(value, "Mode : ", x$mod,"\n")
  if (length(as.character(x$method)) == 1){
    value <- paste0(value, "method : ", as.character(x$method),"\n")
  }
  else{
    value <- paste0(value, "method : ", attr(x$method,"JLTYPE"),"\n")
  }
  value <- paste0(value, "Number of inferences : ", length(x$minimum),"\n")
  value <- paste0(value, "Convergence : ", length(which(x$converged)),"/",
                  length(x$minimum),"\n")
  value <- paste0(value,"Mean of minimizer :\n")
  value <- paste0(value, paste(colMeans(x$minimizer), collapse = "    "),"\n")
  value <- paste0(value, "Best minimizer :\n")
  value <- paste0(value,
                  paste(x$minimizer[which(x$minimum == min(x$minimum))[1],],
                        collapse = "    "))
  value <- paste0(value, "\nReached minima : from ", min(x$minimum), " to ",
                  max(x$minimum))
  value <- paste0(value, "\nMedian time by inference ", median(x$time_run))
  cat(value)
  value
}

#' @title Summary of flimo results
#'
#' @description Most important information about inference results.
#'
#' @param object Object of class flimo_result from any mode/method algorithm
#' of the Jflimo package.
#' @param ... optional args for generic method summary
#'
#' @return List containing most important information about argument
#' of class flimo_result.

#' @export
#'

summary.flimo_result <- function(object, ...){
  value <- list()
  value$Mode <- object$mode

  if (length(as.character(object$method)) == 1){
    value$method <- as.character(object$method)
  }
  else{
    value$method <- attr(object$method,"JLTYPE")
  }
  value$number_inferences <- length(object$minimum)

  if (!is.null(attr(object$method,"JLTYPE"))){
    value$number_converged <- length(which(object$converged))
  }
  else {
    value$number_converged <- length(which(object$converged))
  }
  value$minimizer <- summary(object$minimizer)
  value$minimum <- summary(object$minimum)
  value$median_time_inference <- median(object$time_run)
  value
}

#_______________________________________________________________________________

#' @title Check if simulator with fixed quantiles is well implemented
#'
#' @description Run simulations to catch random variations.
#' Warning : does not check it formally.
#' Warning : does not check if quantiles are used several times.
#'
#' @param simulatorQ Function of type simulatorQ(Theta, quantiles)
#' where Theta is the parameter set for the simulations and
#' quantiles are drawn in U(0,1).
#' @param ndraw Integer. Number of random variables to draw
#' for one simulation of the model.
#' @param Theta_lower 1D numeric array. Lower bounds of Theta parameters.
#' @param Theta_upper 1D numeric array. Upper bounds of Theta parameters.
#' @param ntheta Integer. Number of Theta parameters to test.
#' @param nruns Integer. For each Theta, number of simulations to run.
#'
#' @return Boolean. True if no random effect was detected, False else.
#'
#' @examples
#' simulatorQ <- function(Theta, quantiles){
#' qpois(quantiles, lambda = Theta)
#'}
#'check_simulator(simulatorQ, 5,
#'Theta_lower = 50, Theta_upper = 150)
#'
#' @export
#'

check_simulator <- function(simulatorQ, ndraw, Theta_lower = 0, Theta_upper = 1,
                            ntheta = 5, nruns = 3){

  if (length(Theta_lower) != length(Theta_upper)){
    warning("Warning : no correct dimension")
    stop()
  }
  if (nruns <= 1){
    message("check_simulator is useless with nruns <= 1. nruns set to 2.")
    nruns <- 2
  }
  for (i in 1:ntheta){
    Theta <- runif(length(Theta_lower))*(Theta_upper-Theta_lower)+Theta_lower
    quantiles_1sim <- runif(ndraw)
    quantiles <- matrix(rep(quantiles_1sim, nruns), nrow = nruns, byrow = TRUE)
    simulations <- tryCatch({
      matrix(simulatorQ(Theta, quantiles), nrow = nruns)
    }, error = function(e){
      message("This simulator is not conceived to work with quantiles")
      return(FALSE)
    })
    tryCatch({if (!all(apply(simulations, 2, function(x)
      length(unique(x)) == 1))){
      return(FALSE)
    }},error = function(e){
      return(FALSE)
    })
  }
  TRUE
}

#_______________________________________________________________________________

#' @title Plot main flimo results
#'
#' @description Shows the plots for most important inference results.
#' Default only shows normalized boxplots for each inferred parameter.
#'
#' @param x Object of class flimo_result.
#' @param y unused generic argument.
#' @param ... optional args for generic method
#' @param bins Integer. Number of bins if hist is True.
#' @param hist Boolean. If True, plots the histogram of each inferred parameter.
#' Default to false.
#' @param par_minimum Boolean.
#' If True, plots each inferred parameter by reached minimum. Default to false.
#' @param pairwise_par Boolean.
#' If True, plots each pairs of inferred parameters.Default to false.
#' @param boxplot Boolean.
#' If True, plots the boxplots of each inferred parameter scaled by their mean.
#' Default to true.
#' @param par_names Vector of names for parameters.
#' Default is "par1", ..., "parn".
#'
#' @return Nothing. Prints the asked ggplot objects.
#'
#' @export


plot.flimo_result <- function(x, y, ...,
                              hist = FALSE,
                              bins = 1+as.integer(nrow(x$minimizer)^(1/3)),
                              par_minimum = FALSE,
                              pairwise_par = FALSE,
                              boxplot = TRUE,
                              par_names = NULL){

  ind <- NULL #to remove R CHECK NOTE
  values <- NULL #to remove R CHECK NOTE
  npar <- ncol(x$minimizer)

  for (par in 1:npar){
    if (hist){ #Histogram for each inferred parameter
      print(ggplot2::ggplot()+
              geom_histogram(aes(x$minimizer[,par]), bins = bins)+
              labs(x = paste0("parameter ", par))+
              ggtitle(paste0("Histogram of inferred parameter ", par)))
    }

    if (par_minimum){#Plot inferred par = f(minimum)
      print(ggplot2::ggplot()+
              geom_point(aes(x$minimizer[,par], x$minimum))+
              scale_y_log10()+
              labs(x = paste0("parameter ", par), y = "reached minimum")+
              ggtitle(paste0("inferred parameter ", par,
                             " by reached minimum")))
    }

    if (pairwise_par){ #Plot par1 = f(par2)
      if (npar > 1 & par < npar){
        for (par2 in (par+1):npar){
          print(ggplot2::ggplot()+
                  geom_point(aes(x$minimizer[,par], x$minimizer[,par2]))+
                  labs(x = paste0("Parameter ", par),
                       y = paste0("Parameter ", par2))+
                  ggtitle(paste0("Parameter ", par2, " by parameter ", par)))
        }
      }
    }
  }
  if (boxplot){ #Normalized boxplots
    aux <- as.data.frame(x$minimizer)
    if (is.null(par_names)) colnames(aux) <- paste0("par", 1:ncol(aux))
    else colnames(aux) <- par_names
    aux <- stack(aux)
    for (par in unique(aux$ind)){
      aux[aux$ind == par,]$values <-
        aux[aux$ind == par,]$values/mean(aux[aux$ind == par,]$values)
    }
    print(ggplot2::ggplot(aux)+
            geom_boxplot(aes(x = ind, y = values))+
            labs(x = "Parameters", y = "Value/mean(Value)")+
            ggtitle(paste0("Normalized boxplots\nfor each inferred parameter")))
  }
}

#_______________________________________________________________________________

#' @title Plot the objective to be minimized using flimo
#'
#' @description Plot of the objective function with one parameter moving
#' (objective = f(theta_index)).
#' You need either to provide data, dsumstats AND simulatorQ
#' OR obj.
#'
#' @param ndraw Integer.Number of random variables to draw
#' for one simulation of the model.
#' @param data 1D array containing the observations.
#' @param dsumstats Function computing the distance
#' between simulations and data of form dsumstats(simulations, data)
#' where simulations : 2D array and data : 1D array.
#' ncol(simulations) = length(data) mandatory.
#' @param simulatorQ Function of type simulatorQ(Theta, quantiles)
#' where Theta is the parameter set for the simulations and
#' quantiles are drawn in U(0,1).
#' @param obj objective function of type objective(Theta, quantiles).
#' Default : directly computed with "dsumstats" and "simulatorQ".
#' @param quantiles 2D array containing values drawn in U(0,1).
#' Row number = number of simulations. Default: simulated within the function.
#' Column number = number of random variables to draw in one simulation.
#' @param index Integer. Index of the moving parameter.
#' @param other_param Other parameters of the model. If NULL : assume 1D-model.
#' If numeric : 2D-model, one curve.
#' If 1D-array and dim2 is True (default) :
#' 2D-model, one curve by value in other_param.
#' If 1D-array and dim2 is False or 2D-array : (n>2)D-model,
#' one curve by row in other_param.
#' If your model has n>2 dimensions,
#' you should define other_param as a matrix even if
#' you have only one parameter set to test
#' (with as.matrix(t(vect_param)) where vect_param is a 1D-array).
#' @param nsim Integer. Number of simulations to run for each step
#' of the optimization algorithm.
#' Computation time grows linearly with this number. Default to 10.
#' @param lower Numeric. Lower value of the plot.
#' @param upper Numeric. Upper value of the plot.
#' @param dim2 Boolean. True if model is 2-dimensional.
#' @param visualize_min Boolean. If True, show explicitly the minimum point.
#' @param plot_legend Boolean. If True (default), plots the legend.
#' @param npoints Integer. Number of points evaluated. Default = 300.
#' @param add_to_plot ggplot object. If not NULL,
#' will add all curves/points on previous plot instead of creating a new one.
#' Does not change title/labels/limits defined in previous plot.
#'
#' @return ggplot object representing the objective function to be minimized.
#'
#' @examples
#'data <- rep(100, 5)
#'
#'dsumstats <- function(simulations, data){
#' mean_simu <- mean(rowMeans(simulations))
#' mean_data <- mean(data)
#' (mean_simu-mean_data)^2
#'}
#'
#'simulatorQ <- function(Theta, quantiles){
#' qpois(quantiles, lambda = Theta)
#'}
#'
#' plot_objective(5, data, dsumstats, simulatorQ,
#' lower = 0, upper = 200)
#'
#' @export

plot_objective <- function(ndraw,
                           data = NULL,
                           dsumstats = NULL,
                           simulatorQ = NULL,
                           obj = NULL,
                           quantiles = NULL,
                           index = NULL,
                           other_param = NULL,
                           nsim = 10,
                           lower = 0,
                           upper = 1,
                           dim2 = TRUE,
                           visualize_min = TRUE,
                           plot_legend = TRUE,
                           npoints = 300,
                           add_to_plot = NULL){

  if (is.null(quantiles)){q <- matrix(runif(ndraw*nsim), nrow = nsim)}
  else {q <- quantiles}
  if (is.null(obj)){intern_obj <- function(Theta) flimobjective(Theta, q, data,
                                                dsumstats, simulatorQ)}
  else {intern_obj <- function(Theta) obj(Theta, q)}

  if (is.null(other_param)){
    return(plot_objective1D(intern_obj, lower = lower, upper = upper,
                            visualize_min = visualize_min,
                            npoints = npoints,
                            add_to_plot = add_to_plot))
  }
  else if (is.null(dim(other_param))){
    if (dim2){
      other_param <- as.matrix(other_param) #convert to length*1 matrix

      return(plot_objectivenD(intern_obj, index, other_param,
                              lower = lower, upper = upper,
                              visualize_min = visualize_min,
                              plot_legend = plot_legend,
                              npoints = npoints,
                              add_to_plot = add_to_plot))
    }
    else{
      other_param <- as.matrix(t(other_param)) #convert to 1*length matrix

      return(plot_objectivenD(intern_obj, index, other_param,
                              lower = lower,
                              upper = upper,
                              visualize_min = visualize_min,
                              plot_legend = plot_legend,
                              npoints = npoints,
                              add_to_plot = add_to_plot))
    }
  }
  else {
    return(plot_objectivenD(intern_obj, index, other_param,
                            lower = lower,
                            upper = upper,
                            visualize_min = visualize_min,
                            plot_legend = plot_legend,
                            npoints = npoints,
                            add_to_plot = add_to_plot))
  }
}

plot_objective1D <- function(obj,
                             lower = 0,
                             upper = 1,
                             visualize_min = TRUE,
                             npoints = 300,
                             add_to_plot = NULL){
  #plot_objective when model has 1 parameter
  #internal function

  x <- seq(lower, upper, length.out = npoints)
  y <- sapply(x, FUN = obj)

  if (!is.null(add_to_plot)){
    p <- add_to_plot
  }
  else{
    p<-ggplot2::ggplot()
  }

  p <- p+geom_line(aes(x,y))

  if (is.null(add_to_plot)){
    p <- p + xlim(lower, upper)+
      ggtitle(paste("Objective function to minimize by parameter 1"))+
      labs(x = paste0("theta_1"), y = "objective value")
  }

  if (visualize_min){
    #approaching minimum with y values
    ymin <- min(y, na.rm = TRUE)
    xmin <- x[which.min(y)]

    p <- p + geom_point(aes(xmin, ymin)) +
      geom_segment(aes(x = xmin, xend = xmin, y = 0, yend = ymin), size = 0.25)
  }

  p
}

plot_objectivenD <- function(obj,
                             index,
                             other_param,
                             lower = 0,
                             upper = 1,
                             visualize_min = TRUE,
                             plot_legend = TRUE,
                             npoints = 300,
                             add_to_plot = NULL){
  #plot_objective when model has n>=2 parameters
  #internal function
  #other_param has to be a 2D-array (conversion done in plot_objective function)

  other_par <- NULL #remove R CHECK NOTE
  xmin <- #remove R CHECK NOTE
    ymin <- #remove R CHECK NOTE

    if (index < 1 | index > nrow(other_param)+1){
      warning("FAILURE : Index not valid")
      stop()
    }

  x <- seq(lower, upper, length.out = npoints)
  obj_values <- NULL

  if (visualize_min){
    minima <- NULL
  }

  aux <- function(x, other_par){
    Theta <- rep(NA, (length(other_par)+1))
    Theta[index] <- x
    if (index > 1){
      Theta[1:(index-1)] <- other_par[1:(index-1)]
    }
    if (index < length(other_par)+1){
      Theta[(index+1):(length(other_par)+1)] <-
        other_par[index:length(other_par)]
    }
    obj(Theta)
  }

  for (rpar in 1:nrow(other_param)){
    par <- other_param[rpar,]

    par_ref <- rep(NA, length(par)+1)
    par_ref[index] <- "x"
    if (index > 1){
      par_ref[1:(index-1)] <- as.character(par[1:(index-1)])
    }
    if (index < length(par)+1){
      par_ref[(index+1):(length(par)+1)] <- as.character(par[index:length(par)])
    }
    par_ref <- paste0(par_ref, collapse = ";")
    obj_values <- rbind(obj_values,
                        data.frame(x = x,
                                   y = sapply(x, FUN = aux, other_par = par),
                                   other_par = as.factor(par_ref)))

    if (visualize_min){
      x <- obj_values[obj_values$other_par == par_ref,]$x
      y <- obj_values[obj_values$other_par == par_ref,]$y
      minima <- rbind(minima, data.frame(xmin = x[which.min(y)],
                                         ymin = min(y, na.rm = TRUE),
                                         other_par = as.factor(par_ref)))
    }
  }

  if (!is.null(add_to_plot)){
    p <- add_to_plot
  }
  else{
    p <- ggplot2::ggplot(obj_values)
  }

  p <- p + geom_line(data = obj_values, aes(x,y, color = other_par))+
    xlim(lower, upper)+
    ggtitle(paste("Objective function to minimize by parameter", index))+
    labs(x = paste0("theta_", index), y = "objective value",
         color = "Param.")

  if (is.null(add_to_plot)){
    p <- p +
      ggtitle(paste("Objective function to minimize by parameter",index))+
      labs(x = paste0("theta_", index), y = "objective value")
  }

  if (visualize_min){
    p <- p + geom_point(data = minima, aes(xmin, ymin, color = other_par))+
      geom_segment(data = minima, aes(x = xmin, xend = xmin,
                                      y = 0, yend = ymin, color = other_par),
                   size = 0.25)
  }
  if (!plot_legend){
    p <- p + theme(legend.position = "none")
  }
  p
}

