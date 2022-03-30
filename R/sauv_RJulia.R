#sauvegarde mode R-Julia


#' #_______________________________________________________________________________
#'
#' #' @title flimoptim_RJulia
#' #'
#' #' @description Computes several parameter inferences with Julia optimizer (imported to R) and IPNewton without Automatic Differentiation or Brent optimization.
#' #'
#' #' @param data 1D array containing the observations.
#' #' @param ndraw Integer. Number of random variables to draw for one trajectory of the model.
#' #' @param sumstats Summary statistics to measure distance between simulations and data.
#' #' R function of type sumstats(simulations, data) where simulations : 2D array and data : 1D array.
#' #' ncol(simulations) = length(data) mandatory.
#' #' @param simulatorQ Simulator of the stochastic process with fixed quantiles (see Readme).
#' #  R function of type simulatorQ(Theta, quantiles)
#' #' Theta is the parameter set for the simulations and quantiles are drawn in U(0,1).
#' #' @param obj Objective function to minimize. Default : is directly computed from sumstats and simulatorQ.
#' #' R function of type objective(Theta, quantiles)
#' #' @param nsim Integer. Number of simulations to run for each step of optimisation algorithm.
#' #' Computation time grows linearly with this number. Default to 1000.
#' #' @param ninfer Integer. Number of independent inferences to run. Default to 1.
#' #' @param lower 1D array. Lower bounds for parameters. Same length as upper.
#' #' @param upper 1D array. Upper bounds for parameters. Same length as lower.
#' #' @param Theta0 1D array. Initial values of the parameters. Default : mean(lower, upper).
#' #' @param randomTheta0 Boolean. If True, Theta0 is randomly drawn between lower and upper bounds.
#' #' @param method String. Allows to choose the optimization method : "Brent", "IPNewton". Default : IPNewton.
#' #' @param maxit Integer. Max number of iterations during optimization. Default to 1000.
#' #' @param time_lim Float. Time limit in second for each inference. Default to no limit. Not available for Brent method.
#' #' @param xtol Float. With IPNewton method : xtol option in Optim.Options. Default to 0.
#' #' @param ftol Float. With IPNewton method : ftol option in Optim.Options. Default to 0.
#' #' @param gtol Float. With IPNewton method : gtol option in Optim.Options. Default to 1e-8.
#' #' @param reltol Float. With Brent method : reltol of Optim.optimize. Default is sqrt(.Machine$double.eps), about 1e-8.
#' #' @param abstol Float. With Brent method : abstol of Optim.optimize. Default is .Machine$double.eps, about 1e-16.
#' #' @param show_trace Boolean. If true, shows standard trace. Default to false.
#' #' @param store_trace Boolean. If true, stores standard trace as an array of strings. Default to false. Not available for R mod.
#' #' @param store_quantiles Boolean. If true, stores every quantiles used for inference, to reproduce the results. Default to false.
#' #'
#' #' @return Object of class flimo_result (list) converted from Julia object containing every information about convergence results.
#'
#'
#' flimoptim_RJulia <- function(data, ndraw, sumstats, simulatorQ,
#'                              obj = NULL,
#'                              nsim = 1000,
#'                              ninfer = 1,
#'                              lower = 0,
#'                              upper = 1,
#'                              Theta0 = (lower+upper)/2,
#'                              randomTheta0 = F,
#'                              method = "IPNewton",
#'                              maxit = 1e3,
#'                              time_lim = NaN,
#'                              xtol = 0,
#'                              ftol = 0,
#'                              gtol = 1e-8,
#'                              reltol = sqrt(.Machine$double.eps),
#'                              abstol = .Machine$double.eps,
#'                              show_trace = F,
#'                              store_trace = F,
#'                              store_quantiles = F){
#'
#'   if (!(method %in% c("Brent", "IPNewton"))){
#'     print("Warning : unknown method... Set to IPNewton.")
#'     method <- "IPNewton"
#'   }
#'
#'   julia_optimize <- juliaFun("Optim.optimize")
#'   julia_Options <- juliaFun("Optim.Options")
#'   Julia_TwiceDifferentiable <- juliaFun("Optim.TwiceDifferentiable")
#'   Julia_TwiceDifferentiableConstraints <- juliaFun("Optim.TwiceDifferentiableConstraints")
#'   julia_IPNewton <- juliaFun("Optim.IPNewton")
#'   optim_result <- list()
#'
#'   minimizer <- matrix(rep(NA, ninfer*length(lower)), nrow = ninfer)
#'   colnames(minimizer) <- paste0("par", 1:length(lower))
#'   minimum <- rep(NA, ninfer)
#'   iterations <- rep(NA, ninfer)
#'   iteration_converged <- rep(NA, ninfer)
#'   f_calls <- rep(NA, ninfer)
#'   converged <- rep(NA, ninfer)
#'   time_run <- rep(NA, ninfer)
#'
#'   if (store_trace){
#'     trace <- rep(list(), ninfer)
#'   }
#'
#'   if (store_quantiles){
#'     all_quantiles <- array(rep(NA, ninfer*nsim*ndraw), dim=c(ninfer, nsim, ndraw))
#'   }
#'
#'   if (method == "IPNewton"){
#'     initial_x <- matrix(rep(NA, ninfer*length(lower)), nrow = ninfer)
#'     x_converged <- rep(NA, ninfer)
#'     x_abschange <- rep(NA, ninfer)
#'     f_converged <- rep(NA, ninfer)
#'     f_abschange <- rep(NA, ninfer)
#'     f_relchange <- rep(NA, ninfer)
#'     f_increased <- rep(NA, ninfer)
#'     g_converged <- rep(NA, ninfer)
#'     g_residual <- rep(NA, ninfer)
#'     g_calls <- rep(NA, ninfer)
#'     h_calls <- rep(NA, ninfer)
#'   }
#'
#'   for (infer in 1:ninfer){
#'     quantiles <- matrix(runif(ndraw*nsim), nrow = nsim)
#'     if (is.null(obj)){
#'       intern_obj <- function(Theta){
#'         return(flimobjective(Theta, quantiles, data, sumstats, simulatorQ))
#'       }
#'     }
#'     else {
#'       intern_obj <- function(Theta){
#'         obj(Theta, quantiles)
#'       }
#'     }
#'     if (method == "Brent"){
#'       start_time <- Sys.time()
#'       opt <- juliaGet(julia_optimize(intern_obj, lower, upper,
#'                                      iterations = as.integer(maxit),
#'                                      store_trace = store_trace,
#'                                      show_trace = show_trace,
#'                                      rel_tol = reltol,
#'                                      abs_tol = abstol))
#'       end_time <- Sys.time()
#'
#'       time_run[infer] <- end_time - start_time
#'
#'       converged[infer] <- opt$converged
#'     }
#'     else {
#'       if (randomTheta0){
#'         Theta0 = runif(length(lower))*(upper-lower)+lower
#'       }
#'
#'       df = Julia_TwiceDifferentiable(intern_obj, as.list(Theta0))
#'       dfc = Julia_TwiceDifferentiableConstraints(as.list(lower), as.list(upper))
#'
#'       opt <- juliaGet(julia_optimize(df, dfc, as.list(Theta0),
#'                                      julia_IPNewton(),
#'                                      julia_Options(iterations=as.integer(maxit),
#'                                                    store_trace = store_trace,
#'                                                    show_trace = show_trace,
#'                                                    time_limit = time_lim,
#'                                                    x_tol = xtol,
#'                                                    f_tol = ftol,
#'                                                    g_tol = gtol,
#'                                                    allow_f_increases = T, successive_f_tol = as.integer(2))))
#'
#'       x_converged[infer] <- opt$x_converged
#'       x_abschange[infer] <- opt$x_abschange
#'       f_converged[infer] <- opt$f_converged
#'       f_abschange[infer] <- opt$f_abschange
#'       f_relchange[infer] <- opt$f_relchange
#'       f_increased[infer] <- opt$f_increased
#'       g_converged[infer] <- opt$g_converged
#'       g_residual[infer] <- opt$g_residual
#'       g_calls[infer] <- opt$g_calls
#'       h_calls[infer] <- opt$h_calls
#'       if (store_trace){
#'         trace[infer] <- opt$trace
#'       }
#'
#'       initial_x[infer,] <- opt$initial_x
#'       time_run[infer] <- opt$time_run
#'       converged[infer] <- opt$x_converged | opt$f_converged | opt$g_converged
#'     }
#'     minimizer[infer,] <- opt$minimizer
#'     minimum[infer] <- opt$minimum
#'     iterations[infer] <- opt$iterations
#'     iteration_converged[infer] <- opt$iteration_converged
#'     f_calls[infer] <- opt$f_calls
#'     if (store_quantiles){
#'       all_quantiles[infer,,] <- quantiles
#'     }
#'   }
#'   optim_result <- opt
#'   #Change non constant values
#'   optim_result$minimizer <- minimizer
#'   optim_result$minimum <- minimum
#'   optim_result$iterations <- iterations
#'   optim_result$iteration_converged <- iteration_converged
#'   optim_result$f_calls <- f_calls
#'   optim_result$mod <- "R-Julia"
#'   optim_result$time_run <- time_run
#'   if (method == "Brent"){
#'     optim_result$converged <- converged
#'   }
#'   else{
#'     optim_result$initial_x <- initial_x
#'     optim_result$converged <- converged
#'     optim_result$x_converged <- x_converged
#'     optim_result$x_abschange <- x_abschange
#'     optim_result$f_converged <- f_converged
#'     optim_result$f_abschange <- f_abschange
#'     optim_result$f_relchange <- f_relchange
#'     optim_result$f_increased <- f_increased
#'     optim_result$g_converged <- g_converged
#'     optim_result$g_residual <- g_residual
#'     optim_result$g_calls <- g_calls
#'     optim_result$h_calls <- h_calls
#'   }
#'   optim_result$ls_success <- NULL
#'   optim_result$AD <- F
#'   names(optim_result)[names(optim_result) == "iteration_converged"] <- "iteration_limit_reached"
#'   if (store_quantiles){
#'     optim_result$quantiles <- all_quantiles
#'   }
#'   if (store_trace){
#'     optim_result$trace <- trace
#'   }
#'
#'   class(optim_result) <- "flimo_result"
#'
#'   optim_result
#' }





# else if (mod == "R-Julia"){
#   flimoptim_RJulia(data, ndraw, sumstats, simulatorQ,
#                    obj = obj,
#                    nsim = nsim,
#                    ninfer = ninfer,
#                    lower = lower,
#                    upper = upper,
#                    Theta0 = Theta0,
#                    randomTheta0 = randomTheta0,
#                    method = method,
#                    maxit = maxit,
#                    time_lim = time_lim,
#                    xtol = xtol,
#                    ftol = ftol,
#                    gtol = gtol,
#                    reltol = reltol,
#                    abstol = abstol,
#                    show_trace = show_trace,
#                    store_trace = store_trace,
#                    store_quantiles = store_quantiles)
# }





#optim_gkR10$minimizer

# system.time(optim_gkRJ500 <- flimoptim(NULL, 8, NULL, NULL,
#                 obj = J_order_gk,
#                  ninfer = 1,
#                 mod = "R-Julia",
#                  nsim = 500, lower = rep(0.01, 4), upper = rep(10, 4)))

# #optim_gkRJ500$minimizer
#
# system.time(optim_gkRJ10 <- flimoptim(NULL, 8, NULL, NULL,
#                                       obj = J_order_gk,
#                                       ninfer = 10,
#                                       mod = "R-Julia",
#                                       nsim = 10, lower = rep(0.01, 4), upper = rep(10, 4)))
#
# #optim_gkRJ10$minimizer
#
#
# #Second mod : R-Julia
#
# system.time(optim1RJ <- flimoptim(data1, ndraw1, sumstats1, simulatorQ1N,
#                                   ninfer = 10, randomTheta0 = T, nsim = nsim1, lower = 1, upper = 1000, mod = "R-Julia"))
# optim1RJ
# summary(optim1RJ)
# attributes(optim1RJ)
#
# plot(optim1RJ)
#
#
# system.time(optim1RJBrent <- flimoptim(data1, ndraw1, sumstats1, simulatorQ1N,
#                                        ninfer = 10, randomTheta0 = T, nsim = nsim1, lower = 1, upper = 1000, mod = "R-Julia", method = "Brent"))
# optim1RJBrent
# summary(optim1RJBrent)
# attributes(optim1RJBrent)
#
# plot(optim1RJBrent)
#
# #Second mod : R-Julia
#
# optim2RJ <- flimoptim(data2, ndraw2, sumstats2, simulatorQ2,
#                       ninfer = 10, nsim = nsim2,
#                       lower = c(-5, 0), upper = c(10, 10),
#                       randomTheta0 = T,
#                       mod = "R-Julia")
# optim2RJ
# summary(optim2RJ)
# plot(optim2RJ, pairwise_par = T, hist = T, par_minimum = T)
