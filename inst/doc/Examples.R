## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(flimo)
library(JuliaConnectoR)
library(ggplot2)
library(gridExtra)
theme_set(theme_bw())

## ----first_run----------------------------------------------------------------
#julia_setup()

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

#just one simulation
simulatorQ1(Theta11, quantiles1[1,])
#Same value :
simulatorQ1(Theta11, quantiles1[1,]) 

#independent values :
simulatorQ1(Theta11, quantiles1[2,])


#Matrix of nsim simulations

simu11 <- simulatorQ1(Theta11, quantiles1)
simu21 <- simulatorQ1(Theta21, quantiles1)

#Sample Comparison : summary statistics

sumstats1 <- function(simulations, data){
  #simulations : 2D array
  #data : 1D array
  mean_simu <- mean(rowMeans(simulations))
  mean_data <- mean(data)
  (mean_simu-mean_data)^2
}


## ----objective_function_1-----------------------------------------------------
#Plot objective function

#Objective by parameter :
plot_objective(ndraw1, nsim1, data1, sumstats1, simulatorQ1, lower = 0, upper = 200)

#Objective with Normal approximation :


plot_objective(ndraw1, nsim1, data1, sumstats1, simulatorQ1N, lower = 0, upper = 200)

#both plots
#We use same quantiles for both

quantiles1 <- matrix(runif(ndraw1*nsim1), nrow = nsim1)

p <- plot_objective(NULL, NULL, data1, sumstats1, simulatorQ1,
                    lower = 0, upper = 200, quantiles = quantiles1)
plot_objective(NULL, NULL,
               data1, sumstats1, simulatorQ1N,
               lower = 0, upper = 200,
               visualize_min = FALSE,
               add_to_plot = p, quantiles = quantiles1)

#Locally, Poisson quantiles and then objective function are constant by pieces
#To compare with normal approximation

p <- plot_objective(NULL, NULL, data1, sumstats1,
                    simulatorQ1, lower = 71, upper = 72,
                    visualize_min = FALSE, quantiles = quantiles1)

plot_objective(ndraw1, nsim1, data1, sumstats1, simulatorQ1N,
               lower = 71, upper = 72, visualize_min = FALSE, add_to_plot = p,
               quantiles = quantiles1)

## ----figure_method------------------------------------------------------------
nsimplot <- 5
quantilesplot <- matrix(runif(ndraw1*nsimplot), nrow = nsimplot)

intern_obj <- function(Theta) flimobjective(Theta, quantilesplot, data1, sumstats1, simulatorQ1)

intern_objN <- function(Theta) flimobjective(Theta, quantilesplot, data1, sumstats1, simulatorQ1N)

intern_objRand <- function(Theta){
  sim <- rpois(nsim1, Theta)
  (mean(sim)-mean(data1))^2
}

x <- c(seq(0, 200, length.out = 1e4), seq(70.5, 72.5, length.out = 1e4))
y <- sapply(x, intern_obj)
yN <- sapply(x, intern_objN)

xR <- seq(0, 200, length.out = 2e3)
yR <- sapply(xR, intern_objRand)

df <- data.frame(x = x, y = y, Method = "Fixed Landscape")
df <- rbind(df, data.frame(x = x, y = yR, Method = "Naive computation"))
df <- rbind(df, data.frame(x = x, y = yN, Method = "Fixed Landscape\nwith Normal Approximation"))

ggplot()+
  geom_line(aes(xR,yR), color = "grey")+
  geom_line(aes(x,y), color = "blue")+
  ggtitle("Inference for a Poisson distribution : value of objective")+
  labs(x = "Theta", y = "J(Theta)")+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))

# ggsave("../../manuscript/Figures/poisson.png", dpi = 350, units = "mm", width = 180, height = 150)

x1 <- 71
x2 <- x1+1
xin <- x[which(x<=x1)[length(which(x<=x1))]]
yax <- max(y[which(x<=x1)[length(which(x<=x1))]], yN[which(x<=x1)[length(which(x<=x1))]])
xax <- x[which(x>=x2)[1]]
yin <- min(y[which(x>=x2)[1]], yN[which(x>=x2)[1]])

ggplot()+
  geom_line(aes(x,y), color = "blue")+
  geom_line(aes(x,yN), color = "red")+
  ggtitle("Inference for a Poisson distribution : value of objective")+
  labs(x = "Theta", y = "J(Theta)")+
  theme(legend.position='none',
        plot.title = element_text(hjust = 0.5))+
  coord_cartesian(xlim = c(xin, xax), ylim = c(yin, yax))
# ggsave("../../manuscript/Figures/poisson_zoom.png", dpi = 350, units = "mm", width = 180, height = 150)


## ----optimization_1-----------------------------------------------------------
#Optimization with normal approximation of Poisson distribution

#First mode : full R
#default mode

system.time(optim1R <- flimoptim(data1, ndraw1, sumstats1, simulatorQ1N,
                 ninfer = 10,
                 nsim = nsim1,
                 lower = 1, upper = 1000,
                 randomTheta0 = TRUE))
optim1R
summary(optim1R)
attributes(optim1R)

plot(optim1R)

#Second mode : full J
#Optimization with Automatic Differentiation
#Warning : you need to translate sumstats and simulatorQ to Julia
#Both of these names for the Julia functions are mandatory !

julia_simulatorQ1N <-"
function simulatorQ(Theta, quantiles)
  quantile.(Normal.(Theta, sqrt.(Theta)), quantiles)
end
"

julia_sumstats1 <-"
function sumstats(simulations, data)
  (mean(mean(simulations, dims = 2))-mean(data))^2
end
"

#Most accurate config for complex problems : IPNewton with AD

if (juliaSetupOk()){
  system.time(optim1JAD <- flimoptim(data1, ndraw1, julia_sumstats1, julia_simulatorQ1N,
                 ninfer = 10,
                 nsim = nsim1,
                 lower = 0,
                 upper = 1000,
                 randomTheta0 = TRUE,
                 mode = "Julia",
                 load_julia = TRUE))
  optim1JAD
  summary(optim1JAD)
  attributes(optim1JAD)

  plot(optim1JAD)

#IPNewton without AD
  system.time(optim1J <- flimoptim(data1, ndraw1, julia_sumstats1, julia_simulatorQ1N,
                 ninfer = 10, nsim = nsim1, lower = 0, upper = 1000,
                 randomTheta0 = TRUE, mode = "Julia", AD = FALSE))
  optim1J
  summary(optim1J)
  attributes(optim1J)

  plot(optim1J)

#Brent
system.time(
  optim1JBrent <- flimoptim(data1,
                            ndraw1,
                            julia_sumstats1,
                            julia_simulatorQ1N,
                            ninfer = 10,
                            nsim = nsim1,
                            lower = 0,
                            upper = 1000,
                            randomTheta0 = TRUE,
                            mode = "Julia",
                            method = "Brent")
)
  optim1JBrent
  summary(optim1JBrent)
  attributes(optim1JBrent)

  plot(optim1JBrent)
}

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


sumstats2 <-function(simulations, data){
  mean_simu <- mean(rowMeans(simulations))
  mean_data <- mean(data)
  sd_simu <-mean(apply(simulations, 1, sd))
  sd_data <- sd(data)
  (mean_simu-mean_data)^2+(sd_simu-sd_data)^2
}

nsim2 <- 10

## ----objective_2--------------------------------------------------------------

plot_objective(ndraw2, nsim2, data2, sumstats2, simulatorQ2, index = 1, other_param = c(1, 2, 10),
                           lower = -5, upper = 10)



## ----optimization_2-----------------------------------------------------------
#Optimization


#First mode : full R
#default mode

optim2R <- flimoptim(data2, ndraw2, sumstats2, simulatorQ2,
                 ninfer = 10, nsim = nsim2,
                 lower = c(-5, 0), upper = c(10, 10),
                 randomTheta0 = TRUE)

optim2R
summary(optim2R)
plot(optim2R, pairwise_par = TRUE, hist = TRUE, par_minimum = TRUE)

#Second mode : full Julia
#Optimization with Automatic Differentiation
#Warning : you need to translate sumstats and simulatorQ to Julia
#Both of these names for the Julia functions are mandatory !

if (juliaSetupOk()){
  julia_simulatorQ2 <-"
function simulatorQ(Theta, quantiles)
  quantile.(Normal.(Theta[1], Theta[2]), quantiles)
end
"

  julia_sumstats2 <-"
 function sumstats(simulations, data)
  (mean(mean(simulations, dims = 2))-mean(data))^2+
  (mean(std(simulations, dims = 2))-std(data))^2
end
"


#Most accurate config for complex problems : IPNewton with AD
optim2JAD <- flimoptim(data2, ndraw2, julia_sumstats2, julia_simulatorQ2,
                 ninfer = 10, nsim = nsim2,
                 lower = c(-5, 0), upper = c(10, 10),
                 randomTheta0 = TRUE,
                 mode = "Julia")
optim2JAD
summary(optim2JAD)

plot(optim2JAD)

#IPNewton without AD
optim2J <- flimoptim(data2, ndraw2, julia_sumstats2, julia_simulatorQ2,
                 ninfer = 10, nsim = nsim2,
                 lower = c(-5, 0), upper = c(10, 10),
                 randomTheta0 = TRUE,
                 mode = "Julia")
  optim2J
  summary(optim2J)
  plot(optim2J)
}



