## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(denim)

## -----------------------------------------------------------------------------
sir_parametric <- denim_dsl({
  S -> I = beta * (I/N) * S * timeStep
  I -> R = d_weibull(scale = r_scale, shape = r_shape)
})

## -----------------------------------------------------------------------------
sir_nonparametric <- denim_dsl({
  S -> I = beta * (I/N) * S * timeStep
  I -> R = nonparametric(dwelltime_dist)
})

## -----------------------------------------------------------------------------
# parameters
mod_params <- list(
  beta = 0.4,
  N = 1000,
  r_scale = 4,
  r_shape = 3
)
# initial population
init_vals <- c(S = 950, I = 50, R = 0)
# simulation duration and timestep
sim_duration <- 30
timestep <- 0.05

## -----------------------------------------------------------------------------
parametric_mod <- sim(sir_parametric,
    initialValues = init_vals,
    parameters = mod_params,
    simulationDuration = sim_duration,
    timeStep = timestep) 

plot(parametric_mod, ylim = c(0, 1000))

## -----------------------------------------------------------------------------
# Compute discrete distribution of dwell-tinme
# dist_func - R distribution function for dwell time (pexp, pgamma, etc.)
# ... - parameters for dist_func
compute_dist <- function(dist_func,..., timestep=0.05, error_tolerance=0.0001){
  maxtime <- timestep
  prev_prob <- 0
  prob_dist <- numeric()
  
  while(TRUE){
     # get current cumulative prob and check whether it is sufficiently close to 1
     temp_prob <-  ifelse(
       dist_func(maxtime, ...) < (1 - error_tolerance), 
       dist_func(maxtime, ...), 
       1);

     # get f(t)
     curr_prob <- temp_prob - prev_prob
     prob_dist <- c(prob_dist, curr_prob)
     
     prev_prob <- temp_prob
     maxtime <- maxtime + timestep
     
     if(temp_prob == 1){
       break
     }
  }
  
  prob_dist
}

## -----------------------------------------------------------------------------
# Compute the discrete distribution
dwelltime_dist <- compute_dist(pweibull, 
                               scale = mod_params$r_scale, shape = mod_params$r_shape,
                               timestep = timestep)

# Compute the discrete distribution
nonparametric_mod <- sim(sir_nonparametric,
    initialValues = init_vals,
    parameters = list(
      beta = mod_params$beta,
      N = mod_params$N,
      dwelltime_dist = dwelltime_dist
    ),
    simulationDuration = sim_duration,
    timeStep = timestep) 
plot(nonparametric_mod, ylim = c(0, 1000))

## ----echo=FALSE---------------------------------------------------------------
first_dist <- compute_dist(pweibull, 
                   scale = 1.5, shape = 4,
                   timestep = timestep)
second_dist <- compute_dist(pweibull, 
                   scale = 3, shape = 3.5,
                   timestep = timestep)
first_dist <- c(rep(0, length(second_dist) - length(first_dist)), first_dist)
multimodal_dist <- first_dist + second_dist

## -----------------------------------------------------------------------------
timestep <- 0.05
plot(seq(0, by = 0.05, length.out = length(multimodal_dist)), 
     multimodal_dist, 
     type = "l", col = "#374F77", lty = 1, lwd = 3,
     xlab = "Length of stay (days)", ylab = "", yaxt = 'n')

## -----------------------------------------------------------------------------
# model parameter
parameters <- list(beta = 0.4, N = 1000,
                   dwelltime_dist = multimodal_dist)
# initial population
init_vals <- c(S = 950, I = 50, R = 0)
# simulation duration and timestep
sim_duration <- 30
timestep <- 0.05

# Run the model with multimodel distribution
nonparametric_mod <- sim(
  sir_nonparametric,
  initialValues = init_vals,
  parameters = parameters,
  simulationDuration = sim_duration,
  timeStep = timestep) 
plot(nonparametric_mod, ylim = c(0, 1000))

