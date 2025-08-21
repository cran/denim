## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(denim)

## -----------------------------------------------------------------------------
sir_model <- denim_dsl({
  S -> I = beta * (I/N) * S
  I -> R = d_exponential(rate = gamma)
})

## -----------------------------------------------------------------------------
sir_model

## ----eval=FALSE---------------------------------------------------------------
# sir_model <- denim_dsl({
#   S -> I = beta*(I/N)*S
#   I -> R = d_exponential(rate = 1/4)
# })

## ----eval=FALSE---------------------------------------------------------------
# sir_model <- denim_dsl({
#   # this is a comment
#   S -> I = beta*(I/N)*S
#   I -> R = d_exponential(rate = 1/4) # this is another comment
# })

## -----------------------------------------------------------------------------
# parameters for the model
parameters <- c(
  beta = 0.4,
  N = 1000,
  gamma = 1/7
)
# initial population for each compartment 
initValues <- c(
  S = 999, 
  I = 50,
  R = 0
)

## -----------------------------------------------------------------------------
mod <- sim(sir_model, 
    parameters = parameters, 
    initialValues = initValues, 
    timeStep = 0.01,
    simulationDuration = 40)

## -----------------------------------------------------------------------------
plot(mod, ylim = c(1, 1000))

## -----------------------------------------------------------------------------
time_varying_mod <- denim_dsl({
  A -> B = 20 * (1+cos(omega * time)) 
})

# parameters for the model
parameters <- c(
  omega = 2*pi/10
)
# initial population for each compartment 
initValues <- c(A = 1000, B = 0)

mod <- sim(time_varying_mod, 
    parameters = parameters, 
    initialValues = initValues, 
    timeStep = 0.01,
    simulationDuration = 40)

plot(mod, ylim = c(0, 1000))

## -----------------------------------------------------------------------------
sir_model_list <- list(
  "S -> I" = "beta * (I/N) * S",
  "I -> R" = d_exponential(rate = "gamma")
)

sir_model_list

## -----------------------------------------------------------------------------
# parameters for the model
parameters <- c(
  beta = 0.4,
  N = 1000,
  gamma = 1/7
)
# initial population for each compartment 
initValues <- c(
  S = 999, 
  I = 50,
  R = 0
)
# run the simulation
mod <- sim(sir_model_list, 
    parameters = parameters, 
    initialValues = initValues, 
    timeStep = 0.01,
    simulationDuration = 40)
# plot output
plot(mod, ylim = c(1, 1000))

## ----message=FALSE, warning=FALSE---------------------------------------------
library(tidyverse)
# configurations for 3 different I->R transitions
model_config <- tibble(
  IR_dists = c(d_gamma, d_weibull, d_lognormal),
  IR_pars = list(c(rate = 0.1, shape = 3), c(scale = 5, shape = 0.3), c(mu = 0.3, sigma = 2))
)

walk2(
  model_config$IR_dists, model_config$IR_pars, \(dist, par){
    transitions <- list(
      "S -> I" = "beta * S * (I / N)",
      # This is not applicable when using denim_dsl()
      "I -> R" = do.call(dist, as.list(par))
    )
    
    # model settings
    denimInitialValues <- c(S = 980, I = 20, R = 0)
    parameters <- c(
      beta = 0.4,
      N = 1000
    )
    
    # compare output 
    mod <- sim(transitions = transitions, 
               initialValues = denimInitialValues, 
               parameters = parameters, 
               simulationDuration = 60, 
               timeStep = 0.05)
    
    plot(mod, ylim = c(0,1000))
})

