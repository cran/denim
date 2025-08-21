## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center"
)
library(DiagrammeR) # for flowchart diagram
library(denim)

## ----echo=FALSE---------------------------------------------------------------
rates <- c(0.5, 1, 1.5)
x <- seq(0, 5, 0.001)
y <- dexp(x = x, rate = rates[1])
y2 <- dexp(x = x, rate = rates[2])
y3 <- dexp(x = x, rate = rates[3])

col_codes <- c("#374F77", "#EF9F32", "#6ECBD7")
plot(x, y, type = "l", col = col_codes[1], lty = 1, lwd = 3,
     xlab = "Length of stay (days)", ylab = "", 
     ylim = c(0, 1.5), yaxt = 'n')
lines(x, y2, col = col_codes[2], lty = 1, lwd = 3)
lines(x, y3, col = col_codes[3], lty = 1, lwd = 3)
legend("right", legend = c(0.5, 1.0, 1.5), 
       col = col_codes, lty = 1, lwd = 3, bty = "\n")

## ----echo=FALSE---------------------------------------------------------------
x <- seq(0, 20, 0.001)
y <- dgamma(x = x, shape = 3, rate = 1/2)

plot(x, y, type = "l", col = col_codes[1], lty = 1, lwd = 3,
     xlab = "Length of stay (days)", ylab = "", yaxt = 'n')

## -----------------------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * I / N",
  "I -> R" = d_gamma(rate = 1/2, 3)
)

## -----------------------------------------------------------------------------
transitions <- denim_dsl({
  S -> I = beta * (I/N) * S
  I -> R = d_gamma(rate = 1/2, shape = 3)
})

## -----------------------------------------------------------------------------
initialValues <- c(
  S = 999, 
  I = 1, 
  R = 0
)

## -----------------------------------------------------------------------------
parameters <- c(
  beta = 1.2,
  N = 1000
)

## -----------------------------------------------------------------------------
simulationDuration <- 30
timeStep <- 0.01

## ----fig.width = 6------------------------------------------------------------
mod <- sim(transitions = transitions, 
           initialValues = initialValues, 
           parameters = parameters, 
           simulationDuration = simulationDuration, 
           timeStep = timeStep)
head(mod)
plot(mod)

## ----echo=FALSE---------------------------------------------------------------
# DiagrammeR::grViz("digraph {
#   graph [layout = dot, rankdir = LR]
#   
#   node [shape = rectangle]
#   
#   S -> I1 [label = '&#946;S(I@_{1} + I@_{2} + I@_{3} + I@_{4})/N']
#   I1 -> R1 [label = '&#947;@_{1}I@_{1}']
#   I2 -> R2 [label = '&#947;@_{2}I@_{2}']
#   I3 -> R3 [label = '&#947;@_{3}I@_{3}']
#   I4 -> R4 [label = 'I@_{4}']
#   I1 -> I2 [label = '(1-&#947;@_{1})I@_{1}']
#   I2 -> I3 [label = '(1-&#947;@_{2})I@_{2}']
#   I3 -> I4 [label = '(1-&#947;@_{3})I@_{3}']
#   }", height = "100%", width = "100%")

## ----eval=FALSE---------------------------------------------------------------
# transitions <- denim_dsl({
#   S -> I = beta * S * I / N
#   I -> R = d_gamma(rate = 1/2, 3, dist_init=TRUE)
# })

## ----fig.width = 5------------------------------------------------------------
transitions <- denim_dsl({
  S -> I = beta * S * I / N 
  S -> V = 5
  0.9 * I -> R = d_gamma(1/3, 2)
  0.1 * I -> D = d_lognormal(2, 0.5)
})

initialValues <- c(
  S = 999, 
  I = 1, 
  R = 0,
  V = 0,
  D = 0
)

parameters <- c(
  beta = 1.2,
  N = 1000
)

simulationDuration <- 20
timeStep <- 0.01

mod <- sim(transitions = transitions, 
           initialValues = initialValues, 
           parameters = parameters, 
           simulationDuration = simulationDuration, 
           timeStep = timeStep)

head(mod)
plot(mod, ylim = c(0, 1000))

## ----fig.width = 5------------------------------------------------------------
initialValues <- c(S = 999, I = 1, R = 0, V = 0, IV = 0, D = 0) 
 
modelStructure <- denim_dsl({ 
    S -> I = beta * S * (I + IV) / N 
    S -> V = d_exponential(0.01) 
    0.1 * I -> D = d_lognormal(2, 0.5) 
    0.9 * I -> R = d_gamma(1/3, 2) 
    V -> IV = 0.2 * beta * V * (I + IV) / N 
    IV -> R = nonparametric(iv_r_dist) 
    IV -> D = d_weibull(scale = 2, shape = 1.5) 
}) 

parameters <- list( 
  beta = 0.9,  
  N = 1000,
  iv_r_dist = c(0, 0.15, 0.15,  0.05, 0.2, 0.2, 0.25)
)


simulation <- sim(transitions = modelStructure,  
           initialValues = initialValues,  
           parameters = parameters,  
           simulationDuration = 30,  
           timeStep = 0.5) 

plot(simulation)
head(simulation)

## ----echo=FALSE---------------------------------------------------------------
# ------ Code to generate the distribution to demonstrate nonparametric --------- 
# Helper function
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

timestep <- 0.05
# create a multimodal 
first_dist <- compute_dist(pweibull, 
                   scale = 5, shape = 5,
                   timestep = timestep)
second_dist <- compute_dist(pweibull, 
                   scale = 2.5, shape = 4,
                   timestep = timestep)
second_dist <- c(second_dist, rep(0, length(first_dist) - length(second_dist)))
ir_dist <- first_dist + second_dist

## -----------------------------------------------------------------------------
timestep <- 0.05
plot(seq(0, by = 0.05, length.out = length(ir_dist)), 
     ir_dist, 
     type = "l", col = "#374F77", lty = 1, lwd = 3,
     xlab = "Length of stay (days)", ylab = "", yaxt = 'n')

## -----------------------------------------------------------------------------
transitions <- denim_dsl({
  S -> I = beta*I/N*S
  I -> R = nonparametric(ir_dist)
  I -> D = d_exponential(d_rate)
})

parameters <- list(
  beta = 0.7,
  ir_dist = ir_dist,
  d_rate = 0.1,
  N = 1000
)

initialValues <- c(
  S = 999,
  I = 1,
  R = 0,
  D = 0
)

simulationDuration <- 30
timeStep <- 0.05

mod <- sim(transitions = transitions, 
           initialValues = initialValues, 
           parameters = parameters, 
           simulationDuration = simulationDuration, 
           timeStep = timeStep)

plot(mod, ylim = c(0, 1000))

