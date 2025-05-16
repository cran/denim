## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center"
)
library(DiagrammeR) # for flowchart diagram
library(denim)

## ----echo=FALSE---------------------------------------------------------------
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle]
  
  S -> I [label = '&#946;SI/N']
  I -> R [label = '&#947;I']
  }",
  width = 300, height = "100%")

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

## ----echo=FALSE---------------------------------------------------------------
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle]
  
  S -> I [label = '&#946;SI/N']
  I -> R [label = 'd_gamma(1/2, 3)']
  }",
  width = 400, height = "100%")

## ----echo=FALSE---------------------------------------------------------------
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle]
  
  S -> I [label = '&#946;SI/N']
  I -> R [label = 'd_gamma(1/2, 3)']
  }",
  width = 400, height = "100%")

## -----------------------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * I / N",
  "I -> R" = d_gamma(rate = 1/2, 3)
)

## -----------------------------------------------------------------------------
initialValues <- c(
  S = 999, 
  I = 1, 
  R = 0
)

## -----------------------------------------------------------------------------
parameters <- c(
  beta = 0.012,
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
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle]
  
  S -> I [label = '&#946;SI/N']
  I -> R [label = '&#947;I']
  }",
  width = 300, height = "100%")

## ----echo=FALSE---------------------------------------------------------------
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle]
  
  S -> I1 [label = '&#946;S(I@_{1} + I@_{2} + I@_{3} + I@_{4})/N']
  I1 -> R1 [label = '&#947;@_{1}I@_{1}']
  I2 -> R2 [label = '&#947;@_{2}I@_{2}']
  I3 -> R3 [label = '&#947;@_{3}I@_{3}']
  I4 -> R4 [label = 'I@_{4}']
  I1 -> I2 [label = '(1-&#947;@_{1})I@_{1}']
  I2 -> I3 [label = '(1-&#947;@_{2})I@_{2}']
  I3 -> I4 [label = '(1-&#947;@_{3})I@_{3}']
  }", height = "100%", width = "100%")

## ----eval=FALSE---------------------------------------------------------------
# transitions <- list(
#   "S -> I" = "beta * S * I / N",
#   "I -> R" = d_gamma(rate = 1/2, 3, dist_init=TRUE)
# )

## ----echo=FALSE---------------------------------------------------------------
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle]
  
  S -> I [label = '&#946;SI/N']
  S -> V [label = '5']
  I -> R [label = '0.9 -> d_gamma(1/3, 2)']
  I -> D [label = '0.1 -> d_lognormal(2, 0.5)']
  }",
  width = 500, height = "100%")

## ----fig.width = 5------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * I / N",
  "S -> V" = 5,
  "0.9 * I -> R" = d_gamma(1/3, 2),
  "0.1 * I -> D" = d_lognormal(2, 0.5)
)

initialValues <- c(
  S = 999, 
  I = 1, 
  R = 0,
  V = 0,
  D = 0
)

parameters <- c(
  beta = 0.12,
  N = 1000
)

simulationDuration <- 10
timeStep <- 0.01

mod <- sim(transitions = transitions, 
           initialValues = initialValues, 
           parameters = parameters, 
           simulationDuration = simulationDuration, 
           timeStep = timeStep)

head(mod)
plot(mod)

## ----echo=FALSE---------------------------------------------------------------
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = LR]
  
  node [shape = rectangle]
  
  S -> I [label = '&#946;S(I + IV)/N']
  S -> V [label = '2']
  I -> D [label = '0.1 -> d_lognormal(2, 0.5)']
  I -> R [label = '0.9 -> d_gamma(1/3, 2)']
  V -> IV [label = '0.1 * &#946;V(I + IV)/N']
  IV -> R [label = 'd_exponential(2)']
  }",
  width = 700, height = "100%")

## ----fig.width = 5------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * (I + IV) / N",
  "S -> V" = 2,
  "0.1 * I -> D" = d_lognormal(2, 0.5),
  "0.9 * I -> R" = d_gamma(rate = 1/3, shape = 2),
  "V -> IV" = "0.1 * beta * V * (I + IV) / N",
  "IV -> R" = d_exponential(2)
)

initialValues <- c(
  S = 999, 
  I = 1, 
  R = 0,
  V = 0,
  IV = 0,
  D = 0
)

parameters <- c(
  beta = 0.12,
  N = 1000
)

simulationDuration <- 10
timeStep <- 0.01

mod <- sim(transitions = transitions, 
           initialValues = initialValues, 
           parameters = parameters, 
           simulationDuration = simulationDuration, 
           timeStep = timeStep)
plot(mod)
mod

