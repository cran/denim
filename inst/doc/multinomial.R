## ----message=FALSE------------------------------------------------------------
library(deSolve)
library(denim)

## -----------------------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * (I / N) * timeStep",
  "0.9 * I -> R" = d_gamma(1/3, 2),
  "0.1 * I -> D" = d_exponential(0.1)
)

## -----------------------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * (I / N) * timeStep",
  "36 * I -> R" = d_gamma(1/3, 2),
  "4 * I -> D" = d_exponential(0.1)
)

## -----------------------------------------------------------------------------
# model in denim
transitions <- list(
  "S -> I" = "beta * S * (I / N) * timeStep",
  "0.9 * I -> R" = d_gamma(1/3, 2),
  "0.1 * I -> D" = d_exponential(0.1)
)

denimInitialValues <- c(
  S = 999, 
  I = 1, 
  R = 0,
  D = 0
)

## -----------------------------------------------------------------------------
# model in deSolve
transition_func <- function(t, state, param){
  with(as.list( c(state, param) ), {

      dS = - beta * S * (IR1 + IR2 + ID)/N
      # apply linear chain trick for I -> R transition
      # 0.9 * to specify prop of I that goes to I->R transition
      dIR1 = 0.9 * beta * S * (IR1 + IR2 + ID)/N - rate*IR1
      dIR2 = rate*IR1 - rate*IR2
      dR = rate*IR2 
      
      # handle I -> D transition
      # 0.1 * to specify prop of I that goes to I->D transition
      dID = 0.1 * beta * S * (IR1 + IR2 + ID)/N - exp_rate*ID
      dD = exp_rate*ID
      list(c(dS, dIR1, dIR2, dID, dR, dD))
  })
}

desolveInitialValues <- c(
  S = 999, 
  # note that internally, denim also allocate initial value based on specified proportion
  IR1 = 0.9,
  IR2 = 0,
  ID = 0.1, 
  R = 0,
  D = 0
)

## -----------------------------------------------------------------------------
parameters <- c(
  beta = 0.2,
  N = 1000,
  rate = 1/3,
  exp_rate = 0.1
)

simulationDuration <- 200
timeStep <- 0.05

## -----------------------------------------------------------------------------
# --- run denim model ---- 
mod <- sim(transitions = transitions, 
           initialValues = denimInitialValues, 
           parameters = parameters, 
           simulationDuration = simulationDuration, 
           timeStep = timeStep)

# run deSolve model
times <- seq(0, simulationDuration)
ode_mod <- ode(y = desolveInitialValues, times = times, parms = parameters, func = transition_func) 
ode_mod <- as.data.frame(ode_mod)
ode_mod$I<- rowSums(ode_mod[, c("IR1", "IR2", "ID")])

## ----echo=FALSE---------------------------------------------------------------
# ---- Plot S compartment
plot(x = mod$Time, y = mod$S,xlab = "Time", ylab = "Count", main="S compartment",
     col = "#4876ff", type="l", lwd=3)
lines(ode_mod$time, ode_mod$S, lwd=3, lty=3)
legend(x = 150, y = 850,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot I compartment
plot(x = mod$Time, y = mod$I, xlab = "Time", ylab = "Count", main="I compartment",
      col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$I, lwd=3, lty=3)
legend(x = 150, y = 25,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot R compartment
plot(x = mod$Time, y = mod$R, xlab = "Time", ylab = "Count", main="R compartment",
      col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$R, lwd=3, lty=3)
legend(x = 150, y = 250,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot D compartment
plot(x = mod$Time, y = mod$D, xlab = "Time", ylab = "Count", main="D compartment",
     col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$D, lwd=3, lty=3)
legend(x = 150, y = 25,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

## -----------------------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * (I / N) * timeStep",
  "I -> R" = d_gamma(rate = 1/3, shape = 2),
  "I -> D" = d_gamma(rate = 1/4, shape = 2)
)

## -----------------------------------------------------------------------------
transitions <- list(
  "S -> I" = "beta * S * (I / N) * timeStep",
  "I -> R" = d_gamma(rate = 1/3, shape = 2),
  "I -> D" = d_gamma(rate = 1/4, shape = 2)
)
denimInitialValues <- c(
  S = 950, 
  I = 50, 
  R = 0, 
  D = 0
)

## -----------------------------------------------------------------------------
transition_func <- function(t, state, param){
  with(as.list( c(state, param) ), {
    dS = - beta * S * (I1 + I2 + IR + ID)/N

    dI1 = beta * S * (I1 + I2 + IR + ID)/N - (rate+d_rate)*I1

    dIR = rate*I1 - (d_rate + rate)*IR
    dID = d_rate*I1 - (d_rate + rate)*ID

    dI2 = d_rate*IR + rate*ID - (d_rate + rate)*I2

    dR = rate*IR + rate*I2
    dD = d_rate*ID + d_rate*I2

    list(c(dS, dI1, dIR, dID, dI2, dR, dD))
  })
}

desolveInitialValues <- c(
  S = 950, 
  I1 = 50,
  IR = 0,
  ID = 0,
  I2 = 0,
  R = 0,
  D = 0
)

## -----------------------------------------------------------------------------
parameters <- c(
  beta = 0.2,
  N = 1000,
  rate = 1/3,
  d_rate = 1/4
)

simulationDuration <- 50
timeStep <- 0.05

## -----------------------------------------------------------------------------
# run denim model
mod <- sim(transitions = transitions, 
           initialValues = denimInitialValues, 
           parameters = parameters, 
           simulationDuration = simulationDuration, 
           timeStep = timeStep)

# run deSolve model
times <- seq(0, simulationDuration)
ode_mod <- ode(y = desolveInitialValues, times = times, parms = parameters, func = transition_func) 
ode_mod <- as.data.frame(ode_mod)
ode_mod$I <- rowSums(ode_mod[,c("I1", "ID", "IR", "I2")])

## ----echo=FALSE---------------------------------------------------------------
# ---- Plot S compartment
plot(x = mod$Time, y = mod$S,xlab = "Time", ylab = "Count", main="S compartment",
     col = "#4876ff", type="l", lwd=3)
lines(ode_mod$time, ode_mod$S, lwd=3, lty=3)
legend(x = 30, y = 925,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot I compartment
plot(x = mod$Time, y = mod$I, xlab = "Time", ylab = "Count", main="I compartment",
      col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$I, lwd=3, lty=3)
legend(x = 30, y = 50,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot R compartment
plot(x = mod$Time, y = mod$R, xlab = "Time", ylab = "Count", main="R compartment",
      col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$R, lwd=3, lty=3)
legend(x = 30, y = 80,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot D compartment
plot(x = mod$Time, y = mod$D, xlab = "Time", ylab = "Count", main="D compartment",
     col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$D, lwd=3, lty=3)
legend(x = 30, y = 60,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

