## ----define variables---------------------------------------------------------
library(denim)
library(deSolve)

# --- Transition def for denim
transitions <- list(
  "S -> I" = d_exponential(0.2),
  "I -> R" = d_gamma(3, 2)
)
parameters <- c(rate = 0.2, scale = 3, shape=2) 
initialValues <- c(S = 999, I = 1, I1 = 1, I2=0, R=0)

# --- Transition def for deSolve
transition_func <- function(t, state, param){
  with(as.list( c(state, param) ), {
      gamma_rate = 1/scale
      dS = -rate*S
      # apply linear chain trick
      dI1 = rate*S - gamma_rate*I1
      dI2 = gamma_rate*I1 - gamma_rate*I2
      dI =  dI1 + dI2
      dR = gamma_rate*I2
      list(c(dS, dI, dI1, dI2, dR))
  })
}

# --- Timestep definition
simulationDuration <- 20 
timestep <- 0.001 # small timestep required for comparison

## -----------------------------------------------------------------------------
denim_start <- Sys.time()
mod <- sim(transitions = transitions, initialValues = initialValues, parameters, simulationDuration = simulationDuration, timeStep = timestep)
denim_end <- Sys.time()

# --- show output
head(mod[mod$Time %in% 1:simulationDuration,])

## -----------------------------------------------------------------------------
times <- seq(0, simulationDuration, timestep)

desolve_start <- Sys.time()
ode_mod <- ode(y = initialValues, times = times, parms = parameters, func = transition_func)
desolve_end <- Sys.time()

# --- show output
ode_mod <- as.data.frame(ode_mod)
head(ode_mod[ode_mod$time %in% 1:simulationDuration, c("time", "S", "I", "R")])

## ----include=FALSE------------------------------------------------------------
time_diff <- round(as.numeric(difftime(denim_end, denim_start, units="sec")) /
                     as.numeric(difftime(desolve_end, desolve_start, units="sec")), digits = 2)

## -----------------------------------------------------------------------------
# increase timestep before plotting
mod <- mod[mod$Time %in% seq(0, simulationDuration, 0.2),]
ode_mod <- ode_mod[ode_mod$time %in% seq(0, simulationDuration, 0.2),]

## ----fig.width=8, fig.height=5------------------------------------------------
# ---- Plot S compartment
plot(x = mod$Time, y = mod$S,xlab = "Time", ylab = "Count", main="S compartment",
     col = "#4876ff", type="l", lwd=3)
lines(ode_mod$time, ode_mod$S, lwd=3, lty=3)
legend(x = 15, y = 900,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot I compartment
plot(x = mod$Time, y = mod$I, xlab = "Time", ylab = "Count", main="I compartment",
      col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$I, lwd=3, lty=3)
legend(x = 15, y = 350,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot R compartment
plot(x = mod$Time, y = mod$R, xlab = "Time", ylab = "Count", main="R compartment",
     col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$R, lwd=3, lty=3)
legend(x = 15, y = 300,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))


