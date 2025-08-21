## -----------------------------------------------------------------------------
library(denim)
library(deSolve)

## -----------------------------------------------------------------------------
# --- Model definition in deSolve
transition_func <- function(t, state, param){
  with(as.list( c(state, param) ), {
      dS = -beta*S*I/N
      dI1 = beta*S*I/N - rate*I1
      dI2 = rate*I1 - rate*I2
      dI =  dI1 + dI2
      dR = rate*I2
      list(c(dS, dI, dI1, dI2, dR))
  })
}

# ---- Model configuration 
parameters <- c(beta = 0.3, rate = 1/3, N = 1000) 
initialValues <- c(S = 999, I = 1, I1 = 1, I2=0, R=0)

# ---- Run simulation
times <- seq(0, 100) # simulation duration
ode_mod <- ode(y = initialValues, times = times, parms = parameters, func = transition_func)

# --- show output
ode_mod <- as.data.frame(ode_mod)
head(ode_mod[-1, c("time", "S", "I", "R")])

## ----eval=FALSE---------------------------------------------------------------
# # --- Model definition in deSolve
# transition_func <- function(t, state, param){
#   with(as.list( c(state, param) ), {
# 
#       # For S -> I transition, since it involves parameters (beta, N),
#       # the best transition to describe this is using a mathematical formula
#       dS = -beta*S*I/N
# 
#       # For I -> R transition, linear chain trick is applied --> implies Erlang distributed dwell time
#       # Hence, we can use d_gamma from denim
#       dI1 = beta*S*I/N - rate*I1
#       dI2 = rate*I1 - rate*I2
#       dI =  dI1 + dI2
#       dR = rate*I2
#       list(c(dS, dI, dI1, dI2, dR))
#   })
# }

## -----------------------------------------------------------------------------
# --- Model definition in denim
transitions <- denim_dsl({
  S -> I = beta * S * I/N
  # shape is 2 from number of I sub compartments
  I -> R = d_gamma(rate = 1/3, shape = 2) 
})

## -----------------------------------------------------------------------------
# remove I1, I2 compartments
denim_initialValues <- c(S = 999, I = 1, R=0)
denim_parameters <- c(beta = 0.3, N = 1000) 

## ----eval=FALSE---------------------------------------------------------------
# transitions <- denim_dsl({
#   S -> I = beta * S * I/N
#   I -> R = d_gamma(rate = 1/3, shape = 2, dist_init = TRUE)
# })

## -----------------------------------------------------------------------------
mod <- sim(transitions = transitions,
             initialValues = denim_initialValues, 
             parameters = denim_parameters,
             simulationDuration = 100,
             timeStep = 0.01)

head(mod[mod$Time %% 1 == 0, ])

## ----echo=FALSE, fig.width=8, fig.height=5------------------------------------
# ---- Plot S compartment
plot(x = mod$Time, y = mod$S,xlab = "Time", ylab = "Count", main="S compartment",
     col = "#4876ff", type="l", lwd=3)
lines(ode_mod$time, ode_mod$S, lwd=3, lty=3)
legend(x = 15, y = 4e5,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot I compartment
plot(x = mod$Time, y = mod$I, xlab = "Time", ylab = "Count", main="I compartment",
      col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$I, lwd=3, lty=3)
legend(x = 15, y = 1e5,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

# ---- Plot R compartment
plot(x = mod$Time, y = mod$R, xlab = "Time", ylab = "Count", main="R compartment",
     col = "#4876ff", type="l", lwd=2)
lines(ode_mod$time, ode_mod$R, lwd=3, lty=3)
legend(x = 15, y = 4e5,legend=c("denim", "deSolve"), col = c("#4876ff", "black"), lty=c(1,3))

