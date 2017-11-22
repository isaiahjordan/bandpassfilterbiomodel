library(deSolve)
library(lme4)
#ampicillin concentration and position
ampconc <- seq(0, 1, by = 0.01)
ampco <- rev((exp(ampconc) / 2.72))
getmax = function(x) {
  parameters = c(
    ampc = 1,
    trate = 1.5,
    degtc = 1,
    synrate = 1.5
  )
  #ode model of tetC levels and tet level in the cell
  state <- c(tet = 0.5,
             tetc = .025)
  
  tetconc <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dtet = (tet) - (trate * tetc) * (tet) / (0.5 + tet)
      dtetc = synrate * ampc - (degtc * tetc)
      list(c(dtet, dtetc))
    })
  }
  out <-
    as.data.frame(ode(
      y = state,
      times = seq(1, 50),
      func = tetconc,
      parms = parameters,
      #method = "lsoda"
    ))
  #return(c(amp = a, tet = max(out$tet)))
  max(out$tet)
}

#results <- as.data.frame(ampco, lapply(ampco, getmax))
#head(results)
#plot(x = results$amp, y = results$tet)
#curve(getmax, from = 1, to = 2, n = 101)
curve(getmax, -pi, 2*pi)