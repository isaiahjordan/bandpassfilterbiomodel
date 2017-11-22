library(deSolve)
library(lme4)
#ampicillin concentration and position
ampconc <- seq(0, 1, by = 0.01)
ampco <- rev((exp(ampconc) / 2.72))
for (a in ampco) {
  parameters = c(
    ampc = a,
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
  
  tail(out)
  plot(x = out$time,
       y = out$tet)
}
