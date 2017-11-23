library(deSolve)
library(lme4)
library(dplyr)
#ampicillin concentration and position
outt<-NULL
ampconc <- seq(0, 1, by = 0.01)
ampco <- rev((exp(ampconc) / 2.72))
getmax = for(a in ampco) {
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
      method = "lsoda"
    ))
  maxt<-c(a,max(out$tet))
  outt<-as.data.frame(rbind(outt,maxt))
}
names(outt)<-c("a","t")
survive<-filter(outt, outt$a < 0.8, outt$t<1)
#results <- as.data.frame(ampco, lapply(ampco, getmax))
#head(results)
#plot(x = results$amp, y = results$tet)
#curve(getmax, from = 1, to = 2, n = 101)
plot(outt)
points(survive, col="green",type="l",lwd=10)