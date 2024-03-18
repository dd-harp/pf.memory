## -----------------------------------------------------------------------------
library(ramp.pf.infection)
library(deSolve)
library(knitr)
library(viridisLite)
library(viridis)

## ----fig.height=4, fig.width=7------------------------------------------------
m = seq(0.1, 24)
N = pmax(24, 4*m)
plot(m, ppois(N,m), type = "l", ylim = c(1-1e-07,1+1e-07))

## ----fig.height=4, fig.width=7, echo=F----------------------------------------
foiP3 = list(hbar = 5/365, 
             agePar = par_type2Age(), 
             seasonPar = par_sinSeason(), 
             trendPar = par_flatTrend())

a2years  = 0:730

plot(a2years, FoI(a2years, foiP3), type = "l", 
     xlab = "a - age (in days)", ylab = expression(FoI[tau](a)))

## -----------------------------------------------------------------------------
solveMMinfty(3/365, foiP3) -> out

## ----fig.height=4, fig.width=7, echo=F----------------------------------------
MoIDistPlot(out$moi, 500)

## ----fig.height=4, fig.width=7------------------------------------------------
solve_dm(3/365, foiP3) -> mt 
plot(mt$time, mt$m, type = "l", 
     xlab = "Age (in Days)", ylab = expression(m[tau](a)))

## ----compute 3 ways-----------------------------------------------------------
aa = 0:1825
v1 = meanMoI(aa, foiP3, hhat=5/365)
v2 = solveMMinfty(5/365, foiP3, Amax=1825)$m
v3 = solve_dm(5/365, foiP3, Amax=1825)[,2]

## ----fig.height=7, fig.width=8------------------------------------------------
par(mfrow = c(2,2))

plot(aa, v1, type = "l", col = "red", lwd=4, 
     xlab = "Age (in Days)", 
     ylab = expression(m[tau](a)), 
     main = "meanMoI")

plot(aa, v2, type = "l", lwd=2, 
     xlab = "Age (in Days)", 
     ylab = expression(m[tau](a)), 
     main = "solveMMinfty")
lines(aa, v2, type = "l", col = "yellow", lwd=2, lty=2)

plot(aa, v3, type = "l", col = "darkblue", lwd=1, 
     xlab = "Age (in Days)", 
     ylab = expression(m[tau](a)), 
     main = "solve_dm")

plot(aa, v1, type = "l", col = "red", lwd=4, 
     xlab = "Age (in Days)", 
     ylab = expression(m[tau](a)), 
     main = "All Three")
lines(aa, v2, type = "l", col = "yellow", lwd = 2)
lines(aa, v3, type = "l", col = "darkblue", lwd = 2, lty =2)

## ----fig.height=4, fig.width=6, echo=F----------------------------------------
foiP3 = list(hbar=1, agePar = par_type2Age(), seasonPar = par_sinSeason(), trendPar = par_flatTrend()) 
tm = 0:1825
Xt= solve_dpda(5/365, foiP3, Amax=5*365)
tPR = truePR(tm, foiP3, hhat=5/365)
plot(tm, tPR, type = "l", ylab = expression(p(a)), xlab = "Age (Days)", lwd=5, col = "red")
lines(Xt$time, 1-exp(-Xt$m), col = "white", lwd=2)
lines(Xt$time, Xt$p, col = "blue", lwd=2, lty=2)

## ----fig.height=5, fig.width=7------------------------------------------------
par(mar = c(7, 4, 2, 2))
solve_dAoI(5/365, foiP3) -> mt 
plot(mt$time, mt$xn1, type = "l", xlab = "Age (in Days)", ylab = expression(list(x, sqrt(x[paste("[2]")]), sqrt(x[paste("[3]")], 3))), ylim = range(mt$xn3^(1/3))) 
lines(mt$time, mt$xn2^(1/2), col = "purple")
lines(mt$time, mt$xn3^(1/3), col = "darkgreen") 

## -----------------------------------------------------------------------------
aa = seq(5, 5*365, by = 5) 
moment1 = momentAoI(aa, foiP3)
moment2 = momentAoI(aa, foiP3, n=2)
moment3 = momentAoI(aa, foiP3, n=3)

## -----------------------------------------------------------------------------
solve_dAoI(5/365, foiP3, Amax = 5*365, dt=5) -> mt 

## ----fig.height=7, fig.width=7------------------------------------------------
par(mfrow = c(2,2), mar = c(7, 4, 2, 2))
# Top Left
plot(mt$time, mt$xn1, type = "l", xlab = "Age (in Days)", ylab = expression(x), lwd=3) 
lines(aa, moment1, col = "yellow", lwd=2, lty=2) 
# Top Right 
plot(mt$time, mt$xn2, type = "l", xlab = "Age (in Days)", ylab = expression(x[paste("[2]")]), lwd=3, col = "purple") 
lines(aa, moment2, col = "yellow", lwd=2, lty=2)
plot(mt$time, mt$xn3, type = "l", lwd=3, col = "darkgreen",
     xlab = "Age (in Days)", ylab = expression(x[paste("[3]")]))
lines(aa, moment3, col = "yellow", lwd=2, lty=2)
par(mar = c(7, 4, 2, 2))
solve_dAoI(5/365, foiP3) -> mt 
plot(mt$time, mt$xn1, type = "l", xlab = "Age (in Days)", ylab = expression(list(x, sqrt(x[paste("[2]")]), sqrt(x[paste("[3]")], 3))), lwd=3, ylim = range(mt$xn3^(1/3))) 
lines(mt$time, mt$xn2^(1/2), col = "purple", lwd=3)
lines(mt$time, mt$xn3^(1/3), col = "darkgreen", lwd=3) 
lines(aa, moment1, col = "yellow", lwd=2, lty=2)
lines(aa, moment2^(1/2), col = "yellow", lwd=2, lty=2)
lines(aa, moment3^(1/3), col = "yellow", lwd=2, lty=2)

## -----------------------------------------------------------------------------
solve_dAoYda(5/365, foiP3, Amax=5*365, dt=5, n=9) -> mt

## -----------------------------------------------------------------------------
moment1y = momentAoY(aa, foiP3, hhat=5/365)

## ----fig.height=4, fig.width=6, echo=F----------------------------------------
plot(aa, moment1y, type = "l", lwd=2, xlab = "Age (in Days)")
with(mt, lines(time, y, col = "darkred", lwd=2))

## ----fig.height=4, fig.width=6, echo=F, eval=F--------------------------------
#  plot(moment1y-mt$y[-1], type = "l", xlab = "Age (in Days)", ylab = "Errors")

