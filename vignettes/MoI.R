## ----suppressMessages=T-------------------------------------------------------
library(ramp.pf.infection)
library(viridisLite)

## -----------------------------------------------------------------------------
library(ramp.pf.infection)
library(deSolve)
library(knitr)
library(viridisLite)
library(viridis)

## ----echo=F-------------------------------------------------------------------
aa = 0:1095
foiP3 = list(hbar = 5/365, 
             agePar = par_type2Age(), 
             seasonPar = par_sinSeason(), 
             trendPar = par_flatTrend())

## ----fig.height=4, fig.width=7, echo=F, eval=F--------------------------------
#  plot(aa, FoI(aa, foiP3), type = "l",
#       xlab = "a - age (in days)", ylab = expression(FoI[tau](a)))

## -----------------------------------------------------------------------------
MMinf <- solveMMinfty(5/365, foiP3, Amax=1095)


## ----fig.height=4, fig.width=7------------------------------------------------
with(MMinf, plot(time, m, type = "l", ylab = expression(m[tau](a)), xlab = "a - cohort age (in days)"))

## -----------------------------------------------------------------------------
hybrid = solve_dm(5/365, foiP3, Amax=1095)

## ----fig.height=4, fig.width=7------------------------------------------------
with(hybrid, plot(time, m, type = "l", ylab = expression(m[tau](a)), xlab = "a - cohort age (in days)"))

## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
moi = meanMoI(aa, foiP3, hhat=5/365)

## ----fig.height=4, fig.width=7------------------------------------------------
plot(aa, moi, type = "l", ylab = expression(m[tau](a)), xlab = "a - cohort age (in days)")

## -----------------------------------------------------------------------------
c(mean(abs(moi - hybrid$m)) < 1e-9,
mean(abs(MMinf$m- hybrid$m)) < 1e-10,
max(abs(moi - hybrid$m)) < 1e-8,
max(abs(MMinf$m- hybrid$m))< 1e-10)

