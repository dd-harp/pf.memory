## ----suppressMessages=T-------------------------------------------------------
library(ramp.pf.infection)
library(viridisLite)

## ----echo=F-------------------------------------------------------------------
foiP3 = list(hbar = 5/365, 
             agePar = par_type2Age(), 
             seasonPar = par_sinSeason(), 
             trendPar = par_flatTrend())

foiP3t = list(hbar = 5/365, 
             agePar = par_flatAge(), 
             seasonPar = par_sinSeason(), 
             trendPar = par_flatTrend())

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,1))
a3years  = 1:1095
plot(a3years, FoI(a3years, foiP3t), type = "l", 
     xlab = "Time or Age (in days)", 
     ylab = expression(z[tau](alpha, a)))
lines(a3years, FoI(a3years, foiP3), col = "darkred", lwd=2)
text(270, 0.02, "Population")
text(100, 0.008, "Cohort", col = "darkred")

## ----zde eg-------------------------------------------------------------------
alpha = 60
a = 6*365
zda(60, 6*365, foiP3) 

## -----------------------------------------------------------------------------
zz = zda(a3years, max(a3years), foiP3)

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,1))
plot(a3years, zz, type = "l", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(z[tau](alpha,a)))

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,1))
clrs12a <- turbo(12) 
zz = zda(a3years, 730, foiP3, tau=90)
plot(a3years, zz, type = "n",
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(z[tau](alpha,a)))
for(i in 1:12){
  birthday =  30*(i-1) +15
  lines(a3years, zda(a3years, 730, foiP3, tau=birthday), col = clrs12a[i]) 
  bp = (0.5 + 2.5*i/12)*365
  points(bp, 0.015, col = clrs12a[i], pch =15)
  text(bp, 0.015, paste(i), col = clrs12a[i], pos=1)
}
  text(650, 0.018, "Birth Month") 

## -----------------------------------------------------------------------------
mm = meanMoI(a3years, foiP3, hhat=5/365)

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,1))
plot(a3years, mm, type = "l",
     xlab = expression(list(a, paste("Host Age (in Days)"))), 
     ylab = expression(m[tau](alpha,a)))

## -----------------------------------------------------------------------------
f_A = dAoI(a3years, max(a3years), foiP3)

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,1))
plot(a3years, f_A, type ="l", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(f[A](alpha, a, tau, h)))

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,1))
clrs12a <- turbo(12) 
FA = dAoI(a3years, 1095, foiP3, tau=70)
plot(a3years, FA, type = "n", ylim = range(FA)*1.15,
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(f[A](alpha, tau)))
for(i in 1:12){
  birthday =  30*(i-1) +15
  lines(a3years, dAoI(a3years, 730, foiP3, tau=birthday), col = clrs12a[i]) 
  bp = (0.5 + 2.5*i/12)*365
  points(bp, 0.01, col = clrs12a[i], pch =15)
  text(bp, 0.01, paste(i), col = clrs12a[i], pos=1)
}
  text(650, 0.0115, "Birth Month") 

## -----------------------------------------------------------------------------
F_A = pAoI(a3years, max(a3years), foiP3)

## -----------------------------------------------------------------------------
F_A_alt = cumsum(f_A)

## ----fig.height=3.5, fig.width=6----------------------------------------------
par(mar = c(5,4,1,1))

plot(a3years, F_A, type = "l", 
     xlab = "Parasite Cohort Age", 
     ylab = expression(1-F[X](alpha, a, tau)), lwd=3)

lines(a3years, F_A_alt, col = "red", lwd=2, lty =2)

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,1))
clrs12a <- turbo(12) 
FA = pAoI(a3years, 1095, foiP3, tau=70)
plot(a3years, FA, type = "n", ylim = range(FA)*1.15,
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(F[A](alpha)))
for(i in 1:12){
  birthday =  30*(i-1) +15
  lines(a3years, pAoI(a3years, 730, foiP3, tau=birthday), col = clrs12a[i]) 
  bp = (0.5 + 2.5*i/12)*365
  points(bp, 0.2, col = clrs12a[i], pch =15)
  text(bp, 0.2, paste(i), col = clrs12a[i], pos=1)
}
text(650, 0.3, "Birth Month") 

## -----------------------------------------------------------------------------
rhx = rAoI(10000, 3*365, foiP3)

## ----fig.height=3.5, fig.width=6----------------------------------------------
par(mar = c(5,4,1,2))
plot(stats::ecdf(rhx), xlim = c(0,1095), cex=0.2, main = "", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(list(F[A](alpha), paste("ecdf"))))
lines(a3years, F_A, col = "red", lty = 2, lwd=2)

## ----fig.height=3.5, fig.width=6, echo=F--------------------------------------
par(mar = c(5,4,1,2))
hist(rhx, breaks = seq(0, 1095, by=15), 
     right=F, probability=T, main = "", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     border = grey(0.5))
lines(a3years, f_A, col="red", lwd=2)

## ----eval=T-------------------------------------------------------------------
moment1 = momentAoI(a3years, foiP3)
moment2 = momentAoI(a3years, foiP3, n=2)
moment3 = momentAoI(a3years, foiP3, n=3)

## ----fig.width=6, fig.height=8, echo=FALSE, eval=T----------------------------
par(mfrow = c(3,1), mar = c(5, 4, 0.5, 2))
plot(a3years, moment1, type = "l", xlab = "a - Host Age (in Days)", 
     lwd=2, ylim = range((moment3)^(1/3)),
     ylab = expression(list(x, sqrt(x[paste("[2]")]), sqrt(x[paste("[3]")], 3))))

lines(a3years, sqrt(moment2), col = "darkgreen")
lines(a3years, (moment3)^(1/3), col = "purple")

plot(a3years, moment2, type = "l", xlab = "a - Host Age (in Days)",
     lwd=2, col = "darkgreen",
     ylab = expression(x[paste("[2]")])) 
mar = c(4, 4, 0.5, 2)
plot(a3years, moment3, type = "l", xlab = "a - Host Age (in Days)",
     lwd=2, col = "purple", 
     ylab = expression(x[paste("[3]")])) 

## ----eval=T-------------------------------------------------------------------
f_Y = dAoY(a3years, 3*365, foiP3)

## ----fig.height=3.5, fig.width=6, echo=F, eval=T------------------------------
par(mar = c(5,4,1,1))
plot(a3years, f_Y, type ="l", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     ylab = expression(list(f[Y](alpha), f[A](alpha))), ylim = range(f_Y), lwd=2)
lines(a3years, f_A, col = grey(0.5))

## ----eval=T-------------------------------------------------------------------
raoy = rAoY(10^5, 3*365, foiP3)

## ----fig.height=4, fig.width=7, eval=T----------------------------------------
hist(raoy, breaks=seq(0, 1095, by = 15), 
     right=F, probability=T, main = "", 
     xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
     border = grey(0.5)) -> out
lines(a3years, f_Y, type = "l", col = "red") 

## ----eval=T-------------------------------------------------------------------
aa = seq(5, 3*365, by = 5) 
moment1y = momentAoY(aa, foiP3)
moment2y = momentAoY(aa, foiP3, n=2)
moment3y = momentAoY(aa, foiP3, n=3)

## ----fig.width=6, fig.height=8, echo=FALSE, eval=T----------------------------
par(mfrow = c(4,1), mar = c(0.5, 4, 0.5, 2))
plot(aa, moment1y, type = "l", xlab = "", ylab = expression(E(Y)), lwd=2, xaxt="n", ylim = range((moment3y)^(1/3)) ) 
lines(aa, sqrt(moment2y), col = "darkgreen")
lines(aa, (moment3y)^(1/3), col = "purple")
plot(aa, moment2y, type = "l", xlab = "", ylab = expression(E(Y^2)), lwd=2, xaxt="n", col = "darkgreen") 
plot(aa, moment3y, type = "l", xlab = "", ylab = expression(E(Y^2)), lwd=2, col = "purple") 
mtext("Age (in Days)", 1, 3)

