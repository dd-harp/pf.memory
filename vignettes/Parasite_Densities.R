## ----suppressMessages=T-------------------------------------------------------
library(ramp.pf.infection)
library(viridisLite)

## -----------------------------------------------------------------------------
base_pars = par_Fmu_base()
base_pars1 = par_Fmu_base(tildeb=11, Sa=0.0045)
chronic_pars = par_Fmu_chronic()

## -----------------------------------------------------------------------------
aalpha = 7:366
mu_alpha =  Fmu(aalpha, 0, base_pars)
mu_alpha1 =  Fmu(aalpha, 0, base_pars1)
mu_alpha2 =  Fmu(aalpha, 0, chronic_pars)

## ----fig.height=4, fig.width=6, echo=F----------------------------------------
plot(aalpha, mu_alpha, type = "l", ylim = c(0, 13), xlab = expression(alpha), ylab = expression (mu(alpha)))
lines(aalpha, mu_alpha1, col = "darkblue")
lines(aalpha, mu_alpha2, col = "darkred")
segments(0,6,365,6, lty=2)

## -----------------------------------------------------------------------------
parD = par_Omega_beta() 
d_Omega(9, 9, 13, parD)

## ----fig.height=8, fig.width=8------------------------------------------------
par(mfrow = c(2,2))
xx = seq(0, 13, by=0.1)
plot(xx, d_Omega(xx, 9, 13, parD), type = "l", 
     xlab = expression(xi), ylab = "d_Omega", main = "d_Omega")
lines(xx, d_Omega(xx, 7,  13, parD), col = "darkred")
pp = seq(0,13, by = 0.01)
plot(pp, p_Omega(pp, 9,  13, parD), type = "l", 
     xlab = expression(xi), ylab = "p_Omega", main = "p_Omega")
lines(xx, p_Omega(xx, 7, 13, parD), col = "darkred")
hist(r_Omega(1000, 9, 13, parD), xlab = expression(xi), 
     main = "r_Omega", xlim = c(0,13))
qq = seq(0.01,.99, length.out=100)
plot(qq, q_Omega(qq, 9, 13, parD), type = "l", 
     xlab = expression(xi), ylab = "q_Omega", main = "q_Omega")

## ----fig.height=4, fig.width=7, echo=F----------------------------------------
foiP3 = list(hbar = 5/365, 
             agePar = par_type2Age(), 
             seasonPar = par_sinSeason(), 
             trendPar = par_flatTrend())

a5years  = 0:(5*365)
clrs = turbo(2)
plot(a5years, FoI(a5years, foiP3, tau=0), type = "l", col = clrs[1],
     xlab = "a - host cohort age (in days)", ylab = expression(FoI[tau](a)))
lines(a5years, FoI(a5years, foiP3, tau=180), col = clrs[2])
#lines(a5years, FoI(a5years, foiP3, tau=270), col = clrs[4])

## -----------------------------------------------------------------------------
mu_alpha =  Fmu(a5years, 0, base_pars)
dAoI(a5years, 5*365, foiP3) -> dalpha

## ----echo=F, fig.height=4, fig.width=7, echo=F--------------------------------
plot(a5years, dalpha, type = "l", 
    xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
    ylab = expression(f[A](alpha)))

## ----echo=F, fig.height=4.5, fig.width=7, echo=F------------------------------
par(mar = c(5,5,3,4))
dd = dalpha*max(mu_alpha)/max(dalpha)
plot(a5years, mu_alpha, type = "l", ylim = c(0,11),
    xlab = expression(list(alpha, paste("Parasite Age (in Days)"))), 
    ylab = expression(list(xi, paste("Parasite Densities"))), 
    main = expression(list(mu(alpha), f[A]())))

lines(a5years, dd, col = grey(0.8), lwd=2)

clrs = rev(turbo(50))
i = 1:50
points(20*i, dd[20*i], col =clrs, pch = 15, cex=0.7) 
tks = pretty(dalpha)
axis(4, tks*max(mu_alpha)/max(dalpha), tks)
mtext(expression(f[A](mu)), 4,3)

## -----------------------------------------------------------------------------
alpha = 60
fA = dAoI(alpha, 5*365, foiP3)
Pmu = fA*d_alpha2density(xx, alpha, 5*365)

## ----echo=F, fig.height=4, fig.width=7, echo=F--------------------------------
plot(xx, Pmu, type = "l", ylim = range(Pmu)*1.2)

## -----------------------------------------------------------------------------
Pa = d_clone_density(xx, 5*365, foiP3)

## ----echo=F, fig.height=4.5, fig.width=7--------------------------------------

parD = par_Omega_beta()
alpha = 20
mu = Fmu(alpha, 0, base_pars)
fA = dAoI(alpha, 5*365, foiP3)
Pmu = fA*d_Omega(xx, mu, 13, parD)

par(mar = c(5,5,3,4))
plot(xx, Pa, type = "l",
     ylab = expression(P[a](xi)),
     xlab = expression(xi),
     main = "Parasite Density Distributions" )
segments(6,0, 6, 0.004, col = grey(0.5), lty=2)
for(i in (1:50)){
  alpha = 1020-20*i
  mu = Fmu(alpha, 0, base_pars)
  fA = dAoI(alpha, 5*365, foiP3)
  lines(xx, fA*d_Omega(xx, mu, 13, parD)*max(Pa)/max(Pmu), col = clrs[51-i])
}
tks = pretty(Pa/40, 6)
axis(4, tks*40, tks)
mtext(expression(P[mu](alpha)), 4, 3)

## -----------------------------------------------------------------------------
integrate(d_clone_density, 0, 13, a=5*365, FoIpar=foiP3)$value 

