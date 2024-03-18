## ----include = FALSE----------------------------------------------------------
library(knitr) 
library(ramp.pf.infection)

## ----plotFoI functions--------------------------------------------------------
plotFoI= function(h, FoIpar, tau=0, Tmax=5*365, clrs=NULL){
  L = length(h)
  if(is.null(clrs)) clrs = rep("black", L)
  tms = 1:Tmax
  tm = tms/365 
  plot(tm, 365*h[1]*FoI(tms, FoIpar, tau), type = "n", 
       xlab = "Host Age (Years)", ylab = "aFoI") 
  for(i in 1:L) lines(tm, 365*h[i]*FoI(tms,FoIpar, tau), col = clrs[i])
}

## ----plotFoI------------------------------------------------------------------
clrs8b <- viridis::viridis(8, option = "B")

FoIfigure = function(hh=c(8, 4, 2, 1, .5, .2, 0.1)/365, 
                     tau=0, clrs = clrs8b, Tmax=5*365){
  par(mfcol = c(1,2), mar =c(4,4,1,1))
  plotFoI(hh, foiP4, tau, Tmax, clrs) 
  plotFoI(hh, foiP3, tau, Tmax, clrs) 
} 

