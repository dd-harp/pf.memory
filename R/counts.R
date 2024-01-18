

pCountsPa = function(a, FoIpar, bins=NULL, dx=0.1,
                     hhat=NULL,tau=0, r=1/200,
                     pMu=par_alpha2mu.0(),
                     pRBC=par_lRBC.0(),
                     pSig=par_sigma.0(),
                     pWda=par_Wda.delta(),
                     pC = par_nbCounts()){
  if(is.null(bins)) bins=c(1:5,13)
  lRBC = log10RBC(a, pRBC)
  meshX = seq(0, lRBC, by=dx)
  mix = 1:length(meshX)
  Pda = dDensityPa(meshX, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda)
  Pda = Pda/sum(Pda)
  cdfP = function(xi, a, Cpar){
    pD = Detect(xi,a,Cpar)
    counts = pCounts(bins, xi, a, Cpar)
    #browser()
  }
  Bx=sapply(meshX, cdfP, a=a, Cpar=pC)
  list(bins=bins,cdf=rowSums(Bx*Pda))
}



dCountsPa = function(a, FoIpar,  bins=NULL, dx=0.1,
                     hhat=NULL,tau=0, r=1/200,
                     pMu=par_alpha2mu.0(),
                     pRBC=par_lRBC.0(),
                     pSig=par_sigma.0(),
                     pWda=par_Wda.delta(),
                     pC = par_nbCounts()){
  DP = DetectPa(a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  PC = pCountsPa(a, FoIpar, bins, dx, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  list(bins=PC$bins, pdf=diff(c(DP, PC$cdf))/(1-DP), detect=DP, cdf=PC$cdf, fullpdf = c(DP,diff(c(DP, PC$cdf))))
}







pCountsBa = function(a, FoIpar, bins=NULL, dx=0.1,
                     hhat=NULL,tau=0, r=1/200,
                     pMu=par_alpha2mu.0(),
                     pRBC=par_lRBC.0(),
                     pSig=par_sigma.0(),
                     pWda=par_Wda.delta(),
                     pC = par_nbCounts()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = DetectPa(a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  out= pCountsPa(a, FoIpar, bins, dx, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  list(bins=out$bins, cdf=(1 - exp(-moi*D*out$cdf))/(1-exp(-moi*D)), DM = 1 - exp(-moi*D))
}




dCountsBa = function(a, FoIpar, bins=NULL, dx=0.1,
                     hhat=NULL,tau=0, r=1/200,
                     pMu=par_alpha2mu.0(),
                     pRBC=par_lRBC.0(),
                     pSig=par_sigma.0(),
                     pWda=par_Wda.delta(),
                     pC = par_nbCounts()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = DetectPa(a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  PC = pCountsBa(a, FoIpar, bins, dx, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  DP = PC$DM
  list(bins=PC$bins, pdf=diff(c(DP, PC$cdf))/(1-DP), detect=DP, cdf=PC$cdf, fullpdf = c(DP, diff(c(DP, PC$cdf))))
}


meanCountsPa = function(a, FoIpar, dx=0.1,
                        hhat=NULL,tau=0, r=1/200,
                        pMu=par_alpha2mu.0(),
                        pRBC=par_lRBC.0(),
                        pSig=par_sigma.0(),
                        pWda=par_Wda.delta(),
                        pC = par_nbCounts()){
  bins = seq(0,7, by=dx)
  PC = dCountsPa(a, FoIpar, bins, dx, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  sum(PC$bins*PC$pdf)
}



meanCountsBa = function(a, FoIpar, dx=0.1,
                        hhat=NULL,tau=0, r=1/200,
                        pMu=par_alpha2mu.0(),
                        pRBC=par_lRBC.0(),
                        pSig=par_sigma.0(),
                        pWda=par_Wda.delta(),
                        pC = par_nbCounts()){
  bins = seq(0,7, by = dx)
  PC = dCountsBa(a, FoIpar, bins, dx, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  sum(PC$bins*PC$pdf)
}
