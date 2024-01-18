

dDensityPmu = function(xi, mu, a=0,
                       pRBC=par_lRBC.0(),
                       pSig=par_sigma.0()){
  # xi is log10 parasite densities
  # lRBC is the upper bound
  lRBC = log10RBC(a, pRBC)
  xi1 = xi/lRBC
  mu1 = mu/lRBC
  dbeta1(xi1, mu1, pSig)/lRBC
}

rDensityPmu = function(n, mu, a=0,
                       pRBC=par_lRBC.0(),
                       pSig=par_sigma.0()){
  # lRBC is the upper bound
  # xi is log10 parasite densities
  lRBC = log10RBC(a, pRBC)
  mu1 = mu/lRBC
  rbeta1(n, mu1, pSig)*lRBC
}

pDensityPmu = function(xi, mu, a=0,
                       pRBC=par_lRBC.0(),
                       pSig=par_sigma.0()){
  # lRBC is the upper bound
  # xi is log10 parasite densities
  lRBC = log10RBC(a, pRBC)
  mu1 = mu/lRBC
  xi1 = xi/lRBC
  pbeta1(xi1, mu1, pSig)
}

qDensityPmu = function(xi, mu, a=0,
                       pRBC=par_lRBC.0(),
                       pSig=par_sigma.0()){
  # lRBC is the upper bound
  # xi is log10 parasite densities
  lRBC = log10RBC(a, pRBC)
  mu1 = mu/lRBC
  qbeta1(xi, mu1, pSig)*lRBC
}





dDensityPalpha = function(x, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu.0(),
                          pRBC=par_lRBC.0(),
                          pSig=par_sigma.0()){
  # lRBC is the upper bound of realistic values
  # x is log10 parasite densities
  mu = alpha2mu(alpha, W, pMu)
  dDensityPmu(x, mu, a, pRBC, pSig)
}

rDensityPalpha = function(n, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu.0(),
                          pRBC=par_lRBC.0(),
                          pSig=par_sigma.0()){
  # lRBC is the upper bound of realistic values
  # x is log10 parasite densities
  mu = alpha2mu(alpha, W, pMu)
  rDensityPmu(n, mu, a, pRBC, pSig)
}

pDensityPalpha = function(x, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu.0(),
                          pRBC=par_lRBC.0(),
                          pSig=par_sigma.0()){
  # lRBC is the upper bound of realistic values
  # x is log10 parasite densities
  mu = alpha2mu(alpha, W, pMu)
  pDensityPmu(x, mu, a, pRBC, pSig)
}

qDensityPalpha = function(x, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu.0(),
                          pRBC=par_lRBC.0(),
                          pSig=par_sigma.0()){
  # lRBC is the upper bound of realistic values
  # x is log10 parasite densities
  mu = alpha2mu(alpha, W, pMu)
  qDensityPmu(x, mu, a, pRBC, pSig)
}

dDensityPa = function(x, a, FoIpar,
                      hhat=NULL,
                      tau=0, r=1/200,
                      pMu=par_alpha2mu.0(),
                      pRBC=par_lRBC.0(),
                      pSig=par_sigma.0(),
                      pWda=par_Wda.0()){
  dB = function(x, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda){
    px = function(alpha, x, W, a, m){
      mu = alpha2mu(alpha, W, pMu)
      y = zda(alpha, a, FoIpar, hhat, tau, r)
      dDensityPmu(x, mu, a, pRBC, pSig)*y/m
    }
    W = Wda(a, FoIpar, hhat, tau, pWda)
    m = meanMoI(a, FoIpar, hhat, tau, r)

    lowlim = 8
    uplim = min(a, 6*365)
    dff = integrate(dAoI, lowlim, a, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r)$value
    integrate(px, lowlim, uplim, x=x, W=W, a=a, m=m)$value/dff
  }
  if(length(x)==1) return(dB(x, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda))
  return(sapply(x, dB, a=a, FoIpar=FoIpar, hhat=hhat,  tau=tau, r=r, pMu=pMu, pRBC=pRBC, pSig=pSig, pWda=pWda))
}











pDensityPa = function(x, a, FoIpar,
                      hhat=NULL,tau=0, r=1/200,
                      pMu=par_alpha2mu.0(),
                      pRBC=par_lRBC.0(),
                      pSig=par_sigma.0(),
                      pWda=par_Wda.delta()){
  pB = function(x, a, FoIpar, hhat,tau, r, pMu, pRBC, pSig, pWda){
    px = function(alpha, x, W, a){
      mu = alpha2mu(alpha, W, pMu)
      A = dAoI(alpha, a, FoIpar, hhat, tau, r)
      dDensityPmu(x, mu, a, pRBC, pSig)*A
    }
    W = Wda(a, FoIpar, hhat, tau, pWda)
    dff = integrate(syda, 7, a, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r)$value
    integrate(px, 7, a, x=x, W=W, a=a)$value/dff
  }
  if(length(x)==1) return(pB(x, a, FoIpar, hhat,tau, r, pMu, pRBC, pSig, pWda))
  return (sapply(x, pB, a=a, FoIpar=FoIpar, h=h, tau=tau, r=r, pMu=pMu, pSig=pSig, pRBC=pRBC, pWda=pWda))
}

#dDensityPaConvolve2 = function(x, a, FoIpar, h, r=1/200, tau=0, pMu=par_alpha2mu.0(), pRBC=par_lRBC.0(), pSig=par_sigma.0(),pWda=pWda0){
#  sapply(x, dDensityPaConvolve2, a=a, FoIpar=FoIpar, h=h, r=r,d=d,pMu=pMu, pRBC=pRBC,pSig=pSig,pWda=pWda)
#}


