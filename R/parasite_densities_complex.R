

Bda = function(meshX, a, FoIpar,
               hhat=NULL, tau=0, r=1/200,
               pMu=par_alpha2mu.0(),
               pRBC=par_lRBC.0(),
               pSig=par_sigma.0(),
               pWda=par_Wda.delta()){
  PDFx = dDensityPa(meshX, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda)
  CDFx = cumsum(PDFx)
  CDFx = CDFx/max(CDFx)
  PDFx = c(CDFx[1], diff(CDFx))
  PDFx = PDFx/sum(PDFx)

  moi = meanMoI(a, FoIpar, hhat, tau, r)
  N = max(4*moi, 10)

  cdflist = list()
  cdflist$cdf = list()
  cdflist$pdf = list()
  cdflist$cdf[[1]] = CDFx
  cdflist$pdf[[1]] = PDFx
  CDFm = dpois(1, moi)*CDFx

  CDFn=CDFx
  PDFn=PDFx
  for(i in 2:N){
    cdfn = cdfConvolve2(meshX, CDFx, CDFn)
    #cdfn2 = cdfConvolve2(meshX, CDFn, CDFx)
    #browser()
    CDFn = cdfn
    PDFn = c(CDFn[1], diff(CDFn))
    PDFn = PDFn/sum(PDFn)
    cdflist$cdf[[i]] = CDFn
    cdflist$pdf[[i]] = PDFn
    CDFm = CDFm + dpois(i, moi)*CDFn
  }
  cdflist$CDFm = CDFm/(1-dpois(0,moi))
  PDFm =c(CDFm[1], diff(CDFm))
  cdflist$PDFm = PDFm/sum(PDFm)
  cdflist
}

dBda = function(meshX, a, FoIpar, hhat=NULL, tau=0, r=1/200, pMu=par_alpha2mu.0(), pRBC=par_lRBC.0(), pSig=par_sigma.0(), pWda=par_Wda.delta()){
  Bda(meshX, a, FoIpar,hhat, tau, r, pMu, pRBC, pSig, pWda)$PDFm
}

pBda = function(meshX, a, FoIpar, hhat=NULL, tau=0, r=1/200, pMu=par_alpha2mu.0(), pRBC=par_lRBC.0(), pSig=par_sigma.0(), pWda=par_Wda.delta()){
  Bda(meshX, a, FoIpar,hhat, tau, r, pMu, pRBC, pSig, pWda)$CDFm
}





rBda = function(N, a, FoIpar,
                hhat=NULL, r=1/200, tau=0, alphamin=7,
                pMu=par_alpha2mu.0(),
                pRBC=par_lRBC.0(),
                pSig=par_sigma.0(),
                pWda=par_Wda.delta()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  W = Wda(a, FoIpar, hhat, tau, pWda)
  # MoI, excluding zeros
  hatm = nzPois(N, moi)
  Ny = sum(hatm)
  hatalpha = rAoI(Ny, a, FoIpar, hhat, tau, r, alphamin)
  # their expected values
  hatmu = alpha2mu(hatalpha, W, pMu)
  hatx = rDensityPmu(Ny,hatmu,a,pRBC,pSig)
  lRBC = 10^hatx
  first = sum(lRBC[1:hatm[1]])
  rest = diff(cumsum(lRBC)[cumsum(hatm)])
  log10(c(first, rest))
}



rRda = function(N, R, a, FoIpar, hhat=NULL, tau=0, r=1/200, alphamin=7,
                pMu=par_alpha2mu.0(),
                pRBC=par_lRBC.0(),
                pSig=par_sigma.0(),
                pWda=par_Wda.delta()){
  W = Wda(a, FoIpar, hhat, tau, pWda)
  # MoI, excluding zeros
  Ny = R*N
  hatalpha = rAoI(Ny, a, FoIpar, hhat, tau, r, alphamin)
  # their expected values
  hatmu = alpha2mu(hatalpha, W, pMu)
  hatx = rDensityPmu(Ny,hatmu,a,pRBC,pSig)
  lRBC = 10^hatx
  matrix(lRBC,nrow=R, ncol=N)
}

