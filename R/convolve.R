
dDensityPaConvolve2 = function(x, a, FoIpar,
                               hhat=NULL, tau=0,  r=1/200,
                               pMu=par_alpha2mu.0(),
                               pRBC=par_lRBC.0(),
                               pSig=par_sigma.0(),
                               pWda=par_Wda.delta()){
  px = function(x, log10B, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda){
    lB2 = log10(10^log10B - 10^x)
    dDensityPa(x,a,FoIpar, hhath, tau, r, pMu, pRBC, pSig, pWda)*dDensityPa(lB2,a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda)
  }
  integrate(px, 0, x, log10B=x, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r, pMu=pMu, lRBC=pRBC, pSig=pSig,pWda=pWda)$value
}



cdfConvolve2b = function(meshX, CDF1, CDF2){
  cX = CDF1*0
  L = length(meshX)
  PDF1 = c(CDF1[1], diff(CDF1))
  PDF1 = PDF1/sum(PDF1)
  B1 = B2 = 10^meshX
  for(i in 1:L){
    for(j in 1:L){
      x = log10(B1[j] + B2)
      ix2 = which(x <= meshX[i])
      if(length(ix2>0)){
        cX[i] = cX[i] + PDF1[j]*CDF2[max(ix2)]
      }
    }
  }
  cX
}

cdfConvolve2a = function(meshX, CDF1, CDF2){
  cX = CDF1*0
  L = length(meshX)
  PDF1 = c(CDF1[1], diff(CDF1))
  PDF1 = PDF1/sum(PDF1)
  B1 = B2 = B = 10^meshX
  for(i in 1:L){
    p1=0
    ix1 = which(B1[i]>B)
    if(length(ix1>0)) p1=PDF1[min(ix1)]
    for(j in 1:L){
      x = log10(B1[j] + B2)
      ix2 = which(x > meshX[i])
      if(length(ix2>0)){
        cX[i] = cX[i] + PDF1[j]*CDF2[min(ix2)]
      }
    }
  }
  cX
}

cdfConvolve2 = function(meshX, CDF1, CDF2){
  CDFa = cdfConvolve2a(meshX, CDF1, CDF2)
  CDFb = cdfConvolve2a(meshX, CDF1, CDF2)
  (CDFa+CDFb)/2
}
