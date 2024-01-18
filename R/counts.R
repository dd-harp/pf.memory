

#' Compute the distribution of parasite counts for simple infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param bins a set of break points for computing counts
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#' @param pC
#'
#' @return a [list]
#' @export
#'
#' @examples
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
  }
  Bx=sapply(meshX, cdfP, a=a, Cpar=pC)
  list(bins=bins,cdf=rowSums(Bx*Pda))
}



#' Compute the density of parasite counts for simple infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param bins a set of break points for computing counts
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#' @param pC
#'
#' @return a [list]
#' @export
#'
#' @examples
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







#' Compute the distribution of parasite counts in complex infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param bins a set of break points for computing counts
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#' @param pC
#'
#' @return a [list]
#' @export
#'
#' @examples
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




#' Compute the density of parasite counts in complex infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param bins a set of break points for computing counts
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#' @param pC
#'
#' @return a [list]
#' @export
#'
#' @examples
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


#' Compute the mean parasite counts in simple infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#' @param pC
#'
#' @return a [list]
#' @export
#'
#' @examples
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



#' Compute the mean parasite counts in complex infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#' @param pC
#'
#' @return a [list]
#' @export
#'
#' @examples
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
