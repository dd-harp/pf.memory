

#' Compute the CDF and PDF of parasite densities in a host cohort
#'
#' @param meshX a mesh over parasite densities
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length meshX
#' @export
#'
Bda = function(meshX, a, FoIpar,
               hhat=NULL, tau=0, r=1/200,
               pMu=par_alpha2mu_base(),
               pRBC=par_lRBC_static(),
               pSig=par_sigma_abc(),
               pWda=par_Wda_none()){
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
  CDFm = stats::dpois(1, moi)*CDFx

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
    CDFm = CDFm + stats::dpois(i, moi)*CDFn
  }
  cdflist$CDFm = CDFm/(1-stats::dpois(0,moi))
  PDFm =c(CDFm[1], diff(CDFm))
  cdflist$PDFm = PDFm/sum(PDFm)
  cdflist
}

#' The density function for parasite densities in a host cohort
#'
#' @param meshX a mesh over parasite densities
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length meshX
#' @export
#'
dBda = function(meshX, a, FoIpar, hhat=NULL, tau=0, r=1/200, pMu=par_alpha2mu_base(), pRBC=par_lRBC_static(), pSig=par_sigma_abc(), pWda=par_Wda_none()){
  Bda(meshX, a, FoIpar,hhat, tau, r, pMu, pRBC, pSig, pWda)$PDFm
}

#' The distribution function for parasite densities in a host cohort
#'
#' @param meshX a mesh over parasite densities
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a numeric vector of length meshX
#' @export
#'
pBda = function(meshX, a, FoIpar, hhat=NULL, tau=0, r=1/200, pMu=par_alpha2mu_base(), pRBC=par_lRBC_static(), pSig=par_sigma_abc(), pWda=par_Wda_none()){
  Bda(meshX, a, FoIpar,hhat, tau, r, pMu, pRBC, pSig, pWda)$CDFm
}

#' Random generation for parasite densities in a host cohort
#'
#' @param N number of observations
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param alphamin the minimum value of alpha allowed
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length N
#' @export
#'
rBda = function(N, a, FoIpar,
                hhat=NULL, r=1/200, tau=0, alphamin=7,
                pMu=par_alpha2mu_base(),
                pRBC=par_lRBC_static(),
                pSig=par_sigma_abc(),
                pWda=par_Wda_none()){
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



#' Random generation for M parasite densities in a host cohort with MoI
#'
#' @param R the number of observations
#' @param M the MoI
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of alpha allowed
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a R by M [matrix]
#' @export
#'
rRda = function(M, R, a, FoIpar, hhat=NULL, tau=0, r=1/200, alphamin=7,
                pMu=par_alpha2mu_base(),
                pRBC=par_lRBC_static(),
                pSig=par_sigma_abc(),
                pWda=par_Wda_none()){
  W = Wda(a, FoIpar, hhat, tau, pWda)
  # MoI, excluding zeros
  Ny = R*M
  hatalpha = rAoI(Ny, a, FoIpar, hhat, tau, r, alphamin)
  # their expected values
  hatmu = alpha2mu(hatalpha, W, pMu)
  hatx = rDensityPmu(Ny,hatmu,a,pRBC,pSig)
  lRBC = 10^hatx
  matrix(lRBC,nrow=R, ncol=M)
}

