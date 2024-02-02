#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
Detect = function(xi, a=0, Cpar=par_nbCounts()){
  UseMethod("Detect", Cpar)
}

#' density function of counts
#'
#' @param hatxi mean log10 parasite counts
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
dCounts = function(hatxi, xi, a=0, Cpar=par_nbCounts()){
  UseMethod("dCounts", Cpar)
}

#' cumulative density function of counts
#'
#' @param hatxi mean log10 parasite densities
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
pCounts = function(hatxi, xi, a=0, Cpar=par_nbCounts()){
  UseMethod("pCounts", Cpar)
}

#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param bins breakpoints for summarizing outputs
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
binnedCounts = function(xi, a, bins, Cpar=par_nbCounts()){
  UseMethod("binnedCounts", Cpar)
}

#' Detection of infection given parasitemia
#'
#' @param q the fraction of blood examined
#' @param sz the negative binomial size parameter
#' @param bvm blood volume as log10 red blood cells
#'
#' @return binary detection result
#' @export
#'
par_nbCounts = function(q=6, sz=0.31, bvm = par_lRBC_static()){
  par = list()
  class(par) <- "nb"
  par$q=q
  par$sz=sz
  par$bvm = bvm
  par
}

#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
Detect.nb = function(xi, a=0, Cpar=par_nbCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  1-stats::dnbinom(0, mu=10^(xi-lRBC+q), size=sz)
})}



#' Detection of infection given parasitemia
#'
#' @param hatxi mean log10 parasite counts
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
dCounts.nb = function(hatxi, xi, a=0, Cpar=par_nbCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  stats::dnbinom(round(10^hatxi), mu=10^(xi-lRBC+q), size=sz)
})}



#' Detection of infection given parasitemia
#'
#' @param hatxi mean log10 parasite counts
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
pCounts.nb = function(hatxi, xi, a=0, Cpar=par_nbCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  pnbinom(10^hatxi, mu=10^(xi-lRBC+q), size=sz)
})}

#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param bins break points for binning data
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
binnedCounts.nb = function(xi, a, bins, Cpar=par_nbCounts()){
  if(is.null(bins)) bins=c(1:5,13)
  p0 = Detect(xi,a,Cpar)
  p1 = pCounts(bins, xi, a, Cpar)
  diff(c(1-p0,p1))/p0
}

#' Detection of infection given parasitemia

#' @param q the fraction of blood volume sampled
#' @param bvm blood volume as log10 red blood cells
#'
#' @return par a [list]
#' @export
par_poisCounts = function(q=6, bvm = par_lRBC_static()){
  par = list()
  class(par) <- "pois"
  par$q=q
  par$bvm = bvm
  par
}

#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
Detect.pois = function(xi, a=0, Cpar=par_poisCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  1-stats::dpois(0, 10^(xi-lRBC+q))
})}



#' Detection of infection given parasitemia
#'
#' @param hatxi mean log10 parasite densities
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
dCounts.pois = function(hatxi, xi, a=0, Cpar=par_poisCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  stats::dpois(round(10^hatxi), 10^(xi-lRBC+q))
})}



#' Detection of infection given parasitemia
#'
#' @param hatxi mean log10 parasite counts
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
pCounts.pois = function(hatxi, xi, a=0, Cpar=par_poisCounts()){with(Cpar,{
  lRBC = log10RBC(a, bvm)
  ppois(10^hatxi, 10^(xi-lRBC+q))
})}



#' Detection of infection given parasitemia
#'
#' @param xi mean log10 parasite densities
#' @param a host cohort age
#' @param bins breakpoints for binning data
#' @param Cpar parameters that define a detection function
#'
#' @return binary detection result
#' @export
binnedCounts.pois = function(xi, a, bins, Cpar=par_poisCounts()){
  if(is.null(bins)) bins=c(1:5,13)
  p0 = Detect(xi,a,Cpar)
  p1 = pCounts(bins, xi, a, Cpar)
  diff(c(1-p0,p1))/p0
}


#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return a [numeric] vector of length(a)
DetectPa = function(a, FoIpar,
                    hhat=NULL,tau=0, r=1/200,
                    pMu=par_alpha2mu_base(),
                    pRBC=par_lRBC_static(),
                    pSig=par_sigma_abc(),
                    pWda=par_Wda_none(),
                    pC = par_nbCounts()){
  pD = function(a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC){
    Dx = function(x, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC){
      dDensityPa(x, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda)*Detect(x,a,pC)
    }
    hatb=log10RBC(a,pRBC)
    stats::integrate(Dx, 0, hatb, a=a, FoIpar=FoIpar, hhat=hhat,tau=tau,r=r,pMu=pMu,pRBC=pRBC, pSig=pSig,pWda=pWda, pC=pC)$value
  }
  if(length(a)==1) return(pD(a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC))
  return (sapply(a, pD, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r, pMu=pMu, pSig=pSig, pRBC=pRBC, pWda=pWda, pC=pC))
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return binary detection result
#' @export
DetectPM = function(a, FoIpar,
                    hhat=NULL,tau=0, r=1/200,
                    pMu=par_alpha2mu_base(),
                    pRBC=par_lRBC_static(),
                    pSig=par_sigma_abc(),
                    pWda=par_Wda_none(),
                    pC = par_nbCounts()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = DetectPa(a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  1 - exp(-moi*D)
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param dx width of the mesh
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return binary detection result
#' @export
DetectBda = function(a, FoIpar, dx=0.1,
                     hhat=NULL,tau=0, r=1/200,
                     pMu=par_alpha2mu_base(),
                     pRBC=par_lRBC_static(),
                     pSig=par_sigma_abc(),
                     pWda=par_Wda_none(),
                     pC = par_nbCounts()){

  lRBC = log10RBC(a, pRBC)
  meshX = seq(0, lRBC, by = dx)
  Bx = Bda(meshX, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda)$PDFm
  Bx = Bx/sum(Bx)
  t(Detect(meshX, a, pC)) -> Bs
  sum(Bx*Bs)
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return detection probability
#' @export
MoIPM = function(a, FoIpar,
                 hhat=NULL,tau=0, r=1/200,
                 pMu=par_alpha2mu_base(),
                 pRBC=par_lRBC_static(),
                 pSig=par_sigma_abc(),
                 pWda=par_Wda_none(),
                 pC = par_nbCounts()){

  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = DetectPa(a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda, pC)
  moi*D/(1-exp(-moi*D))
}
