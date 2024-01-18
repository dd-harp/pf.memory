
#' The density function for parasite densities in a simple malaria infection
#'
#' @param xi a vector of quantiles
#' @param mu the expected value for log10 parasite densities
#' @param a host cohort age
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
#' @examples
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

#' Random generation of parasite densities from a simple malaria infection
#'
#' @param n the number of observations
#' @param mu the expected value for log10 parasite densities
#' @param a host cohort age
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
#' @examples
rDensityPmu = function(n, mu, a=0,
                       pRBC=par_lRBC.0(),
                       pSig=par_sigma.0()){
  # lRBC is the upper bound
  # xi is log10 parasite densities
  lRBC = log10RBC(a, pRBC)
  mu1 = mu/lRBC
  rbeta1(n, mu1, pSig)*lRBC
}

#' The density function for parasite densities as a function of the mean
#'
#' @param xi a vector of probabilities
#' @param mu the expected value for log10 parasite densities
#' @param a host cohort age
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
#' @examples
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

#' The quantile function for parasite densities in a simple malaria infection
#'
#' @param xi a vector of quantiles
#' @param mu the expected value for log10 parasite densities
#' @param a host cohort age
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
#' @examples
qDensityPmu = function(xi, mu, a=0,
                       pRBC=par_lRBC.0(),
                       pSig=par_sigma.0()){
  # lRBC is the upper bound
  # xi is log10 parasite densities
  lRBC = log10RBC(a, pRBC)
  mu1 = mu/lRBC
  qbeta1(xi, mu1, pSig)*lRBC
}





#' The density function for parasite densities in a simple malaria infection of age alpha
#'
#' @param x a vector of quantiles
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
#' @examples
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

#' The distribution function for parasite densities in a simple malaria infection of age alpha
#'
#' @param n the number of observations
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to compute [sigma]
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
#' @examples
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

#' The distribution function for parasite densities in a simple malaria infection of age alpha
#'
#' @param x a vector of probabilities
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
#' @examples
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

#' The quantile function for parasite densities in a simple malaria infection of age alpha
#'
#' @param x a vector of quantiles
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
#' @examples
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

#' The density function for parasite densities from a host cohort population of age a
#'
#' @param x
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
#' @examples
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

#' The distribution function for parasite densities from a host cohort population of age a
#'
#' @param x a vector of probabilities
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param pMu parameters to compute [alpha2mu]
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
#' @examples
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


