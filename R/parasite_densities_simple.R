
#' The density function for parasite densities in a simple malaria infection
#'
#' @param xi a vector of quantiles
#' @param mu the expected value for log10 parasite densities
#' @param a host cohort age
#' @param pRBC parameters to compute [log10RBC]
#' @param pSig parameters to dispatch [sigma_mu]
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
dDensityPmu = function(xi, mu, a=0,
                       pRBC=par_lRBC_static(),
                       pSig=par_sigma_abc()){
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
#' @param pSig parameters to dispatch [sigma_mu]
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
rDensityPmu = function(n, mu, a=0,
                       pRBC=par_lRBC_static(),
                       pSig=par_sigma_abc()){
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
#' @param pSig parameters to dispatch [sigma_mu]
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
pDensityPmu = function(xi, mu, a=0,
                       pRBC=par_lRBC_static(),
                       pSig=par_sigma_abc()){
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
#' @param pSig parameters to dispatch [sigma_mu]
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
qDensityPmu = function(xi, mu, a=0,
                       pRBC=par_lRBC_static(),
                       pSig=par_sigma_abc()){
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
#' @param pSig parameters to dispatch [sigma_mu]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
dDensityPalpha = function(x, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu_base(),
                          pRBC=par_lRBC_static(),
                          pSig=par_sigma_abc()){
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
#' @param pSig parameters to compute [sigma_mu]
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
rDensityPalpha = function(n, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu_base(),
                          pRBC=par_lRBC_static(),
                          pSig=par_sigma_abc()){
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
#' @param pSig parameters to dispatch [sigma_mu]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
pDensityPalpha = function(x, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu_base(),
                          pRBC=par_lRBC_static(),
                          pSig=par_sigma_abc()){
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
#' @param pSig parameters to dispatch [sigma_mu]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
qDensityPalpha = function(x, alpha,
                          W=0, a=0,
                          pMu=par_alpha2mu_base(),
                          pRBC=par_lRBC_static(),
                          pSig=par_sigma_abc()){
  # lRBC is the upper bound of realistic values
  # x is log10 parasite densities
  mu = alpha2mu(alpha, W, pMu)
  qDensityPmu(x, mu, a, pRBC, pSig)
}

#' The density function for parasite densities from a host cohort population of age a
#'
#' @param x host cohort age
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
#' @return a [numeric] vector of length(x)
#' @export
#'
dDensityPa = function(x, a, FoIpar,
                      hhat=NULL,
                      tau=0, r=1/200,
                      pMu=par_alpha2mu_base(),
                      pRBC=par_lRBC_static(),
                      pSig=par_sigma_abc(),
                      pWda=par_Wda_none()){
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
    dff = stats::integrate(dAoI, lowlim, a, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r)$value
    stats::integrate(px, lowlim, uplim, x=x, W=W, a=a, m=m)$value/dff
  }
  if(length(x)==1) return(dB(x, a, FoIpar, hhat, tau, r, pMu, pRBC, pSig, pWda))
  return(sapply(x, dB, a=a, FoIpar=FoIpar, hhat=hhat,  tau=tau, r=r, pMu=pMu, pRBC=pRBC, pSig=pSig, pWda=pWda))
}
