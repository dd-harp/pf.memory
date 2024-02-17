
#' Compute \eqn{P_\tau(a |h)}
#'
#' @description
#' Compute \deqn{P_\tau(a| h) \sim f_P(\xi; a, \tau |h ) = \int_0^a \Omega(\xi|F_\mu(\alpha)) \; f_A(\alpha; a, \tau | h) d\alpha}
#'
#' @param xi log10
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length(x)
#' @export
d_Pdensity = function(xi, a, FoIpar, tau=0,
                      hhat=1, r=1/200,
                      par_RBC = par_lRBC_static(),
                      par_Fmu=par_Fmu_base(),
                      par_mu2dens = par_mu2dens_beta(),
                      pWda=par_Wda_none()){

  simpleDensities = function(x, a, FoIpar, tau, hhat, r, lRBC, m, W, par_Fmu, par_mu2dens){
    Pa = function(alpha, x, W, a, m){
      mu = alpha2mu(alpha, W, par_Fmu)
      y = zda(alpha, a, FoIpar, tau, hhat, r)
      d_mu2density(x, mu, par_mu2dens, lRBC)*y/m
    }
    lowlim = 8
    uplim = min(a, 6*365)
    dff = stats::integrate(dAoI, lowlim, a, a=a, FoIpar=FoIpar,  tau=tau, hhat=hhat, r=r)$value
    stats::integrate(Pa, lowlim, uplim, x=x, W=W, a=a, m=m)$value/dff
  }

  m = meanMoI(a, FoIpar, tau, hhat,  r)
  W = Wda(a, FoIpar, tau, hhat, pWda)
  lRBC = log10RBC(a, par_RBC)
  if(length(xi)==1) return(simpleDensities(xi, a, FoIpar, tau, hhat, r,
                                           lRBC, m, W, par_Fmu, par_mu2dens))
  return(sapply(xi, simpleDensities, a=a, FoIpar=FoIpar,  tau=tau, hhat=hhat, r=r,
                lRBC=lRBC, m=m, W=W, par_Fmu=par_Fmu, par_mu2dens=par_mu2dens))
}

#' Compute the moments of P_density
#'
#' @description
#' Compute \deqn{P_\tau(a| h) \sim f_P(\xi; a, \tau |h ) = \int_0^a \Omega(\xi|F_\mu(\alpha)) \; f_A(\alpha; a, \tau | h) d\alpha}
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param n the moment to compute
#' @param dt the mesh size over xi
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length(x)
#' @export
moment_Pdensity = function(a, FoIpar, tau=0, n=1, dt=0.1,
                      hhat=1, r=1/200,
                      par_RBC = par_lRBC_static(),
                      par_Fmu=par_Fmu_base(),
                      par_mu2dens = par_mu2dens_beta(),
                      pWda=par_Wda_none()){

  lRBC = log10RBC(a, par_RBC)
  xi = seq(0, lRBC, by = dt)
  Pd = d_Pdensity(xi, a, FoIpar, tau, hhat, r,
                  par_RBC,  par_Fmu, par_mu2dens, pWda)
  sum(xi^n*Pd*dt)
}
