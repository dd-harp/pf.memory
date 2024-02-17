
#' Compute \eqn{B_\tau(a | h)}
#'
#' @description
#' This function first computes the PDF for simple infections, [d_Pdensity], and
#' then computes the CDF for the full distribution, summed from MoI = 1 up to a large number
#' using convolutions.
#'
#'
#' @param meshX a mesh over parasite densities
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param RBC_par parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length meshX
#' @export
#'
Bdensity = function(meshX,
                    a, FoIpar, tau=0,
                    hhat=1, r=1/200,
                    RBC_par = par_lRBC_static(),
                    Fmu_par = par_Fmu_base(),
                    mu_par = par_mu2dens_beta(),
                    pWda = par_Wda_none()){
  PDFx = d_Pdensity(meshX, a, FoIpar, tau, hhat, r, RBC_par, Fmu_par, mu_par, pWda)
  CDFx = cumsum(PDFx)
  CDFx = CDFx/max(CDFx)
  PDFx = c(CDFx[1], diff(CDFx))
  PDFx = PDFx/sum(PDFx)

  moi = meanMoI(a, FoIpar, tau, hhat, r)
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



#' Compute the distribution function for \eqn{B_\tau(a)}
#'
#' @description
#' Call [Bdensity] and return the PDF:
#' \deqn{f_B(\xi; a, \tau |h) = \log_{10} \left( \sum_{M_\tau(a|h)>0} 10^{P_\tau(a |h)}\right)}
#'
#' @inheritParams Bdensity
#'
#' @return a [numeric] vector of length meshX
#' @export
#'
d_Bdensity = function(meshX,
                      a, FoIpar, tau=0, hhat=1, r=1/200,
                      RBC_par = par_lRBC_static(),
                      Fmu_par = par_Fmu_base(),
                      mu_par  = par_mu2dens_beta(),
                      pWda    = par_Wda_none()){
  Bdensity(meshX, a, FoIpar, tau, hhat, r, RBC_par, Fmu_par, mu_par, pWda)$PDFm
}


#' Call [Bdensity] and return the density vector
#'
#' @inheritParams Bdensity
#'
#' @return a numeric vector of length meshX
#' @export
#'
p_Bdensity = function(meshX, a, FoIpar, tau=0, hhat=1, r=1/200,
                      RBC_par = par_lRBC_static(),
                      Fmu_par = par_Fmu_base(),
                      mu_par = par_mu2dens_beta(),
                      pWda=par_Wda_none()){
  Bdensity(meshX, a, FoIpar, tau, hhat, r, RBC_par, Fmu_par, mu_par, pWda)$CDFm
}

#' Random generation for parasite densities in a host cohort
#'
#' @param N number of observations
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of alpha allowed
#' @param RBC_par parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] vector of length N
#' @export
#'
r_Bdensity = function(N, a, FoIpar, tau=0,
                hhat=1, r=1/200, alphamin=7,
                RBC_par = par_lRBC_static(),
                Fmu_par= par_Fmu_base(),
                mu_par = par_mu2dens_beta(),
                pWda=par_Wda_none()){
  moi = meanMoI(a, FoIpar, tau, hhat, r)
  W = Wda(a, FoIpar, tau, hhat, pWda)

  hatm = nzPois(N, moi)
  Ny = sum(hatm)
  hatalpha = rAoI(Ny, a, FoIpar, tau, hhat, r, alphamin)
  # their expected values
  hatmu = alpha2mu(hatalpha, W, Fmu_par)
  bvm = log10RBC(a, RBC_par)
  hatx = r_mu2density(Ny, hatmu, mu_par, bvm)
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
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of alpha allowed
#' @param RBC_par parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a R by M [matrix]
#' @export
#'
rRda = function(M, R, a, FoIpar, tau=0, hhat=1, r=1/200, alphamin=7,
                RBC_par = par_lRBC_static(),
                Fmu_par = par_Fmu_base(),
                mu_par = par_mu2dens_beta(),
                pWda=par_Wda_none()){
  W = Wda(a, FoIpar, tau, hhat, pWda)
  # MoI, excluding zeros
  Ny = R*M
  hatalpha = rAoI(Ny, a, FoIpar, tau, hhat, r, alphamin)
  # their expected values
  hatmu = alpha2mu(hatalpha, W, Fmu_par)
  bvm = log10RBC(a, RBC_par)
  hatx = r_mu2density(Ny, hatmu, mu_par, bvm)
  lRBC = 10^hatx
  matrix(lRBC,nrow=R, ncol=M)
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
moment_Bdensity = function(a, FoIpar, tau=0, n=1, dt=0.1,
                           hhat=1, r=1/200,
                           par_RBC = par_lRBC_static(),
                           par_Fmu=par_Fmu_base(),
                           par_mu2dens = par_mu2dens_beta(),
                           pWda=par_Wda_none()){

  lRBC = log10RBC(a, par_RBC)
  xi = seq(0, lRBC, by = dt)
  Bd = d_Bdensity(xi, a, FoIpar, tau, hhat, r,
                  par_RBC,  par_Fmu, par_mu2dens, pWda)
  sum(xi^n*Bd)
}
