

#' Compute the distribution of parasite counts for simple infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param bins a set of break points for computing counts
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param RBC_par parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param Cpar parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
pCountsPa = function(a, FoIpar, bins=c(1:5, 13), dx=0.1,
                     hhat=1,tau=0, r=1/200,
                     RBC_par = par_lRBC_static(),
                     Fmu_par = par_Fmu_base(),
                     mu_par = par_mu2dens_beta(),
                     pWda=par_Wda_none(),
                     Cpar = par_nbCounts()){
  lRBC = log10RBC(a, RBC_par)
  meshX = seq(0, lRBC, by=dx)
  mix = 1:length(meshX)
  Pda = d_Pdensity(meshX, a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, mu_par, pWda)
  Pda = Pda/sum(Pda)
  cdfP = function(xi, a, Cpar){
    pD = Detect(xi,a,Cpar)
    counts = pCounts(bins, xi, a, Cpar)
  }
  Bx=sapply(meshX, cdfP, a=a, Cpar=Cpar)
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
#' @param RBC_par parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param Cpar parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
dCountsPa = function(a, FoIpar,  bins=1, dx=0.1,
                     hhat=1,tau=0, r=1/200,
                     RBC_par = par_lRBC_static(),
                     Fmu_par=par_Fmu_base(),
                     mu_par = par_mu2dens_beta(),
                     pWda=par_Wda_none(),
                     Cpar = par_nbCounts()){
  DP = DetectPa(a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, mu_par, pWda, Cpar)
  PC = pCountsPa(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, mu_par, pWda, Cpar)
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
#' @param RBC_par parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param Cpar parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
pCountsBa = function(a, FoIpar, bins=1, dx=0.1,
                     hhat=1,tau=0, r=1/200,
                     RBC_par = par_lRBC_static(),
                     Fmu_par=par_Fmu_base(),
                     mu_par = par_mu2dens_beta(),
                     pWda=par_Wda_none(),
                     Cpar = par_nbCounts()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = DetectPa(a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, mu_par, pWda, Cpar)
  out= pCountsPa(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, mu_par, pWda, Cpar)
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
#' @param RBC_par parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param Cpar parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
dCountsBa = function(a, FoIpar, bins=1, dx=0.1,
                     hhat=1,tau=0, r=1/200,
                     RBC_par = par_lRBC_static(),
                     Fmu_par=par_Fmu_base(),
                     mu_par = par_mu2dens_beta(),
                     pWda=par_Wda_none(),
                     Cpar = par_nbCounts()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = DetectPa(a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, mu_par, pWda, Cpar)
  PC = pCountsBa(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, mu_par, pWda, Cpar)
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
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param Cpar parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
meanCountsPa = function(a, FoIpar, dx=0.1,
                        hhat=1,tau=0, r=1/200,
                        Fmu_par=par_Fmu_base(),
                        mu_par = par_mu2dens_beta(),
                        pWda=par_Wda_none(),
                        Cpar = par_nbCounts()){
  bins = seq(0,7, by=dx)
  PC = dCountsPa(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, mu_par, pWda, Cpar)
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
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param mu_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param Cpar parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
meanCountsBa = function(a, FoIpar, dx=0.1,
                        hhat=1,tau=0, r=1/200,
                        Fmu_par=par_Fmu_base(),
                        mu_par = par_mu2dens_beta(),
                        pWda=par_Wda_none(),
                        Cpar = par_nbCounts()){
  bins = seq(0,7, by = dx)
  PC = dCountsBa(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, mu_par, pWda, Cpar)
  sum(PC$bins*PC$pdf)
}
