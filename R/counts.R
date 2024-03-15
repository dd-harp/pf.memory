
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
#' @param Fmu_par parameters to compute [Fmu]
#' @param Omega_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
p_clone_counts = function(a, FoIpar, bins=c(1:5, 13), dx=0.1,
                          hhat=1,tau=0, r=1/200,
                          RBC_par = par_lRBC_static(),
                          Fmu_par = par_Fmu_base(),
                          Omega_par = par_Omega_beta(),
                          pWda=par_Wda_none(),
                          par_sample = par_nb()){
  bvm = log10RBC(a, RBC_par)
  meshX = seq(0, bvm, by=dx)
  mix = 1:length(meshX)
  Pda = d_clone_density(meshX, a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, Omega_par, pWda)
  Pda = Pda/sum(Pda)
  cdfP = function(xi, bvm, par_sample){
    pD = d_detect(xi, bvm, par_sample)
    counts = p_counts(bins, xi, bvm, par_sample)
  }
  Bx=sapply(meshX, cdfP, bvm=bvm, par_sample=par_sample)
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
#' @param Fmu_par parameters to compute [Fmu]
#' @param Omega_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
d_clone_counts = function(a, FoIpar,  bins=1, dx=0.1,
                          hhat=1,tau=0, r=1/200,
                          RBC_par = par_lRBC_static(),
                          Fmu_par=par_Fmu_base(),
                          Omega_par = par_Omega_beta(),
                          pWda=par_Wda_none(),
                          par_sample = par_nb()){
  DP = d_clone_detect(a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, Omega_par, pWda, par_sample)
  PC = p_clone_counts(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, Omega_par, pWda, par_sample)
  list(bins=PC$bins, pdf=diff(c(DP, PC$cdf))/(1-DP), detect=DP, cdf=PC$cdf, fullpdf = c(DP,diff(c(DP, PC$cdf))))
}

#' Compute the mean parasite counts in simple infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param Fmu_par parameters to compute [Fmu]
#' @param Omega_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
mean_clone_counts = function(a, FoIpar, dx=0.1,
                             hhat=1,tau=0, r=1/200,
                             Fmu_par=par_Fmu_base(),
                             Omega_par = par_Omega_beta(),
                             pWda=par_Wda_none(),
                             par_sample = par_nb()){
  bins = seq(0,7, by=dx)
  PC = d_clone_counts(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, Omega_par, pWda, par_sample)
  sum(PC$bins*PC$pdf)
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
#' @param Fmu_par parameters to compute [Fmu]
#' @param Omega_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
p_parasite_counts = function(a, FoIpar, bins=1, dx=0.1,
                     hhat=1,tau=0, r=1/200,
                     RBC_par = par_lRBC_static(),
                     Fmu_par=par_Fmu_base(),
                     Omega_par = par_Omega_beta(),
                     pWda=par_Wda_none(),
                     par_sample = par_nb()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = d_parasite_detect(a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, Omega_par, pWda, par_sample)
  out= p_clone_counts(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, Omega_par, pWda, par_sample)
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
#' @param Fmu_par parameters to compute [Fmu]
#' @param Omega_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
d_parasite_counts = function(a, FoIpar, bins=1, dx=0.1,
                     hhat=1,tau=0, r=1/200,
                     RBC_par = par_lRBC_static(),
                     Fmu_par=par_Fmu_base(),
                     Omega_par = par_Omega_beta(),
                     pWda=par_Wda_none(),
                     par_sample = par_nb()){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  D = d_parasite_detect(a, FoIpar, hhat, tau, r, RBC_par, Fmu_par, Omega_par, pWda, par_sample)
  PC = p_parasite_counts(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, Omega_par, pWda, par_sample)
  DP = PC$DM
  list(bins=PC$bins, pdf=diff(c(DP, PC$cdf))/(1-DP), detect=DP, cdf=PC$cdf, fullpdf = c(DP, diff(c(DP, PC$cdf))))
}

#' Compute the mean parasite counts in complex infections
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param dx the width of the mesh for computing the CDF of parasite densities
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param Fmu_par parameters to compute [Fmu]
#' @param Omega_par parameters defining parasite densities as a function of the mu
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return a [list]
#' @export
#'
mean_parasite_counts = function(a, FoIpar, dx=0.1,
                        hhat=1,tau=0, r=1/200,
                        Fmu_par=par_Fmu_base(),
                        Omega_par = par_Omega_beta(),
                        pWda=par_Wda_none(),
                        par_sample = par_nb()){
  bins = seq(0,7, by = dx)
  PC = d_parasite_counts(a, FoIpar, bins, dx, hhat, tau, r, Fmu_par, Omega_par, pWda, par_sample)
  sum(PC$bins*PC$pdf)
}
