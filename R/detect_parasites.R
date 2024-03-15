#' Detection of infection given parasitemia
#'
#' @param xi_mesh a mesh over \eqn{\log_{10}} parasite densities
#' @param xi_density the density distribution for parasites over [xi_mesh]
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param par_sample parameters that define a detection function
#'
#' @return the fraction of infected hosts that would test positive
#' @export
d_detect_mesh = function(xi_mesh, xi_density, bvm=13, par_sample = par_nb()){
  d_detect(xi_mesh, bvm, par_sample) -> pr_pos
  sum(pr_pos*xi_density)
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param dx width of the mesh
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_parasite_detect = function(a, FoIpar, dx=0.1, tau=0,
                               hhat=1,r=1/200,
                               par_RBC = par_lRBC_static(),
                               Fmu_par = par_Fmu_base(),
                               par_Omega = par_Omega_beta(),
                               pWda=par_Wda_none(),
                               par_sample = par_nb()){

  lRBC = log10RBC(a, par_RBC)
  xi_mesh = seq(0, lRBC, by = dx)
  Bx = d_parasite_density(xi_mesh, a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_Omega, pWda)
  d_detect_mesh(xi_mesh, Bx, lRBC, par_sample)
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return binary detection result
#' @export
d_parasite_detect_moi = function(a, FoIpar,tau=0,
                                  hhat=1,r=1/200,
                                  par_RBC = par_lRBC_static(),
                                  Fmu_par = par_Fmu_base(),
                                  par_Omega = par_Omega_beta(),
                                  pWda = par_Wda_none(),
                                  par_sample = par_nb()){
  moi = meanMoI(a, FoIpar, tau, hhat, r)
  D = cohort_detect_clone(a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_Omega, pWda, par_sample)
  1 - exp(-moi*D)
}


#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#' @param par_sample parameters that define a detection function
#'
#' @return detection probability
#' @export
d_moi_count = function(a, FoIpar,tau=0,
                       hhat=1,r=1/200,
                       par_RBC = par_lRBC_static(),
                       Fmu_par = par_Fmu_base(),
                       par_Omega = par_Omega_beta(),
                       pWda=par_Wda_none(),
                       par_sample = par_nb()){

  moi = meanMoI(a, FoIpar, tau, hhat, r)
  D = detect_clone_by_cohort(a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_Omega, pWda, par_sample)
  moi*D/(1-exp(-moi*D))
}
