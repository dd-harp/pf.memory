

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return a [numeric] vector of length(a)
DetectPa = function(a, FoIpar, tau=0,
                    hhat=1, r=1/200,
                    par_RBC = par_lRBC_static(),
                    Fmu_par=par_Fmu_base(),
                    par_mu2dens = par_mu2dens_beta(),
                    pWda=par_Wda_none(),
                    pC = par_nbCounts()){
  pD = function(a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_mu2dens, pWda, pC){
    Dx = function(x, a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, p, pWda, pC){
      d_Pdensity(x, a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_mu2dens, pWda)*Detect(x,a,pC)
    }
    hatb=log10RBC(a,par_mu2dens$pRBC)
    stats::integrate(Dx, 0, hatb, a=a, FoIpar=FoIpar, hhat=hhat,tau=tau,r=r,par_RBC, Fmu_par=par_RBC, Fmu_par,par_mu2dens=par_mu2dens, pWda=pWda, pC=pC)$value
  }
  if(length(a)==1) return(pD(a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_mu2dens, pWda, pC))
  return (sapply(a, pD, FoIpar=FoIpar, hhat=tau, hhat=tau, r=r, par_RBC, Fmu_par=par_RBC, Fmu_par, par_mu2dens=par_mu2dens, pWda=pWda, pC=pC))
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return binary detection result
#' @export
DetectPM = function(a, FoIpar,tau=0,
                    hhat=1,r=1/200,
                    par_RBC = par_lRBC_static(),
                    Fmu_par=par_Fmu_base(),
                    par_mu2dens = par_mu2dens_beta(),
                    pWda=par_Wda_none(),
                    pC = par_nbCounts()){
  moi = meanMoI(a, FoIpar, tau, hhat, r)
  D = DetectPa(a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_mu2dens, pWda, pC)
  1 - exp(-moi*D)
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
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return binary detection result
#' @export
DetectBda = function(a, FoIpar, dx=0.1,tau=0,
                     hhat=1,r=1/200,
                     par_RBC = par_lRBC_static(),
                     Fmu_par = par_Fmu_base(),
                     par_mu2dens = par_mu2dens_beta(),
                     pWda=par_Wda_none(),
                     pC = par_nbCounts()){

  lRBC = log10RBC(a, par_mu2dens$pRBC)
  meshX = seq(0, lRBC, by = dx)
  Bx = d_Bdensity(meshX, a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_mu2dens, pWda)$PDFm
  t(Detect(meshX, a, pC)) -> Bs
  sum(Bx*Bs)
}

#' Detection of infection given parasitemia
#'
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param tau the cohort birthday
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param Fmu_par parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [d_mu2density]
#' @param pWda parameters to dispatch [Wda]
#' @param pC parameters that define a detection function
#'
#' @return detection probability
#' @export
MoIPM = function(a, FoIpar,tau=0,
                 hhat=1,r=1/200,
                 par_RBC = par_lRBC_static(),
                 Fmu_par = par_Fmu_base(),
                 par_mu2dens = par_mu2dens_beta(),
                 pWda=par_Wda_none(),
                 pC = par_nbCounts()){

  moi = meanMoI(a, FoIpar, tau, hhat, r)
  D = DetectPa(a, FoIpar, tau, hhat, r, par_RBC, Fmu_par, par_mu2dens, pWda, pC)
  moi*D/(1-exp(-moi*D))
}
