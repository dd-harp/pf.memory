
#' The density function for the sum of two infections
#'
#' @param x host cohort age
#' @param a host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#' @param pWda parameters to dispatch [Wda]
#'
#' @return a [numeric] value
#' @export
#'
dDensityPaConvolve2 = function(x, a, FoIpar,
                               hhat=NULL, tau=0,  r=1/200,
                               par_RBC = par_lRBC_static(),
                               par_Fmu=par_Fmu_base(),
                               par_Omega = par_Omega_beta(),
                               pWda=par_Wda_none()){
  px = function(x, log10B, a, FoIpar, hhat, tau, r,
                par_RBC, par_Fmu, par_Omega, pWda){
    lB2 = log10(10^log10B - 10^x)
    d_clone_density(x,a,FoIpar, hhat, tau, r, par_RBC, par_Fmu, par_Omega, pWda)*d_clone_density(lB2,a, FoIpar, hhat, tau, r, par_RBC, par_Fmu, par_Omega, pWda)
  }
  stats::integrate(px, 0, x, log10B=x, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r, par_RBC=par_RBC, par_Fmu=par_Fmu, par_Omega=par_Omega, pWda=pWda)$value
}



#' The density function for the sum of two infections, method a
#'
#' @param meshX a mesh of host cohort ages
#' @param CDF1 the CDF for a complex distribution
#' @param CDF2 the CDF for a complex distribution
#'
#' @return a [numeric] vector of length(meshX)
#' @export
#'
cdfConvolve2b = function(meshX, CDF1, CDF2){
  cX = CDF1*0
  L = length(meshX)
  PDF1 = c(CDF1[1], diff(CDF1))
  PDF1 = PDF1/sum(PDF1)
  B1 = B2 = 10^meshX
  for(i in 1:L){
    for(j in 1:L){
      x = log10(B1[j] + B2)
      ix2 = which(x <= meshX[i])
      if(length(ix2>0)){
        cX[i] = cX[i] + PDF1[j]*CDF2[max(ix2)]
      }
    }
  }
  cX
}

#' The density function for the sum of two infections, method b
#'
#' @param meshX a mesh of host cohort ages
#' @param CDF1 the CDF for a complex distribution
#' @param CDF2 the CDF for a complex distribution
#'
#' @return a [numeric] vector of length(meshX)
#' @export
#'
cdfConvolve2a = function(meshX, CDF1, CDF2){
  cX = CDF1*0
  L = length(meshX)
  PDF1 = c(CDF1[1], diff(CDF1))
  PDF1 = PDF1/sum(PDF1)
  B1 = B2 = B = 10^meshX
  for(i in 1:L){
    p1=0
    ix1 = which(B1[i]>B)
    if(length(ix1>0)) p1=PDF1[min(ix1)]
    for(j in 1:L){
      x = log10(B1[j] + B2)
      ix2 = which(x > meshX[i])
      if(length(ix2>0)){
        cX[i] = cX[i] + PDF1[j]*CDF2[min(ix2)]
      }
    }
  }
  cX
}

#' The density function for the sum of two infections
#'
#' @param meshX a mesh of host cohort ages
#' @param CDF1 the CDF for a complex distribution
#' @param CDF2 the CDF for a complex distribution
#'
#' @return a [numeric] vector of length(meshX)
#' @export
#'
cdfConvolve2 = function(meshX, CDF1, CDF2){
  CDFa = cdfConvolve2a(meshX, CDF1, CDF2)
  CDFb = cdfConvolve2a(meshX, CDF1, CDF2)
  (CDFa+CDFb)/2
}
