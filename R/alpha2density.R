
#' The density function for parasite densities in a simple malaria infection of age alpha
#'
#' @param xi a vector of quantiles
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [d_mu2density]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
d_alpha2density = function(xi, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu = par_Fmu_base(),
                           par_mu2dens = par_mu2dens_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = alpha2mu(alpha, W, par_Fmu)
  d_mu2density(xi, mu, par_mu2dens, bvm)
}

#' The distribution function for parasite densities in a simple malaria infection of age alpha
#'
#' @param n the number of observations
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [r_mu2density]
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
r_alpha2density = function(n, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu=par_Fmu_base(),
                           par_mu2dens = par_mu2dens_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = alpha2mu(alpha, W, par_Fmu)
  r_mu2density(n, mu, par_mu2dens, bvm)
}

#' The distribution function for parasite densities in a simple malaria infection of age alpha
#'
#' @param x a vector of probabilities
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [p_mu2density]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
p_alpha2density = function(x, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu = par_Fmu_base(),
                           par_mu2dens = par_mu2dens_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = alpha2mu(alpha, W, par_Fmu)
  p_mu2density(x, mu, par_mu2dens, bvm)
}

#' The quantile function for parasite densities in a simple malaria infection of age alpha
#'
#' @param x a vector of quantiles
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [alpha2mu]
#' @param par_mu2dens parameters to compute [q_mu2density]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
q_alpha2density = function(x, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu = par_Fmu_base(),
                           par_mu2dens = par_mu2dens_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = alpha2mu(alpha, W, par_Fmu)
  q_mu2density(x, mu, par_mu2dens, bvm)
}

