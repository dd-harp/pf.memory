
#' The density function for parasite densities in a simple malaria infection of age alpha
#'
#' @param xi a vector of quantiles
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [Fmu]
#' @param par_Omega parameters to compute [d_Omega]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
d_alpha2density = function(xi, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu = par_Fmu_base(),
                           par_Omega = par_Omega_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = Fmu(alpha, W, par_Fmu)
  d_Omega(xi, mu, bvm, par_Omega)
}

#' The distribution function for parasite densities in a simple malaria infection of age alpha
#'
#' @param n the number of observations
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [Fmu]
#' @param par_Omega parameters to compute [r_Omega]
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
r_alpha2density = function(n, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu=par_Fmu_base(),
                           par_Omega = par_Omega_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = Fmu(alpha, W, par_Fmu)
  r_Omega(n, mu, bvm, par_Omega)
}

#' The distribution function for parasite densities in a simple malaria infection of age alpha
#'
#' @param x a vector of probabilities
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [Fmu]
#' @param par_Omega parameters to compute [p_Omega]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
p_alpha2density = function(x, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu = par_Fmu_base(),
                           par_Omega = par_Omega_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = Fmu(alpha, W, par_Fmu)
  p_Omega(x, mu, bvm, par_Omega)
}

#' The quantile function for parasite densities in a simple malaria infection of age alpha
#'
#' @param x a vector of quantiles
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param a host cohort age
#' @param par_RBC parameters to compute [log10RBC]
#' @param par_Fmu parameters to compute [Fmu]
#' @param par_Omega parameters to compute [q_Omega]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
q_alpha2density = function(x, alpha, W=0, a=0,
                           par_RBC = par_lRBC_static(),
                           par_Fmu = par_Fmu_base(),
                           par_Omega = par_Omega_beta()){
  bvm = log10RBC(a, par_RBC)
  mu = Fmu(alpha, W, par_Fmu)
  q_Omega(x, mu, bvm, par_Omega)
}

