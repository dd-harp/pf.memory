
#' The density function for parasite densities in a simple malaria infection
#'
#' @param xi vector of quantiles for \eqn{\log_{10}} parasite densities
#' @param mu the expected value for \eqn{\log_{10}} parasite densities
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param par_Omega parameters to compute parasite densities as a function of mu
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
d_Omega = function(xi, mu, bvm=13, par_Omega=par_Omega_beta()){
  UseMethod("d_Omega", par_Omega)
}

#' Random generation of parasite densities from a simple malaria infection
#'
#' @param n the number of observations
#' @param mu the expected value for \eqn{\log_{10}} parasite densities
#' @param bvm blood volume as \eqn{\log_{10}} red blood cells
#' @param par_Omega parameters defining parasite densities as a function of the mu
#'
#' @return a [numeric] vector of length(n)
#' @export
#'
r_Omega = function(n, mu, bvm=13, par_Omega=par_Omega_beta()){
  UseMethod("r_Omega", par_Omega)
}

#' The density function for parasite densities as a function of the mean
#'
#' @param xi a vector of probabilities for \eqn{log_{10}} parasite densities
#' @param mu the expected value for \eqn{log_{10}} parasite densities
#' @param bvm blood volume as \eqn{log_{10}} red blood cells
#' @param par_Omega parameters defining parasite densities as a function of the mu
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
p_Omega = function(xi, mu, bvm=13, par_Omega=par_Omega_beta()){
  UseMethod("p_Omega", par_Omega)
}

#' The quantile function for parasite densities in a simple malaria infection
#'
#' @param xi a vector of quantiles
#' @param mu the expected value for \eqn{log_{10}} parasite densities
#' @param bvm blood volume as \eqn{log_{10}} red blood cells
#' @param par_Omega parameters defining parasite densities as a function of the mu
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
q_Omega = function(xi, mu, bvm=13, par_Omega=par_Omega_beta()){
  UseMethod("q_Omega", par_Omega)
}
