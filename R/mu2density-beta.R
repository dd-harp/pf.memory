
#' Modified beta distribution, density function
#'
#' @description
#' The beta distribution, parameterized by the mean and variance, modified to return
#' a number between 0 and bvm (the number of red blood cells)
#'
#'
#' @inheritParams d_mu2density
#'
#' @return a [numeric] vector of length(xi)
#' @export
d_mu2density.beta = function(xi, mu, par_mu2dens=par_mu2dens_beta(), bvm=13){
  with(par_mu2dens,{
    xi1 = xi/bvm
    mu1 = mu/bvm
    dbeta1(xi1, mu1, pSig)/bvm
})}

#' Modified beta distribution, random numbers
#'
#' @description
#' The beta distribution, parameterized by the mean and variance, modified to return
#' a number between 0 and bvm (the number of red blood cells)
#'
#' @inheritParams r_mu2density
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
r_mu2density.beta = function(n, mu, par_mu2dens=par_mu2dens_beta(), bvm=13){
  with(par_mu2dens,{
    mu1 = mu/bvm
    rbeta1(n, mu1, pSig)*bvm
})}


#' Modified beta distribution, distribution function
#'
#' @inheritParams p_mu2density
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
p_mu2density.beta = function(xi, mu, par_mu2dens=par_mu2dens_beta(), bvm=13){
  with(par_mu2dens,{
    # bvm is the upper bound
    # xi is log10 parasite densities
    mu1 = mu/bvm
         xi1 = xi/bvm
       pbeta1(xi1, mu1, pSig)
})}

#' Modified beta distribution, distribution function
#'
#' @inheritParams q_mu2density
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
q_mu2density.beta = function(xi, mu, par_mu2dens=par_mu2dens_beta(), bvm=13){
  with(par_mu2dens,{
    mu1 = mu/bvm
    qbeta1(xi, mu1, pSig)*bvm
})}

#' The quantile function for parasite densities in a simple malaria infection
#'
#' @param par_sigma parameters to compute sigma_mu
#'
#' @return a compound [list]
#' @export
par_mu2dens_beta = function(par_sigma = par_sigma_abc()){
  par = list()
  class(par) <- "beta"
  par$pSig = par_sigma
  return(par)
}
