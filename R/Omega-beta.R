
#' Modified beta distribution, density function
#'
#' @description
#' The beta distribution, parameterized by the mean and variance, modified to return
#' a number between 0 and bvm (the number of red blood cells)
#'
#'
#' @inheritParams d_Omega
#'
#' @return a [numeric] vector of length(xi)
#' @export
d_Omega.beta = function(xi, mu, bvm=13, par_Omega=par_Omega_beta()){
  with(par_Omega,{
    xi1 = xi/bvm
    mu1 = mu/bvm
    dbeta1(xi1, mu1, pSig)/bvm
})}


#' Density function for the beta distribution, an alternative parameterization
#' @description
#' The beta distribution is parameterized using the mean and a function `sigma_mu` that computes the variance
#' as a function of the mean
#'
#' @param x a vector of quantiles
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma_mu]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
dbeta1 = function(x, mu,
                  pSig=par_sigma_abc()){

  var = sigma_mu(mu, pSig)
  stats::dbeta(x, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

#' Modified beta distribution, random numbers
#'
#' @description
#' The beta distribution, parameterized by the mean and variance, modified to return
#' a number between 0 and bvm (the number of red blood cells)
#'
#' @inheritParams r_Omega
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
r_Omega.beta = function(n, mu, bvm=13, par_Omega=par_Omega_beta()){
  with(par_Omega,{
    mu1 = mu/bvm
    rbeta1(n, mu1, pSig)*bvm
})}


#' The random generation function for the beta distribution, an alternative parameterization
#' Title
#'
#' @param n number of observations
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma_mu]
#'
#' @return a [numeric] vector of length n
#' @export
#'
rbeta1 = function(n, mu,
                  pSig=par_sigma_abc()){

  var = sigma_mu(mu, pSig)
  stats::rbeta(n, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}


#' Modified beta distribution, distribution function
#'
#' @inheritParams p_Omega
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
p_Omega.beta = function(xi, mu, bvm=13, par_Omega=par_Omega_beta()){
  with(par_Omega,{
    # bvm is the upper bound
    # xi is log10 parasite densities
    mu1 = mu/bvm
         xi1 = xi/bvm
       pbeta1(xi1, mu1, pSig)
})}


#' Disribution function for the beta distribution, an alternative parameterization
#'
#' @param p a vector of probabilities
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma_mu]
#'
#' @return a [numeric] vector of length(p)
#' @export
#'
pbeta1 = function(p, mu,
                  pSig=par_sigma_abc()){

  var = sigma_mu(mu, pSig)
  stats::pbeta(p, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

#' Modified beta distribution, distribution function
#'
#' @inheritParams q_Omega
#'
#' @return a [numeric] vector of length(xi)
#' @export
#'
q_Omega.beta = function(xi, mu, bvm=13, par_Omega=par_Omega_beta()){
  with(par_Omega,{
    mu1 = mu/bvm
    qbeta1(xi, mu1, pSig)*bvm
})}


#' The quantile function for the beta distribution, an alternative parameterization
#'
#' @param x a vector of quantiles
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma_mu]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
qbeta1 = function(x, mu,
                  pSig=par_sigma_abc()){

  var = sigma_mu(mu, pSig)
  stats::qbeta(x, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

#' The quantile function for parasite densities in a simple malaria infection
#'
#' @param par_sigma parameters to compute sigma_mu
#'
#' @return a compound [list]
#' @export
par_Omega_beta = function(par_sigma = par_sigma_abc()){
  par = list()
  class(par) <- "beta"
  par$pSig = par_sigma
  return(par)
}

#' A function to compute the variance of the beta distrution as a function of the mean.
#'
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param par parameters to dispatch and configure the instances
#'
#' @return a [numeric] vector of length(mu)
#' @export
#'
sigma_mu = function(mu, par){
  UseMethod("sigma_mu", par)
}

#' A function that returns constrained values of the variance for the beta distrution as a function of the mean.
#'
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param par parameters to dispatch and configure the instances
#'
#' @return a [numeric] vector of length(mu)
#' @export
#'
sigma_mu.abc = function(mu, par){with(par,{
  pmin(abs(cc)*mu^(1+abs(bb))*(1-mu)^(1+abs(aa)), mu*(1-mu))
})}

#' Parameters to configure [sigma_mu.abc]
#'
#' @param aa a shape parameter
#' @param bb a shape parameter
#' @param cc a shape parameter
#'
#' @return a [list]
#' @export
#'
par_sigma_abc = function(aa=3.11, bb=2.14, cc=1.74){
  par = list()
  class(par) <- "abc"
  par$aa=aa
  par$bb=bb
  par$cc=cc
  par
}
