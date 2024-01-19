

#' Density function for the beta distribution, an alternative parameterization
#' @description
#' The beta distribution is parameterized using the mean and a function `sigma` that computes the variance
#' as a function of the mean
#'
#' @param x a vector of quantiles
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
dbeta1 = function(x, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  stats::dbeta(x, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

#' Disribution function for the beta distribution, an alternative parameterization
#'
#' @param p a vector of probabilities
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma]
#'
#' @return a [numeric] vector of length(p)
#' @export
#'
pbeta1 = function(p, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  stats::pbeta(p, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

#' The quantile function for the beta distribution, an alternative parameterization
#'
#' @param x a vector of quantiles
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma]
#'
#' @return a [numeric] vector of length(x)
#' @export
#'
qbeta1 = function(x, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  stats::qbeta(x, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}


#' The random generation function for the beta distribution, an alternative parameterization
#' Title
#'
#' @param n number of observations
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param pSig parameters to dispatch the S3 function [sigma]
#'
#' @return a [numeric] vector of length n
#' @export
#'
rbeta1 = function(n, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  stats::rbeta(n, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}


#' A function to compute the variance of the beta distrution as a function of the mean.
#'
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param par parameters to dispatch and configure the instances
#'
#' @return a [numeric] vector of length(mu)
#' @export
#'
sigma = function(mu, par){
  UseMethod("sigma", par)
}

#' A function that returns constrained values of the variance for the beta distrution as a function of the mean.
#'
#' @param mu the mean value for the distribution (0 <= mu <= 1)
#' @param par parameters to dispatch and configure the instances
#'
#' @return a [numeric] vector of length(mu)
#' @export
#'
sigma.0 = function(mu, par){with(par,{
  pmin(abs(cc)*mu^(1+abs(bb))*(1-mu)^(1+abs(aa)), mu*(1-mu))
})}

#' Parameters to configure [sigma.0]
#'
#' @param aa a shape parameter
#' @param bb a shape parameter
#' @param cc a shape parameter
#'
#' @return a [list]
#' @export
#'
par_sigma.0 = function(aa=3.11, bb=2.14, cc=1.74){
  par = list()
  class(par) <- "0"
  par$aa=aa
  par$bb=bb
  par$cc=cc
  par
}
