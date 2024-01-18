

#' Title
#'
#' @param x
#' @param mu
#' @param pSig
#'
#' @return
#' @export
#'
#' @examples
dbeta1 = function(x, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  dbeta(x, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

pbeta1 = function(x, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  pbeta(x, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

qbeta1 = function(x, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  qbeta(x, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}

rbeta1 = function(n, mu,
                  pSig=par_sigma.0()){

  var = sigma(mu, pSig)
  rbeta(n, mu*(mu*(1-mu)/var-1), (1-mu)*(mu*(1-mu)/var-1))
}


sigma = function(mu, par){
  UseMethod("sigma", par)
}

sigma.0 = function(mu, par){with(par,{
  pmin(abs(cc)*mu^(1+abs(bb))*(1-mu)^(1+abs(aa)), mu*(1-mu))
})}

par_sigma.0 = function(aa=3.11, bb=2.14, cc=1.74){
  par = list()
  class(par) <- "0"
  par$aa=aa
  par$bb=bb
  par$cc=cc
  par
}
