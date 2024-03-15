
#' Parameters for Poisson sampling
#' @param s sample volume as log10 RBCs
#' @return par a [list]
#' @export
par_pois = function(s=6){
  par = list()
  class(par) <- "pois"
  par$s=s
  par
}

#' Detection of infection given parasitemia
#'
#' @inheritParams d_detect
#'
#' @return binary detection result
#' @export
d_detect.pois = function(xi, bvm=13, pars=par_pois()){with(pars,{
  1-stats::dpois(0, 10^(xi+s-bvm))
})}

#' Detection of infection given parasitemia
#'
#' @inheritParams d_counts
#'
#' @return binary detection result
#' @export
d_counts.pois = function(x, xi, bvm=13, pars=par_pois()){with(pars,{
  stats::dpois(x, 10^(xi+s-bvm))
})}

#' Detection of infection given parasitemia
#'
#' @inheritParams p_counts
#'
#' @return binary detection result
#' @export
p_counts.pois = function(x, xi, bvm=13, pars=par_pois()){with(pars,{
  q = 10^(s-bvm)
  stats::ppois(x, q*10^xi)
})}

