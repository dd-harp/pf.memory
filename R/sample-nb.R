
#' Parameters for negative binomial sampling
#'
#' @param sz the negative binomial size parameter
#' @param s sample volume as log10 RBCs
#'
#' @return binary detection result
#' @export
#'
par_nb = function(sz=0.31, s=6){
  par = list()
  class(par) <- "nb"
  par$sz=sz
  par$s=6
  par
}

#' Detection of infection given parasitemia
#'
#' @inheritParams d_detect
#'
#' @return binary detection result
#' @export
d_detect.nb = function(xi, bvm=13, pars=par_nb()){with(pars,{
  1-stats::dnbinom(0, mu=10^(xi+s-bvm), size=sz)
})}

#' Negative binomial PDF for raw parasite counts
#'
#' @inheritParams d_counts
#'
#' @return binary detection result
#' @export
d_counts.nb = function(x, xi, bvm=13, pars=par_nb()){with(pars,{
  stats::dnbinom(x, mu=10^(xi+s-bvm), size=sz)
})}

#' Detection of infection given parasitemia
#'
#' @inheritParams p_counts
#'
#' @return binary detection result
#' @export
p_counts.nb = function(x, xi, bvm=13, pars=par_nb()){with(pars,{
  stats::pnbinom(x, mu=10^(xi+s-bvm), size=sz)
})}

