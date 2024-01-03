#' Add a trend to the FoI trace function
#'
#' @param t the time
#' @param par a [list]
#'
#' @return [numeric]
#' @export
trendFoI = function(t, par){
  UseMethod("trendFoI", par)
}

#' The "no trend" function
#'
#' @inheritParams seasonalFoI
#'
#' @return [numeric]
#' @export
trendFoI.flat= function(t, par){return(1)}

#' Return a list that dispatches trendFoI.flat
#'
#' @return a [list]
#' @export
par_flatTrend = function(){
  pars = list()
  class(pars) <- "flat"
  pars
}
