#' Functions to modify the FoI by age
#'
#' @description This method dispatches on `par`.
#' @param a host age
#' @param par a [list]
#'
#' @return [numeric]
#' @export
ageFoI = function(a, par){
  UseMethod("ageFoI", par)
}

#' A function that does not modify the FoI by age
#'
#' @inheritParams ageFoI
#'
#' @return [numeric]
#' @export
ageFoI.flat = function(a, par){
  1+0*a
}

#' Make a parameter list to dispatch ageFoI.flat
#'
#' @return [list]
#' @export
par_flatAge = function(){
  par = list()
  class(par) <- "flat"
  return(par)
}

#' A function that does not modify the FoI by age
#'
#' @inheritParams ageFoI
#'
#' @return [numeric]
#' @export
ageFoI.type2 = function(a, par){with(par,{
  aa = (a+shft)/365
  A*aa/(B+aa)
})}

#' Make a parameter list to dispatch ageFoI.type2
#'
#' @param shift a shape parameter
#' @param A a shape parameter
#' @param B the age weight for adults
#' @return [list]
#' @export
par_type2Age = function(shift=0.1, A=1.5, B=2){
  ageFoIpar0 = list()
  class(ageFoIpar0) <- "type2"
  ageFoIpar0$shft=shift
  ageFoIpar0$A = 1.5
  ageFoIpar0$B = 2
  ageFoIpar0
}
