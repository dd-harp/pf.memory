#' Compute mean, expected parasite densities `mu` as a function of the age of infection `alpha`
#'
#' @param alpha the age of a parasite infection
#' @param W the immune tracking variables
#' @param par_Fmu a [list] that defines a model
#'
#' @return mean log10 parasite densities, a [numeric] vector of length(alpha)
#' @export
#'
Fmu = function(alpha, W, par_Fmu){
  UseMethod("Fmu", par_Fmu)
}

#' Compute mean, expected parasite densities `mu` for multiple values of the age of infection `alpha`
#'
#' @param alpha a parasite's age of infection
#' @param W immune tracking variables
#' @param par_Fmu a [list] that defines a model
#'
#' @return mean log10 parasite densities, a [numeric] vector of length(alpha)
#' @export
#'
sFmu = function(alpha, W, par_Fmu){
  sapply(alpha, Fmu, W=W, par_Fmu=par_Fmu)
}


