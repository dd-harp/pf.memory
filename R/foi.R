
#' Force of Infection, \eqn{h_\tau(a)}
#'
#' @description
#' Given a birth day, \eqn{\tau,} age is related
#' to time by the formula \eqn{t = a+\tau}.
#' The FoI is computed as \deqn{\bar h \; \omega(a) S(a+\tau) T(a+\tau)}
#'
#' @param a host age
#' @param FoIpar a compound [list] to compute \eqn{h_\tau(a)}
#' @param tau the host cohort's birthday
#'
#' @return [numeric]
#' @export
FoI = function(a, FoIpar, tau=0){
  with(FoIpar,{
    hbar*ageFoI(a,agePar)*
    seasonalFoI(a+tau,seasonPar)*
    trendFoI(a+tau,trendPar)
})}


