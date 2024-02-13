
#' The mean MoI in a host cohort of age \eqn{a}
#'
#' @description
#' The mean multiplicity of infection (MoI) is \deqn{m_\tau(a) = \int_0^a z_\tau(\alpha, a) d \alpha}
#'
#' @inheritParams FoI
#' @param hhat a local scaling parameter for the FoI
#' @param r the clearance rate for a simple infection
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
meanMoI = function(a, FoIpar, tau=0, hhat=1, r=1/200){
  moif = function(a, FoIpar, tau, hhat,r){
    stats::integrate(zda, 0, a, a=a, FoIpar=FoIpar, tau=tau, hhat=hhat, r=r)$value
  }
  if(length(a)==1){return(moif(a,FoIpar,tau,hhat,r))} else{
    (return(sapply(a,moif,FoIpar=FoIpar,tau=tau,hhat=hhat,r=r)))}
}

#' Derivatives for the queuing model \eqn{M/M/\infty}
#'
#' @description
#'
#' This queuing model \eqn{M/M/\infty} tracks the MoI in a cohort of humans
#' as it ages. It assumes a time- and age-dependent hazard rate for infection,
#' called the force of infection (FoI, \eqn{h_\tau(a)}). Infections do not affect
#' each other, and each one clears independently at the rate \eqn{r}.
#'
#' Let \eqn{\zeta_i} the fraction of the population with
#' MoI = i, then
#' \deqn{\frac{d\zeta_0}{da}= -h_\tau(a) \zeta_0 + r \zeta_1}
#' and for \eqn{i\geq 1}
#' \deqn{\frac{d\zeta_i}{da}= h_\tau(a) \left( \zeta_{i-1} - \zeta_i \right)  - ri \zeta_i + r(i+1)\zeta_{i+1}}
#'
#' This function computes the derivatives in a form that can be used by [deSolve::ode].
#'
#' @param a the host age
#' @param M the state variables
#' @param pars the parameters
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#'
#' @return the derivatives as a [list]
#' @seealso [solveMMinfty]
#' @export
#'
dMoIda = function(a, M, pars, FoIpar){with(as.list(c(M,pars)),{
  foi = h*FoI(a, FoIpar, tau)
  i = 1:N
  m = i-1
  dM = 0*M-(foi + r*m)*M
  dM[-N] = dM[-N] + r*m[-1]*M[-1]
  dM[-1] = dM[-1] + foi*M[-N]
  list(c(dM))
})}

#' Solve the queuing model \eqn{M/M/\infty}
#'
#' @description
#'
#' A wrapper to solve the queuing model \eqn{M/M/\infty} (see [dMoIda]).
#'
#' The function automatically sets the maximum MoI to be computed, and it
#' sets initial conditions. The equations are solved using [deSolve::ode] and returned at
#' regular intervals dt from age 0 up to Amax (in days).
#'
#' @param h the force of infection
#' @param FoIpar \eqn{h_\tau(a)}, a [list] formatted to compute [FoI]
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return a [list] with the orbits by name
#' @seealso [dMoIda]
#' @export
#'
solveMMinfty = function(h, FoIpar, r=1/200, tau=0, Amax=730, dt=1){
  tms = seq(0, Amax, by = dt)
  N = round(max(4*h/r,20))
  prms = c(h=h,r=r,N=N,tau=tau)
  inits = rep(0,N)
  inits[1]=1
  out = deSolve::ode(inits, times=tms, dMoIda, prms, FoIpar=FoIpar)
  time = out[,1]; moi = out[,-1]
  m = moi %*% c(0:(N-1))
  list(time=time, moi=moi, m=m)
}

#' Plot the output of [solveMMinfty]
#'
#' @param moi the mean moi
#' @param t the time
#' @param clr1 the color
#' @export
MoIDistPlot = function(moi, t, clr1 = "red"){
  N = dim(moi)[2]-2
  mm = 1:N -1
  plot(mm, moi[t,1:N +1], type="h", xlab = "MoI", ylab = expression(M[tau](a)), main = paste ("Age = ", t, "Days"))
}


