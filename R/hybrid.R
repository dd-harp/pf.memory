
#' Compute the derivatives for MoI using a hybrid model
#'
#' @param a the host age
#' @param M the state variables
#' @param p the parameters
#' @param FoIpar parameters that define an FoI function
#'
#' @return the derivatives as a [list]
#' @export
dmda = function(a,M,p,FoIpar){with(as.list(c(M,p)),{
  foi = h*FoI(a,FoIpar,tau)
  dm = foi - r*m
  list(c(dm))
})}

#' Solve the hybrid model for the MoI
#'
#' @param h the average, annual force of infection
#' @param FoIpar a FoI trace function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return a [matrix] with the orbits
#' @export
#'
solve_dm = function(h, FoIpar, r=1/200, tau=0, Amax=730, dt=1){
  tms = seq(0, Amax, by = dt)
  prms = c(h=h,r=r,tau=tau)
  inits = c(m=0)
  data.frame(deSolve::ode(inits, times=tms, dmda, prms, FoIpar=FoIpar))
}



#' Compute the derivatives for MoI and true PR
#'
#' @param a the host age
#' @param M state variables
#' @param par the model parameters
#' @param FoIpar parameters that define an FoI function
#'
#' @return the derivatives, as a [list]
#' @export
#'
dpda = function(a,M,par,FoIpar){with(as.list(c(M,par)),{
  foi = h*FoI(a,FoIpar,tau)
  R = function(m){ ifelse(m==0,r,r*m/(exp(m)-1))}
  dp  = foi*(1-p) - R(m)*p
  dm  = foi - r*m
  list(c(dp, dm))
})}

#' Solve a system of differential equations to compute the true PR and the MoI
#'
#' @param h a scaling factor on the FoI
#' @param FoIpar parameters that define an FoI function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return a [data.frame] with the orbits
#' @export
#'
solve_dpda = function(h, FoIpar, r=1/200, tau=0, Amax=730, dt=1){
  tms = seq(0, Amax, by = dt)
  prms = c(h=h, r=r, tau=tau)
  inits = c(p=0, m=0)
  data.frame(deSolve::ode(inits, times=tms, dpda, prms, FoIpar=FoIpar))
}


#' Compute the derivatives of the AoI moments dynamically
#' @description The
#'
#' @param a the host age
#' @param M the state variables
#' @param p the parameters
#' @param FoIpar parameters that define an FoI function
#'
#' @return the derivatives of the MoI and the AoI as a [list]
#' @export
#'
dAoIda = function(a,M,p,FoIpar){with(as.list(c(M,p)),{
  foi = h*FoI(a,FoIpar,tau)
  m0 = pmax(m,1e-6)
  x1 = M[2]
  xn = M[1+1:N]
  dm = foi - r*m
  dx1 = 1 - foi*x1/m0
  dxn = (2:N)*xn[-N] - foi*xn[-1]/m0
  list(c(dm, dx1, dxn))
})}

#' Solve the system of differential equations to compute the moments of the AoI over time.
#'
#' @param h the force of infection
#' @param FoIpar a FoI trace function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#' @param N The total number of moments to compute
#'
#' @return a data.frame with the orbits
#' @export
#'
solve_dAoI = function(h, FoIpar, r=1/200, tau=0, Amax=730, dt=1, N=3){
  stopifnot(N>2)
  tms = seq(0, Amax, by = dt)
  prms = c(h=h, r=r, tau=tau, N=N)
  inits = c(m=0, xn = rep(0,N))
  data.frame(deSolve::ode(inits, times=tms, dAoIda, prms, FoIpar=FoIpar))
}



#' Compute the derivatives of the approximate moments of the AoY dynamically
#' @description The
#'
#' @param a the host age
#' @param vars the state variables
#' @param pars the parameters
#' @param FoIpar parameters that define an FoI function
#'
#' @return the derivatives of the MoI and the AoI as a [list]
#' @export
#'
dAoYda = function(a, vars, pars, FoIpar){with(as.list(c(vars,pars)),{

  m0 = function(m){pmax(m,1e-7)}
  p = 1-exp(-m0(m))

  foi = h*FoI(a,FoIpar,tau)

  F2rm = function(m, n){
    r*(sum(dpois(2:n, m)/c(2:n)))
  }

  dm  = foi - r*m
  dx  = 1 - foi*x/m0(m)
  dy  = 1 - foi*y/p + F2rm(m, n)*x

  list(c(dm, dx, dy))
})}

#' Solve the system of differential equations to compute the approximate moments of the AoY over time.
#'
#' @param h the force of infection
#' @param FoIpar a FoI trace function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Amax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#' @param n The number of terms to use in \eqn{\phi(r,m)}
#'
#' @return a [data.frame] describing the orbits
#' @export
solve_dAoYda = function(h, FoIpar, r=1/200, tau=0, Amax=730, dt=1, n=8){
  tms = seq(0, Amax, by = dt)
  prms = c(h=h, r=r, tau=tau, n=n)
  inits = c(m=1e-8, x=0, y=0)
  data.frame(deSolve::ode(inits, times=tms, dAoYda, prms, FoIpar=FoIpar))
}


