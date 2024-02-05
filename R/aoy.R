
#' The distribution function for the age of the youngest infection (AoY)
#'
#' @param alpha the age of infection
#' @param a the age of a cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
pAoY = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  X = 1-exp(-moi)
  cdf = dAoI(alpha, a, FoIpar, hhat, tau, r)
  (1-exp(-moi*cdf))/X
}

#' Alternative method for computing the distribution function for the age of the youngest infection (AoY)
#' @description
#' This method computes the distribution function for the AoY by summing over AoI distributions
#'
#' @param alpha the age of infection
#' @param a the age of a cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
pAoY_long= function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  py = dAoI(alpha, a, FoIpar, hhat, tau, r)
  ix = max(15, 5*moi)
  aoy=py*0
  for(m in 1:ix){
    aoy = aoy+stats::dpois(m,moi)*(1-(1-py)^m)
  }
  aoy = aoy/truePRa(a,FoIpar,hhat,tau,r)
  return(aoy)
}


#' The density function for the age of the youngest infection (AoY)
#'
#' @param alpha the age of infection
#' @param a the age of a cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length(alpha)
#' @export
#'
dAoY = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  # The function call
  dAoYcompute = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
    moi = meanMoI(a,FoIpar,hhat,tau,r)
    moi*dAoI(alpha,a,FoIpar,hhat,tau,r)*exp(-moi*pAoI(alpha,a,FoIpar,hhat,tau,r))/truePRa(a,FoIpar,hhat,tau,r)
  }
  # Use sapply to call dAoYcompute multiple times
  if(length(alpha)==1) return(dAoYcompute(alpha, a,FoIpar,hhat,tau,r))
  return(sapply(alpha, dAoYcompute, a=a, FoIpar=FoIpar,hhat=hhat,tau=tau,r=r))
}


#' The random generation function for the age of the youngest infection (AoY)
#'
#' @param N the number of observations
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of the AoI to return
#'
#' @return a [numeric] vector of length(N)
#' @export
#'
rAoY = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200, alphamin=0){
  minit = function(i, alphas, iix, jix){
    min(alphas[iix[i]:jix[i]])
  }
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  hatm = nzPois(N, moi)
  Ny = sum(hatm)
  hatalpha = rAoI(Ny, a, FoIpar, hhat, tau, r, alphamin)
  jix = cumsum(hatm)
  iix = c(1,jix+1)
  sapply(1:N, minit, alphas=hatalpha, iix=iix, jix=jix)
}




#' Compute the moments for the AoI density function for a cohort of age a
#'
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param n the moment desired
#'
#' @return a [numeric] vector of length(a)
#' @export
#'
momentAoY = function(a, FoIpar, hhat=5/365, tau=0, r=1/200, n=1){

  ffAoYda = function(a, FoIpar, hhat, tau, r, n){
    ff = function(alpha, a, FoIpar, hhat, tau, r, n){
      alpha^n*dAoY(alpha, a, FoIpar, hhat, tau, r)
    }
    stats::integrate(ff, 0, a,a=a,FoIpar=FoIpar,hhat=hhat,r=r,tau=tau,n=n)$value
  }

  if(length(a)==1){return(ffAoYda(a, FoIpar, hhat, tau, r, n))} else{
    sapply(a, ffAoYda, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r, n=n)}
}

#' The distribution function for the youngest age of N infections
#'
#' @param N the MoI
#' @param a the age of a cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length a + 1
#' @export
#'
pAoYN = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  alpha = 0:a
  py = pAoI(alpha, a, FoIpar, hhat, tau, r)
  1-(1-py)^N
}

#' The density function for the youngest age of N infections
#'
#' @param N the number of infections
#' @param a the age of a cohort
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return a [numeric] vector of length a + 1
#' @export
#'
dAoYN = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  cdf = pAoYN(N, a, FoIpar, hhat, tau, r)
  pdf = diff(cdf)
  pdf/sum(pdf)
}

#' The random generation function for the age of the youngest of N infections
#'
#' @param R the number of observations
#' @param N the number of infections
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#' @param alphamin the minimum value of the AoI to return
#'
#' @return a [numeric] vector of length R
#' @export
#'
rAoYN = function(R, N, a, FoIpar, hhat=NULL, tau=0, r=1/200, alphamin=0){
  matrix(rAoI(R*N,a,FoIpar,hhat,tau,r,alphamin), nrow=N, ncol=R)
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
dAoYda = function(a,vars, pars,FoIpar){with(as.list(c(vars,pars)),{

  m0 = function(m){pmax(m,1e-7)}
  p = 1-exp(-m0(m))

  foi = FoI(a,FoIpar,tau,h)

  F2rm = function(m, n){
    r*(sum(dpois(2:n, m)/c(2:n)))
  }

  dm  = foi - r*m
  dx  = 1 - foi*x/m0(m)
  dy  = 1 - foi*y/p + F2rm(m, n)*x/p

  list(c(dm, dx, dy))
})}

#' Solve the system of differential equations to compute the approximate moments of the AoY over time.
#'
#' @param h the force of infection
#' @param FoIpar a FoI trace function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Tmax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#' @param n The number of terms to use in F(r,m)
#'
#' @return a [data.frame] describing the orbits
#' @export
#'
solve_dAoYda = function(h, FoIpar, r=1/200, tau=0, Tmax=730, dt=1, n=8){
  tms = seq(1, Tmax, by = dt)
  prms = c(h=h, r=r, tau=tau, n=n)
  inits = c(m=1e-7, x=0, y=0)
  data.frame(deSolve::ode(inits, times=tms, dAoYda, prms, FoIpar=FoIpar))
}

