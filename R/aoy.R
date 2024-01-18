
#' The distribution
#'
#' @param alpha the age of infection
#' @param a the age of a cohort
#' @param FoIpar
#' @param hhat
#' @param tau the cohort birthday
#' @param r the clearance rate of a simple infection
#'
#' @return
#' @export
#'
#' @examples
pAoY = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  X = 1-exp(-moi)
  cdf = dAoI(alpha, a, FoIpar, hhat, tau, r)
  (1-exp(-moi*cdf))/X
}

pAoY_long= function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  moi = meanMoI(a, FoIpar, hhat, tau, r)
  py = dAoI(alpha, a, FoIpar, hhat, tau, r)
  ix = max(15, 5*moi)
  aoy=py*0
  for(m in 1:ix){
    aoy = aoy+dpois(m,moi)*(1-(1-py)^m)
  }
  aoy = aoy/truePRa(a,FoIpar,hhat,tau,r)
  return(aoy)
}



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





momentAoY = function(a, FoIpar, hhat=5/365, tau=0, r=1/200, n=1){

  ffAoYda = function(a, FoIpar, hhat, tau, r, n){
    ff = function(alpha, a, FoIpar, hhat, tau, r, n){
      alpha^n*dAoY(alpha, a, FoIpar, hhat, tau, r)
    }
    integrate(ff, 0, a,a=a,FoIpar=FoIpar,hhat=hhat,r=r,tau=tau,n=n)$value
  }

  if(length(a)==1){return(ffAoYda(a, FoIpar, hhat, tau, r, n))} else{
    sapply(a, ffAoYda, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r, n=n)}
}









pAoYN = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  alpha = 0:a
  py = pAoI(alpha, a, FoIpar, hhat, tau, r)
  1-(1-py)^N
}

dAoYN = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  cdf = pAoYN(N, a, FoIpar, hhat, tau, r)
  pdf = diff(cdf)
  pdf/sum(pdf)
}

rAoYN = function(R, N, a, FoIpar, hhat=NULL, tau=0, r=1/200, alphamin=0){
  matrix(rAoI(R*N,a,FoIpar,hhat,tau,r,alphamin), nrow=N, ncol=R)
}











dAoYda = function(a,M,par,FoIpar){with(as.list(c(M,par)),{
  foi = FoI(a,FoIpar,tau,h)

  R = function(m){ ifelse(m==0,r,r*m/(exp(m)-1))}
  frM = function(m){
    m^sqrt(2)/4 + m^3/3/factorial(3) + m^4/4/factorial(4) + m^5/5/factorial(5)
  }
  m0 = function(m){pmax(m,1e-6)}
  p0 = function(p){pmax(p,1e-6)}

  dp  = foi*(1-p) - R(m)*p
  dm  = foi - r*m
  dy  = 1 - foi/p0(p)*y + R(m)*frM(m)*y
  dx  = 1 - foi/m0(m)*x

  list(c(dp, dm, dy,dx))
})}

solve_dAoYda = function(h, FoIpar, r=1/200, tau=0, Tmax=730, dt=1){
  tms = seq(0, Tmax, by = dt)
  prms = c(h=h,r=r,FoI=FoI,tau=tau)
  inits = c(p=0,m=0,y=0,x=0)
  data.frame(lsode(inits, times=tms, dAoYda, prms, FoIpar=FoIpar))
}

