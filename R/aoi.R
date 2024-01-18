


dAoI = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  dAoIcompute = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
    zda(alpha, a, FoIpar, hhat, tau, r)/meanMoI(a, FoIpar, hhat, tau, r)
  }

  if(length(alpha)==1){
    return(dAoIcompute(alpha, a, FoIpar, hhat, tau, r))
  }else{
    return(sapply(alpha, dAoIcompute, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r))
  }
}





pAoI = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  pAoIfunction = function(alpha, a, FoIpar, hhat, tau, r){
    integrate(dAoI,0,alpha,a=a,FoIpar=FoIpar,hhat=hhat,tau=tau,r=r)$value
  }
  if(length(alpha)==1) {return(pAoIfunction(alpha, a, FoIpar, hhat, tau, r))} else{
    return(sapply(alpha, pAoIfunction, a=a, FoIpar=FoIpar, hhat=hhat, tau=tau, r=r))}
}









rAoI = function(N, a, FoIpar, hhat=NULL, tau=0, r=1/200, alphamin=0){
  if(N==0) return(-1)
  alpha = alphamin:a
  scdf = pAoI(alpha, a, FoIpar, hhat, tau, r)
  pdf = diff(scdf)
  pdf / sum(pdf)
  sample(alpha[-length(alpha)], N, replace=T, prob=pdf)
}




momentAoI = function(a, FoIpar, hhat=5/365, tau=0, r=1/200, n=1){

  fAda = function(a, FoIpar, hhat=5/365, tau=0, r=1/200, n=1){
    ff = function(alpha, a, FoIpar, hhat, tau, r, n){
      alpha^n*zda(alpha, a, FoIpar, hhat, tau, r)
    }
    m =  meanMoI(a,FoIpar,hhat,tau,r)
    integrate(ff, 0, a,a=a,FoIpar=FoIpar,hhat=hhat,r=r,tau=tau,n=n)$value/m
  }

  if(length(a)==1){return(fAda(a, FoIpar, hhat, tau, r, n))} else{
    sapply(a, fAda,FoIpar=FoIpar,hhat=hhat,tau=tau, r=r, n=n)}
}









dAoIda = function(a,M,p,FoIpar){with(as.list(c(M,p)),{
  foi = FoI(a,FoIpar,tau,h)
  m0 = pmax(m,1e-6)
  x1 = M[2]
  xn = M[1+1:N]
  dm = foi - r*m
  dx1 = 1 - foi*x1/m0
  dxn = (2:N)*xn[-N] - foi*xn[-1]/m0
  list(c(dm, dx1, dxn))
})}

solve_dAoI = function(h, FoIpar, r=1/200, tau=0, Tmax=730, dt=1, N=3){
  tms = seq(0, Tmax, by = dt)
  prms = c(h=h,r=r,tau=tau,N=N)
  inits = c(m=0, xn = rep(0,N))
  data.frame(lsode(inits, times=tms, dAoIda, prms, FoIpar=FoIpar))
}









