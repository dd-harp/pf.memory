
truePRa = function(a, FoIpar, hhat=NULL, tau=0, r=1/200){
  1-exp(-meanMoI(a, FoIpar, hhat, tau, r))
}



dpda = function(a,M,par,FoIpar){with(as.list(c(M,par)),{
  foi = FoI(a,FoIpar,tau,h)
  R = function(m){ ifelse(m==0,r,r*m/(exp(m)-1))}
  dp  = foi*(1-p) - R(m)*p
  dm  = foi - r*m
  list(c(dp, dm))
})}

solve_dpda = function(h, FoIpar, r=1/200, tau=0, Tmax=730, dt=1){
  tms = seq(0, Tmax, by = dt)
  prms = c(h=h,r=r,FoI=FoI,tau=tau)
  inits = c(p=0,m=0)
  data.frame(lsode(inits, times=tms, dpda, prms, FoIpar=FoIpar))
}
