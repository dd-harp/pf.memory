
meanMoI = function(a, FoIpar, hhat=NULL, tau=0, r=1/200){
  moif = function(a, FoIpar, hhat, tau,r){
    integrate(zda, 0, a, a=a, FoIpar=FoIpar,hhat=hhat,tau=tau,r=r)$value
  }
  if(length(a)==1){return(moif(a,FoIpar,hhat,tau,r))} else{
    (return(sapply(a,moif,FoIpar=FoIpar,hhat=hhat,tau=tau,r=r)))}
}

dMoIda = function(a,M,p,FoIpar){with(as.list(c(M,p)),{
  foi = FoI(a,FoIpar,tau,h)
  i = 1:N
  m = i-1
  dM = 0*M-(foi + r*m)*M
  dM[-N] = dM[-N] + r*m[-1]*M[-1]
  dM[-1] = dM[-1] + foi*M[-N]
  list(c(dM))
})}

solveMMinfty = function(h,FoIpar,r=1/200,tau=0,Tmax=730, dt=1){
  tms = seq(0, Tmax, by = dt)
  N = round(max(4*h/r,20))
  prms = c(h=h,r=r,N=N,tau=tau)
  inits = rep(0,N)
  inits[1]=1
  lsode(inits, times=tms, dMoIda, prms, FoIpar=FoIpar) -> out
  time = out[,1]; moi = out[,-1]
  m = moi %*% c(0:(N-1))
  list(time=time, moi=moi, m=m)
}



MoIDistPlot = function(moi, t, clr1 = "red", withm = FALSE){
  N = dim(moi)[2]-2
  mm = 1:N -1
  plot(mm, moi[t,1:N +1], type="h", xlab = "MoI", ylab = expression(M[tau](a)), main = paste ("Age = ", t, "Days"))
  if(withm == TRUE){
    m = out$m[t]
    points(mm, dpois(mm, m), col = clr1)
  }
}



dmda = function(a,M,p,FoIpar){with(as.list(c(M,p)),{
  foi = FoI(a,FoIpar,tau,h)
  dm = foi - r*m
  list(c(dm))
})}

solve_dm = function(h, FoIpar, r=1/200, tau=0, Tmax=730, dt=1){
  tms = seq(0, Tmax, by = dt)
  prms = c(h=h,r=r,tau=tau)
  inits = c(m=0)
  data.frame(lsode(inits, times=tms, dmda, prms, FoIpar=FoIpar))
}

