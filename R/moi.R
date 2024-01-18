
#' Compute the mean MoI directly
#'
#' @param a the host cohort age
#' @param FoIpar parameters that define an FoI function
#' @param hhat a local scaling parameter for the FoI
#' @param tau the cohort birthday
#' @param r the clearance rate for a simple infection
#'
#' @return
#' @export
#'
#' @examples
meanMoI = function(a, FoIpar, hhat=NULL, tau=0, r=1/200){
  moif = function(a, FoIpar, hhat, tau,r){
    integrate(zda, 0, a, a=a, FoIpar=FoIpar,hhat=hhat,tau=tau,r=r)$value
  }
  if(length(a)==1){return(moif(a,FoIpar,hhat,tau,r))} else{
    (return(sapply(a,moif,FoIpar=FoIpar,hhat=hhat,tau=tau,r=r)))}
}

#' Compute the first derivatives for the queuing model M/M/infinity
#'
#' @param a the host age
#' @param M the state variables
#' @param p the parameters
#' @param FoIpar parameters that define an FoI function
#'
#' @return
#' @export
#'
#' @examples
dMoIda = function(a,M,p,FoIpar){with(as.list(c(M,p)),{
  foi = FoI(a,FoIpar,tau,h)
  i = 1:N
  m = i-1
  dM = 0*M-(foi + r*m)*M
  dM[-N] = dM[-N] + r*m[-1]*M[-1]
  dM[-1] = dM[-1] + foi*M[-N]
  list(c(dM))
})}

#' Solve the queuing model M/M/infinity
#'
#' @param h the force of infection
#' @param FoIpar a FoI trace function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Tmax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return
#' @export
#'
#' @examples
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



#' Plot the distribution of the MoI at time t
#'
#' @param moi the mean moi
#' @param t the time
#' @param clr1 the color
#' @param withm
#'
#' @return
#' @export
#'
#' @examples
MoIDistPlot = function(moi, t, clr1 = "red", withm = FALSE){
  N = dim(moi)[2]-2
  mm = 1:N -1
  plot(mm, moi[t,1:N +1], type="h", xlab = "MoI", ylab = expression(M[tau](a)), main = paste ("Age = ", t, "Days"))
  if(withm == TRUE){
    m = out$m[t]
    points(mm, dpois(mm, m), col = clr1)
  }
}

#' Compute the derivatives for MoI using a hybrid model
#'
#' @param a the host age
#' @param M the state variables
#' @param p the parameters
#' @param FoIpar parameters that define an FoI function
#'
#' @return
#' @export
#'
#' @examples
dmda = function(a,M,p,FoIpar){with(as.list(c(M,p)),{
  foi = FoI(a,FoIpar,tau,h)
  dm = foi - r*m
  list(c(dm))
})}

#' Solve the hybrid model for the MoI
#'
#' @param h the force of infection
#' @param FoIpar a FoI trace function
#' @param r the clearance rate for a simple infection
#' @param tau the cohort birthday
#' @param Tmax The maximum runtime (in days)
#' @param dt The output frequency (in days)
#'
#' @return a [matrix] with the orbits
#' @export
#'
#' @examples
solve_dm = function(h, FoIpar, r=1/200, tau=0, Tmax=730, dt=1){
  tms = seq(0, Tmax, by = dt)
  prms = c(h=h,r=r,tau=tau)
  inits = c(m=0)
  data.frame(lsode(inits, times=tms, dmda, prms, FoIpar=FoIpar))
}

