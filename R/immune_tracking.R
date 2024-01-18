

Wda = function(a, FoIpar, hhat=NULL, tau=0, par=par_Wda.0()){
  UseMethod("Wda", par)
}

Wda.0 = function(a, FoIpar, hhat=NULL, tau=0,  par=par_Wda.0()){
  0*a
}

par_Wda.0 = function(){
  par = list()
  class(par) <- "0"
  return(par)
}

par_Wda.delta = function(delta=0.001){
  par = list()
  class(par) <- "delta"
  par$delta=delta
  return(par)
}

Wda.delta = function(a, FoIpar, hhat=NULL, tau=0, par=par_Wda.0()){with(par,{
  Wd = function(a,FoIpar,hhat,tau,delta){
    ff = function(s,a,FoIpar,hhat,tau,delta)  FoI(a-s,FoIpar,tau,hhat)*exp(-delta*(a-s))
    integrate(ff,0,a,a=a,FoIpar=FoIpar,hhat=hhat,tau=tau,delta=delta)$value
  }
  if(length(a)==1){
    return(Wd(a,FoIpar,hhat,tau,delta))
  } else {
    return(sapply(a,Wd,FoIpar=FoIpar,hhat=hhat,tau=tau,delta=delta))
  }
})}
