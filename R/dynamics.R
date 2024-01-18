
zda = function(alpha, a, FoIpar, hhat=NULL, tau=0, r=1/200){
  FoI(a-alpha,FoIpar,tau,hhat)*exp(-r*alpha)
}
