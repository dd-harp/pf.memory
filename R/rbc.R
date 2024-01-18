log10RBC = function(a,par){
  UseMethod("log10RBC", par)
}

log10RBC.0 = function(a,par){par$lRBCmax }

par_lRBC.0 = function(lRBCmax=13){
  par = list()
  class(par) <- "0"
  par$lRBCmax = lRBCmax
  par
}

