alpha2mu = function(alpha, W, par){
  UseMethod("alpha2mu", par)
}

salpha2mu = function(alpha, W, par){
  sapply(alpha, alpha2mu, W=W, par=par)
}

alpha2mu.0 = function(alpha, W, par){with(par,{
  B = tildel + (tildeb-tildel)*exp(-(Sa*(alpha-D)))
  ix = which(alpha<D)
  if(length(ix>0)) B[ix] = tildel+(tildeb-tildel)*alpha[ix]/D
  #  ix = which(alpha<liver)
  #  if(length(ix>0)) B[ix] = sqrt(-1)
  B
})}

par_alpha2mu.0 = function(D=20, liver=7, tildeb=10.3, tildel=2, Sa=0.0033, family ="0"){
  par = list()
  class(par) <- "0"
  par$D=D
  par$liver=liver
  par$tildeb=tildeb
  par$tildel=tildel
  par$Sa=Sa
  par
}



alpha2mu.W = function(alpha, W, par){with(par,{
  B = tildel + (tildeb-tildel)*exp(-Sa*(alpha-D)-Sw*W)
  ix = which(alpha<=D)
  if(length(ix>0)) B[ix] = tildel+(tildeb-tildel)*exp(-Sw*W)*alpha[ix]/D
  # ix = which(alpha<=liver)
  #  if(length(ix>0)) B[ix] = sqrt(-1)
  B
})}

par_alpha2mu.W = function(D=20, liver=7, tildeb=10.3, tildel=2, Sa=0.0033, Sw=0.001){
  par = list()
  class(par) <- "W"
  par$D=D
  par$liver=liver
  par$tildeb=tildeb
  par$tildel=tildel
  par$Sa=Sa
  par$Sw=Sw
  par
}
