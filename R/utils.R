
nzPois = function(N, moi){
  n = 3*N/(1-dpois(0, moi))
  X = rpois(n, moi)
  X = X[which(X>0)]
  if(length(X)<N) print("error")
  X[1:N]
}
