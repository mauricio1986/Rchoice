#Jacobian for ordered model
jacobian <- function(kappa){
  ## ij elemenet is dkappa_i/d alpha_j  
  k <- length(kappa)
  etheta <- exp(kappa)
  mat <- matrix(0 ,k, k)
  for (i in 1:k) mat[i:k, i] <- etheta[i]
  mat
}

make.time <- function(object){
  et <- object$time[3]
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  paste(h, "h:", m, "m:", s, "s", sep="")
}