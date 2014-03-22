## logLik function for standard poisson model

#lnpoisson<-function(theta, y, X,
#                    weights = NULL, ...)
#{
#  if (is.null(weights)) weights <- rep(1, nrow(X))
#  index <- tcrossprod(X, t(theta))
#  mu    <- pmax(exp(index), .Machine$double.eps) #This avoid NaN
#  pi    <- dpois(y, lambda = mu)
#  if (all(pi > 0)) ll <- sum(weights * log(pi)) 
#  else  ll <-NA
#  
#  #Gradient
#  lambda <- y - mu
#  G      <- as.vector(lambda) * X  
#  if(all(pi > 0))  attr(ll, 'gradient') <- weights * G 
#  else attr(ll, 'gradient') <- rep(NA, length(theta))
#  
#  #Hessian
#  lambda2 <- as.vector((-1) * mu)
#  H       <- crossprod(lambda2 * X, X)
#  attr(ll,'hessian') <- H
#  
#  if (all(pi > 0)) attr(ll, 'hessian') <- H 
#  else attr(ll, 'hessian') <- matrix(NA, nrow = length(theta), ncol = length(theta))
#  
#  attr(ll,'probabilities') <- pi
#  ll
#}

lnpoisson<-function(theta, y, X,
                    weights = NULL, ...)
{
  if (is.null(weights)) weights <- rep(1, nrow(X))
  index <- tcrossprod(X, t(theta))
  mu    <- pmax(exp(index), .Machine$double.eps) # This avoid NA on log?
  pi    <- dpois(y, mu)
  pi    <- ifelse(pi == 0, .Machine$double.eps, pi) # Avoiding -Inf in log(pi)
  ll    <- sum(weights * log(pi))
  
  #Gradient
  lambda <- y - mu
  G      <- as.vector(lambda) * X  
  attr(ll, 'gradient') <- weights * G 
  
  #Hessian
  lambda2 <- as.vector((-1) * mu)
  H       <- crossprod(lambda2 * X, X)
  attr(ll,'hessian') <- H
  
  attr(ll,'probabilities') <- pi
  ll
}