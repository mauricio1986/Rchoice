## logLik function for standard poisson model

lnpoisson<-function(theta, y, X,
                    weights = NULL, ...)
{
  if (is.null(weights)) weights <- 1
  index <- tcrossprod(X, t(theta))
  mu    <- pmin(exp(index), 700)
  pi    <- dpois(y, mu)
  pi    <- ifelse(pi == 0, .Machine$double.eps, pi) 
  ll    <- sum(weights * log(pi))
  
  G      <- as.vector(y - mu) * X 
  colnames(G) <- names(theta)
  attr(ll, 'gradient') <- weights * G 
  
  H       <- crossprod(drop(-mu) * X, X)
  colnames(H) <- rownames(H) <- names(theta)
  attr(ll,'hessian') <- H
  
  attr(ll,'probabilities') <- pi
  ll
}