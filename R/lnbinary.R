lnbinary<-function(theta, y, X, link, 
                   weights = NULL, ...)
{
  if (is.null(weights)) weights <- rep(1, nrow(X))
  index <- tcrossprod(X,t(theta))
  qi    <- 2 * y - 1
  
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  qfun <- switch(link,
                 "probit" = qnorm,
                 "logit"  = qlogis)
  ddfun <- switch(link,
                  "logit" = function(x) (1 - 2 * plogis(x)) * plogis(x) * (1 - plogis(x)),
                  "probit"= function(x) -x * dnorm(x))
  
  thresh <- -qfun(.Machine$double.eps)
  eta    <- pmin(pmax(index, -thresh), thresh)
  pi     <- pfun(qi * eta)
  ll     <- sum(weights * log(pi))
  
  #Gradient
  dens   <- pmax(dfun(qi * index), .Machine$double.eps)
  lambda <- qi * (dens / pi)
  G <- as.vector(lambda) * X
  colnames(G) <- names(theta)
  attr(ll,'gradient') <- weights * G
  
  #Hessian
  ddens   <- ddfun(qi * index)
  lambda2 <- ((ddens / pi)-(dens / pi) ^ 2) * qi^2
  H <- crossprod(as.vector(lambda2) * X, X)
  colnames(H) <- rownames(H) <- names(theta)
  attr(ll,'hessian') <- H
  
  attr(ll,'probabilities') <- pi
  ll
}