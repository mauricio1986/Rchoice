lnbinary<-function(theta, y, X, link, 
                   weights = NULL, ...)
{
  if (is.null(weights)) weights <- 1
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis
                 )
  dfun <- switch(link,
                  "probit" = dnorm,
                  "logit"  = dlogis
                 )
  ddfun <- switch(link,
                  "logit" = function(x) (1 - 2 * pfun(x)) * pfun(x) * (1 - pfun(x)),
                  "probit"= function(x) -x * dnorm(x)
                 )  
  mill  <- function(x) dfun(x) / pfun(x)
  millh <- function(x) ddfun(x) / pfun(x) - (dfun(x) / pfun(x))^2
  
  index  <- tcrossprod(X, t(theta))
  q      <- 2 * y - 1 
  pi     <- pfun(q * index)
  pi     <- ifelse(pi <= 0, .Machine$double.eps, pi) 
  ll     <- sum(weights * log(pi))
  
  ## Gradient
  G <- as.vector(q * mill(q * index)) * X
  colnames(G) <- names(theta)
  attr(ll,'gradient') <- weights * G
  
  ## Hessian
  H <- crossprod(as.vector(millh(q * index)) * X, X)
  colnames(H) <- rownames(H) <- names(theta)
  attr(ll,'hessian') <- H
  
  attr(ll,'probabilities') <- pi
  ll
}