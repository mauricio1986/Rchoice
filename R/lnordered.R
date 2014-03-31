## Log-Likelihood standard ordered model.
## Based on porl from MASS package

lnordered<-function(theta, y, X, link,
                    weights = NULL, ... )
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
  ## it is assumed that y is a "factor" class
  m       <- sort(unique(y))
  m       <- unclass(m)
  J       <- length(levels(y))
  ## Just J-2 alpha are estimated, so the gradient is just for 2:(J-1) cuts
  m       <- as.matrix(m[2:(J - 1)], nrow = (J - 2))               
  y       <- unclass(y)
  delta   <- t(kronecker(m, t(rep(1, nrow(X)))))    
  ## indicator for alternative and previous alternative
  deltaj  <- delta == y
  deltak  <- delta == y - 1
  
  alpha    <- theta[1:(J - 2)]
  kappa    <- c(-Inf, cumsum(c(0, exp(alpha))) , +Inf)
  beta     <- theta[-c(1:(J - 2))]
  index    <- tcrossprod(X, t(beta))
  eta1     <- kappa[y + 1] - index
  eta2     <- kappa[y] - index
  pi       <- pfun(eta1) - pfun(eta2)
  pi       <- ifelse(pi <= 0, .Machine$double.eps, pi) 
  ll       <- sum(weights * log(pi))
  
  ## Gradient of parameters
  phi1     <- dfun(eta1) ; phi2 <- dfun(eta2)
  lambda   <- drop((phi2 - phi1) / pi)
  gbeta    <- as.vector(lambda) * X
  
  ## Gradient of kappa
  lamkap   <- (deltaj * drop(phi1) - deltak * drop(phi2)) / as.vector(pi)
  gkappa   <- lamkap %*% jacobian(alpha)
  G        <- cbind(gkappa, gbeta)
  colnames(G) <- names(theta)
  attr(ll,'gradient') <- weights * G  
  
  attr(ll,'probabilities') <- pi
  ll
}