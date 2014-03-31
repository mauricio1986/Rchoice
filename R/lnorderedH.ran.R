# Log-likelihood of Hierarchical Ordered with Random Parameters

lnorderedH.ran<-function(theta, y, X, S, ranp, R, correlation, link,
                         weights = NULL, haltons = NULL,
                         seed = 123, make.estb = FALSE,
                         ... )
{
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis
  )
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis
  )
  N    <- nrow(X)
  K    <- ncol(X)
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  Xa   <- X[, Vara, drop = F]                        
  Xc   <- X[, Varc, drop = F]
  M    <- ncol(S)
  
  m       <- unclass(sort(unique(y)))
  J       <- length(levels(y))
  m       <- as.matrix(m[2:(J - 1)], nrow = (J - 2))               
  y       <- unclass(y)
  delta   <- t(kronecker(m, t(rep(1, nrow(X)))))    
  deltaj  <- delta == y
  deltak  <- delta == y-1
  
  alpha    <- theta[1:(J-2)]
  kappa    <- c(-Inf, cumsum(c(0, exp(alpha))), +Inf)
  gamma    <- as.matrix(theta[((J - 2) + 1):((J - 2) + Kc)])
  betasr   <- theta[-c(1:((J - 2) + Kc))]
  beta.bar <- as.matrix(betasr[c(1:Ka)])
  phis     <- betasr[(Ka + 1):(Ka + Ka * M)]         
  Phi      <- matrix(phis, nrow = Ka)                   
  colnames(Phi) <-colnames(S)
  rownames(Phi) <-colnames(Xa)
  sigma    <- betasr[-c(1:(Ka + Ka * M))]          
  
  
  set.seed(seed)
  Omega <- make.draws(R * N, Ka, haltons) 
  
  ZB <- as.vector(crossprod(t(Xc), gamma)) 
  XB  <- matrix(NA, N, R)
  XaS <- matrix(NA, N, length(phis))
  colnames(XaS) <- names(phis)
  Br <- array(NA, dim = c(N, R, Ka))
  for (i in 1:N){
    beta.r    <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                            omega = Omega[,((i - 1) * R + 1):(i * R) , drop = FALSE], 
                            correlation = correlation, Pi = Phi, S = S[i, , drop = FALSE])
    XB[i, ]   <- crossprod(t(Xa[i, , drop = FALSE]), beta.r$br)
    XaS[i, ]  <- kronecker(t(S[i, , drop = FALSE]), t(Xa[i, , drop = FALSE]))
    Br[i, , ] <- t(beta.r$br)
  }
  
  index  <- ZB + XB
  eta1   <- kappa[y + 1] - index
  eta2   <- kappa[y] - index
  Pir    <- pfun(eta1) - pfun(eta2)
  Pir    <- ifelse(Pir <= 0, .Machine$double.eps, Pir)
  Pi     <- rowSums(Pir) / R
  lls    <- sum(weights * log(Pi))
  
  phi1     <- dfun(eta1) ; phi2 <- dfun(eta2)
  lambda   <- (phi2 - phi1) / Pir
  Qir      <- Pir / (Pi * R)
  eta      <- Qir * lambda
  
  gkappa <- vector(mode = "list", length = (J - 2))
  for (j in 1:(J - 2)){
    gkappa[[j]] <- matrix(NA, N, R)
    gkappa[[j]] <- drop(deltaj[ ,j]) * (phi1 / Pir) - drop(deltak[, j]) * (phi2 / Pir)
  }
  gkappa <- lapply(gkappa, function(x) x * Qir)
  gkappa <- lapply(gkappa, function(x) apply(x, 1, sum))
  gkappa <- Reduce(cbind, gkappa) %*% jacobian(alpha)
  
  
  dUdb   <- matrix(NA, N, Ka)
  dUdphi <- matrix(NA, N, length(phis))
  if (correlation){
    dUds <- matrix(NA, N, (0.5 * Ka * (Ka + 1))) 
  } else {
    dUds <- matrix(NA, N, Ka)  
  }
  for(i in 1:N){
    beta.r      <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                              omega = Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                              correlation = correlation, Pi = Phi, S = S[i, , drop = FALSE])
    dUdb[i, ]   <- tcrossprod(eta[i, ], beta.r$d.mu) 
    dUdphi[i, ] <- tcrossprod(eta[i, ], beta.r$d.pis)
    dUds[i, ]   <- tcrossprod(eta[i, ], beta.r$d.sigma)
  }
  
  if (correlation) {
    vecX <- c()
    for (i in 1:Ka){
      vecX <- c(vecX, i:Ka)
    }
    Xac <- Xa[ , vecX]
  } else {
    Xac <- Xa  
  }
  
  gbarfi  <- Xc  * rowSums(eta)
  gbarmi  <- Xa  * dUdb
  gbarphi <- XaS * dUdphi
  gbarvi  <- Xac * dUds
  
  gbari <- cbind(gkappa, gbarfi , gbarmi, gbarphi, gbarvi)
  colnames(gbari) <- names(theta)
  attr(lls, 'gradient') <- weights * gbari
  
  if (make.estb){ 
    b.ran  <- array(NA, dim = c(N, R, Ka))
    b.ran2 <- array(NA, dim = c(N, R, Ka))
    for(j in 1:Ka){
      b.ran[, ,j]  <- Br[, , j] * Qir
      b.ran2[, ,j] <- (Br[, , j]^2) * Qir
    }
    b.ran  <- apply(b.ran, c(1, 3), sum)
    b.ran2 <- apply(b.ran2, c(1, 3), sum)
    colnames(b.ran) <- colnames(Xa)
    
    sd.ran <- sqrt(b.ran2 - b.ran^2)
    colnames(sd.ran) <- colnames(Xa)
    attr(lls,'b.random') <- b.ran
    attr(lls,'sd.random') <- sd.ran
    attr(lls,'probabilities') <- Pi
  } 
  lls
}
