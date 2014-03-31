# Log-likelihood of Hierarchical Binary with Random Parameters

lnlbinaryH.ran<-function(theta, y, X, S, ranp, R, link, seed, correlation, 
                   weights = NULL, haltons = NULL, make.estb = FALSE,
                   ...)
{
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis
  )
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis
  )
  mill <- function(x) dfun(x) / pfun(x)
  
  N    <- nrow(X)
  K    <- ncol(X)
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  M    <- ncol(S)
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  Xa   <- X[ , Vara, drop = F]                        
  Xc   <- X[ , Varc, drop = F]                        

  gamma    <- as.matrix(theta[1:Kc])               
  beta.bar <- as.matrix(theta[(Kc + 1):(Kc + Ka)])  
  phis     <- theta[(Kc + Ka + 1):(Kc + Ka + Ka * M)]         
  Phi      <- matrix(phis, nrow = Ka)                   
  colnames(Phi) <- colnames(S)
  rownames(Phi) <- colnames(Xa)
  sigma    <- theta[-c(1:(Kc + Ka + Ka * M))]           
  
  set.seed(seed)
  Omega <- make.draws(R * N, Ka, haltons) #Ka * N*R
  
  ZB <- as.vector(crossprod(t(Xc), gamma)) 
  XB  <- matrix(NA, N, R)  
  XaS <- matrix(NA, N, length(phis))
  colnames(XaS) <- names(phis)
  Br <- array(NA, dim = c(N, R, Ka))
  for (i in 1:N){
    beta.r   <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                          omega = Omega[ , ((i - 1) * R + 1):(i * R) , drop = FALSE], 
                          correlation = correlation, Pi = Phi, S = S[i , , drop = FALSE])
    XB[i, ]   <- crossprod(t(Xa[i, , drop = FALSE]), beta.r$br)
    XaS[i, ]  <- kronecker(t(S[i, , drop = FALSE]), t(Xa[i, , drop = FALSE]))
    Br[i, , ] <- t(beta.r$br) 
  }
  
  q <- 2 * y - 1
  index <- ZB + XB
  Pir   <- pfun(q * index)
  Pir   <- ifelse(Pir <= 0, .Machine$double.eps, Pir) 
  Pi    <- rowSums(Pir) / R
  lls   <- sum(weights * log(Pi))
  
  lambda <- q * mill(q * index)
  Qir    <- Pir / (Pi * R) 
  eta    <- Qir * lambda
  
  dUdb   <- matrix(NA, N, Ka)
  dUdphi <- matrix(NA, N, length(phis))
  if(correlation){
    dUds <- matrix(NA, N, (0.5 * Ka * (Ka + 1))) 
  }else{
    dUds <- matrix(NA, N, Ka)  
  }
  for(i in 1:N){
    beta.r      <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, 
                              omega = Omega[, ((i - 1) * R + 1):(i * R), drop = FALSE], 
                              correlation = correlation, Pi = Phi, S = S[i, , drop = FALSE])
    dUdb[i, ]   <- tcrossprod(eta[i, ], beta.r$d.mu)  
    dUdphi[i, ] <- tcrossprod(eta[i, ], beta.r$d.pis)
    dUds[i, ]   <- tcrossprod(eta[i, ], beta.r$d.sigma)
  }
  
  if(correlation){
    vecX <- c()
    for (i in 1:Ka){
      vecX <- c(vecX, i:Ka)
    }
    Xac <- Xa[ , vecX]
  }else{
    Xac <- Xa  
  }
  
  gbarfi  <- Xc  * rowSums(eta)
  gbarmi  <- Xa  * dUdb
  gbarphi <- XaS * dUdphi
  gbarvi  <- Xac * dUds
  
  gbari <- cbind(gbarfi, gbarmi, gbarphi, gbarvi)
  colnames(gbari) <- names(theta)
  attr(lls, 'gradient') <- weights * gbari
  
  if (make.estb){
    b.ran  <- array(NA, dim = c(N, R, Ka))
    b.ran2 <- array(NA, dim = c(N, R, Ka))
    for(j in 1:Ka){
      b.ran[, , j]  <- Br[, , j] * Qir
      b.ran2[, , j] <- (Br[, , j]^2) * Qir
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