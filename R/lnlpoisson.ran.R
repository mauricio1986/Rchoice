## Log-Likelihood function for poisson model with random parameters

lnlpoisson.ran <- function(theta, y, X, ranp, R, correlation, link,
                        weights = NULL, haltons = NULL,
                        seed = 123, make.estb = FALSE,
                        ... )
{
  N    <- nrow(X)
  K    <- ncol(X)
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  Xa   <- X[ , Vara, drop = F]                        
  Xc   <- X[ , Varc, drop = F]                        
  
  gamma    <- as.matrix(theta[1:Kc])               
  beta.bar <- as.matrix(theta[(Kc + 1):(Kc + Ka)])  
  sigma    <- theta[-c(1:(Kc + Ka))]                 
  
  set.seed(seed)
  Omega    <- make.draws(R * N, Ka, haltons) 
  
  ZB <- as.vector(crossprod(t(Xc), gamma))
  
  XB <- matrix(NA, N, R)  
  Br <- array(NA, dim = c(N, R, Ka))
  for (i in 1:N){
    beta.r    <- Make.rcoef(beta.bar, sigma, ranp, Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE], 
                          correlation, Pi = NULL, S = NULL)
    XB[i, ]   <- crossprod(t(Xa[i, , drop = FALSE]), beta.r$br)
    Br[i, , ] <- t(beta.r$br) 
  }
  
  index <- ZB + XB
  # Note that exp(709) = Inf ==> gradient = Inf
  mu    <- pmin(exp(index), 700)
  Pir   <- dpois(y, lambda = mu)
  Pir   <- ifelse(Pir <= 0, .Machine$double.eps, Pir)
  Pi    <- rowSums(Pir) / R
  lls   <- sum(weights * log(Pi))
  
  lambda <- y - mu
  Qir    <- Pir / (Pi * R) 
  eta    <- Qir * lambda            
  
  dUdb   <- matrix(NA, N, Ka)
  if(correlation){
    dUds <- matrix(NA, N, (0.5 * Ka * (Ka + 1))) 
  } else {
    dUds <- matrix(NA, N, Ka)  
  }
  for (i in 1:N){
    beta.r  <- Make.rcoef(beta.bar, sigma, ranp, Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE], 
                          correlation, Pi = NULL, S = NULL)
    dUdb[i,] <- tcrossprod(eta[i, ], beta.r$d.mu)  
    dUds[i,] <- tcrossprod(eta[i, ], beta.r$d.sigma)
  }
  
  if (correlation){
    vecX <- c()
    for (i in 1:Ka){
      vecX <- c(vecX, i:Ka)
    }
    Xac <- Xa[ , vecX]
  } else {
    Xac <- Xa  
  }
  
  gbarfi <- Xc  * rowSums(eta)
  gbarmi <- Xa  * dUdb
  gbarvi <- Xac * dUds
  
  gbari  <- cbind(gbarfi, gbarmi, gbarvi)
  colnames(gbari) <- names(theta)
  attr(lls, 'gradient') <- weights * gbari
  
  if (make.estb){
    b.ran  <- array(NA, dim = c(N, R, Ka))
    b.ran2 <- array(NA, dim = c(N, R, Ka))
    for(j in 1:Ka){
      b.ran[, , j]  <- Br[, ,j ] * Qir
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