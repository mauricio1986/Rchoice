#########################################
# Standard log likelihood functions
# for binary, poisson and ordered model
#########################################

lnpoisson<-function(theta, y, X,
                    weights = NULL, ...)
{
  if (is.null(weights)) weights <- rep(1, nrow(X))
  index <- tcrossprod(X, t(theta))
  mu    <- pmax(exp(index), .Machine$double.eps)
  pi    <- dpois(y, lambda = mu)
  ll    <- sum(weights * log(pi))
  
  #Gradient
  lambda <- y-mu
  G      <- as.vector(lambda) * X
  attr(ll,'gradient') <- weights * G
  
  #Hessian
  lambda2 <- as.vector((-1) * mu)
  H       <- crossprod(lambda2 * X, X)
  attr(ll,'hessian') <- H
  
  attr(ll,'probabilities') <- pi
  ll
}


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
  attr(ll,'gradient') <- weights * G
  
  #Hessian
  ddens   <- ddfun(qi * index)
  lambda2 <- ((ddens / pi)-(dens / pi) ^ 2) * qi^2
  H <- crossprod(as.vector(lambda2) * X, X)
  attr(ll,'hessian') <- H
  
  attr(ll,'probabilities') <- pi
  ll
}

#Jacobian for ordered model
jacobian <- function(kappa){
  ## ij elemenet is dkappa_i/d alpha_j  
  k <- length(kappa)
  etheta <- exp(kappa)
  mat <- matrix(0 ,k, k)
  for (i in 1:k) mat[i:k, i] <- etheta[i]
  mat
}


lnordered<-function(theta, y, X, link,
                    weights = NULL, ... )
{
  if (is.null(weights)) weights <- rep(1, nrow(X))
  pfun <- switch(link,
               "ordered probit" = pnorm,
               "ordered logit"  = plogis)
  dfun <- switch(link,
               "ordered probit" = dnorm,
               "ordered logit"  = dlogis)
  qfun <- switch(link,
               "ordered probit" = qnorm,
               "ordered logit"  = qlogis)
  # it is assumed that y is a "factor" class
  m       <- sort(unique(y))
  m       <- unclass(m)
  J       <- length(levels(y))
  #Just J-2 alpha are estimated, so the gradient is just for 2:(J-1) cuts
  m       <- as.matrix(m[2:(J - 1)], nrow = (J - 2))               
  y       <- unclass(y)
  delta<-t(kronecker(m, t(rep(1,nrow(X)))))    
  #indicator for alternative and previous alternative
  deltaj<- delta == y
  deltak<- delta == y - 1
  
  alpha    <- theta[1:(J - 2)]
  kappa    <- c(-Inf, cumsum(c(0, exp(alpha))) , +Inf)
  beta     <- theta[-c(1:(J - 2))]
  index    <- tcrossprod(X, t(beta))
  #thresh   <- -qfun(.Machine$double.eps)
  thresh   <- 100          # As polr of MASS
  eta1     <- pmin(kappa[y + 1] - index, thresh)
  eta2     <- pmax(kappa[y] - index, -thresh)
  pi       <- pfun(eta1) - pfun(eta2)
  if (all(pi > 0)) ll <- sum(weights * log(pi)) else ll <- Inf

  
  #Gradient
  phi1     <- dfun(eta1) ; phi2 <- dfun(eta2)
  lambda   <- drop((phi2 - phi1) / pi)
  gbeta    <- as.vector(lambda) * X
  
  #Gradient of kappa
  lamkap   <- (deltaj * drop(phi1) - deltak * drop(phi2)) / as.vector(pi)         #N*(j-2)
  gkappa   <- lamkap %*% jacobian(alpha)
  G        <- cbind(gkappa, gbeta)
  
  if (all(pi > 0))  attr(ll,'gradient') <- weights * G else attr(ll,'gradient') <- rep(NA, length(theta)) 
  
  attr(ll,'probabilities') <- pi
  ll
}

#########################################
# Simulated maximum likelihood functions
# for binary, poisson and ordered model
#########################################

make.prob<-function(link)
{
  switch(link,
         "probit" = {
           pfun <- function(y, index){
                    thresh <- -qnorm(.Machine$double.eps)
                    eta    <- pmin(pmax(index, -thresh), thresh)
                    pnorm((2 * y - 1) * eta)
                  }
           dfun <- function(y, index){
                    dens  <- pmax(dnorm((2 * y - 1) * index), .Machine$double.eps)
                    (2 * y - 1) * (dens / pfun(y ,index))
                  }           
         },
         "poisson" = {
           pfun <- function(y, index){
                    mu<-pmax(exp(index), .Machine$double.eps)
                    dpois(y, lambda = mu)
                  }
           dfun <- function(y, index){
                    mu<-pmax(exp(index), .Machine$double.eps)             
                    y-mu       
                  }
         },
         "logit" = {
           pfun <- function(y,index){
                    thresh <- -qlogis(.Machine$double.eps)
                    eta    <- pmin(pmax(index, -thresh), thresh)
                    plogis((2 * y - 1) * eta)        
                  }
           dfun <- function(y,index){
                    dens <- pmax(dlogis((2 * y - 1) * index), .Machine$double.eps)
                    (2 * y - 1) * (dens / pfun(y,index))   
                  }
         },
         stop(gettextf("%s link not recognised", sQuote(link)), 
              domain = NA)
   )
list(pfun=pfun, dfun=dfun)  
}

#SLL for Random Parameters
lnl.ran<-function(theta, y, X, ranp, R, correlation, link,
                  weights = NULL, haltons = NULL,
                  seed = 123, make.estb = FALSE,
                  ... )
{
  #Get Variables
  N    <- nrow(X)
  K    <- ncol(X)
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  Xa   <- X[ , Vara, drop = F]                        
  Xc   <- X[ , Varc, drop = F]                        
  
  #Get vector of parameters
  gamma    <- as.matrix(theta[1:Kc])               
  beta.bar <- as.matrix(theta[(Kc + 1):(Kc + Ka)])  
  sigma    <- theta[-c(1:(Kc + Ka))]                 
  
  #Make Random Draws
  set.seed(seed)
  Omega    <- make.draws(R*N, Ka, haltons) 
  
  #Fixed Part of index
  ZB <- as.vector(crossprod(t(Xc), gamma)) #N*1
  
  #Random Part of index
  XB <- matrix(NA, N, R)  
  #Matrix for Est.E(b_i)
  Br <- array(NA, dim = c(N, R, Ka))
  for (i in 1:N){
    omega   <- Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE] #Ka*R for individual i
    beta.r  <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, correlation = correlation,
                       Pi = NULL, S = NULL)
    XB[i, ]  <- crossprod(t(Xa[i, , drop = FALSE]), beta.r$br)
    Br[i, , ] <- t(beta.r$br) #Now is R*Ka
  }
  
  #Get probability & Lls
  index <- ZB + XB
  Pir   <- make.prob(link)$pfun(y,index)
  Pi    <- apply(Pir, 1, mean)
  lls   <- sum(weights * log(Pi))
 
  #Make gradient
  lambda <- make.prob(link)$dfun(y,index)
  Qir    <- Pir / (Pi * R) 
  eta    <- Qir * lambda            
  
  dUdb<-  matrix(NA, N, Ka)
  if(correlation){
    dUds <- matrix(NA, N, (0.5 * Ka * (Ka + 1))) 
  }else{
    dUds <- matrix(NA, N, Ka)  
  }
  for(i in 1:N){
    omega  <- Omega[ , ((i - 1) * R + 1):(i * R), drop = FALSE]
    beta.r <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, correlation = correlation, 
                         Pi = NULL, S = NULL)
    dUdb[i,] <- tcrossprod(eta[i, ], beta.r$d.mu)  
    dUds[i,] <- tcrossprod(eta[i, ], beta.r$d.sigma)
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
  
  gbarfi <- Xc  * as.vector(apply(Qir * lambda, 1, sum))
  gbarmi <- Xa  * dUdb
  gbarvi <- Xac * dUds
  
  gbari  <- cbind(gbarfi, gbarmi, gbarvi)
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

#LL for hierarchical random parameter model
lnlH.ran<-function(theta, y, X, S, ranp, R, link, seed, correlation, 
                   weights = NULL, haltons = NULL, make.estb = FALSE,
                   ...)
{
  #Get Variables and Global Parameters
  N    <- nrow(X)
  K    <- ncol(X)
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  M    <- ncol(S)
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  Xa   <- X[ , Vara, drop = F]                        
  Xc   <- X[ , Varc, drop = F]                        
  
  #Get estimated coefficients
  gamma    <- as.matrix(theta[1:Kc])               
  beta.bar <- as.matrix(theta[(Kc + 1):(Kc + Ka)])  
  phis     <- theta[(Kc + Ka + 1):(Kc + Ka + Ka * M)]         
  Phi      <- matrix(phis, nrow = Ka)                   
  colnames(Phi) <- colnames(S)
  rownames(Phi) <- colnames(Xa)
  sigma    <- theta[-c(1:(Kc + Ka + Ka * M))]           
  
  #Make Random Draws
  set.seed(seed)
  Omega <- make.draws(R * N, Ka, haltons) #Ka * N*R
  
  #Fixed part of index
  ZB <- as.vector(crossprod(t(Xc), gamma)) 
  
  XB  <- matrix(NA, N, R)  
  XaS <- matrix(NA, N, length(phis))
  colnames(XaS) <- names(phis)
  Br <- array(NA, dim = c(N, R, Ka))
  for (i in 1:N){
    omega   <- Omega[ , ((i - 1) * R + 1):(i * R) , drop = FALSE] #Ka*R for individual i
    Si      <- S[i , , drop = FALSE]
    beta.r  <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, 
                          correlation = correlation, Pi = Phi, S = Si)
    XB[i,]  <- crossprod(t(Xa[i, , drop = FALSE]), beta.r$br)
    XaS[i,] <- kronecker(t(S[i, , drop = FALSE]), t(Xa[i, , drop = FALSE]))
    Br[i,,] <- t(beta.r$br) 
  }
  
  index <- ZB + XB
  Pir   <- make.prob(link)$pfun(y,index)
  Pi    <- apply(Pir, 1, mean) 
  lls   <- sum(weights * log(Pi))
  
  #Make gradient
  lambda <- make.prob(link)$dfun(y,index)
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
    omega<-Omega[,((i - 1) * R + 1):(i * R) , drop = FALSE] #Ka*R for individual i
    Si<-S[i,,drop=FALSE]
    beta.r      <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, 
                              correlation = correlation, Pi = Phi, S = Si)
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
  
  gbarfi  <- Xc  * as.vector(apply(eta, 1, sum))
  gbarmi  <- Xa  * dUdb
  gbarphi <- XaS * dUdphi
  gbarvi  <- Xac * dUds
  
  gbari <- cbind(gbarfi, gbarmi, gbarphi, gbarvi)
  colnames(gbari)<-names(theta)
  attr(lls, 'gradient') <-weights * gbari
  
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
    
    sd.ran <- sqrt(b.ran2-b.ran^2)
    colnames(sd.ran) <- colnames(Xa)
    attr(lls,'b.random') <- b.ran
    attr(lls,'sd.random') <- sd.ran
    attr(lls,'probabilities') <- Pi
  }
  
  lls
}


lnordered.ran<-function(theta, y, X, ranp, R, correlation, link,
                        weights = NULL, haltons = NULL,
                        seed = 123, make.estb = FALSE,
                        ... )
{
  pfun <- switch(link,
               "ordered probit" = pnorm,
               "ordered logit"  = plogis)
  dfun <- switch(link,
               "ordered probit" = dnorm,
               "ordered logit"  = dlogis)
  #Get Variables and Global Parameters
  N    <- nrow(X)
  K    <- ncol(X)
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  Xa   <- X[, Vara, drop = F]                        
  Xc   <- X[, Varc, drop = F] 
  
  #For gradient
  m       <- unclass(sort(unique(y)))
  J       <- length(levels(y))
  m       <- as.matrix(m[2:(J - 1)], nrow = (J - 2))               
  y       <- unclass(y)
  delta   <- t(kronecker(m, t(rep(1, nrow(X)))))    
  deltaj  <- delta == y
  deltak  <- delta == y-1
  
  #Get estimated coefficients
  alpha    <- theta[1:(J - 2)]
  kappa    <- c(-Inf, cumsum(c(0, exp(alpha))), +Inf)
  gamma    <- as.matrix(theta[((J - 2) + 1):((J - 2) + Kc)])
  betasr   <- theta[-c(1:((J - 2) + Kc))]
  beta.bar <- as.matrix(betasr[c(1:Ka)])  
  sigma    <- betasr[-c(1:Ka)]            
  
  #Make Random Draws
  set.seed(seed)
  Omega <- make.draws(R * N, Ka, haltons) #Ka * N*R
  
  #Fixed part of index
  ZB    <- as.vector(crossprod(t(Xc), gamma)) #N*1
  
  #Random Part of index
  XB    <- matrix(NA, N, R)  #N*R
  #Matrix for Est.E(b_i)
  Br    <- array(NA, dim = c(N, R, Ka))
  for (i in 1:N){
    omega    <-Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE] #Ka*R for individual i
    beta.r   <-Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, 
                          correlation = correlation, Pi = NULL, S = NULL)
    XB[i, ]  <-crossprod(t(Xa[i, , drop = FALSE]), beta.r$br)
    Br[i, , ] <-t(beta.r$br) #Now is R*Ka
  }
  
  #Get probability
  index    <- ZB + XB
  #thresh   <- -qfun(.Machine$double.eps)
  #thresh  <- 100   # as polr in MASS
  #eta1     <- pmin(kappa[y+1] - index, thresh)
  #eta2     <- pmax(kappa[y]   - index, -thresh)
  eta1     <- kappa[y + 1] - index
  eta2     <- kappa[y]   - index
  Pir      <- pfun(eta1) - pfun(eta2)
  Pir      <- ifelse(Pir == 0, 0.000001, Pir)
  Pi       <- apply(Pir, 1, mean)
  
  if (all(Pi > 0)) lls <- sum(weights * log(Pi)) else lls <- Inf
  
  #Gradient
  phi1     <- dfun(eta1) ; phi2 <- dfun(eta2)
  lambda   <- (phi2 - phi1) / Pir
  Qir      <- Pir / (Pi * R)
  eta      <- Qir * lambda
  
  #Kappa part
  gkappa   <- vector(mode = "list", length = (J - 2))
  for (j in 1:(J - 2)){
    gkappa[[j]] <- matrix(NA, N, R)
    gkappa[[j]] <- drop(deltaj[ , j]) * (phi1 / Pir) - drop(deltak[ ,j]) * (phi2 / Pir)
  }
  gkappa <- lapply(gkappa, function(x) x * Qir)
  gkappa <- lapply(gkappa, function(x) apply(x, 1, sum))
  gkappa <- Reduce(cbind, gkappa) %*% jacobian(alpha)


  dUdb <- matrix(NA, N, Ka)
  if(correlation){
    dUds <- matrix(NA, N, (0.5 * Ka * (Ka + 1))) 
  }else{
    dUds <- matrix(NA, N, Ka)  
  }
  for(i in 1:N){
    omega     <- Omega[ , ((i - 1) * R + 1):(i * R) , drop = FALSE]
    beta.r    <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, 
                           correlation = correlation, Pi = NULL, S = NULL)
    dUdb[i, ] <- tcrossprod(eta[i, ], beta.r$d.mu)  
    dUds[i, ] <- tcrossprod(eta[i, ], beta.r$d.sigma)
  }
  
  if(correlation){
    vecX <- c()
    for (i in 1:Ka){
      vecX <- c(vecX, i:Ka)
    }
    Xac <- Xa[ ,vecX]
  }else{
    Xac <- Xa  
  }
  
  gbarfi <- Xc  * as.vector(apply(Qir * lambda, 1, sum))
  gbarmi <- Xa  * dUdb
  gbarvi <- Xac * dUds
  
  gbari  <- cbind(gkappa, gbarfi , gbarmi , gbarvi)
  if(all(Pi > 0))  attr(lls, 'gradient') <- weights * gbari else attr(lls, 'gradient') <- rep(NA, length(theta))
  
  if (make.estb){
    b.ran  <- array(NA, dim = c(N, R, Ka))
    b.ran2 <- array(NA, dim = c(N, R, Ka))
    for(j in 1:Ka){
      b.ran[, , j]  < Br[, , j] * Qir
      b.ran2[, , j] <- (Br[, , j]^2) * Qir
    }
    b.ran  <- apply(b.ran,  c(1,3), sum)
    b.ran2 <- apply(b.ran2, c(1,3), sum)
    colnames(b.ran) <- colnames(Xa)
    
    sd.ran <- sqrt(b.ran2 - b.ran^2)
    colnames(sd.ran) <- colnames(Xa)
    attr(lls,'b.random') <- b.ran
    attr(lls,'sd.random') <- sd.ran
    attr(lls,'probabilities') <- Pi
  } 
  
  lls
}


lnorderedH.ran<-function(theta, y, X, S, ranp, R, correlation, link,
                         weights = NULL, haltons = NULL,
                         seed = 123, make.estb = FALSE,
                         ... )
{
  pfun <- switch(link,
               "ordered probit" = pnorm,
               "ordered logit"  = plogis)
  dfun <- switch(link,
               "ordered probit" = dnorm,
               "ordered logit"  = dlogis)
  #qfun <- switch(link,
  #             "ordered probit" = qnorm,
  #             "ordered logit"  = qlogis)
  
  #Get Variables and Global Parameters
  N    <- nrow(X)
  K    <- ncol(X)
  Vara <- sort(match(names(ranp), colnames(X)))
  Varc <- (1:K)[- Vara]
  Ka   <- length(Vara)
  Kc   <- length(Varc)
  Xa   <- X[, Vara, drop = F]                        
  Xc   <- X[, Varc, drop = F]
  M    <- ncol(S)
  
  #For gradient
  m       <- unclass(sort(unique(y)))
  J       <- length(levels(y))
  m       <- as.matrix(m[2:(J - 1)], nrow = (J - 2))               
  y       <- unclass(y)
  delta   <- t(kronecker(m, t(rep(1, nrow(X)))))    
  deltaj  <- delta == y
  deltak  <- delta == y-1
  
  #Get estimated coefficients
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
         
  
  #Make Random Draws
  set.seed(seed)
  Omega <- make.draws(R * N, Ka, haltons) #Ka * N*R
  
  #Fixed part of index
  ZB <- as.vector(crossprod(t(Xc), gamma)) #N*1
  
  #Random Part of index
  XB  <- matrix(NA, N, R)  #N*R
  XaS <- matrix(NA, N, length(phis))
  colnames(XaS) <- names(phis)
  Br <- array(NA, dim = c(N, R, Ka))
  for (i in 1:N){
    omega<-Omega[,((i - 1) * R + 1):(i * R) , drop = FALSE] #Ka*R for individual i
    Si        <- S[i, , drop = FALSE]
    beta.r    <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, 
                            correlation = correlation, Pi = Phi, S = Si)
    XB[i, ]   <- crossprod(t(Xa[i, , drop = FALSE]), beta.r$br)
    XaS[i, ]  <- kronecker(t(S[i, , drop = FALSE]), t(Xa[i, , drop = FALSE]))
    Br[i, , ] <- t(beta.r$br) #Now is R*Ka
  }
  
  #Get probability
  index  <- ZB+XB
  #thresh <- -qfun(.Machine$double.eps)
  #thresh<- 100
  #eta1   <- pmin(kappa[y + 1] - index, thresh)
  #eta2   <- pmax(kappa[y] - index, -thresh)
  eta1    <- kappa[y + 1] - index
  eta2    <- kappa[y] - index
  Pir    <- pfun(eta1) - pfun(eta2)
  Pir    <- ifelse(Pir == 0, 0.00000001, Pir)
  Pi     <- apply(Pir, 1, mean)
  if (all(Pi > 0)) lls <- sum(weights * log(Pi)) else lls <- Inf

  
  #Gradient
  phi1     <- dfun(eta1) ; phi2 <- dfun(eta2)
  lambda   <- (phi2 - phi1) / Pir
  Qir      <- Pir / (Pi * R)
  eta      <- Qir * lambda
  
  #Kappa part
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
  if(correlation){
    dUds <- matrix(NA, N, (0.5 * Ka * (Ka + 1))) 
  }else{
    dUds <- matrix(NA, N, Ka)  
  }
  for(i in 1:N){
    omega       <- Omega[, ((i - 1) * R + 1):(i * R) , drop = FALSE]
    Si          <- S[i, , drop = FALSE]
    beta.r      <- Make.rcoef(beta = beta.bar, sigma = sigma, ranp = ranp, omega = omega, 
                              correlation = correlation, Pi = Phi, S = Si)
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
  }else{
    Xac <- Xa  
  }
  
  gbarfi  <- Xc  * as.vector(apply(Qir * lambda, 1, sum))
  gbarmi  <- Xa  * dUdb
  gbarphi <- XaS * dUdphi
  gbarvi  <- Xac * dUds
  
  gbari <- cbind(gkappa, gbarfi , gbarmi, gbarphi, gbarvi)
  if(all(Pi > 0))  attr(lls, 'gradient') <- weights * gbari else attr(lls, 'gradient') <- rep(NA, length(theta))
  
  if (make.estb){ 
    b.ran  <- array(NA, dim = c(N, R, Ka))
    b.ran2 <- array(NA, dim = c(N, R, Ka))
    for(j in 1:Ka){
      b.ran[, ,j]  <-Br[, , j] * Qir
      b.ran2[, ,j] <-(Br[, , j]^2) * Qir
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


