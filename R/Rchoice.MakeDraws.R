##############################
# Make Pseudo and Halton Draws
# Codes from Yves Croissant
##############################

## Halton function: copied from mlogit package
halton <- function(prime = 3, length = 100, drop = 10){
  halt <- 0
  t <- 0
  while(length(halt) < length + drop){
    t <- t + 1
    halt <- c(halt, rep(halt, prime - 1) + rep(seq(1, prime - 1, 1) / prime ^ t, each = length(halt)))
  }
  halt[(drop + 1):(length + drop)]
}

## Make draw function. Modified from mlogit package
make.draws <- function(R, Ka, haltons){
  # Create the matrix of random numbers
  if (!is.null(haltons)){
    length.haltons <- rep(R,Ka)
    prime <- c(3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
               47, 53, 59, 61, 71, 73, 79, 83, 89, 97, 101, 103,
               107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167,
               173, 179, 181, 191, 193, 197, 199)
    drop.haltons <- rep(100,Ka)
    if (!is.na(haltons) && !is.null(haltons$prime)){
      if (length(haltons$prime) != Ka){
        stop("wrong number of prime numbers indicated")
      }
      else{
        prime <- haltons$prime
      }
      if (!is.na(haltons) && !is.null(haltons$drop)){
        if (!length(haltons$drop) %in% c(1,Ka)) stop("wrong number of drop indicated")
        if (length(haltons$drop) == 1){
          drop.haltons <- rep(haltons$drop, Ka)
        }
        else{
          drop.haltons <- haltons$drop
        }
      }
    }
    random.nb <- numeric(0)
    i <- 0
    for (i in 1:Ka){
      random.nb <- cbind(random.nb,qnorm(halton(prime[i],R,drop.haltons[i])))
    }
  }
  else{
    random.nb <- matrix(rnorm(R*Ka), ncol = Ka, nrow = R)
  }
  t(random.nb)  
}

## Make Lower Triangular
makeL <- function(x){
  K <- (-1+sqrt(1 + 8 * length(x)))/2
  mat <- matrix(0, K, K)
  mat[lower.tri(mat, diag = TRUE)] <- x
  mat
}


## Make Random Coefficients
Make.rcoef<-function(beta, sigma, ranp, omega, correlation, Pi=NULL, S=NULL){
  names.r    <- names(ranp)
  censored   <- names.r[ranp == "cn"]
  lognormal  <- names.r[ranp == "ln"]
  normal     <- names.r[ranp == "n"]
  uniform    <- names.r[ranp == "u"]
  triangular <- names.r[ranp == "t"]
  
  Ka    <- nrow(omega) 
  R     <- ncol(omega)  
  
  br    <- matrix(NA, Ka, R) 
  d.mu  <- d.sigma <- br   
  
  beta  <- drop(beta)
  sigma <- drop(sigma)
  names(beta) <- names(sigma) <- names.r
  
  rownames(br) <- rownames(d.mu) <- rownames(d.sigma) <- rownames(omega) <- names.r
  if (!is.null(Pi)){
    rownames(Pi)    <- names.r
    M               <- ncol(Pi)
    d.pis           <- matrix(NA, M*Ka, R)
    names.pi        <- rep(names.r, M)
    rownames(d.pis) <- names.pi
  }
  
  
  if (correlation){
    rownames(omega) <- NULL
    L <- makeL(sigma)
    L.omega <- tcrossprod(L, t(omega))
    rownames(L.omega) <- names.r
    if (!is.null(Pi)){
       br <- drop(beta) + drop(Pi %*% t(S)) + L.omega #Ka*R
       d.pis[,] <- 1
    }else{
      br <- drop(beta) + L.omega
    }
    d.mu <- matrix(1, Ka, R)
    d.sigma <- omega[rep(1:Ka, Ka:1), ]   #d U/ d omega
    for (i in 1:Ka){
      sigi <- i + cumsum(c(0, (Ka-1):1))[1:i]
      if (ranp[i] == "cn"){
        br[i, ]  <- pmax(br[i, ], 0)
        d.mu[i,] <- as.numeric(br[i, ]>0)
        d.sigma[sigi,] <- as.numeric(br[i, ] > 0) * d.sigma[sigi, ]
      }
      if (ranp[i] == "ln"){
        br[i, ] <- exp(br[i, ])
        d.mu[i, ] <- br[i, ]
        d.sigma[sigi, ] <- br[i, ] * d.sigma[sigi, ]
        if (!is.null(Pi)) d.pis[i, ] <- br[i, ]
      }
    }
  ### no correlation ##
  } else {
    if(length(censored) > 0){
      sel <- censored
      if (!is.null(Pi)){
        br[sel, ] <- pmax(beta[sel] + drop(Pi[sel, , drop = F] %*% t(S)) + sigma[sel] * omega[sel, , drop = F], 0)
      }else{
        br[sel, ] <- pmax(beta[sel] + sigma[sel]*omega[sel, , drop = F], 0)
      }
      d.mu[sel, ] <- as.numeric(br[sel, ] > 0)
      d.sigma[sel, ] <- d.mu[sel, ] * omega[sel, ]
      if (!is.null(Pi)) {
        therows <- which(rownames(d.pis) == sel)
        d.pis[therows, ] <- 1
      }
    }
    if(length(lognormal) > 0){
      sel <- lognormal
      if(!is.null(Pi)){
        br[sel,] <-exp(beta[sel] + drop(Pi[sel,, drop = F] %*% t(S)) + sigma[sel] * omega[sel,, drop = F])
      }else {
        br[sel,]<- exp(beta[sel] + sigma[sel] * omega[sel,, drop = F])
      } 
      d.mu[sel,]    <- br[sel,, drop = F]
      d.sigma[sel,] <- d.mu[sel,] * omega[sel,]
      if (!is.null(Pi)) {
        therows <- which(rownames(d.pis) == sel)
        d.pis[therows,] <- br[sel,, drop = F]
      }
    }
    if(length(normal) > 0){
      sel <- normal
      if(!is.null(Pi)){
        br[sel,] <- beta[sel] + drop(Pi[sel,, drop = F] %*% t(S)) + sigma[sel]*omega[sel,, drop = F]
      } else {
        br[sel,] <-beta[sel] + sigma[sel] * omega[sel,, drop = F]
      }
      d.mu[sel,]  <- 1
      d.sigma[sel,] <- omega[sel,, drop = F]
      if (!is.null(Pi)) {
        therows <- which(rownames(d.pis) == sel)
        d.pis[therows,] <- 1
      }
    }
    if (length(uniform) > 0){
      sel <- uniform
      etauni <- pnorm(omega[sel,, drop = F])
      if(!is.null(Pi)){
        br[sel,] <- beta[sel] - sigma[sel] + drop(Pi[sel,, drop = F] %*% t(S)) + 2 * omega[sel,,drop = F] * sigma[sel]
      }else{
        br[sel,] <- beta[sel] - sigma[sel] + 2 * omega[sel,, drop = F]*sigma[sel]
      }
      d.mu[sel,] <- 1
      d.sigma[sel,] <- 2 * etauni - 1
      if (!is.null(Pi)) {
        therows <- which(rownames(d.pis) == sel)
        d.pis[therows,] <- 1
      }
    }
    if (length(triangular) > 0){
      sel <- triangular
      eta05 <- pnorm(omega[sel,, drop = F]) <= 0.5
      d.mu[sel,] <- 1
      d.sigma[sel,] <- eta05 * (sqrt(2 * pnorm(omega[sel,, drop = F])) - 1) + 
        !eta05 * (1 - sqrt(2 * (1 - pnorm(omega[sel,, drop = F]))))
      if (!is.null(Pi)){
        br[sel,] <- beta[sel] + drop(Pi[sel,, drop = F] %*% t(S)) + sigma[sel] * d.sigma[sel,]
      }else{
        br[sel,] <- beta[sel] + sigma[sel] * d.sigma[sel,]
      }
      if(!is.null(Pi)) {
        therows <- which(rownames(d.pis) == sel)
        d.pis[therows,] <- 1
      }
    }    
  }
  if(!is.null(Pi)){
    list(br = br, d.mu = d.mu, d.sigma = d.sigma, d.pis = d.pis)
  }else{
    list(br = br, d.mu = d.mu, d.sigma = d.sigma)
  }
}






