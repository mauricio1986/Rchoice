#Jacobian for ordered model
jacobian <- function(kappa){
  ## ij elemenet is dkappa_i/d alpha_j  
  k <- length(kappa)
  etheta <- exp(kappa)
  mat <- matrix(0 ,k, k)
  for (i in 1:k) mat[i:k, i] <- etheta[i]
  mat
}

# Compute time
make.time <- function(object){
  et <- object$time[3]
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  paste(h, "h:", m, "m:", s, "s", sep="")
}

#' @export
ordinal <- function(link = c('probit', 'logit')){
  link <- match.arg(link)
  list(family = 'ordinal', link = link)
}

# Added in version 0.2
repRows <- function(x, n){
  matrix(rep(x, each = n), nrow = n)
}

# Added in version 0.2
make.add <- function(row, col, Ka){
  sa <- makeL(1:rep(0.5 * Ka * (Ka + 1)))
  for (k in row:col){ 
    cb <- sa[k, row]
    form <- paste(paste("x",  cb:(cb + (Ka - k)), sep = ""), paste("x", cb, sep = ""), sep = "*")
  }
  form
}

# Added in version 0.3-3 for marginal effects
make.inter.num  <- function(the.var, beta, beta.hat, X){
  b             <-  if (the.var %in% beta)  beta.hat[beta   %in% the.var] else 0
  # Interactions
  all.interact  <- beta[grepl(":", beta, fixed = TRUE)]
  interactions  <- all.interact[grepl(the.var, all.interact, fixed = TRUE)]
  xinteractions <- gsub(paste(":", the.var, "|", the.var, ":", sep = ""), "", interactions)
  if (length(xinteractions) > 0) {
    bint <- X[, xinteractions, drop = F] %*% beta.hat[interactions] 
  } else {
    bint <- 0
  }
  # Quadratics
  quadratics <- beta.hat[grepl(paste(the.var, "^2", sep = ""), beta, fixed = TRUE)]
  if (length(quadratics) > 0){
    bsq <- 2 * quadratics * X[, the.var, drop = TRUE]
  } else {
    bsq <- 0
  }
  bk <- b + bint + bsq
  return(as.vector(bk))
}

make.inter.factor <- function(the.var, beta, levs){
  lev           <- levs
  the.vars      <- paste0("factor(",the.var,")",lev, sep = "")
  names.b       <- unlist(beta[beta %in% the.vars])
  all.interact  <- beta[grepl(":", beta)]
  interactions  <- all.interact[unlist(lapply(lev, function(x)
    grep(paste0("factor\\(", the.var, "\\)", x), all.interact)))]
  if (length(interactions) > 0){
    temp <- paste0("factor\\(", the.var, "\\)", lev, ":", "|", ":", 
                   "factor\\(", the.var, "\\)", lev, collapse = "|") 
    xinteractions <- gsub(temp, 
                          "", interactions)
  } else {
    interactions  <- NULL
    xinteractions <- NULL
  }
  names       <- c(names.b, interactions)
  names.inter <- c(names.b, xinteractions) 
  return(list(names = names, names.inte = names.inter))
}
