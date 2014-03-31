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

is.hierarchical <- function(object) {
  ifelse(length(object)[2] == 2, TRUE, FALSE)
}

has.intercept <- function(object, ...) {
  UseMethod("has.intercept")
}

has.intercept.default <- function(object, ...) {
  has.intercept(formula(object), ...)
}

has.intercept.formula <- function(object, ...) {
  attr(terms(object), "intercept") == 1L
}

has.intercept.Formula <- function(object, rhs = NULL, ...) {
  if(is.null(rhs)) rhs <- 1:length(attr(object, "rhs"))
  sapply(rhs, function(x) has.intercept(formula(object, lhs = 0, rhs = x)))
}

