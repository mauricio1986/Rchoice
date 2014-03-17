############################
# S3 method for Rchoice
###########################

#' @rdname Rchoice
#' @S3method plot Rchoice
#' @export
terms.Rchoice <- function(x, ...){
  terms(x$formula)
}

#' @rdname Rchoice
#' @S3method model.matrix Rchoice 
#' @export
model.matrix.Rchoice <- function(object, ...){
  model.matrix(object$formula, object$mf)
}


#' @rdname Rchoice
#' @export
vcov.Rchoice <- function(object,...)
{
# FIXME: See what happens when tyring linear hyp with kappas  
  H<-object$logLik$hessian
  if(object$link=="ordered probit" || object$link=="ordered logit"){
    bhat  <- coef(object)
    J     <- length(attr(object$coefficients, "kappas"))
    ahat  <- bhat[1:J]
    A     <- diag(length(bhat))
    z.id  <- seq(1,J,1)
    Jacob <- jacobian(ahat)
    A[z.id, z.id] <- Jacob
    vcov <- A %*% solve(-H) %*% t(A)
    rownames(vcov) <- colnames(vcov) <- names(bhat)
  } else {
  vcov<-(solve(-H))
  rownames(vcov) <- colnames(vcov) <- names(coef(object))
  }
  return(vcov)
}


#' @rdname Rchoice
#' @S3method coef Rchoice
coef.Rchoice <- function(object, ...){
  result <- object$coefficients
  return(result)
}


#' @rdname Rchoice
#' @export nObs.Rchoice 
nObs.Rchoice <- function( x, ... ) {
  return(x$logLik$nobs )
}

#' @rdname Rchoice
#' @S3method fitted Rchoice
#' @export
fitted.Rchoice <- function(object, ...){
  result <- object$probabilities
  return(result)
}

#' @rdname Rchoice
#' @S3method df.residual Rchoice
#' @export
df.residual.Rchoice <- function(object, ...){
  n <- length(residuals(object))
  K <- length(coef(object))
  return(n-K)
}

#' @rdname Rchoice
#' @S3method update Rchoice
#' @export
update.Rchoice <- function (object, new, ...){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(new))
    call$formula <- update(formula(object), new)
  if(length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    ## do these individually to allow NULL to remove entries.
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if(any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}

#' Akaike's Information Criterion
#' 
#' Calculate Akaike's information Criterion (AIC) or the Bayesian
#' information Criterion (BIC) for a model of object of class.
#' 
#' @param object a fitted model of class \code{Rchoice}
#' @param ... additional arguments to be passed to or from other functions
#' @param k a numeric value, use as penalty coefficient for number of parameters
#' in the fitted model
#' @return a numeric value with the corresponding AIC or BIC value.
#' @seealso \code{\link[Rchoice]{Rchoice}}
#' @method AIC Rchoice
#' @export AIC.Rchoice
#' @examples
#' ## Probit model
#' data("Workmroz")
#' probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  data = Workmroz , link="probit")
#' summary(probit)
#' 
#' AIC(probit)
#' BIC(probit)
AIC.Rchoice <- function(object, ..., k = 2) {
  return(- 2 * object$logLik$maximum[[1]] + k* length(coef(object)) )
}

#' @rdname AIC.Rchoice 
#' @export BIC.Rchoice
BIC.Rchoice <- function( object, ...) {
  return(AIC(object, k = log(nObs(object))) )
}

#' @rdname Rchoice
#' @S3method logLik Rchoice
#' @export
logLik.Rchoice <- function(object,...){
  structure(-object$logLik$maximum[[1]], df = length(object$coefficients),
            nobs = nObs(object), class = "logLik")
}


#' Bread for sandwiches
#' 
#' Computes the bread of the sandwich covariance matrix
#' 
#' @param x a fitted model of class \code{Rchoice}
#' @param ... Other arguments when \code{bread} is applied to another
#' class object
#' @return the covariance matrix times observations
#' @references Zeileis A (2006), Object-oriented Computation of Sandwich 
#' Estimators. Journal of Statistical Software, 16(9), 1--16.
#' @method bread Rchoice
#' @export bread.Rchoice
#' @examples
#' ## Probit model
#' data("Workmroz")
#' probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  data = Workmroz , link="probit")
#' summary(probit)
#' 
#' library(sandwich)
#' bread(probit) 
bread.Rchoice <- function( x, ... ) {
  return( vcov( x ) * nObs(x))
}

#' Gradient for observations
#' 
#' It extracts the gradient for each observations evaluated at the estimated parameters
#' 
#' @param x a fitted model of class \code{Rchoice}
#' @param ... Other arguments when \code{estfun} is applied to another
#' class object
#' @return the gradient matrix of dimension n times k 
#' @references Zeileis A (2006), Object-oriented Computation of Sandwich 
#' Estimators. Journal of Statistical Software, 16(9), 1--16.
#' @method estfun Rchoice
#' @export estfun.Rchoice
#' @examples
#' ## Probit model
#' data("Workmroz")
#' probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  data = Workmroz , link="probit")
#' summary(probit)
#' 
#' estfun(probit) 
estfun.Rchoice <- function( x, ... ) {
  return(x$logLik$gradientObs )
}


#' @S3method print Rchoice
#' @export
print.Rchoice <- function(x, digits = max(3,getOption("digits")-3),
                          width = getOption("width"),...)
{
  cat("\nCall:\n", deparse(x$call),"\n\n", sep="")
  
  cat("\nCoefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname Rchoice
#' @S3method summary Rchoice
#' @method summary Rchoice
#' @export
summary.Rchoice <- function (object,...){
  if(object$link=="probit" || object$link=="logit" || object$link=="poisson" ){
    b <- object$coefficients
    std.err <- sqrt(diag(vcov(object)))
  }else{
    b<-attr(object$coefficients, "fixed")
    std.err <- sqrt(diag(vcov(object)))[names(attr(object$coefficients, "fixed"))]
  }
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "t-value", "Pr(>|t|)")
  
  et <- object$time[3]
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  
  result <- structure(
    list(
      CoefTable     = CoefTable,
      link          = object$link,
      mf            = object$mf,
      logLik        = object$logLik,
      formula       = object$formula,
      time          = object$time,
      freq          = object$freq,
      call          = object$call,
      draws         = object$draws,
      R.model       = object$R.model,
      R             = object$R,
      tstr          = paste(h, "h:", m, "m:", s, "s", sep="")),
    class = c("summary.Rchoice","Rchoice")
  )
  
  if(object$link =="ordered probit" || object$link =="ordered logit") result$kappa<-attr(object$coefficients,"kappas")
  return(result)
}


##' @S3method print summary.Rchoice
print.summary.Rchoice<- function(x,digits = max(3, getOption("digits") - 2),
                                 width = getOption("width"),
                                 ...)
{
  cat(paste("\nModel:",x$link))
  cat(paste("\nModel estimated on:",format(Sys.time(), "%a %b %d %X %Y"),"\n"))
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  if(!(x$link=="poisson")){
    cat("\nFrecuencies of categories:\n")
    print(prop.table(x$freq),digits = digits)
  }
  
  cat(paste("The estimation took:", x$tstr,"\n"))
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable,digits=digits)
  
  if(x$link=="ordered probit" || x$link == "ordered logit"){
    cat("\nThresholds:\n")
    print(x$kappa)
  }
  
  
  cat(paste("\nOptimization of log-likelihood by",x$logLik$type))
  cat(paste("\nLL:",round(x$logLik$maximum,digits)))
  cat(paste("\nNumber of Individuals:",x$logLik$nobs))
  cat(paste("\nNumber of Iterations:",x$logLik$iterations))
  cat(paste("\nExit of MLE:",x$logLik$message))
  
  
  if(x$R.model){
    if(is.null(x$draws)){
      cat(paste("\nSimulation based on",x$R,"pseudo-random draws"))
    }else {
      cat(paste("\nSimulation based on",x$R,"Halton draws"))
    }
  }
  invisible(x)
}


#' Plot random parameters
#' 
#' Plot the conditional expectation of random parameters estimated by \code{Rchoice}.
#' 
#' @param x a object of class \code{Rchoice},
#' @param par a string giving the name of the variable with random parameter,
#' @param type a string indicating the type of distribution: it can be a \code{histogram} or a \code{density} of
#' the conditional expectation of the random coefficients,
#' @param ind a boolean. If \code{TRUE}, a 95% interval of conditional distribution for each individual is plotted. 
#' As default, the conditional expectation of \code{par} for the first 10 individual is plotted,
#' @param id only relevant if \code{ind} is not \code{NULL}. This is a vector indicating the position of the individual
#' for whom the user want to plot the conditional coefficients, 
#' @param bin bin of histrogram,
#' @param adjust  bandwidth for the kernel density,
#' @param ... further arguments to be passed to \code{qplot} or \code{plotCI}, 
#' @return a plot with the distribution or a confident interval of the conditional random coefficients.
#' @references
#' \itemize{
#' \item Greene, W. H. (2003). Econometric analysis. Pearson Education India.
#' \item Train, K. (2009). Discrete choice methods with simulation. Cambridge university press.
#' }
#' @seealso \code{\link[Rchoice]{Rchoice}}, \code{\link[ggplot2]{ggplot2}} 
#' @S3method plot Rchoice
#' @method plot Rchoice
#' @export
plot.Rchoice<-function(x, par=NULL, ind=FALSE, id=NULL, type = c("density", "histogram"), bin = 1 , adjust = 1,...){
  if(!x$R.model) stop("the plot method is only relevant for random parameters")
  if (is.null(par)) stop("Must specified the name of the random parameters")
  type <- match.arg(type)
  
  if(!ind){
    ylab<-switch(type,
               "density"   = "Density",
               "histogram" = "Frequency")
  
    rpar<-x$b.random[,par]
    ggplot2::qplot(as.vector(rpar), geom=type,
        main=paste("Conditional Distribution: ", par), xlab=expression(E(hat(beta[i]))), 
        ylab=ylab, binwidth = bin, adjust=adjust)
  }
  else{
    if(is.null(id)) id<-seq(1,10,1) 
    f.bran<-x$b.random[,par]
    f.sran<-x$sd.random[,par]
    lower<-f.bran-2*f.sran
    upper<-f.bran+2*f.sran
    plotrix::plotCI(id,f.bran[id],ui=upper[id],li=lower[id], col="red",
                     xlab="Individual", ylab=expression(E(hat(beta[i]))),
                     lty = 2,main=paste("Conditional Distribution: ", par),
                     pch=21)
  } 
}


#' Covariance and Correlation matrix of random parameters
#' 
#' Computes the Variance-Covariance matrix and the Correlation matrix of the random parameters
#' 
#' @param x a object of class \code{Rchoice},
#' @param ... further arguments
#' @return a matrix with the variance of the random parameters if model is fitted with random coefficients or the correlation matrix if argument 
#' \code{correlation = TRUE} in the fitted model.
#' @references
#' \itemize{
#' \item Greene, W. H. (2003). Econometric analysis. Pearson Education India.
#' \item Train, K. (2009). Discrete choice methods with simulation. Cambridge university press.
#' }
#' @seealso \code{\link[Rchoice]{Rchoice}}
#' @export
cov.Rchoice <- function(x){
  if (is.null(x$ranp)) stop('cov.Rchoice only relevant for random coefficient model')
  K<-length(x$ranp)
  nr<-names(x$ranp)
  if (x$correlation){
    Ktot <- length(x$coefficients)
    v    <- x$coefficients[(Ktot - 0.5 * K * (K + 1) + 1) : Ktot]
    V    <- tcrossprod(makeL(v))
    colnames(V) <- rownames(V) <- nr
    sv <- sqrt(diag(V))
  } else{
    Ktot <-length(x$coefficients)
    sv   <-tail(x$coefficients, K)
    V    <-matrix(0, K, K)
    diag(V) <- sv^2
    colnames(V)<-rownames(V)<-nr
  }
  V
}

#' @rdname cov.Rchoice
#' @export
cor.Rchoice <- function(x){
  if (!x$correlation) stop('cor.Rchoice only relevant for correlated random coefficient')
  V   <- cov.Rchoice(x)
  nr  <- names(x$ranp)
  D   <- diag(sqrt(diag(V)))
  Rho <- solve(D) %*% V %*% solve(D)
  colnames(Rho) <- rownames(Rho) <- nr
  Rho
}

#' @rdname Rchoice
#' @export getSummary.Rchoice
getSummary.Rchoice<-function (obj, alpha = 0.05, ...)
{
  smry <- summary(obj)
  coef <- smry$CoefTable
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  nObs(obj)
  ll <- logLik(obj)
  sumstat <- c(logLik = ll, deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}



