############################
# S3 method for Rchoice
###########################

#' @rdname Rchoice
#' @method terms Rchoice
#' @export
terms.Rchoice <- function(x, ...){
  terms(x$formula)
}

#' @rdname Rchoice
#' @method model.matrix Rchoice
#' @export
model.matrix.Rchoice <- function(object, ...){
  X <- model.matrix(object$formula, object$mf)
  if (has.intercept(object$formula, rhs = 1)) {
    namesX <- colnames(X)
    namesX[1L] <- "constant"
    colnames(X) <- namesX
  }
  X
}

#' vcov method for Rchoice objects
#' 
#' The \code{vcov} method for \code{Rchoice} objects extracts the covariance matrix of the coefficients or the random parameters. It also allows to get the standard errors for the variance-covariance matrix of the random parameters
#' 
#' @param object a fitted model of class \code{Rchoice},
#' @param what indicates which covariance matrix has to be extracted. The default is \code{coefficient}. In this case the \code{vcov} behaves as usual. If \code{what = "ranp"} the covariance matrix of the random parameters is returned as default, 
#' @param type if the model is estimated with random parameters, then this argument indicates what matrix should be returned. If \code{type = "cov"}, then the covariance matrix of the random parameters is returned; if \code{type = "cor"} then the correlation matrix of the random parameters is returned; if \code{type = "sd"} then the standard deviation of the random parameters is returned,
#' @param se if \code{TRUE} and \code{type = "cov"} then the standard error of the covariance matrix of the random parameters is returned; if \code{TRUE} and \code{type = "sd"} the standard error of the standard deviation of the random parameter is returned. This argument if valid only if the model is estimated using correlated random parameters,
#' @param digits number of digits,
#' @param x a fitted model of class \code{Rchoice},
#' @param sd if \code{TRUE}, then the standard deviation of the random parameters are returned,
#' @param ... further arguments
#' @details This new interface replaces the \code{cor.Rchoice}, \code{cov.Rchoice} and \code{se.cov.Rchoice} functions which are deprecated.
#' @seealso \code{\link[Rchoice]{Rchoice}} for the estimation of discrete choice models with random parameters.
#' @method vcov Rchoice
#' @import stats
#' @export
vcov.Rchoice <- function(object, what = c('coefficient', 'ranp'), type = c('cov', 'cor', 'sd'), 
                         se = FALSE, digits = max(3, getOption("digits") - 2), ...)
{
  what <- match.arg(what)
  type <- match.arg(type)
  if (what == 'coefficient') {
    H <- object$logLik$hessian
    if (object$family == "ordinal"){
      bhat  <- coef(object)
      ahat  <- attr(object$coefficients, "alphas")
      J     <- length(ahat)
      A     <- diag(length(bhat))
      z.id  <- seq(1, J, 1)
      Jacob <- jacobian(ahat)
      A[z.id, z.id] <- Jacob
      result <- A %*% solve(-H) %*% t(A)
      rownames(result) <- colnames(result) <- names(bhat)
    } else {
      result <- (solve(-H))
      rownames(result) <- colnames(result) <- names(coef(object))
    }
    return(result)
  }
  if (what == 'ranp') {
    if (se) {
      if (type == 'cov') se.cov.Rchoice(object, sd = FALSE, digits = digits)
      if (type == 'sd')  se.cov.Rchoice(object, sd = TRUE, digits = digits)
      if (type == 'cor') stop("standard error for correlation coefficients not implemented yet")
    } else {
      if (type == 'cov') print(cov.Rchoice(object)) 
      if (type == 'cor') print(cor.Rchoice(object))
      if (type == 'sd')  print(sqrt(diag(cov.Rchoice(object))))
    }
  }
}

#' @rdname Rchoice
#' @method coef Rchoice
#' @export
coef.Rchoice <- function(object, ...){
  result <- object$coefficients
  return(result)
}


#nObs.Rchoice <- function(x, ... ) {
#  return(x$logLik$nobs)
#}

#' @rdname Rchoice
#' @method fitted Rchoice
#' @export
fitted.Rchoice <- function(object, ...){
  result <- object$probabilities
  return(result)
}

#' @rdname Rchoice
#' @method residuals Rchoice
#' @export
residuals.Rchoice <- function(object, ...){
  result <- object$residuals
  result
}


#' @rdname Rchoice
#' @method df.residual Rchoice
#' @import stats
#' @export
df.residual.Rchoice <- function(object, ...){
  n <- length(residuals(object))
  K <- length(coef(object))
  return(n - K)
}

#' @rdname Rchoice
#' @method update Rchoice
#' @import stats
#' @export
update.Rchoice <- function(object, new, ...){
  call <- object$call
  if (is.null(call))
    stop("need an object with call component")
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(new))
    call$formula <- update(formula(object), new)
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, parent.frame())
}



#' @rdname Rchoice
#' @method logLik Rchoice
#' @export
logLik.Rchoice <- function(object,...){
  structure(object$logLik$maximum[[1]], df = length(object$coefficients),
            nobs = object$logLik$nobs, class = "logLik")
}

#' Bread for sandwiches
#' 
#' Computes the ``bread'' of the sandwich covariance matrix for a model of class \code{Rchoice}
#' 
#' @param x a fitted model of class \code{Rchoice},
#' @param ... Other arguments when \code{bread} is applied to another class object.
#' @return the covariance matrix times observations
#' @references Zeileis A (2006), Object-oriented Computation of Sandwich 
#' Estimators. Journal of Statistical Software, 16(9), 1--16.
#' @details For more information see \code{\link[sandwich]{bread}} from the package \pkg{sandwich}.
#' @method bread Rchoice
#' @import stats
#' @export bread.Rchoice
#' @examples
#' ## Probit model
#' data("Workmroz")
#' probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  
#'                   data = Workmroz , family = binomial('probit'))
#' summary(probit)
#' 
#' library("sandwich")
#' bread(probit) 
bread.Rchoice <- function(x, ... ) {
  return( vcov( x ) * x$logLik$nobs)
}

#' Gradient for observations
#' 
#' It extracts the gradient for each observations evaluated at the estimated parameters for a model of class \code{Rchoice}
#' 
#' @param x a fitted model of class \code{Rchoice},
#' @param ... Other arguments when \code{estfun} is applied to another
#' class object
#' @return the gradient matrix of dimension n times k 
#' @references Zeileis A (2006), Object-oriented Computation of Sandwich 
#' Estimators. Journal of Statistical Software, 16(9), 1--16.
#' @details For more information see \code{\link[sandwich]{estfun}} from package \pkg{sandwich}.
#' @method estfun Rchoice
#' @export estfun.Rchoice
estfun.Rchoice <- function(x, ... ) {
  return(x$logLik$gradientObs )
}

#' @rdname Rchoice
#' @method print Rchoice
#' @import stats
#' @export
print.Rchoice <- function(x, digits = max(3,getOption("digits") - 3),
                          width = getOption("width"), ...)
{
  cat("\nCall:\n", deparse(x$call),"\n\n", sep = "")
  
  cat("\nCoefficients:\n")
  print.default(format(coef(x), digits = digits), print.gap = 2,
                quote = FALSE)
  cat("\n")
  invisible(x)
}

#' @rdname Rchoice
#' @method summary Rchoice
#' @import stats
#' @export
summary.Rchoice <- function(object, ...){
  b <- object$coefficients
  std.err <- sqrt(diag(vcov(object)))
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  CoefTable <- cbind(b, std.err, z, p)
  colnames(CoefTable) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)")
  object$CoefTable    <- CoefTable
  
  class(object) <- c("summary.Rchoice","Rchoice")
  return(object)
}

#' @rdname Rchoice
#' @method print summary.Rchoice
#' @import stats
#' @export
print.summary.Rchoice <- function(x, digits = max(3, getOption("digits") - 3),
                                 width = getOption("width"),
                                 ...)
{
  cat(paste("\nModel:", x$family))
  cat(paste("\nModel estimated on:", format(Sys.time(), "%a %b %d %X %Y"), "\n"))
  cat("\nCall:\n")
  cat(paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  if (!(x$family == "poisson")) {
    cat("\nFrequencies of categories:\n")
    print(prop.table(x$freq), digits = digits)
  }
  
  cat(paste("The estimation took:", make.time(x) ,"\n"))
  
  cat("\nCoefficients:\n")
  printCoefmat(x$CoefTable, digits = digits)
  
  cat(paste("\nOptimization of log-likelihood by", x$logLik$type))
  cat(paste("\nLog Likelihood:", signif(x$logLik$maximum, digits)))
  cat(paste("\nNumber of observations:", x$logLik$nobs))
  cat(paste("\nNumber of iterations:" , x$logLik$iterations))
  cat(paste("\nExit of MLE:", x$logLik$message))
  
  if (x$R.model) {
    if (is.null(x$draws)) {
      cat(paste("\nSimulation based on", x$R, "pseudo-random draws"))
    }else {
      cat(paste("\nSimulation based on", x$R, "Halton draws"))
    }
  }
  invisible(x)
}

#' Get Model Summaries for use with "mtable" for object of class Rchoice
#' 
#' A generic function to collect coefficients and summary statistics from a \code{Rchoice} object. It is used in \code{mtable}
#' 
#' @param obj a \code{Rchoice} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @import stats
#' @importFrom memisc getSummary
#' @export 
getSummary.Rchoice <- function(obj, alpha = 0.05, ...){
  smry <- summary(obj)
  coef <- smry$CoefTable
  lower <- coef[, 1] - coef[, 2] * qnorm(alpha/2)
  upper <- coef[, 1] + coef[, 2] * qnorm(alpha/2)
  coef <- cbind(coef, lower, upper)
  colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  N <-  obj$logLik$nobs
  ll <- logLik(obj)
  sumstat <- c(logLik = ll, deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = N, 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA,
               Nagelkerke = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts,
       xlevels = NULL, call = obj$call)
}

##=================================
## Method for Random Parameters
##=================================

#' Plot the distribution of conditional expectation for random parameters.
#' 
#' Plot the distribution of the conditional expectation of the random parameters or compensating variations for objects of class \code{Rchoice}. 
#' 
#' @param x a object of class \code{Rchoice},
#' @param par a string giving the name of the variable with random parameter,
#' @param effect a string indicating what should be plotted: the conditional expectation of the individual coefficients "\code{ce}", or the conditional expectation of the individual compensating variations "\code{cv}",
#' @param wrt a string indicating repect to which variable should be computed the compensating variation,
#' @param type a string indicating the type of distribution: it can be a \code{histogram} or a \code{density} of
#' the conditional expectation,
#' @param ind a boolean. If \code{TRUE}, a 95% interval of conditional distribution for each individual is plotted. 
#' As default, the conditional expectation of \code{par} for the first 10 individual is plotted,
#' @param id only relevant if \code{ind} is not \code{NULL}. This is a vector indicating the individuals for which the confidence intervals are plotted, 
#' @param main an overall title for the plot,
#' @param xlab  a title for the x axis,
#' @param ylab a title for the y axis,
#' @param adjust  bandwidth for the kernel density,
#' @param breaks number of breaks for the histrogram if \code{type = "histogram"},
#' @param col color for the graph, 
#' @param ... further arguments. Ignored.
#' @references
#' \itemize{
#' \item Greene, W. H. (2012). Econometric analysis, Seventh Edition. Pearson Hall.
#' \item Train, K. (2009). Discrete choice methods with simulation. Cambridge university press.
#' }
#' @seealso \code{\link[Rchoice]{Rchoice}} for the estimation of different discrete choice models with individual parameters.
#' @method plot Rchoice
#' @author Mauricio Sarrias
#' @export
#' @import graphics
#' @import stats
#' @examples
#' \donttest{
#' # Poisson with random parameters
#' data("Articles")
#' poisson.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
#'                        data = Articles,  family = poisson,
#'                        ranp = c(kid5 = "n", phd = "n", ment = "n"), 
#'                        R = 10)
#' 
#' ## Plot the distribution of the conditional mean for ment
#' plot(poisson.ran, par = "ment", type = "density")
#' 
#' ## Plot the conditional mean for the first 20 individuals
#' plot(poisson.ran, par = "ment", ind =  TRUE, id = 1:20, col = "blue")
#' 
#' ## Plot the compensating variation with respect to fem
#' plot(poisson.ran, par = "ment", effect = "cv", wrt = "fem", type = "histogram")
#' }
#' @importFrom plotrix plotCI
plot.Rchoice <- function(x, par = NULL, effect = c("ce", "cv"), wrt = NULL,
                         type = c("density", "histogram"), adjust = 1, 
                         main = NULL, col = "indianred1", breaks = 10, ylab = NULL,
                         xlab = NULL, ind = FALSE, id = NULL, ...){
  if (!x$R.model) stop("the plot method is only relevant for random parameters")
  if (is.null(par)) stop("Must specified the name of the random parameters")
  type <- match.arg(type)
  effect <- match.arg(effect)
  if (is.null(xlab)) {
    xlab <- switch(effect,
                   "cv" = expression(E(hat(cv[i]))),
                   "ce"  = expression(E(hat(beta[i]))))
  }
  if (!ind) {
    if (is.null(main)) main <- paste("Conditional Distribution for", par)
    if (is.null(ylab)) {
      ylab <- switch(type,
                     "density"   = "Density",
                     "histogram" = "Frequency")
    }
    rpar <- effect.Rchoice(x, par,  effect = effect, wrt =  wrt)$mean
    if (type == "density") {
      pdens <- density(rpar, adjust = adjust)
      plot(pdens, ylab = ylab, xlab = xlab, main = main, col =  col)
      has.pos <- any(pdens$x > 0)
      if (has.pos) {
        x1 <- min(which(pdens$x >= 0))  
        x2 <- max(which(pdens$x <  max(pdens$x)))
        with(pdens, polygon(x = c(x[c(x1, x1:x2, x2)]), y = c(0, y[x1:x2], 0), 
                            col = col, border = NA))
      }
    } else {
      minb <- round(min(rpar), 2)
      maxb <- round(max(rpar), 2)
      hist(rpar, xlab = xlab, main = main, col = col, breaks = breaks, 
           xaxs = "i", yaxs = "i", las = 1, xaxt = 'n', ylab = ylab)
      axis(1, at = seq(minb, maxb, (maxb - minb) * .05))
    }
  }
  else{
    if (is.null(main)) main <- paste("95% Probability Intervals for ", par)
    if (is.null(id)) id <- seq(1, 10, 1)
    if (is.null(ylab)) ylab <- "Individuals"
    f.bran <- effect.Rchoice(x, par,  effect = effect, wrt =  wrt)$mean
    f.sran <- effect.Rchoice(x, par,  effect = effect, wrt =  wrt)$sd.est
    lower <- f.bran - qnorm(0.975) * f.sran
    upper <- f.bran + qnorm(0.975) * f.sran
    plotrix::plotCI(as.numeric(id), f.bran[id], ui = upper[id], li = lower[id],
                    xlab = ylab, ylab = xlab,
                    lty = 2, main = main,
                    pch = 21, col = col)
  } 
}


#' @rdname vcov.Rchoice
#' @export
cov.Rchoice <- function(x){
  if (!inherits(x, "Rchoice")) stop("not a \"Rchoice\" object")
  if (is.null(x$ranp)) stop('\"cov.Rchoice\"  only relevant for random coefficient model')
  beta.hat <- coef(x)
  K <- length(x$ranp)
  nr <- names(x$ranp)
  if (x$correlation) {
    names.stds <- c()
    for (i in 1:K) names.stds <- c(names.stds, paste('sd', nr[i], nr[i:K], sep = '.'))
    v    <- beta.hat[names.stds]
    V    <- tcrossprod(makeL(v))
    colnames(V) <- rownames(V) <- nr
  } else{
    names.stds <- paste("sd", nr, sep = ".")
    sv   <- beta.hat[names.stds]
    V    <- matrix(0, K, K)
    diag(V) <- sv ^ 2
    colnames(V) <- rownames(V) <- nr
  }
  V
}

#' @rdname vcov.Rchoice
#' @export
cor.Rchoice <- function(x){
  if (!x$correlation) stop('\"cor.Rchoice\"  only relevant for correlated random coefficient model')
  V   <- cov.Rchoice(x)
  nr  <- names(x$ranp)
  D   <- diag(sqrt(diag(V)))
  Rho <- solve(D) %*% V %*% solve(D)
  colnames(Rho) <- rownames(Rho) <- nr
  Rho
}


#' Get the conditional individual coefficients
#' 
#' This a helper function to obtain the individuals' conditional estimate of the random parameters or compensating variations.
#' @param object an object of class \code{Rchoice},
#' @param par a string giving the name of the variable with random parameter,
#' @param effect a string indicating what should be computed: the conditional expectation of the individual coefficients "\code{ce}", or the conditional expectation of the individual compensating variations "\code{cv}",
#' @param wrt a string indicating respect to which variable the compensating variation should be computed,
#' @param ... further arguments. Ignored.
#' 
#' @return A named list where ``mean'' contains the individuals' conditional mean for the random parameter or compensating variation, and where `sd.est' contains their standard errors. 
#' @export
#' @references
#' \itemize{
#' \item Greene, W. H. (2012). Econometric Analysis, Seventh Edition. Pearson Hall.
#' \item Train, K. (2009). Discrete Choice Methods with Simulation. Cambridge university press.
#' }
#' @seealso \code{\link[Rchoice]{Rchoice}} for the estimation of different discrete choice models with individual parameters.
#' @import stats
#' @examples
#' \donttest{
#' # Poisson with random parameters
#' data("Articles")
#' poisson.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
#'                        data = Articles,  family = poisson,
#'                        ranp = c(kid5 = "n", phd = "n", ment = "n"), 
#'                        R = 10)
#' 
#' ## Get the individuals' conditional mean and their standard errors for ment                      
#' bi.ment <- effect(poisson.ran, par = "ment", effect = "ce")
#' summary(bi.ment$mean)
#' summary(bi.ment$sd.est)
#' }
effect.Rchoice <- function(object, par = NULL, effect = c("cv", "ce"), wrt = NULL, ... ){
  if (!inherits(object, "Rchoice")) stop("not a \"Rchoice\" object")
  type <- match.arg(effect)
  ranp <- object$ranp
  if (!is.null(par) && !(par %in% names(ranp))) stop("This parameter is not random: ", par)
  bi   <- object$bi
  Qir  <- object$Qir
  R    <- ncol(Qir)
  N    <- nrow(Qir)
  K <- dim(bi)[[3]]
  
  mean <- mean.sq <- array(NA, dim = c(N, R, K))
  for (j in 1:K) {
    if (type == "cv") {
      # Check if wrt is fixed or random
      if (is.null(wrt)) stop("you need to specify wrt")
      is.ran <- any(names(ranp) %in% wrt)
      gamma <- if (is.ran) bi[, , wrt] else coef(object)[wrt]
      mean[, , j]    <- (bi[, , j] / gamma)   * Qir 
      mean.sq[, , j] <- ((bi[, , j] / gamma) ^ 2 ) * Qir 
    } else {
      mean[, , j]    <- bi[, , j] * Qir 
      mean.sq[, , j] <- (bi[, , j] ^ 2) * Qir
    }  
  }
  mean    <- apply(mean,  c(1,3), sum)
  mean.sq <- apply(mean.sq, c(1,3), sum)
  sd.est  <- suppressWarnings(sqrt(mean.sq - mean ^ 2))
  colnames(mean) <- colnames(mean.sq) <- colnames(sd.est) <- dimnames(bi)[[3]]
  rownames(mean) <- rownames(mean.sq) <- rownames(sd.est) <- dimnames(bi)[[1]]
  if (!is.null(par)) {
    mean   <- mean[, par]
    sd.est <- sd.est[, par]
  }
  effe <- list(mean = mean,
              sd.est = sd.est)
  return(effe)
}

#' @rdname vcov.Rchoice
#' @importFrom msm deltamethod
#' @import stats
#' @export
se.cov.Rchoice <- function(x, sd =  FALSE, digits = max(3, getOption("digits") - 2)){
  if (!inherits(x, "Rchoice")) stop("not a \"Rchoice\" object")
  if (!x$correlation) stop('se.cov.Rchoice only relevant for correlated random coefficient')
  beta.hat <- x$coefficients
  Ka <- length(x$ranp)
  nr <- names(x$ranp)
  names.stds <- c()
  for (i in 1:Ka) names.stds <- c(names.stds, paste('sd', nr[i], nr[i:Ka], sep = '.'))
  stds.hat <- beta.hat[names.stds]
  sel.vcov <- vcov(x)[names.stds, names.stds]
  form <- c()
  if (sd) {
    for (i in 1:Ka) {
      k <- i
      if (i == 1) {
        form <- paste("~ sqrt(", c(form, paste(paste("x",  i, sep = ""), paste("x", k, sep = ""), sep = "*")), ")")
      } else {
        temp <- paste(paste("x",  i, sep = ""), paste("x", k, sep = ""), sep = "*")
        j <- 2
        while (j <= i) {
          temp <- paste(temp, make.add(row = j, col = k, Ka = Ka)[1], sep = "+") 
          j <- j + 1
        }
        form <- c(form, paste("~ sqrt(", temp, ")"))
      }
    }
    b <- sqrt(diag(cov.Rchoice(x)))
    names(b) <- colnames(cov.Rchoice(x))
  } else {
    for (i in 1:Ka) { 
      if (i == 1) {
        form <- paste("~", c(form, paste(paste("x",  i:Ka, sep = ""), paste("x", i, sep = ""), sep = "*")))
      } else {
        temp <- paste(paste("x",  i:Ka, sep = ""), paste("x", i, sep = ""), sep = "*")
        j <- 2
        while (j <= i) {
          temp <- paste(temp, make.add(row = j, col = i, Ka = Ka), sep = "+") 
          j <- j + 1
        }
        form <- c(form, paste("~", temp))
      }
    }
    names.vcov <- c()
    for (i in 1:Ka) names.vcov <- c(names.vcov, paste('v', nr[i], nr[i:Ka], sep = '.'))
    b <- drop(cov.Rchoice(x)[lower.tri(cov.Rchoice(x), diag = TRUE)])
    names(b) <- names.vcov
  }
  std.err <- c()
  for (i in 1:length(form)) {
    std.err <- c(std.err, msm::deltamethod(as.formula(form[i]), stds.hat, sel.vcov, ses =  TRUE))
  }  
  z <- b / std.err
  p <- 2 * (1 - pnorm(abs(z)))
  tableChol <- cbind(b, std.err, z, p)
  if (!sd) cat(paste("\nElements of the variance-covariance matrix \n\n"))
  else cat(paste("\nStandard deviations of the random parameters \n\n"))
  colnames(tableChol) <- c("Estimate", "Std. Error", "z-value", "Pr(>|z|)") 
  printCoefmat(tableChol, digits =  digits)
}

