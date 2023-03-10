#' Estimate Instrumental Variable Probit model by Maximum Likelihood.
#' 
#' Estimation of Probit model with one endogenous and continuous variable by Maximum Likelihood. 
#' 
#' @name ivpml
#' @param formula a symbolic description of the model of the form \code{y ~ x | z} where \code{y} is the binary dependent variable, \code{x} includes the exogenous and the endogenous continuous variable, and \code{z} is the complete set of instruments.  
#' @param data the data of class \code{data.frame}.
#' @param messages if \code{TRUE}, then additional messages for the estimation procedure are displayed. 
#' @param ... arguments passed to \code{maxLik}.
#' @param x,object an object of class \code{ivpml}.
#' @param digits the number of digits.
#' @param eigentol the standard errors are only calculated if the ratio of the smallest and largest eigenvalue of the Hessian matrix is less than \code{eigentol}.  Otherwise the Hessian is treated as singular. 
#' @param newdata optionally, a data frame in which to look for variables with which to predict.
#' @param type the type of prediction required. The default, \code{type = xb}, is on the linear prediction. If \code{type = pr}, the predicted probabilities of a positive outcome is returned. Finally, if \code{type = stdp} the standard errors of the linear predictions for each individual is returned.
#' @param asf  if \code{TRUE}, the average structural function is used. This option is not allowed with \code{xb} or \code{stdp}.
#' @author Mauricio Sarrias.
#' @details 
#' 
#' The IV probit for cross-sectional data has the following structure:
#'
#' \deqn{
#' y_{1i}^*  = x_i^\top\beta + \gamma y_{2i}+ \epsilon_i,
#' }
#' with
#' \deqn{
#'  y_{2i}  = z_i^\top\delta +  \upsilon_i,   
#' }
#' where \eqn{y_{1i}^*} is the latent (unobserved) dependent variable for individual \eqn{i = 1,...,N}; 
#' \eqn{y_{2i}} is the endogenous continuous variable; \eqn{z_i} is the vector of exogenous variables 
#' which also includes the instruments for \eqn{y_{2i}}; 
#' and \eqn{(\epsilon, \upsilon)} are normal jointly distributed.
#'
#' The model is estimated using the \code{maxLik} function from \code{\link[maxLik]{maxLik}} package using 
#' analytic gradient.
#' 
#' @references
#' Greene, W. H. (2012). Econometric Analysis. 7 edition. Prentice Hall.
#' @examples
#' \donttest{ 
#' # Data
#' library("AER")
#' data("PSID1976")
#' PSID1976$lfp  <- as.numeric(PSID1976$participation == "yes")
#' PSID1976$kids <- with(PSID1976, factor((youngkids + oldkids) > 0,
#'                                       levels = c(FALSE, TRUE), 
#'                                       labels = c("no", "yes")))
#'                                       
#' # IV probit model by MLE
#' # (nwincome is endogenous and heducation is the additional instrument)
#' PSID1976$nwincome <- with(PSID1976, (fincome - hours * wage)/1000)
#' fiml.probit <- ivpml(lfp ~  education + experience + I(experience^2) + age + 
#'                             youngkids + oldkids + nwincome |
#'                             education + experience + I(experience^2) + age + 
#'                             youngkids + oldkids + heducation, 
#'                      data = PSID1976)
#' summary(fiml.probit)
#' }
#' @import Formula maxLik stats
#' @keywords models
#' @export
ivpml <- function(formula, data, messages = TRUE, ...){
  callT  <- match.call(expand.dots = TRUE)
  callF  <- match.call(expand.dots = FALSE)
  nframe <- length(sys.calls())
  
  # ============================-
  # 1. Model frame ----
  # ============================-
  mf         <- callT
  m          <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf         <- mf[c(1, m)]
  f          <- Formula(formula)
  if (length(f)[2L] < 2L) stop("No instruments in the second part of formula")
  mf$formula <- f
  mf[[1]]    <- as.name("model.frame")
  mf         <- eval(mf, parent.frame())
  
  # ============================-
  # 2. Variables----
  # ============================-
  y1    <- model.response(mf)
  y.var <- f[[2]]
  X     <- model.matrix(f, data = mf, rhs = 1)
  Z     <- model.matrix(f, data = mf, rhs = 2)
  y2    <- X[, !(colnames(X) %in% colnames(Z)), drop = FALSE]
  
  if (ncol(Z) < ncol(X)) stop("Model underidentified: check your variables in formula")
  if (messages && (ncol(Z) == ncol(X))) cat("\nEstimating a just identified model....\n")
  if (messages && (ncol(Z) > ncol(X)))  cat("\nEstimating an overidentified model....\n")
  end.var     <- colnames(y2)
  instruments <- colnames(Z)
  if (length(end.var) > 1L) stop("ivpml only works with one endogenous variable")
  
  # Check dependent variable
  if (!all(y1 %in% c( 0, 1, TRUE, FALSE))){
    stop( "all dependent variables must be either 0, 1, TRUE, or FALSE")
  }
  if (!is.numeric(y1)) y1 <- as.numeric(y1) 
  
  # ============================-
  # 3. Initial values----
  # ============================-
  # These are different as those in STATA for the first equation
  if (messages) cat("\nObtaining starting values from probit and linear model...\n")
  probit       <- glm.fit(X, y1,  family = binomial("probit"))
  linear       <- lm.fit(Z, y2)
  sigma2       <- sum(resid(linear) ^ 2) / linear$df.residual
  lnsigma      <- log(sqrt(sigma2))
  rho          <- cor(resid(linear), resid(probit))
  athrho       <- 0.5 * log((1 + rho)/(1 - rho))
  start        <- c(coef(probit), coef(linear), lnsigma, athrho)
  nam_prob     <- paste(y.var, colnames(X), sep = ":")
  nam_lin      <- paste(end.var, colnames(Z), sep = ":")
  names(start) <- c(nam_prob, nam_lin,  "lnsigma", "atanhrho")
  
  # ============================-
  # 4. Optimization----
  # ============================-
  if (is.null(callT$method))  callT$method   <- 'nr'
  opt <- callT
  m   <- match(c("print.level", "ftol", "tol", "reltol",
                 "gradtol", "steptol", "lambdatol", "qrtol",
                 "iterlim", "fixed", "activePar", "method", "control", "constraints", "gradient", "hessian"),
               names(opt), 0L)
  opt        <- opt[c(1L, m)]
  opt$start  <- start
  opt[[1]]   <- as.name('maxLik')
  opt$logLik <- as.name('lnbinary_iv')
  opt[c('y1', 'y2', 'X', 'Z')] <- list(as.name('y1'), as.name('y2'), as.name('X'), as.name('Z'))
  out <- eval(opt, sys.frame(which = nframe))
  
  # ============================-
  # 5. Save results----
  # ============================-
  out$end.var     <- end.var
  out$y.var       <- y.var
  out$instruments <- instruments
  out$formula     <- f
  out$mf          <- mf
  out$call        <- callT
  #class(out)      <- c("ivpml", "maxLik", class(out))
  class(out)      <- c("ivpml", class(out)) # JSS reviewer
  return(out)
}



############################---
# S3 method for ivpml class ----
#############################---

#' @rdname ivpml
#' @method terms ivpml
#' @export
terms.ivpml <- function(x, ...){
  formula(x$formula)
}

#' @rdname ivpml
#' @method model.matrix ivpml
#' @import stats
#' @export
model.matrix.ivpml <- function(object, ...){
  f   <- object$formula
  mf  <- object$mf
  X   <- model.matrix(f, data = mf, rhs = 1)
  Z   <- model.matrix(f, data = mf, rhs = 2)
  y2  <- X[, !(colnames(X) %in% colnames(Z)), drop = FALSE]
  out <- list(X = X, Z = Z, y2 = y2)
  return(out)
}

#' @rdname ivpml
#' @method estfun ivpml
#' @importFrom sandwich estfun
#' @export estfun.ivpml
estfun.ivpml <- function(x, ...){
  class(x) <- c("maxLik", "maxim")
  estfun(x, ...)
}

#' @rdname ivpml
#' @method bread ivpml
#' @importFrom sandwich bread
#' @export bread.ivpml
bread.ivpml <- function(x, ...){
  class(x) <- c("maxLik", "maxim")
  bread(x, ...)
}

# #' @rdname ivpml
# #' @import stats
# #' @method AIC ivpml
# #' @export
# AIC.ivpml <- function(object, k = 2, ...){
#   -2*logLik(object) + k * length(coef(object))
# }
# 
# #' @rdname ivpml
# #' @import stats
# #' @method BIC ivpml
# #' @export
# BIC.ivpml <- function(object, ...){
#   AIC(object, k = log(nrow(object$gradientObs)), ...)
# }


#' @rdname ivpml
#' @method vcov ivpml
#' @import stats
#' @export 
vcov.ivpml <- function(object, ...){
  class(object) <- c("maxLik", "maxim")
  vcov(object, ...)
}

#' @rdname ivpml
#' @import stats
df.residual.ivpml <- function(object, ...){
  return(row(object$gradientObs) - length(coef(object)))
}

#' @rdname ivpml
#' @export
coef.ivpml <- function(object, ...){
  class(object) <- c("maxLik", "maxim")
  coef(object, ...)
}


#' @rdname ivpml
#' @export 
logLik.ivpml <- function(object, ...){
  structure(object$maximum, df = length(coef(object)), nobs = nrow(object$gradientObs), class = "logLik")
}

#' @rdname ivpml
#' @method print ivpml
#' @import stats
#' @export 
print.ivpml <- function(x, ...){
  cat("Maximum Likelihood estimation\n")
  cat(maximType(x), ", ", nIter(x), " iterations\n", sep = "")
  cat("Return code ", returnCode(x), ": ", returnMessage(x), 
      "\n", sep = "")
  if (!is.null(x$estimate)) {
    cat("Log-Likelihood:", x$maximum)
    cat(" (", sum(activePar(x)), " free parameter(s))\n", 
        sep = "")
    cat("Estimate(s):", x$estimate, "\n")
  }
}

#' @rdname ivpml
#' @method summary ivpml
#' @import stats
#' @importFrom miscTools stdEr
#' @export
summary.ivpml <- function(object, eigentol = 1e-12, ...){
  result <- object$maxim
  nParam <- length(coef(object))
  activePar <- activePar(object)
  if ((object$code < 100) & !is.null(coef(object))) {
    t <- coef(object)/stdEr(object, eigentol = eigentol)
    p <- 2 * pnorm(-abs(t))
    t[!activePar(object)] <- NA
    p[!activePar(object)] <- NA
    results <- cbind(Estimate = coef(object), `Std. error` = stdEr(object, 
                                                                   eigentol = eigentol), `z value` = t, `Pr(> z)` = p)
  }
  else {
    results <- NULL
  }
  summary <- list(maximType = object$type, iterations = object$iterations, 
                  returnCode = object$code, returnMessage = object$message, 
                  loglik = object$maximum, estimate = results, fixed = !activePar, 
                  NActivePar = sum(activePar), constraints = object$constraints, 
                  end.var = object$end.var, instruments = object$instruments)
  class(summary) <- "summary.ivpml"
  summary
}

#' @rdname ivpml
#' @method print summary.ivpml
#' @import stats
#' @export
print.summary.ivpml <- function(x, digits = max(3, getOption("digits") - 2),
                                ...){
  cat("--------------------------------------------\n")
  cat("Maximum Likelihood estimation of IV Probit model \n")
  cat(maximType(x), ", ", nIter(x), " iterations\n", sep = "")
  cat("Return code ", returnCode(x), ": ", returnMessage(x), 
      "\n", sep = "")
  if (!is.null(x$estimate)) {
    cat("Log-Likelihood:", x$loglik, "\n")
    cat(x$NActivePar, " free parameters\n")
    cat("Estimates:\n")
    printCoefmat(x$estimate, digits = digits)
    cat("\nInstrumented:", x$end.var)
    cat("\nInstruments:",  x$instruments, "\n")
    cat("Wald test of exogeneity (corr = 0): chi2", round(x$estimate["atanhrho",3]^2, 2),
        "with 1 df. Prob > chi2 = ", round(pchisq(x$estimate["atanhrho",3]^2,1, lower.tail =  FALSE), 4), "\n")
  }
  
  if (!is.null(x$constraints)) {
    cat("\nWarning: constrained likelihood estimation.", 
        "Inference is probably wrong\n")
    cat("Constrained optimization based on", x$constraints$type, 
        "\n")
    if (!is.null(x$constraints$code)) 
      cat("Return code:", x$constraints$code, "\n")
    if (!is.null(x$constraints$message)) 
      cat(x$constraints$message, "\n")
    cat(x$constraints$outer.iterations, " outer iterations, barrier value", 
        x$constraints$barrier.value, "\n")
  }
  cat("--------------------------------------------\n")
}

############################-
# Effects and other functions ----
#############################-

#' @rdname ivpml
#' @method predict ivpml
#' @export
predict.ivpml <- function(object, newdata = NULL, 
                          type = c("xb", "pr", "stdp"), 
                          asf = TRUE, 
                          ...){
  # xb: linear prediction excluding endogeneity; the default
  # pr: probability of a positive outcome
  # stdp: standard error of the linear prediction
  # asf: average structural function; the default. This option
  #     is not allowed with xb or stdp
  type  <- match.arg(type)
  mf    <- if (is.null(newdata)) object$mf else newdata
  f     <- object$formula
  X     <- model.matrix(f, data = mf, rhs = 1)
  Z     <- model.matrix(f, data = mf, rhs = 2)
  y2    <- X[, !(colnames(X) %in% colnames(Z)), drop = FALSE]
  K        <- ncol(X)
  P        <- ncol(Z)
  param    <- coef(object)
  beta     <- param[1L:K]
  delta    <- param[(K + 1L):(K + P)]
  lnsigma  <- param[(K + P + 1L)]
  atanhrho <- tail(param, n = 1L)
  sigma    <- exp(lnsigma)
  rho      <- tanh(atanhrho)
  index1   <- crossprod(t(X), beta)
  index2   <- crossprod(t(Z), delta)
  mi       <- (index1 + (rho / sigma) * (y2 - index2)) / sqrt(1 - rho ^ 2)
  if (type == "pr"){
    if (asf){
      out <- as.vector(pnorm(mi))
    } else {
      out <- as.vector(pnorm(index1))
    }
  }
  if (type == "xb"){
     out <- as.vector(index1)
  }
  if (type == "stdp"){
     out <- c()
     for (i in 1:nrow(X)){
       out <- c(out, sqrt(X[i, , drop = F]%*% vcov(object)[1:K, 1:K] %*% t(X[i, , drop = F])))
     }
  }
  return(out)
}


#' Get Model Summaries for use with "mtable" for objects of class ivpml
#' 
#' A generic function to collect coefficients and summary statistics from a \code{ivpml} object. It is used in \code{mtable}
#' 
#' @param obj a \code{ivpml} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @import stats
#' @importFrom memisc getSummary
#' @export 
getSummary.ivpml <- function(obj, alpha = 0.05, ...){
  s       <- summary(obj)$estimate
  f       <- obj$formula
  mf      <- obj$mf
  X       <- model.matrix(f, data = mf, rhs = 1)
  Z       <- model.matrix(f, data = mf, rhs = 2)
  end.var <- obj$end.var
  y.var   <- obj$y.var
  cf.eq1  <- s[rownames(s) %in% c(paste(y.var, colnames(X), sep = ":"), "lnsigma", "atanhrho"), ]
  cf.eq2  <- s[rownames(s) %in%  paste(end.var, colnames(Z), sep = ":"), ]
  cval    <- qnorm(1 - alpha/2)
  cf.eq1  <- cbind(cf.eq1, cf.eq1[, 1] - cval * cf.eq1[, 2], cf.eq1[, 1] + cval * cf.eq1[, 2])
  cf.eq2  <- cbind(cf.eq2, cf.eq2[, 1] - cval * cf.eq2[, 2], cf.eq2[, 1] + cval * cf.eq2[, 2])
  rownames(cf.eq1) <- c(colnames(X), "lnsigma", "atanhrho")
  rownames(cf.eq2) <- colnames(Z)
  all.vars    <- unique(c(colnames(X), colnames(Z), "lnsigma", "atanhrho"))
  # Make Table
  coef        <- array(dim = c(length(all.vars), 6, 2), 
                       dimnames = list(all.vars, c("est", "se", "stat", "p", "lwr", "upr"), c(y.var, end.var)))
  coef[rownames(cf.eq1),,1]    <- cf.eq1
  coef[rownames(cf.eq2),,2]    <- cf.eq2
  
  # Statistics
  sumstat <- c(logLik = logLik(obj), deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = nrow(obj$gradientObs), 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts, xlevels = obj$xlevels, call = obj$call)
}



#### Effects

#' Get Model Summaries for use with "mtable" for objects of class effect.ivpml
#' 
#' A generic function to collect coefficients and summary statistics from a \code{effect.ivpml} object. It is used in \code{mtable}
#' 
#' @param obj an \code{effect.ivpml} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @import stats
#' @importFrom memisc getSummary
#' @export
getSummary.effect.ivpml <- function(obj, alpha = 0.05, ...){
  cf             <- summary(obj)$CoefTable
  cval           <- qnorm(1 - alpha/2)
  coef           <- cbind(cf, cf[, 1] - cval * cf[, 2], cf[, 1] + cval * cf[, 2])
  dim(coef)      <- c(dim(coef)[1], dim(coef)[2], 1)
  dimnames(coef) <- list(rownames(cf), c("est", "se", "stat", "p", "lwr", "upr"), all.vars(obj$formula)[1])
  #colnames(coef) <- c("est", "se", "stat", "p", "lwr", "upr")
  # Statistics
  sumstat <- c(logLik = obj$maximum, deviance = NA, AIC = NA, BIC = NA, N = nrow(obj$gradientObs), 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA)
  list(coef = coef, sumstat = sumstat, contrasts = NULL, xlevels = NULL, call = obj$call)
}


mdydx.ivpml <- function(coeff, object, asf){
  f     <- object$formula
  mf    <- object$mf
  X     <- model.matrix(f, data = mf, rhs = 1)
  Z     <- model.matrix(f, data = mf, rhs = 2)
  y2    <- X[, !(colnames(X) %in% colnames(Z)), drop = FALSE]
  bhat  <- coeff
  sigma <- exp(bhat["lnsigma"])
  rho   <- tanh(bhat["atanhrho"])
  r     <- (rho / sigma)
  sr    <- sqrt(1 - rho^2)
  K     <- ncol(X)
  P     <- ncol(Z)
  beta  <- bhat[1L:K]
  delta <- bhat[(K + 1L):(K + P)]
  
  # Obtain classes    
  all.var  <- all.vars(formula(f, rhs = 1))[-1L]
  classes  <- rep("numeric", length(all.var))
  class.mf <- attributes(terms(mf))[["dataClasses"]][-1L]
  classes[paste0("factor(", all.var, ")") %in% names(class.mf)]  <- class.mf[names(class.mf) %in% paste0("factor(", all.var, ")")]
  names(beta) <- colnames(X)
  
  ## Compute marginal effects
  mes      <- c()
  mes.name <- c()
  for (k in 1:length(all.var)){
    if (classes[k] == "numeric") {
      xb       <- crossprod(t(X), beta)
      upsilon  <- y2 - crossprod(t(Z), delta)
      mi       <- (xb + r * upsilon) / sr
      bk       <- make.inter.num(all.var[k], names(beta), beta, X)
      if (asf) {
        me <- dnorm(mi) * bk / sr
      } else {
        me <- dnorm(xb) * bk
      }
      mes     <- cbind(mes, me)
      mes.name <- c(mes.name, all.var[k])
    }
    if (classes[k] == "factor"){
      levs <- attributes(mf[, paste0("factor(", all.var[k], ")")])$levels
      levs <- levs[-1L]
      ## Make P0
      beta.temp  <- beta
      vb <- make.inter.factor(all.var[k], names(beta), levs)
      if (any(vb$names %in% names(beta)))  beta.temp[names(beta)   %in% vb$names] <- 0
      xb       <- crossprod(t(X), beta.temp)
      upsilon  <- y2 - crossprod(t(Z), delta)
      mi       <- (xb + r * upsilon) / sr
      p0       <- pnorm(mi)
      for (j in 1:length(levs)){
        ## Make P1
        Xtemp  <- X
        if (any(vb$names %in% names(beta)))   Xtemp[, names(beta)  %in% vb$names] <- 0 
        vbj <- make.inter.factor(all.var[k], names(beta), levs[j])
        if (any(vbj$names %in% names(beta)))  Xtemp[, names(beta) %in% vbj$names] <- X[, vbj$names.inte] 
        if (vbj$names[1] %in% names(beta))    Xtemp[, names(beta) %in% vbj$names[1]] <- 1
        xbt <- crossprod(t(Xtemp), beta)
        mit <- (xbt + r * upsilon) / sr
        p1  <- pnorm(mit)
        me  <- p1 - p0
        mes <- cbind(mes, me)
        mes.name <- c(mes.name, paste0("factor(",all.var[k],")",levs[j], sep = ""))
      }
    }
  }
  colnames(mes) <- mes.name 
  mes <- colMeans(mes)
  return(mes)
}


#### Likelihood function ----
#' @importFrom utils tail
lnbinary_iv <- function(param, y1, y2, X, Z, gradient = TRUE, hessian = TRUE){
  # Likelihood function for binary Probit IV model with cross-sectional 
  #   data and one continuous endogenous variable
  F        <- pnorm
  f        <- dnorm
  #ff       <- function(x) -x * dnorm(x)
  K        <- ncol(X)
  P        <- ncol(Z)
  beta     <- param[1L:K]
  delta    <- param[(K + 1L):(K + P)]
  lnsigma  <- param[(K + P + 1L)]
  atanhrho <- tail(param, n = 1L)
  sigma    <- exp(lnsigma)
  rho      <- tanh(atanhrho)
  index1   <- crossprod(t(X), beta)
  index2   <- crossprod(t(Z), delta)
  q        <- 2 * y1 - 1
  ai       <- as.numeric(q * (index1 + (rho / sigma) * (y2 - index2)) / sqrt(1 - rho ^ 2))
  P1       <- F(ai)
  bi       <- as.numeric((y2 - index2) / sigma)
  P2       <- (1 / sigma) * f(bi)
  Pi       <- pmax(P1 * P2, .Machine$double.eps)
  Li       <- log(Pi)
  
  if (gradient){
    m        <- function(x) f(x) / F(x)
    gb       <-  X * (m(ai) * q  / sqrt(1 - rho^2))
    gd       <-  Z * (- m(ai) * q * (rho / sigma) / sqrt(1 - rho^2) + (bi / sigma))
    glnsigma <- -(m(ai) * q * rho / sqrt(1 - rho ^ 2)) * bi + bi^2 - 1
    gathrho  <- m(ai) * q * (index1 * rho + bi) / sqrt(sech(atanhrho)^2) 
    attr(Li,'gradient') <- cbind(gb, gd, glnsigma, gathrho)
  }
  if (hessian){
    h <- function(x) - x * m(x) - m(x)^2
    H <- matrix(0, K + P + 2, K + P + 2)
    cqr <- q / sqrt(1 - rho^2)
    bb <- h(ai) * cqr^2
    H[1:K, 1:K] <- crossprod(bb * X, X)
    bd <- - h(ai) * cqr^2 * (rho / sigma)
    H[1:K, (K + 1):(K + P)] <- crossprod(bd * X, Z)
    bls <- -h(ai) * cqr^2 * rho * bi
    H[1:K, (K + P + 1)] <- colSums(bls * X)
    bt <- h(ai) * cqr * (q * (index1 * rho + bi) / sqrt(sech(atanhrho)^2))
    H[1:K, (K + P + 2)] <- colSums(as.numeric(bt) * X)
    H[(K + 1):(K + P), 1:K] <- t(H[1:K, (K + 1):(K + P)])
    cqr2 <- q  * (rho / sigma) / sqrt(1 - rho^2)
    dd <- h(ai) * cqr2^2 - (1 / sigma^2)
    H[(K + 1):(K + P), (K + 1):(K + P)] <- crossprod(dd * Z, Z)
    dls <- (bi / sigma)* (h(ai)* (q * rho / sqrt(1 - rho^2))^2 -2)      
    H[(K + 1):(K + P), (K + P + 1)] <- colSums(dls * Z)
    dt <- -h(ai) * cqr2 * (q * (index1 * rho + bi) / sqrt(sech(atanhrho)^2))
    H[(K + 1):(K + P), (K + P + 2)] <- colSums(as.numeric(dt) * Z)
    H[(K + P + 1), 1:K]<- t(H[1:K, (K + P + 1)])
    H[(K + P + 1), (K + 1):(K + P)] <- t(H[(K + 1):(K + P), (K + P + 1)]) 
    H[K + P  + 1, K + P + 1] <- sum(h(ai) * (q * rho * bi / sqrt(1 - rho^2))^2 + m(ai) * (q * rho * bi / sqrt(1 - rho^2)) - 2 * bi^2)
    H[K + P  + 1, K + P + 2] <- sum(- bi * (h(ai) * (q * rho / sqrt(1 - rho^2)) *  (q * (index1 * rho + bi) / sqrt(sech(atanhrho)^2)) + m(ai) * (q / sqrt(sech(atanhrho)^2))))
    H[K + P + 2, 1:K] <- t(H[1:K, (K + P + 2)])
    H[K + P + 2, (K + 1):(K + P)] <- t(H[(K + 1):(K + P), (K + P + 2)])
    H[K + P + 2, K + P + 1] <- H[K + P  + 1, K + P + 2]
    H[K + P + 2, K + P + 2] <- sum(h(ai) * (q * (index1 * rho + bi) / sqrt(sech(atanhrho)^2))^2 + q * m(ai) * (index1 + bi * rho)/sqrt(sech(atanhrho)^2))
    attr(Li, 'hessian') <- H
  }
  return(Li)
}

sech <- function(x) 1 / cosh(x)
