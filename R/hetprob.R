#' Estimate heterokedastic binary (Probit or Logit) model.
#' 
#' Estimation of binary dependent variables, either probit or logit, with heteroskedastic error terms for cross-sectional dataset. 
#'
#' @name hetprob
#' @param x,object an object of class \code{hetprob}.
#' @param formula a symbolic description of the model of the form \code{y ~ x | z} where \code{y} is the binary dependent variable and \code{x} and \code{z} are regressors variables for the mean of the model and lnsigma.
#' @param data the data of class \code{data.frame}.
#' @param link the assumption of the distribution of the error term. 
#' @param Hes logical. Should the analytic Hessian to be used? \code{TRUE} as default.
#' @param digits the number of digits.
#' @param eigentol the standard errors are only calculated if the ratio of the smallest and largest eigenvalue of the Hessian matrix is less than \code{eigentol}.  Otherwise the Hessian is treated as singular. 
#' @param k a numeric value, use as penalty coefficient for number of parameters in the fitted model
#' @param newdata optionally, a data frame in which to look for variables with which to predict.
#' @param type the type of prediction required. The default, \code{type = xb}, is on the linear prediction without the variance. If \code{type = pr}, the predicted probabilities of a positive outcome is returned. Finally, if \code{type = sigma} the predictions of \eqn{\sigma} for each individual is returned.
#' @param ... arguments passed to \code{maxLik}. 
#' 
#' @details 
#' 
#' The heterokedastic binary model for cross-sectional data has the following structure:
#'
#' \deqn{
#' y_i^*  = x_i^\top\beta + \epsilon_i,
#' }
#' with
#' \deqn{
#' var(\epsilon_i|x_i, z_i)  = \sigma_i^2 = \left[\exp\left(z_i^\top\delta\right)\right]^2,   
#' }
#' where \eqn{y_i^*} is the latent (unobserved) dependent variable for individual \eqn{i = 1,...,N}; 
#' \eqn{x_i} is a \eqn{K\times 1} vector of independent variables determining the latent variable \eqn{y_i^*} (\code{x} variables in \code{formula}); 
#' and \eqn{\epsilon_i} is the error term distributed either normally or logistically with \eqn{E(\epsilon_i|z_i, x_i) = 0} 
#' and heterokedastic variance \eqn{var(\epsilon_i|x_i, z_i)  = \sigma_i^2, \forall i = 1,...,N}. 
#' The variance for each individual is modeled parametrically assuming that it depends on a \eqn{P\times 1} 
#' vector observed variables \eqn{z_i} (\code{z} in \code{formula}), whereas \eqn{\delta} is the vector of parameters associated with each variable. 
#' It is important to emphasize that \eqn{z_i} does not include a constant, otherwise the parameters are not identified.
#'
#' The models are estimated using the \code{maxLik} function from \code{\link[maxLik]{maxLik}} package using both 
#' analytic gradient and hessian  (if \code{Hess = TRUE}). In particular, the log-likelihood function is:
#' 
#' \deqn{\log L(\theta) = \sum_i^n\log \left\lbrace \left[1- F\left(\frac{x_i^\top\beta}{\exp(z_i^\top\delta)}\right)\right]^{1-y_i}\left[F\left(\frac{x_i^\top\beta}{\exp(z_i^\top\delta)}\right)\right]^{y_i}\right\rbrace.}
#' 
#' @return An object of class ``\code{hetprob}'', a list elements:
#' \item{logLik0}{logLik for the homokedastic model,}
#' \item{f1}{the formula,}
#' \item{mf}{the model framed used,}
#' \item{call}{the matched call.} 
#' @author Mauricio Sarrias. 
#' @import Formula maxLik stats
#' @examples
#' \dontrun{
#' # Estimate a heterokedastic probit and logit model
#' data("Health")
#' 
#' het.probit <- hetprob(working ~ factor(female) + factor(year) + educ + age + I(age^2) | 
#'                                 factor(female) + age + I(age^2), 
#'                      data = Health, 
#'                      link = "probit")
#' summary(het.probit)
#'
#' het.logit <- hetprob(working ~ factor(female) + factor(year) + educ + age + I(age^2) | 
#'                                factor(female) + age + I(age^2), 
#'                     data = Health, 
#'                     link = "logit")
#' summary(het.logit)
#' }
#' @export
hetprob <- function(formula, 
                    data, 
                    link = c("probit", "logit"), 
                    Hes =  TRUE, ...){
  callT  <- match.call(expand.dots = TRUE)
  callF  <- match.call(expand.dots = FALSE)
  nframe <- length(sys.calls())
  link   <- match.arg(link)
  
  ##############################
  ## 1. Model frame 
  ##############################
  mf         <- callT
  m          <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf         <- mf[c(1L, m)]
  f1         <- Formula(formula)
  #Check if there is a second part
  if (length(f1)[2L] < 2L) f1 <- as.Formula(formula(f1), formula(f1, lhs = 0L))
  mf$formula <- f1
  mf[[1]]    <- as.name("model.frame")
  mf         <- eval(mf, parent.frame())
  
  ##############################
  ## 2. Variables
  ##############################
  y <- model.response(mf)
  X <- model.matrix(f1, data = mf, rhs = 1)
  Z <- model.matrix(f1, data = mf, rhs = 2)
  
  ### Drop intercept on Z if included
  zint <- match("(Intercept)", dimnames(Z)[[2]], nomatch = 0)
  if (zint > 0) Z <- Z[, - zint, drop = FALSE]
  
  ##############################
  ## 3. Initial values
  ##############################
  aux1    <- glm(y ~ X - 1, data = mf, family = binomial(link))
  logLik0 <- logLik(aux1)
  betas   <- coef(aux1)
  #aux2   <- lm(residuals(aux1) ~ Z - 1) 
  #gammas <- exp(coef(aux2))
  gammas  <- rep(0, ncol(Z))
  theta   <- c(betas, gammas)
  names(theta) <- c(colnames(X), paste0("het.", colnames(Z)))

  ##############################
  ## 4. Optimization
  ##############################
  if (is.null(callT$method))  callT$method   <- 'nr'
  opt <- callT
  m   <- match(c("print.level", "ftol", "tol", "reltol",
                 "gradtol", "steptol", "lambdatol", "qrtol",
                 "iterlim", "fixed", "activePar", "method", "control", "constraints"),
               names(opt), 0L)
  opt        <- opt[c(1L, m)]
  opt$start  <- theta
  opt[[1]]   <- as.name('maxLik')
  opt$logLik <- as.name('lnbinary_het')
  opt$Hes    <- as.name('Hes')
  opt$link   <- as.name('link')
  opt[c('y', 'X', 'Z')] <- list(as.name('y'), as.name('X'), as.name('Z'))
  out <- eval(opt, sys.frame(which = nframe))
  
  ##############################
  ## 5. Return results
  ##############################
  out$logLik0     <- logLik0 
  out$formula     <- f1
  out$mf          <- mf
  out$call        <- callT
  class(out)      <- c("hetprob", "maxLik", class(out))
  return(out)
}

## Log-likelihood function ====
lnbinary_het <- function(theta, y, X, Z, 
                         Hes =  TRUE, 
                         link = c("probit", "logit")){
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  dfun <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  K     <- ncol(X)
  beta  <- theta[1:K]
  gamma <- theta[-c(1:K)]
  J     <- length(gamma)
  index <- tcrossprod(X, t(beta)) # n x 1
  q     <- 2 * y - 1
  het   <- exp(-tcrossprod(Z, t(gamma)))
  ai    <- q * index * het
  pi    <- pmax(pfun(ai), .Machine$double.eps)
  ll    <- sum(log(pi))
  
  ## Gradient
  mill  <- switch(link,
                  "probit" = function(x) dfun(x) / pmax(pfun(x), .Machine$double.eps),
                  "logit"  = function(x) 1 - pmax(pfun(x), .Machine$double.eps))
  a_theta     <- drop(q * het) * cbind(X, as.vector(-index) * Z)
  G           <- as.vector(mill(ai)) * a_theta
  colnames(G) <- names(theta)
  attr(ll,'gradient') <- G
  
  ## Hessian
  if (Hes) {
    h      <- switch(link,
                     "probit" = function(x) -x * mill(x) - mill(x) ^ 2,
                     "logit"  = function(x) -pfun(x) * (1 - pfun(x)))
    Htemp <- matrix(0, K + J, K + J)
    Htemp[1:K, -c(1:K)]   <- crossprod(as.vector(mill(ai) * q * het) * X, -Z)
    Htemp[-(1:K), 1:K]    <- crossprod(as.vector(mill(ai) * q * het) * -Z, X)
    Htemp[-(1:K), -(1:K)] <- crossprod(as.vector(mill(ai) * q * index * het) * Z, Z)
    H <- crossprod(as.vector(h(ai)) * a_theta, a_theta) + Htemp
    attr(ll,'hessian') <- H
  }
  return(ll)
}


############################
# S3 method for hetprob class
#############################

#' @rdname hetprob
#' @method terms hetprob
#' @export
terms.hetprob <- function(x, ...){
  formula(x$formula)
}

#' @rdname hetprob
#' @method model.matrix hetprob
#' @import stats
#' @export
model.matrix.hetprob <- function(object, ...){
  X <- model.matrix(object$formula, data = object$mf, rhs = 1)
  Z <- model.matrix(object$formula, data = object$mf, rhs = 2)
  ### Drop intercept on Z if included
  zint <- match("(Intercept)", dimnames(Z)[[2]], nomatch = 0)
  if (zint > 0) Z <- Z[, - zint, drop = FALSE]
  out <- list(X = X, Z = Z)
  return(out)
}

#' @rdname hetprob
#' @method estfun hetprob
#' @importFrom sandwich estfun
#' @export estfun.hetprob
estfun.hetprob <- function(x, ...){
  class(x) <- c("maxLik", "maxim")
  estfun(x, ...)
}

#' @rdname hetprob
#' @method bread hetprob
#' @importFrom sandwich bread
#' @export bread.hetprob
bread.hetprob <- function(x, ...){
  class(x) <- c("maxLik", "maxim")
  bread(x, ...)
}

#' @rdname hetprob
#' @import stats
#' @method AIC hetprob
#' @export
AIC.hetprob <- function(object, k = 2, ...){
  -2*logLik(object) + k * length(coef(object))
}

#' @rdname hetprob
#' @import stats
#' @method BIC hetprob
#' @export
BIC.hetprob <- function(object, ...){
  AIC(object, k = log(nrow(object$gradientObs)), ...)
}


#nObs.hetprob <- function(x, ...){
#  return(nrow(x$gradientObs))
#}

#' @rdname hetprob
#' @method vcov hetprob
#' @import stats
#' @export 
vcov.hetprob <- function(object, eigentol = 1e-12, ...){
  class(object) <- c("maxLik", "maxim")
  #vcov(object, eigentol = 1e-12, ...)
  vcov(object, eigentol = eigentol, ...)
}

#' @rdname hetprob
#' @import stats
df.residual.hetprob <- function(object, ...){
  return(nrow(object$gradientObs) - length(coef(object)))
}

#' @rdname hetprob
#' @export
coef.hetprob <- function(object, ...){
  class(object) <- c("maxLik", "maxim")
  coef(object, ...)
}

#' @rdname hetprob
#' @export 
logLik.hetprob <- function(object, ...){
  structure(object$maximum, df = length(coef(object)), nobs = nrow(object$gradientObs), class = "logLik")
}

#' @rdname hetprob
#' @method print hetprob
#' @import stats
#' @export 
print.hetprob <- function(x, ...){
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

#' @rdname hetprob
#' @method summary hetprob
#' @import stats
#' @importFrom miscTools stdEr
#' @export
summary.hetprob <- function(object, eigentol = 1e-12, ...){
  result    <- object$maxim
  nParam    <- length(coef(object))
  activePar <- activePar(object)
  if ((object$code < 100) & !is.null(coef(object))) {
    K <- ncol(model.matrix(object)$X)
    t <- coef(object)/stdEr(object, eigentol = eigentol)
    p <- 2 * pnorm(-abs(t))
    t[!activePar(object)] <- NA
    p[!activePar(object)] <- NA
    results.mean    <- cbind(Estimate = coef(object)[1:K], `Std. error` = stdEr(object, 
                                                                   eigentol = eigentol)[1:K], `z value` = t[1:K], `Pr(> z)` = p[1:K])
    results.lnsigma <- cbind(Estimate = coef(object)[-c(1:K)], `Std. error` = stdEr(object, 
                                                                   eigentol = eigentol)[-c(1:K)], `z value` = t[-c(1:K)], `Pr(> z)` = p[-c(1:K)])
  }
  else {
    results <- NULL
  }
  summary <- list(maximType = object$type, iterations = object$iterations, 
                  returnCode = object$code, returnMessage = object$message, 
                  loglik = object$maximum, results.mean = results.mean, results.lnsigma = results.lnsigma, fixed = !activePar, 
                  NActivePar = sum(activePar), constraints = object$constraints, logLik0 = object$logLik0, 
                  logLik1 = logLik(object))
  class(summary) <- "summary.hetprob"
  summary
}

#' @rdname hetprob
#' @method print summary.hetprob
#' @import stats
#' @export
print.summary.hetprob <- function(x, 
                                  digits = max(3, getOption("digits") - 2),
                                ...){
  cat("------------------------------------------------------------------\n")
  cat("Maximum Likelihood estimation of Heteroskedastic Binary model \n")
  cat(maximType(x), ", ", nIter(x), " iterations\n", sep = "")
  cat("Return code ", returnCode(x), ": ", returnMessage(x), 
      "\n", sep = "")
  if (!is.null(x$results.mean)) {
    cat("Log-Likelihood:", x$loglik, "\n")
    cat(x$NActivePar, " free parameters\n")
    cat(paste("\nEstimates for the mean:\n"), sep = "")
    printCoefmat(x$results.mean, digits = digits)
    cat(paste("\nEstimates for lnsigma:\n"), sep = "")
    printCoefmat(x$results.lnsigma, digits = digits)
    cat("\nLR test of lnsigma = 0: chi2", round(2*(x$logLik1 - x$logLik0), 2),
        "with", attributes(x$logLik1)[["df"]] - attributes(x$logLik0)[["df"]], "df. Prob > chi2 = ", round(pchisq(2*(x$logLik1 - x$logLik0),
                                                  attributes(x$logLik1)[["df"]] - attributes(x$logLik0)[["df"]], 
                                                  lower.tail =  FALSE), 4), "\n")
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
  cat("-------------------------------------------------------------------\n")
}

############################
# Effects and other functions
#############################

#' @rdname hetprob
#' @method predict hetprob
#' @export
predict.hetprob <- function(object, newdata = NULL, 
                            type = c("xb", "pr", "sigma"), 
                            ...){
  # xb: linear prediction xb
  # pr: probability of a positive outcome
  # sigma: sigma_i = exp(zi gamma)
  type <- match.arg(type)
  mf   <- if (is.null(newdata)) object$mf else newdata
  X    <- model.matrix(object$formula, data = mf, rhs = 1)
  Z    <- model.matrix(object$formula, data = mf, rhs = 2)
  ### Drop intercept on Z if included
  zint <- match("(Intercept)", dimnames(Z)[[2]], nomatch = 0)
  if (zint > 0) Z <- Z[, - zint, drop = FALSE]
  K    <- ncol(X)
  call <- object$call
  link <- call[[match("link", names(call))]]
  pfun <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  theta.hat <- coef(object)
  beta.hat  <- theta.hat[1:K]
  delta.hat <- theta.hat[-c(1:K)]
  xb        <- crossprod(t(X), beta.hat)
  sigma     <- exp(crossprod(t(Z), delta.hat))
  pr        <- pfun(xb / sigma)
  if (type == "pr")    out <- as.vector(pr)
  if (type == "xb")    out <- as.vector(xb)
  if (type == "sigma") out <- as.vector(sigma)
  return(out)
}

#' Get average marginal effects for heterokedastic binary models
#' 
#' Obtain the average marginal effects from \code{hetprob} class model.
#' @param object an object of class \code{hetprob} and \code{effect.hetprob} for \code{summary} and \code{print} method. 
#' @param vcov an estimate of the asymptotic variance-covariance matrix of the parameters for a \code{hetprob} object.
#' @param digits the number of digits.
#' @param ... further arguments.Ignored.
#' @param x an object of class \code{effect.hetprob}.
#' @return An object of class \code{effect.heprob}. 
#' @details 
#' This function allows to obtain the average marginal effects (not the marginal effects at the mean). The standard errors are computed using Delta Method. 
#' @examples
#' \dontrun{ 
#' # Average marginal effects
#' data("Health")
#' het.probit <- hetprob(working ~ factor(female) + factor(year) + educ + age + I(age^2) | 
#'                                 factor(female) + age + I(age^2), 
#'                      data = Health, 
#'                      link = "probit")
#' summary(het.probit)
#' eff <- effect.hetprob(het.probit)
#' summary(eff)
#' }
#' @import stats
#' @importFrom numDeriv jacobian
#' @export 
effect.hetprob <- function(object,
                           vcov = NULL, 
                           digits = max(3, getOption("digits") - 2),
                           ...){
  if (!inherits(object, "hetprob")) stop("not a \"hetprob\" object")
  # Variance covariance matrix
  if (is.null(vcov)){
    V <- vcov(object)
  } else {
    V <- vcov
    n.param <- length(coef(object))
    if (dim(V)[1L] != n.param | dim(V)[2L] != n.param)  stop("dim of vcov are not the same as the estimated parameters")
  } 
  
  # Make effects
  me <- mdydx.hetprob(coeff = coef(object), object) 
  
  # Make Jacobian (use numerical jacobian)
  jac <- numDeriv::jacobian(mdydx.hetprob, coef(object), object = object)
  
  # Save results
  se <- sqrt(diag(jac %*% V %*% t(jac))) 
  z  <-  me / se 
  p  <- 2 * pnorm(-abs(z))
  results            <- cbind(`dydx` = me, `Std. error` = se, `z value` = z, `Pr(> z)` = p)
  object$margins     <- results
  class(object)      <- c("effect.hetprob")
  return(object)
}

#' @rdname effect.hetprob
#' @method summary effect.hetprob
#' @import stats
#' @export
summary.effect.hetprob <- function(object, ...){
  CoefTable      <- object$margins
  summary        <- list(CoefTable = CoefTable)
  class(summary) <- "summary.effect.hetprob"
  summary
}

#' @rdname effect.hetprob
#' @method print summary.effect.hetprob
#' @import stats
#' @export
print.summary.effect.hetprob <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("------------------------------------------------------", fill = TRUE)
  cat("Marginal effects for the heteroskedastic binary model:\n")
  cat("------------------------------------------------------",fill = TRUE)
  printCoefmat(x$CoefTable, digits = digits)
  cat("\nNote: Marginal effects computed as the average for each individual", fill = TRUE)
}

mdydx.hetprob <- function(coeff, object){
  # Three scenarios:
  # 1: the variable appears in x
  # 2: the variable appears in z
  # 3: the variable appears in both
  X         <- model.matrix(object)$X
  Z         <- model.matrix(object)$Z
  K         <- ncol(X)
  beta      <- colnames(X)
  delta     <- colnames(Z)
  theta.hat <- coeff
  beta.hat  <- theta.hat[1:K]
  delta.hat <- theta.hat[-c(1:K)]
  call      <- object$call
  link      <- call[[match("link", names(call))]]
  dfun      <- switch(link,
                 "probit" = dnorm,
                 "logit"  = dlogis)
  pfun      <- switch(link,
                 "probit" = pnorm,
                 "logit"  = plogis)
  
  ## Makes classes of parameters
  all.vars <- all.vars(object$formula)[-1L]
  classes  <- rep("numeric", length(all.vars))
  class.mf <- attributes(terms(object$mf))[["dataClasses"]][-1L]
  classes[paste0("factor(", all.vars, ")") %in% names(class.mf)]  <- class.mf[names(class.mf) %in% paste0("factor(", all.vars, ")")]
  names(delta.hat) <- colnames(Z)
  
  ## Compute marginal effects
  mes <- c()
  mes.name <- c()
  for (k in 1:length(all.vars)){
    if (classes[k] == "numeric"){
      xb     <- crossprod(t(X), beta.hat)
      sigma  <- exp(crossprod(t(Z), delta.hat))
      dens   <- dfun( xb / sigma) 
      # check if continuous appears in interaction
      betak  <- make.inter.num(the.var = all.vars[k], beta, beta.hat, X)
      deltak <- make.inter.num(the.var = all.vars[k], delta, delta.hat, Z)
      res    <- dens * (betak- xb * deltak) / sigma  
      mes    <- cbind(mes, res)
      mes.name <- c(mes.name, all.vars[k])
    }
    if (classes[k] == "factor"){
      levs <- attributes(object$mf[, paste0("factor(", all.vars[k], ")")])$levels
      levs <- levs[-1L]
      ## Make P0
      beta.temp  <- beta.hat
      delta.temp <- delta.hat
      vb <- make.inter.factor(all.vars[k], beta, levs)
      vd <- make.inter.factor(all.vars[k], delta, levs)
      if (any(vb$names %in% beta))    beta.temp[beta   %in% vb$names] <- 0 
      if (any(vd$names %in% delta))  delta.temp[delta  %in% vd$names] <- 0 
      p0 <- pfun(crossprod(t(X), beta.temp) / exp(crossprod(t(Z), delta.temp)))
      for (j in 1:length(levs)){
        ## Make P1
        Xtemp  <- X
        Ztemp  <- Z
        if (any(vb$names %in% beta))   Xtemp[, beta   %in% vb$names] <- 0 
        if (any(vd$names %in% delta))  Ztemp[, delta  %in% vd$names] <- 0 
        vbj <- make.inter.factor(all.vars[k], beta, levs[j])
        vdj <- make.inter.factor(all.vars[k], delta, levs[j])
        if (any(vbj$names %in% beta))  Xtemp[, beta   %in% vbj$names] <- X[, vbj$names.inte] 
        if (any(vdj$names %in% delta)) Ztemp[, delta  %in% vdj$names] <- Z[, vdj$names.inte]
        
        if (vbj$names[1] %in% beta)   Xtemp[, beta    %in% vbj$names[1]] <- 1
        if (vdj$names[1] %in% delta)  Ztemp[, delta   %in% vdj$names[1]] <- 1
        p1     <- pfun(crossprod(t(Xtemp), beta.hat) / exp(crossprod(t(Ztemp), delta.hat)))
        res    <- p1 - p0
        mes    <- cbind(mes, res)
        mes.name <- c(mes.name, paste0("factor(",all.vars[k],")",levs[j], sep = ""))
      }
    }
  }
  colnames(mes) <- mes.name 
  mes <- colMeans(mes)
  return(mes)
}

#' Get Model Summaries for use with "mtable" for objects of class hetprob
#' 
#' A generic function to collect coefficients and summary statistics from a \code{hetprob} object. It is used in \code{mtable}
#' 
#' @param obj a \code{hetprob} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @import stats
#' @importFrom memisc getSummary
#' @export 
getSummary.hetprob <- function(obj, alpha = 0.05, ...){
  s           <- summary(obj)
  cf.mean     <- s$results.mean
  cf.lnsigma  <- s$results.lnsigma
  cval        <- qnorm(1 - alpha/2)
  cf.mean     <- cbind(cf.mean,    cf.mean[, 1] - cval * cf.mean[, 2], cf.mean[, 1] + cval * cf.mean[, 2])
  cf.lnsigma  <- cbind(cf.lnsigma, cf.lnsigma[, 1] - cval * cf.lnsigma[, 2], cf.lnsigma[, 1] + cval * cf.lnsigma[, 2])
  rownames(cf.lnsigma) <- colnames(model.matrix(obj)$Z)
  all.vars    <- unique(c(colnames(model.matrix(obj)$X), colnames(model.matrix(obj)$Z)))
  # Make Table
  coef        <- array(dim = c(length(all.vars), 6, 2), 
                        dimnames = list(all.vars, c("est", "se", "stat", "p", "lwr", "upr"), c("mean", "lnsigma")))
  coef[rownames(cf.mean),,1]    <- cf.mean
  coef[rownames(cf.lnsigma),,2] <- cf.lnsigma
  
  # Statistics
  sumstat <- c(logLik = logLik(obj), deviance = NA, AIC = AIC(obj), BIC = BIC(obj), N = nrow(obj$gradientObs), 
               LR = NA, df = NA, p = NA, Aldrich.Nelson = NA, McFadden = NA, Cox.Snell = NA)
  list(coef = coef, sumstat = sumstat, contrasts = obj$contrasts, xlevels = obj$xlevels, call = obj$call)
}


#' Get Model Summaries for use with "mtable" for objects of class effect.hetprob
#' 
#' A generic function to collect coefficients and summary statistics from a \code{effect.hetprob} object. It is used in \code{mtable}
#' 
#' @param obj an \code{effect.hetprob} object,
#' @param alpha level of the confidence intervals,
#' @param ... further arguments,
#' 
#' @details For more details see package \pkg{memisc}.
#' @import stats
#' @importFrom memisc getSummary
#' @export
getSummary.effect.hetprob <- function(obj, alpha = 0.05, ...){
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

