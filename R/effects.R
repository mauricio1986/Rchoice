### Effect methods

#' Get average marginal effects for heterokedastic binary models and IV probit models
#' 
#' Obtain the average marginal effects from \code{hetprob} or \code{ivpml} class models.
#' @param object an object of class \code{hetprob} or \code{ivpml}.
#' @param ... Additional arguments to be passed.
#' 
#' @return Estimates of the average marginal effects computed as the average for each individual.
#' @examples
#' \donttest{
#' # Data
#' library("AER")
#' data("PSID1976")
#' PSID1976$lfp  <- as.numeric(PSID1976$participation == "yes")
#' PSID1976$kids <- with(PSID1976, factor((youngkids + oldkids) > 0,
#'                                       levels = c(FALSE, TRUE), 
#'                                       labels = c("no", "yes")))
#' PSID1976$finc <-  PSID1976$fincome / 10000
#'                                  
#' # Average marginal effects for heteroskedastic Probit model
#' labor_het <- hetprob(lfp ~  age + I(age^2) + finc + education + factor(kids) | 
#'                             factor(kids) + finc,              
#'                      data = PSID1976,                        
#'                      link = "probit")
#' eff_labor_het <- effect(labor_het)
#' summary(eff_labor_het)
#' 
#' # Average marginal effects for IV probit model 
#' # (nwincome is endogenous and heducation is the additional instrument)
#' PSID1976$nwincome <- with(PSID1976, (fincome - hours * wage)/1000)
#' fiml.probit <- ivpml(lfp ~  education + experience + I(experience^2) + age + 
#'                             youngkids + oldkids + nwincome |
#'                             education + experience + I(experience^2) + age + 
#'                             youngkids + oldkids + heducation, 
#'                      data = PSID1976)
#' summary(effect(fiml.probit))
#' summary(effect(fiml.probit, asf = FALSE))
#' }
#' @export
effect <- function(object, ...){
  UseMethod("effect", object)
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
#' @author Mauricio Sarrias. 
#' @examples
#' \donttest{
#' # Data
#' library("AER")
#' data("PSID1976")
#' PSID1976$lfp  <- as.numeric(PSID1976$participation == "yes")
#' PSID1976$kids <- with(PSID1976, factor((youngkids + oldkids) > 0,
#'                                       levels = c(FALSE, TRUE), 
#'                                       labels = c("no", "yes")))
#'                                       PSID1976$finc <-  PSID1976$fincome / 10000
#'                                  
#' # Average marginal effects for heteroskedastic Probit model
#' labor_het <- hetprob(lfp ~  age + I(age^2) + finc + education + factor(kids) | 
#'                             factor(kids) + finc,              
#'                      data = PSID1976,                        
#'                      link = "probit")
#' eff_labor_het <- effect(labor_het)
#' summary(eff_labor_het)
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
#' @method print effect.hetprob
#' @import stats
#' @export 
print.effect.hetprob <- function(x, ...){
  cat("The marginal effects are:\n")
  cat("Estimate(s):", x$margins[, 1], "\n")
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


#' Get average marginal effects for IV Probit model.
#' 
#' Obtain the average marginal effects from \code{ivpml} class model.
#' @param object an object of class \code{ivpml} and \code{effect.ivpml} for \code{summary} and \code{print} method. 
#' @param vcov an estimate of the asymptotic variance-covariance matrix of the parameters for a \code{ivpml} object.
#' @param asf  if \code{TRUE}, the average structural function is used. 
#' @param digits the number of digits.
#' @param ... further arguments.Ignored.
#' @param x an object of class \code{effect.ivpml}.
#' @return An object of class \code{effect.ivpml}. 
#' @details 
#' This function allows to obtain the average marginal effects (not the marginal effects at the mean). The standard errors are computed using Delta Method.
#' @author Mauricio Sarrias.
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
#' # Average marginal effects for IV probit model 
#' # (nwincome is endogenous and heducation is the additional instrument)
#' PSID1976$nwincome <- with(PSID1976, (fincome - hours * wage)/1000)
#' fiml.probit <- ivpml(lfp ~  education + experience + I(experience^2) + age + 
#'                             youngkids + oldkids + nwincome |
#'                             education + experience + I(experience^2) + age + 
#'                             youngkids + oldkids + heducation, 
#'                      data = PSID1976)
#' summary(effect(fiml.probit))
#' summary(effect(fiml.probit, asf = FALSE))
#' } 
#' @import stats
#' @importFrom numDeriv jacobian
#' @export 
effect.ivpml <- function(object,
                         vcov = NULL, 
                         asf = TRUE,
                         digits = max(3, getOption("digits") - 2), 
                         ...){
  if (!inherits(object, "ivpml")) stop("not a \"ivpml\" object")
  # Variance covariance matrix
  if (is.null(vcov)){
    V <- vcov(object)
  } else {
    V <- vcov
    n.param <- length(coef(object))
    if (dim(V)[1L] != n.param | dim(V)[2L] != n.param)  stop("dim of vcov are not the same as the estimated parameters")
  } 
  
  # Make effects
  me <- mdydx.ivpml(coeff = coef(object), object, asf) 
  
  # Make Jacobian (use numerical jacobian from numDeriv package)
  jac <- numDeriv::jacobian(mdydx.ivpml, coef(object), object = object, asf = asf)
  
  # Print results
  se <- sqrt(diag(jac %*% V %*% t(jac))) 
  z  <-  me / se 
  p  <- 2 * pnorm(-abs(z))
  results            <- cbind(`dydx` = me, `Std. error` = se, `z value` = z, `Pr(> z)` = p)
  object$margins     <- results
  class(object)      <- c("effect.ivpml")
  return(object)
}

#' @rdname effect.ivpml
#' @method summary effect.ivpml
#' @import stats
#' @export
summary.effect.ivpml <- function(object, ...){
  CoefTable      <- object$margins
  summary        <- list(CoefTable = CoefTable)
  class(summary) <- "summary.effect.ivpml"
  summary
}

#' @rdname effect.ivpml
#' @method print effect.ivpml
#' @import stats
#' @export 
print.effect.ivpml <- function(x, ...){
  cat("The marginal effects are:\n")
  cat("Estimate(s):", x$margins[, 1], "\n")
}

#' @rdname effect.ivpml
#' @method print summary.effect.ivpml
#' @import stats
#' @export
print.summary.effect.ivpml <- function(x, digits = max(3, getOption("digits") - 3), ...){
  cat("------------------------------------------------------", fill = TRUE)
  cat("Marginal effects for the IV Probit model:\n")
  cat("------------------------------------------------------",fill = TRUE)
  printCoefmat(x$CoefTable, digits = digits)
  cat("\nNote: Marginal effects computed as the average for each individual", fill = TRUE)
}

