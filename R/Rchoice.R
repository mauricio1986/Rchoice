#' Estimate discrete choice model with random parameters
#' 
#' Estimation of discrete choice models such as binary (logit and probit), 
#' poisson and ordered (logit and probit) model with random coefficients for cross-section data by simulated maximum likelihood
#' 
#' @param x,object,obj and object of class \code{Rchoice},
#' @param formula a symbolic description of the model to be estimated,
#' @param new an updated formula for the update method,
#' @param data the data,
#' @param subset an optional vector specifying a subset of observations,
#' @param weights an optional vector of weigths,
#' @param na.action a function wich indicated what should happen when the data
#' contains \code{NA}'s
#' @param start a vector of starting values,
#' @param link a string indicating the model to be estimated: "\code{probit}", "\code{logit}" , "\code{poisson}", 
#' "\code{ordered probit}" or "\code{ordered logit}"
#' @param ranp a named vector whose names are the random parameters and values the distribution:
#' "\code{n}" for normal, "\code{ln}" for log-normal, "\code{cn}" for truncated normal, "\code{u}" for uniform, "\code{t}" for triangular,
#' @param R the number of draws of pseudo-random numbers if \code{ranp} is not \code{NULL}.
#' @param haltons only relevant if \code{ranp} is not \code{NULL}. If not \code{NULL}, halton sequence is used
#' instead of pseudo-random numbers. If \code{haltons=NA}, some default values are used for
#' the prime of the sequence and for the number of element droped. Otherwise, \code{haltons} should
#' be a list with elements \code{prime} and \code{drop}.
#' @param seed ,
#' @param correlation only relevant if \code{ranp} is not \code{NULL}. If true, the correlation between random
#' parameters is taken into account,
#' @param alpha significance value for \code{getSummary},
#' @param ... further arguments passed to \code{maxLik}.
#' @export
#' @details
#' \itemize{
#'  \item The model is estimated using the \code{maxLik} function of \code{\link[maxLik]{maxLik}} package.
#'  \item If \code{ranp} is not NULL, the random parameter (random coefficient) model is estimated. The model 
#'  specified by \code{link} is estimated. The \code{link} model is estimated using Simulated Maximum Likelihood (SMLE)
#'  where the probabilities are simulated using \code{R} pseudo-draws if \code{halton=NULL} or \code{R} halton
#'  draws if \code{halton=NA}. The user can also specified the primes and the number of dropped elements for the halton draws.
#'  For example, if the model consists of two random parameters, the user can specify \code{haltons = list("prime"=c(2,3), "drop" = c(11,11))}. 
#'  The functions to create the random parameters and the random draws are those crated by Yves Croissant for the package \code{\link[mlogit]{mlogit}}.
#'  \item A random parameter hierarchical model can be estimated by including heterogeneity in the mean of the 
#'  random parameters. \pkg{Rchoice} manages the variables in the hierarchical model by 
#'  the \code{formula} object: all the hierarchical variables are included after the \code{|} symbol. See examples below
#' }
#' @return An object of class ``\code{Rchoice}'', a list elements:
#' \item{coefficients}{the named vector of coefficients,}
#' \item{link}{the type of model fitted,}
#' \item{logLik}{a set of values of the maximum likelihood procedure,}   
#' \item{mf}{the model framed used,} 
#' \item{formula}{the formula (a Formula object),}
#' \item{time}{\code{proc.time()} minus the start time,}
#' \item{freq}{frequency of dependent variable,}
#' \item{draws}{type of draws used,}
#' \item{R.model}{\code{TRUE} if a random parameter model is fitted,}
#' \item{R}{number of draws used,}
#' \item{b.random}{matrix of conditional expectation of random parameters,}
#' \item{sd.random}{matrix of standard deviation of conditional expectation of random parameters,}
#' \item{ranp}{vector indicating the variables with random parameters and their distribution,}
#' \item{probabilities}{the fitted probabilities for each individuals,}
#' \item{residuals}{the residuals,}
#' \item{call}{the matched call.}   
#' @author Mauricio Sarrias
#' @seealso \code{\link[mlogit]{mlogit}}, \code{\link[maxLik]{maxLik}}
#' @examples
#' ## Probit model
#' data("Workmroz")
#' probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  data = Workmroz , link="probit")
#' summary(probit)
#' 
#' ## Poisson model
#' data("Articles")
#' poisson <- Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, link = "poisson")
#' summary(poisson)
#' 
#' ## Ordered probit model
#' data("Health")
#' oprobit<-Rchoice(newhsat ~ age + educ + hhinc + married + hhkids, 
#' data = Health, link = "ordered probit", subset = year == 1988)
#' summary(oprobit)
#' 
#'  \dontrun{
#' ## Hierarchical Logit Random Parameter Model 
#' Hran.logit<-Rchoice(lfp ~ k618 + lwg + wc + inc + k5 | age + wc + hc, 
#' ranp = c(inc = "t", k5 = "n"), 
#' link = "logit", data = Workmroz)
#' summary(Hran.logit)
#' }
#' 
#' \dontrun{
#' ## Hierarchical Poisson model with correlated random parameters
#' poissonH.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment | fem, data = Articles,
#' ranp = c(kid5="n", phd = "n", ment = "n"), link = "poisson", correlation =  TRUE)
#' summary(poissonH.ran)
#' }
#' 
#' \dontrun{
#' ## Ordered Probit model with random parameters
#' oprobit.ran<-Rchoice(newhsat ~ age + educ + hhinc + married + hhkids, 
#'                     data = Health, link = "ordered probit", subset = year == 1988, 
#'                     ranp = c(age = "n", hhinc = "n"), print.level=1, 
#'                     start = rep(0,11))
#' summary(oprobit.ran)
#' }
#' @references
#' Greene, W. H. (2012). Econometric analysis. 7 edition. Prentice Hall.
#' 
#' Train, K. (2009). Discretechoice methods with simulation. Cambridge university press.

Rchoice<-function(formula, data , subset , weights , na.action ,
                  start = NULL , link = c("probit","logit","poisson", "ordered probit", "ordered logit") ,
                  ranp = NULL , R = 40 , haltons = NA , seed = 123 , correlation = FALSE,
                  ...)
{
  start.time <- proc.time()
  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE)
  
  #I use Formula package
  formula <- callF$formula <- Formula::as.Formula(formula)
  nframe  <- length(sys.calls())
  
  #Check Model
  link <- match.arg(link)
  R.model <- !is.null(ranp)
  #Hierarchical model?
  Hier <- ifelse(length(formula)[2] == 2, TRUE, FALSE)
  
  #Checking specification of hierarchical model
  if (Hier) if (!R.model) stop('Hierarchical model needs ranp to be specified')
  
  if (R.model){
    if (is.null(callT$method))      callT$method      <- 'bfgs'
    if (is.null(callT$print.level)) callT$print.level <- 0
    if (is.null(callT$iterlim))     callT$iterlim     <- 2000
    if (is.null(callT$tol))         callT$tol         <- 1E-06
    if (is.null(callT$steptol))     callT$steptol     <- 1E-10
    if (is.null(callT$ftol))        callT$ftol        <- 1E-08
  }else{
    if (is.null(callT$method) && (link == "ordered probit" || link == "ordered logit")) callT$method <- 'bfgs'
    if (is.null(callT$method)) callT$method <- 'nr'
  }
  
  ####################
  # Model Frame
  ####################
  mf<-callT
  m <- match(c("formula", "data", "subset", "na.action", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ##########################################
  # Extract the elements of the model
  ##########################################
  
  if (Hier){
    cons <- attr(terms(formula), "intercept")
    if (cons == 1){
      warning("Model assumes no constant in S variables...updating formula", call. = FALSE, immediate. = TRUE)
      formula <- update(formula, ~. | ~. -1)
    }
    S <- model.matrix(formula, mf , rhs = 2) 
  }
  y <- model.response(mf)
  X <- model.matrix(formula, mf, rhs = 1)
  weights <- model.weights(mf)
  freq <- table(y)
  if (link == "ordered probit" || link == "ordered logit") {
    y <- as.factor(y)
    J <- length(levels(y))
    if (J < 3) stop("The alternatives must be >=3 in y")
  }
  
  ########################################################
  #   Initial Values
  ########################################################
  
  # Names of thresholds if ordered model
  names.kappa <- c()
  if (link == "ordered probit" || link == "ordered logit")  names.kappa <- paste('kappa', 1:(J - 2) , sep = ':')

  # Names for models with random parameters
  names.random <- c()
  if (R.model){
    #Check distribution of random parameters
    ndist <- ranp[! (ranp %in% c("cn", "ln", "n", "u", "t"))]
    if (length(ndist) > 0){
      udstr <- paste("unknown distribution", paste(unique(ndist), collapse = ", "))
      stop(udstr)
    }
    Vara <- sort(match(names(ranp), colnames(X))) 
    Varc <- (1:ncol(X))[- Vara]
    Xa   <- X[ , Vara, drop = F]                        
    Xc   <- X[ , Varc, drop = F]  
    colnamesX <- colnames(X)
    names.f <- colnamesX[Varc] #Names for fixed parameters
    names.b <- paste('mean', colnamesX[Vara], sep = ':')
    if (!correlation) {
      names.sd <- paste('sd', colnamesX[Vara], sep = ':')
    } else {
      names.sd <- c()
      Ka <- length(ranp)
      for (i in 1:Ka){
        names.sd <- c(names.sd,
                      paste('sd', names(ranp)[i], names(ranp)[i:Ka], sep = ':')
        )
      }
    }
    names.phi <- c()
    if (Hier) names.phi <- c(outer(names(ranp), colnames(S) , FUN = paste, sep = ":"))
    names.random <- c(names.b, names.phi, names.sd)
  }  else {
    names.f <- colnames(X) 
  }

all.names   <- c(names.kappa, names.f, names.random)
# As if not random parameters and start is null
if (is.null(start)){
  if (link == "ordered probit" || link == "ordered logit"){
    if (R.model)  theta <- coef(lm(unclass(y)~Xc-1+Xa)) else theta <- coef(lm(unclass(y)~X-1))
     z <- as.integer(table(unclass(y)))
     z <- (cumsum(z)/sum(z))[1:(J-1)]
     start.kappa <- switch(link,
                   "ordered probit" = qnorm(z),
                   "ordered logit"  = qlogis(z))
    theta       <- c(log(diff(start.kappa)), theta)
  } else {
    if (R.model)  theta <- coef(lm(y~Xc-1+Xa)) else theta <- coef(lm(y~X-1))
  }
  names(theta) <- c(names.kappa, names.f)
}

  
# Initial value if random parameters and start is null
 if (R.model & is.null(start)) {
   callst        <- callT
   callst$start  <- theta
   callst$method <- 'bfgs'
   callst$iterlim <-  callst$tol <- callst$steptol <- callst$ftol <- NULL
   callst$haltons <- callst$ranp <- callst$R <- callst$data <- callst$formula <- NULL 
   callst$print.level <- 0
   callst[[1]]   <- as.name('maxLik')
   Xst <- cbind(Xc, Xa) 
   callst$X <- Xst ; callst$y <- y
   callst$link  <- link
   if (link == "poisson")                                    callst$logLik <- as.name('lnpoisson')
   if (link == "probit" || link == "logit")                  callst$logLik <- as.name('lnbinary')
   if (link == "ordered probit" || link == "ordered logit")  callst$logLik <- as.name('lnordered')
   start.fixed  <- coef(eval(callst, sys.frame(which=nframe)))
   if (is.null(start.fixed)) stop("attempt to find suitable starting values failed")
   start.random <- rep(0, length(c(names.phi, names.sd)))
   theta        <- c(start.fixed, start.random)
   names(theta) <- all.names
 }
 
# Initial value start is not null
if (!is.null(start)){
  theta <- start
  if (length(start) != length(all.names)) stop('Incorrect Number of Initial Parameters')
  names(theta) <- all.names
}
  
  cat("\nStarting values Parameters:\n")
  print(theta)
   
 #######################################################################
 # Estimate the model using maxLik and passing the correct arguments
 #######################################################################
 opt <- callT
 opt$start <- theta
 
 # Maximization control arguments
 m <- match(c('method', 'print.level', 'iterlim',
              'start','tol', 'ftol', 'steptol'),
               names(opt), 0L)
 opt <- opt[c(1L, m)]
 
 #Optimization code name
 opt[[1]] <- as.name('maxLik')
 
 #Variables
 opt[c('X', 'y')] <- list(as.name('X'), as.name('y'))
 
 #Weights
 if (is.null(weights)) weights <- rep(1, nrow(X))
 opt$weights <- weights
 
 #Link
 opt$link  <- link
 
 #Standard Models
 if (link == "poisson") opt$logLik <- as.name('lnpoisson')
 if (link == "probit" || link == "logit")  opt$logLik <- as.name('lnbinary')
 if (link == "ordered probit" || link == "ordered logit")  opt$logLik <- as.name('lnordered')
 
 #Arguments for random parameters
 if (R.model) {
   opt[c('R', 'seed', 'ranp', 'correlation','haltons')] <-
     list(as.name('R'), as.name('seed'), as.name('ranp'), as.name('correlation'), as.name('haltons'))
   if (Hier) {
     if (link == "ordered probit" || link == "ordered logit") opt$logLik <- as.name('lnorderedH.ran') else opt$logLik <- as.name('lnlH.ran') 
     opt$S <- as.name('S')
   } else {
     if (link == "ordered probit" || link == "ordered logit") opt$logLik <- as.name('lnordered.ran') else opt$logLik <- as.name('lnl.ran')
   }
 }
  
  #Optimizing the ML
  x <- eval(opt, sys.frame(which = nframe))

 ###########################
 # Put results in form
 ###########################
 
 #Get probability, and conditional beta
 if (!is.null(ranp)){
   opt$steptol <- opt$logLik <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol<-opt$ftol <- NULL
   names(opt)[[2]] <- 'theta'
   betahat <- coef(x)
   if (Hier) {
     if (link == "ordered probit" || link == "ordered logit") opt[[1]] <- as.name('lnorderedH.ran') else opt[[1]] <- as.name('lnlH.ran') 
   }else{
     if (link == "ordered probit" || link == "ordered logit") opt[[1]] <- as.name('lnordered.ran') else opt[[1]] <- as.name('lnl.ran')
   } 
   betahat <- ifelse(names(betahat) %in% names.sd, abs(betahat), betahat)
   names(betahat) <- names(coef(x))
   opt[[2]] <-betahat
   opt$make.estb <- TRUE
   again <- eval(opt, sys.frame(which = nframe))
   x$probabilities <- attr(again, 'probabilities')
   x$b.random      <- attr(again, 'b.random')
   x$sd.random     <- attr(again, 'sd.random')
   x$estimate      <- betahat
 } else {
   opt$steptol <- opt$logLik <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol<-opt$ftol <- NULL
   names(opt)[[2]] <- 'theta'
   betahat  <- coef(x)
   opt[[2]] <- betahat
   if (link == "poisson")                                    opt[[1]] <- as.name('lnpoisson')
   if (link == "probit" || link == "logit")                  opt[[1]] <- as.name('lnbinary')
   if (link == "ordered probit" || link == "ordered logit")  opt[[1]] <- as.name('lnordered')
   again <- eval(opt, sys.frame(which = nframe))
   x$probabilities <- drop(attr(again, 'probabilities'))
 }
 
 #Ordered Model
 if (link == "ordered probit" || link == "ordered logit"){
   J <- length(levels(y))
   kappas <- cumsum(c(exp(x$estimate[1:(J - 2)])))
   names(kappas) <- names.kappa
   attr(x$estimate, "kappas") <- kappas
   attr(x$estimate, "fixed" ) <- x$estimate[-c(1:(J - 2))]
 }
 
 resid <- drop(unclass(y) - x$probabilities)
 
 
 
logLik<-structure(list(
                    maximum     = logLik(x),
                    gradient    = x$gradient,
                    nobs        = nObs(x),
                    gradientObs = x$gradientObs,
                    hessian     = hessian(x),
                    iterations  = nIter(x),
                    type        = maximType(x),
                    code        = returnCode(x),
                    nparam      = nParam(x),
                    message     = returnMessage(x)),
                    class = "logLik"
                   )
 
 
  result<-structure(list(
                    coefficients  = x$estimate,
                    link          = link,
                    logLik        = logLik,
                    mf            = mf,
                    formula       = formula,
                    time          = proc.time()-start.time,
                    freq          = freq,
                    draws         = haltons,
                    R.model       = R.model,
                    R             = R,
                    b.random      = x$b.random,
                    sd.random     = x$sd.random,
                    ranp          = ranp,
                    probabilities = x$probabilities,
                    residuals     = resid,
                    correlation   = correlation,
                    call          = callT),
                    class="Rchoice"
                 )
 result
}