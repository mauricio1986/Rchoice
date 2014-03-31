\name{Rchoice}
\alias{Rchoice}
\alias{coef.Rchoice}
\alias{df.residual.Rchoice}
\alias{fitted.Rchoice}
\alias{getSummary.Rchoice}
\alias{logLik.Rchoice}
\alias{model.matrix.Rchoice}
\alias{nObs.Rchoice}
\alias{ordinal}
\alias{print.Rchoice}
\alias{summary.Rchoice}
\alias{terms.Rchoice}
\alias{update.Rchoice}
\alias{vcov.Rchoice}
\title{Estimate discrete choice model with random parameters}
\usage{
Rchoice(formula, data, subset, weights, na.action, family, start = NULL,
  ranp = NULL, R = 40, haltons = NA, seed = 10, correlation = FALSE,
  ...)

\method{terms}{Rchoice}(x, ...)

\method{model.matrix}{Rchoice}(object, ...)

\method{vcov}{Rchoice}(object, ...)

\method{coef}{Rchoice}(object, ...)

\method{nObs}{Rchoice}(x, ...)

\method{fitted}{Rchoice}(object, ...)

\method{df.residual}{Rchoice}(object, ...)

\method{update}{Rchoice}(object, new, ...)

\method{logLik}{Rchoice}(object, ...)

\method{print}{Rchoice}(x, digits = max(3, getOption("digits") - 3),
  width = getOption("width"), ...)

\method{summary}{Rchoice}(object, ...)

getSummary.Rchoice(obj, alpha = 0.05, ...)
}
\arguments{
  \item{x,object,obj}{and object of class \code{Rchoice},}

  \item{formula}{a symbolic description of the model to be
  estimated,}

  \item{new}{an updated formula for the update method,}

  \item{data}{the data,}

  \item{subset}{an optional vector specifying a subset of
  observations,}

  \item{weights}{an optional vector of weigths,}

  \item{na.action}{a function wich indicated what should
  happen when the data contains \code{NA}'s}

  \item{start}{a vector of starting values,}

  \item{family}{the distribution to be used,}

  \item{ranp}{a named vector whose names are the random
  parameters and values the distribution: "\code{n}" for
  normal, "\code{ln}" for log-normal, "\code{cn}" for
  truncated normal, "\code{u}" for uniform, "\code{t}" for
  triangular,}

  \item{R}{the number of draws of pseudo-random numbers if
  \code{ranp} is not \code{NULL}.}

  \item{haltons}{only relevant if \code{ranp} is not
  \code{NULL}. If not \code{NULL}, halton sequence is used
  instead of pseudo-random numbers. If \code{haltons=NA},
  some default values are used for the prime of the
  sequence and for the number of element droped. Otherwise,
  \code{haltons} should be a list with elements
  \code{prime} and \code{drop}.}

  \item{seed}{,}

  \item{correlation}{only relevant if \code{ranp} is not
  \code{NULL}. If true, the correlation between random
  parameters is taken into account,}

  \item{alpha}{significance value for \code{getSummary},}

  \item{digits}{number of digits,}

  \item{width}{width,}

  \item{...}{further arguments passed to \code{maxLik}.}
}
\value{
An object of class ``\code{Rchoice}'', a list elements:
\item{coefficients}{the named vector of coefficients,}
\item{family}{type of model,} \item{link}{distribution of
the errors,} \item{logLik}{a set of values of the maximum
likelihood procedure,} \item{mf}{the model framed used,}
\item{formula}{the formula (a Formula object),}
\item{time}{\code{proc.time()} minus the start time,}
\item{freq}{frequency of dependent variable,}
\item{draws}{type of draws used,}
\item{R.model}{\code{TRUE} if a random parameter model is
fitted,} \item{R}{number of draws used,}
\item{b.random}{matrix of conditional expectation of random
parameters,} \item{sd.random}{matrix of standard deviation
of conditional expectation of random parameters,}
\item{ranp}{vector indicating the variables with random
parameters and their distribution,}
\item{probabilities}{the fitted probabilities for each
individuals,} \item{residuals}{the residuals,}
\item{call}{the matched call.}
}
\description{
Estimation of discrete choice models such as Binary (logit
and probit), Poisson and Ordered (logit and probit) model
with random coefficients for cross-section data by
simulated maximum likelihood
}
\details{
The models are estimated using the \code{maxLik} function
of \code{\link[maxLik]{maxLik}} package.

If \code{ranp} is not NULL, the random parameter (random
coefficient) model is estimated.  A random parameter model
or random coefficient models permits regression parameter
to vary across individuals according to some distribution.
A fully parametric random parameter model specifies the
latent variable \eqn{y^{*}} conditional on regressors
\eqn{x} and given parameters \eqn{\beta_i} to have
conditional density \eqn{f(y|x, \beta_i)} where
\eqn{\beta_i} are iid with density
\eqn{g(\beta_i|\theta_i)}. The density is assumed a priori
by the user by the argument \code{ranp}. If the parameters
are assumed to be normally distributed \eqn{\beta_i ~
N(\beta, \Sigma)}, then the random parameter are
constructed as: \deqn{\beta_{ir}=\beta+L\omega_{ir}} where
\eqn{LL'=\Sigma} and \eqn{\omega_{ir}} is the {r}-th draw
from standard normal distribution for individual \eqn{i}.

Once the model is specified by the argument \code{family},
the model is estimated using Simulated Maximum Likelihood
(SMLE). The probabilities, given by \eqn{f(y|x, \beta_i)},
are simulated using \code{R} pseudo-draws if
\code{halton=NULL} or \code{R} halton draws if \code{halton
= NA}. The user can also specified the primes and the
number of dropped elements for the halton draws.  For
example, if the model consists of two random parameters,
the user can specify \code{haltons = list("prime" = c(2,
3), "drop" = c(11, 11))}.

A random parameter hierarchical model can be estimated by
including heterogeneity in the mean of the random
parameters: \deqn{\beta_{ir}=\beta+\pi's_i+L\omega_{ir}}
\pkg{Rchoice} manages the variables in the hierarchical
model by the \code{formula} object: all the hierarchical
variables (\eqn{s_i}) are included after the \code{|}
symbol. See examples below
}
\examples{
## Probit model
data("Workmroz")
probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,
                 data = Workmroz, family = binomial('probit'))
summary(probit)

## Poisson model
data("Articles")
poisson <- Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, family = poisson)
summary(poisson)

## Ordered probit model
data("Health")
oprobit <- Rchoice(newhsat ~ age + educ + hhinc + married + hhkids,
data = Health, family = ordinal('probit'), subset = year == 1988)
summary(oprobit)

 \dontrun{
## Hierarchical Logit Random Parameter Model
Hran.logit<-Rchoice(lfp ~ k618 + lwg + wc + inc + k5 | age + wc + hc,
ranp = c(inc = "t", k5 = "n"),
family = binomial('logit'), data = Workmroz)
summary(Hran.logit)
}

\dontrun{
## Hierarchical Poisson model with correlated random parameters
poissonH.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment | fem, data = Articles,
ranp = c(kid5="n", phd = "n", ment = "n"), family = poisson, correlation =  TRUE)
summary(poissonH.ran)
}

\dontrun{
## Ordered Probit model with random parameters
oprobit.ran <- Rchoice(newhsat ~ age + educ + hhinc + married + hhkids,
                    data = Health, family = ordinal('probit'), subset = year == 1988,
                    ranp = c(age = "n", hhinc = "n"), print.level=1,
                    start = rep(0,11))
summary(oprobit.ran)
}
}
\author{
Mauricio Sarrias \email{msarrias86@gmail.com}
}
\references{
Greene, W. H. (2012). Econometric analysis. 7 edition.
Prentice Hall.

Train, K. (2009). Discretechoice methods with simulation.
Cambridge university press.
}
\seealso{
\code{\link[mlogit]{mlogit}}, \code{\link[maxLik]{maxLik}}
}
\keyword{package}

