\name{Rchoice}
\alias{Rchoice}
\alias{coef.Rchoice}
\alias{df.residual.Rchoice}
\alias{fitted.Rchoice}
\alias{getSummary.Rchoice}
\alias{logLik.Rchoice}
\alias{model.matrix.Rchoice}
\alias{nObs.Rchoice}
\alias{summary.Rchoice}
\alias{terms.Rchoice}
\alias{update.Rchoice}
\alias{vcov.Rchoice}
\title{Estimate discrete choice model with random parameters}
\usage{
Rchoice(formula, data, subset, weights, na.action, start = NULL,
  link = c("probit", "logit", "poisson", "ordered probit", "ordered logit"),
  ranp = NULL, R = 40, haltons = NA, seed = 123, correlation = FALSE,
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

  \item{link}{a string indicating the model to be
  estimated: "\code{probit}", "\code{logit}" ,
  "\code{poisson}", "\code{ordered probit}" or
  "\code{ordered logit}"}

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

  \item{...}{further arguments passed to \code{maxLik}.}
}
\value{
An object of class ``\code{Rchoice}'', a list elements:
\item{coefficients}{the named vector of coefficients,}
\item{link}{the type of model fitted,} \item{logLik}{a set
of values of the maximum likelihood procedure,}
\item{mf}{the model framed used,} \item{formula}{the
formula (a Formula object),} \item{time}{\code{proc.time()}
minus the start time,} \item{freq}{frequency of dependent
variable,} \item{draws}{type of draws used,}
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
\itemize{ \item The model is estimated using the
\code{maxLik} function of \code{\link[maxLik]{maxLik}}
package.  \item If \code{ranp} is not NULL, the random
parameter (random coefficient) model is estimated. The
model specified by \code{link} is estimated. The
\code{link} model is estimated using Simulated Maximum
Likelihood (SMLE) where the probabilities are simulated
using \code{R} pseudo-draws if \code{halton=NULL} or
\code{R} halton draws if \code{halton=NA}. The user can
also specified the primes and the number of dropped
elements for the halton draws.  For example, if the model
consists of two random parameters, the user can specify
\code{haltons = list("prime"=c(2,3), "drop" = c(11,11))}.
The functions to create the random parameters and the
random draws are those crated by Yves Croissant for the
package \code{\link[mlogit]{mlogit}}.  \item A random
parameter hierarchical model can be estimated by including
heterogeneity in the mean of the random parameters.
\pkg{Rchoice} manages the variables in the hierarchical
model by the \code{formula} object: all the hierarchical
variables are included after the \code{|} symbol. See
examples below }
}
\examples{
## Probit model
data("Workmroz")
probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  data = Workmroz , link="probit")
summary(probit)

## Poisson model
data("Articles")
poisson <- Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, link = "poisson")
summary(poisson)

## Ordered probit model
data("Health")
oprobit<-Rchoice(newhsat ~ age + educ + hhinc + married + hhkids,
data = Health, link = "ordered probit", subset = year == 1988)
summary(oprobit)

 \dontrun{
## Hierarchical Logit Random Parameter Model
Hran.logit<-Rchoice(lfp ~ k618 + lwg + wc + inc + k5 | age + wc + hc,
ranp = c(inc = "t", k5 = "n"),
link = "logit", data = Workmroz)
summary(Hran.logit)
}

\dontrun{
## Hierarchical Poisson model with correlated random parameters
poissonH.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment | fem, data = Articles,
ranp = c(kid5="n", phd = "n", ment = "n"), link = "poisson", correlation =  TRUE)
summary(poissonH.ran)
}

\dontrun{
## Ordered Probit model with random parameters
oprobit.ran<-Rchoice(newhsat ~ age + educ + hhinc + married + hhkids,
                    data = Health, link = "ordered probit", subset = year == 1988,
                    ranp = c(age = "n", hhinc = "n"), print.level=1,
                    start = rep(0,11))
summary(oprobit.ran)
}
}
\author{
Mauricio Sarrias
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

