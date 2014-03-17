\name{estfun.Rchoice}
\alias{estfun.Rchoice}
\title{Gradient for observations}
\usage{
\method{estfun}{Rchoice}(x, ...)
}
\arguments{
  \item{x}{a fitted model of class \code{Rchoice}}

  \item{...}{Other arguments when \code{estfun} is applied
  to another class object}
}
\value{
the gradient matrix of dimension n times k
}
\description{
It extracts the gradient for each observations evaluated at
the estimated parameters
}
\examples{
## Probit model
data("Workmroz")
probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,  data = Workmroz , link="probit")
summary(probit)

estfun(probit)
}
\references{
Zeileis A (2006), Object-oriented Computation of Sandwich
Estimators. Journal of Statistical Software, 16(9), 1--16.
}

