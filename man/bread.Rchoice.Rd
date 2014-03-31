\name{bread.Rchoice}
\alias{bread.Rchoice}
\title{Bread for sandwiches}
\usage{
\method{bread}{Rchoice}(x, ...)
}
\arguments{
  \item{x}{a fitted model of class \code{Rchoice}}

  \item{...}{Other arguments when \code{bread} is applied
  to another class object}
}
\value{
the covariance matrix times observations
}
\description{
Computes the bread of the sandwich covariance matrix
}
\examples{
## Probit model
data("Workmroz")
probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,
                  data = Workmroz , family = binomial('probit'))
summary(probit)

library(sandwich)
bread(probit)
}
\references{
Zeileis A (2006), Object-oriented Computation of Sandwich
Estimators. Journal of Statistical Software, 16(9), 1--16.
}

