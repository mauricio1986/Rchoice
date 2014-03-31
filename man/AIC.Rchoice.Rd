\name{AIC.Rchoice}
\alias{AIC.Rchoice}
\alias{BIC.Rchoice}
\title{Akaike's Information Criterion}
\usage{
\method{AIC}{Rchoice}(object, ..., k = 2)

\method{BIC}{Rchoice}(object, ...)
}
\arguments{
  \item{object}{a fitted model of class \code{Rchoice}}

  \item{...}{additional arguments to be passed to or from
  other functions}

  \item{k}{a numeric value, use as penalty coefficient for
  number of parameters in the fitted model}
}
\value{
a numeric value with the corresponding AIC or BIC value.
}
\description{
Calculate Akaike's information Criterion (AIC) or the
Bayesian information Criterion (BIC) for a model of object
of class.
}
\examples{
## Probit model
data("Workmroz")
probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc,
                 data = Workmroz , family = binomial('probit'))
summary(probit)

AIC(probit)
BIC(probit)
}
\seealso{
\code{\link[Rchoice]{Rchoice}}
}

