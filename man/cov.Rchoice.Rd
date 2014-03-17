\name{cov.Rchoice}
\alias{cor.Rchoice}
\alias{cov.Rchoice}
\title{Covariance and Correlation matrix of random parameters}
\usage{
cov.Rchoice(x)

cor.Rchoice(x)
}
\arguments{
  \item{x}{a object of class \code{Rchoice},}

  \item{...}{further arguments}
}
\value{
a matrix with the variance of the random parameters if
model is fitted with random coefficients or the correlation
matrix if argument \code{correlation = TRUE} in the fitted
model.
}
\description{
Computes the Variance-Covariance matrix and the Correlation
matrix of the random parameters
}
\references{
\itemize{ \item Greene, W. H. (2003). Econometric analysis.
Pearson Education India. \item Train, K. (2009). Discrete
choice methods with simulation. Cambridge university press.
}
}
\seealso{
\code{\link[Rchoice]{Rchoice}}
}

