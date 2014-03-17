\name{plot.Rchoice}
\alias{plot.Rchoice}
\title{Plot random parameters}
\usage{
\method{plot}{Rchoice}(x, par = NULL, ind = FALSE, id = NULL,
  type = c("density", "histogram"), bin = 1, adjust = 1, ...)
}
\arguments{
  \item{x}{a object of class \code{Rchoice},}

  \item{par}{a string giving the name of the variable with
  random parameter,}

  \item{type}{a string indicating the type of distribution:
  it can be a \code{histogram} or a \code{density} of the
  conditional expectation of the random coefficients,}

  \item{ind}{a boolean. If \code{TRUE}, a 95% interval of
  conditional distribution for each individual is plotted.
  As default, the conditional expectation of \code{par} for
  the first 10 individual is plotted,}

  \item{id}{only relevant if \code{ind} is not \code{NULL}.
  This is a vector indicating the position of the
  individual for whom the user want to plot the conditional
  coefficients,}

  \item{bin}{bin of histrogram,}

  \item{adjust}{bandwidth for the kernel density,}

  \item{...}{further arguments to be passed to \code{qplot}
  or \code{plotCI},}
}
\value{
a plot with the distribution or a confident interval of the
conditional random coefficients.
}
\description{
Plot the conditional expectation of random parameters
estimated by \code{Rchoice}.
}
\references{
\itemize{ \item Greene, W. H. (2003). Econometric analysis.
Pearson Education India. \item Train, K. (2009). Discrete
choice methods with simulation. Cambridge university press.
}
}
\seealso{
\code{\link[Rchoice]{Rchoice}},
\code{\link[ggplot2]{ggplot2}}
}

