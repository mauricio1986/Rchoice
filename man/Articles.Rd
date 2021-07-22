\name{Articles}
\alias{Articles}
\docType{data}
\title{
Doctoral Publications
}
\description{
Data from research by Long(1990) that analizes the scientist's level of publications.
}
\usage{data(Articles)}
\format{
  A data frame with 915 observations on the following 6 variables:
  \describe{
    \item{\code{art}}{Articles during last 3 years of Ph.D.,}
    \item{\code{fem}}{1 if female scientist; else 0,}
    \item{\code{mar}}{1 if married; else 0,}
    \item{\code{kid5}}{Number of children 5 or younger,}
    \item{\code{phd}}{Prestige of Ph.D. department,}
    \item{\code{ment}}{Articles by mentor during last 3 years,}
  }
}

\source{
\itemize{
\item Long, J. S. (1990). The origins of sex differences in science. Social Forces, 68(4), 1297-1316.
\item Long, J. S. (1997). Regression models for categorical and limited dependent variables (Vol. 7). Sage.
\item Long, J. S., & Freese, J. (2006). Regression models for categorical and limited dependent variables using Stata. Stata Press, College Station, TX.
}
}
\examples{
data(Articles)
}
\keyword{datasets}
