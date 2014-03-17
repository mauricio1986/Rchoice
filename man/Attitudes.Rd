\name{Attitudes}
\alias{Attitudes}
\docType{data}
\title{
Attituded toward working mothers
}
\description{
In 1997 and 1989, the General Social Survey asked respondents to
evaluate the following statement: "A working mother can establish
 just as warm and secure a relationship with her children
 as a mother who does not work".
}
\usage{data(Attitudes)}
\format{
  A data frame with 2293 observations on the following 10 variables.
  \describe{
    \item{\code{warm}}{1 = Strongly disagree, 2 = disagree, 3 = agree, 4 = strongly agree}
    \item{\code{yr89}}{survey year: 1 =  1989; 0 = 1977}
    \item{\code{male}}{1 = male; 0 = female}
    \item{\code{white}}{1 = white; 0 = nonwhite}
    \item{\code{age}}{age in years}
    \item{\code{ed}}{years of education}
    \item{\code{prst}}{occupational prestige}
  }
}

\source{
\itemize{
 \item Clogg, C. C., & Shihadeh, E. S. (1994). Statistical models for ordinal variables. Thousand Oaks, CA: Sage Publications.
\item Long, J. S. (1997). Regression models for categorical and limited dependent variables (Vol. 7). Sage.
\item Long, J. S., & Freese, J. (2006). Regression models for categorical and limited dependent variables using Stata. Stata Press, College Station, TX.
}
}

\examples{
data(Attitudes)
}
\keyword{datasets}
