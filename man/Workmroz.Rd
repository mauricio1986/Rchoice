\name{Workmroz}
\alias{Workmroz}
\docType{data}
\title{
Labor Force Participation
}
\description{
Data extracted by Mroz(1987) from the 197 Panel Study of Income Dynacmis. The sample consists of 753 white, married women between the ages of 30 and 60.
}
\usage{data(Workmroz)}
\format{
  A data frame with 753 observations on the following 9 variables.
  \describe{
    \item{\code{lfp}}{1 if wife is in the paid labor force; else 0}
    \item{\code{k5}}{Number of children ages 5 and younger}
    \item{\code{k618}}{Number of children ages 6 to 18}
    \item{\code{age}}{Wife's age in years}
    \item{\code{wc}}{1 if wife attended college; else 0}
    \item{\code{hc}}{1 if husband attended college; else 0}
    \item{\code{lwg}}{Log of wife's estimated wage rate}
    \item{\code{inc}}{Family income excluding wife's wage}
    \item{\code{linc}}{Log of Family income excluding wife's wage}
  }
}
\source{
Mroz, T. A. (1987). The sensitivity of an empirical model of married women's hours of work to economic and statistical assumptions. Econometrica, 55(4), 765-799
}
\examples{
data(Workmroz)
}
\keyword{datasets}
