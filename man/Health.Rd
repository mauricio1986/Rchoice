\name{Health}
\alias{Health}
\docType{data}
\title{
German Health Care Data
}
\description{
German Health Care Data, unbalanced panel. 
}
\usage{data(Health)}
\format{
  A data frame with 27326 observations on the following 27 variables.
  \describe{
    \item{\code{id}}{person identification number}
    \item{\code{female}}{female =1, male =0}
    \item{\code{year}}{calendar year of the observation}
    \item{\code{age}}{age in years}
    \item{\code{hsat}}{health satisfaction, 0 (low),...,10 (high)}
    \item{\code{handdum}}{handicapped = 1, 0 otherwise}
    \item{\code{handper}}{degree of handicap in percent; 0,100}
    \item{\code{hhinc}}{household nomial monthly net income in German marks/10,000}
    \item{\code{hhkids}}{children under age 16 in the household = 1; otherwise = 0}
    \item{\code{educ}}{years of schooling}
    \item{\code{married}}{married =1, otherwise = 0}
    \item{\code{haupts}}{highest schooling degree is Hauptschul degree = 1; otherwise = 0}
    \item{\code{reals}}{highest schooling degree is Realschul degree = 1, otherwise = 0}
    \item{\code{fachhs}}{highest schooling degree is Polytechical degree = 1; otherwise = 0}
    \item{\code{abitur}}{highest schooling degree is Abitur = 1; otherwise = 0}
    \item{\code{univ}}{highest schooling degree is university degree =1; otherwise = 0}
    \item{\code{working}}{ employed =1; otherwise = 0}
    \item{\code{bluec}}{blue-collar employee = 1; otherwise = 0}
    \item{\code{whitec}}{white-collar employeee =1; otherwise = 0}
    \item{\code{self}}{self-employed = 1; otherwise = 0}
    \item{\code{beamt}}{civil servant = 1; otherwise = 0}
    \item{\code{docvis}}{number of doctor visits in last three months}
    \item{\code{hospvis}}{number of hospital visits in last calendar year}
    \item{\code{public}}{insured in public health =1; otherwise = 0}
    \item{\code{addon}}{insured by add-on insurance =1; otherwise = 0}
    \item{\code{hsat2}}{40 observations on hsat recorded between 6 and 7 were changed to 7}
    \item{\code{newhsat}}{recording of hsat, (0-2) = 0, (3-5)=1, (6-8)=2, (9)=3 (10)=4}
  }
}

\source{
Riphahn, R. T., Wambach, A., & Million, A. (2003). Incentive effects in the demand for health care: a bivariate panel count data estimation. Journal of applied econometrics, 18(4), 387-405.
}
\references{
Greene, W. H. (2003). Econometric analysis. Pearson Education India.
}
\examples{
data(Health)
}
\keyword{datasets}
