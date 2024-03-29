% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/effects.R
\name{effect}
\alias{effect}
\title{Get average marginal effects for heterokedastic binary models and IV probit models}
\usage{
effect(object, ...)
}
\arguments{
\item{object}{an object of class \code{hetprob} or \code{ivpml}.}

\item{...}{Additional arguments to be passed.}
}
\value{
Estimates of the average marginal effects computed as the average for each individual.
}
\description{
Obtain the average marginal effects from \code{hetprob} or \code{ivpml} class models.
}
\examples{
\donttest{
# Data
library("AER")
data("PSID1976")
PSID1976$lfp  <- as.numeric(PSID1976$participation == "yes")
PSID1976$kids <- with(PSID1976, factor((youngkids + oldkids) > 0,
                                      levels = c(FALSE, TRUE), 
                                      labels = c("no", "yes")))
PSID1976$finc <-  PSID1976$fincome / 10000
                                 
# Average marginal effects for heteroskedastic Probit model
labor_het <- hetprob(lfp ~  age + I(age^2) + finc + education + factor(kids) | 
                            factor(kids) + finc,              
                     data = PSID1976,                        
                     link = "probit")
eff_labor_het <- effect(labor_het)
summary(eff_labor_het)

# Average marginal effects for IV probit model 
# (nwincome is endogenous and heducation is the additional instrument)
PSID1976$nwincome <- with(PSID1976, (fincome - hours * wage)/1000)
fiml.probit <- ivpml(lfp ~  education + experience + I(experience^2) + age + 
                            youngkids + oldkids + nwincome |
                            education + experience + I(experience^2) + age + 
                            youngkids + oldkids + heducation, 
                     data = PSID1976)
summary(effect(fiml.probit))
summary(effect(fiml.probit, asf = FALSE))
}
}
