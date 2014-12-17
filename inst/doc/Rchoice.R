## ----echo = FALSE, include =  FALSE--------------------------------------
# set global chunk options
library(knitr)
opts_chunk$set(tidy = FALSE, fig.width = 7, fig.height = 5, comment = "")

## ----message = FALSE-----------------------------------------------------
library("Rchoice")

## ------------------------------------------------------------------------
data("Articles")
head(Articles,3)

## ----eval=FALSE----------------------------------------------------------
#  help(Articles)

## ------------------------------------------------------------------------
poisson <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
                   data = Articles, 
                   family = poisson)
summary(poisson)

## ----message=FALSE-------------------------------------------------------
library(sandwich)
library(lmtest)
coeftest(poisson, vcov = sandwich)

## ------------------------------------------------------------------------
vcov.stata <- vcovHC(poisson, type = "HC0") * nObs(poisson)/(nObs(poisson)-1)
coeftest(poisson, vcov = vcov.stata)

## ----message=FALSE, warning = FALSE--------------------------------------
library(car)
deltaMethod(poisson, "phd/ment")

## ------------------------------------------------------------------------
data("Workmroz")
probit <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + linc,
                  data = Workmroz,
                  family = binomial('probit'))
summary(probit)

## ------------------------------------------------------------------------
data("Health")
Health$linc <- log(Health$hhinc)
ologit <- Rchoice(newhsat ~ age + educ + married + hhkids + linc,
                   data = Health[1:2000, ],
                   family = ordinal('logit'))
summary(ologit)

## ------------------------------------------------------------------------
poisson.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
                       data = Articles,  
                       family = poisson,
                       ranp = c(kid5 = "n", phd = "n", ment = "n"))

## ------------------------------------------------------------------------
summary(poisson.ran)

## ------------------------------------------------------------------------
pnorm(coef(poisson.ran)["mean.kid5"]/coef(poisson.ran)["sd.kid5"])

## ------------------------------------------------------------------------
waldtest(poisson.ran, poisson)
lrtest(poisson.ran, poisson)

## ----eval = FALSE--------------------------------------------------------
#  poisson.ran2 <- update(poisson.ran,
#                         ranp = c(kid5 = "u", phd = "t" , ment = "cn"),
#                         R = 10)

## ----echo=TRUE, message=FALSE--------------------------------------------
library(memisc)

## ----eval = FALSE--------------------------------------------------------
#  mtable("model 1"= poisson.ran, "model 2" = poisson.ran2,
#         summary.stats = c("N", "Log-likelihood", "BIC", "AIC"))

## ------------------------------------------------------------------------
poissonc.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
                        data = Articles, 
                        ranp = c(kid5 = "n", phd = "n", ment = "n"), 
                        family = poisson, 
                        correlation =  TRUE,
                        R = 10)
summary(poissonc.ran)

## ------------------------------------------------------------------------
cov.Rchoice(poissonc.ran)
cor.Rchoice(poissonc.ran)

## ------------------------------------------------------------------------
se.cov.Rchoice(poissonc.ran)

## ------------------------------------------------------------------------
se.cov.Rchoice(poissonc.ran, sd = TRUE)

## ------------------------------------------------------------------------
data('Unions', package = 'pglm')
Unions$lwage <- log(Unions$wage)

## ------------------------------------------------------------------------
union.ran <- Rchoice(union ~ age + exper + rural + lwage,
                     data = Unions[1:2000, ],
                     family = binomial('probit'),
                     ranp = c(constant = "n", lwage = "t"),
                     R = 10,
                     panel = TRUE,
                     index = "id",
                     print.init = TRUE)

## ------------------------------------------------------------------------
summary(union.ran)

## ----eval = FALSE--------------------------------------------------------
#  oprobit.ran <- Rchoice(newhsat ~ age + educ + married + hhkids + linc,
#                        data = Health[1:2000, ],
#                        family = ordinal('probit'),
#                        ranp = c(constant = "n", hhkids = "n", linc = "n"),
#                        panel = TRUE,
#                        index = "id",
#                        R = 100,
#                        print.init = TRUE)
#  summary(oprobit.ran)

## ----eval = FALSE--------------------------------------------------------
#  poissonH.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment | fem + phd,
#                          data = Articles,
#                          ranp = c(kid5 = "n", phd = "n", ment = "n"),
#                          mvar = list(phd = c("fem"), ment = c("fem", "phd")),
#                          family = poisson,
#                          R = 10)

## ----eval = FALSE--------------------------------------------------------
#  summary(poissonH.ran)

## ----eval = FALSE--------------------------------------------------------
#  lrtest(poissonH.ran, poisson.ran)

## ----eval =  FALSE-------------------------------------------------------
#  plot(union.ran, par = "lwage", type = "density")

## ----plot1, echo = FALSE-------------------------------------------------
plot(union.ran, par = "lwage", type = "density")

## ----eval =  FALSE-------------------------------------------------------
#  plot(union.ran, par = "lwage", ind =  TRUE, id = 1:20, col = "blue")

## ----plot2, echo = FALSE-------------------------------------------------
plot(union.ran, par = "lwage", ind =  TRUE, id = 1:20, col = "blue")

## ------------------------------------------------------------------------
bi.wage <- effect.Rchoice(union.ran, par = "lwage", effect = "ce")

## ------------------------------------------------------------------------
summary(bi.wage$mean)
summary(bi.wage$sd.est)

