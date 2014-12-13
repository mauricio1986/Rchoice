### R code from vignette source 'Rchoice.Rnw'

###################################################
### code chunk number 1: Rchoice.Rnw:416-417
###################################################
library("Rchoice")


###################################################
### code chunk number 2: Rchoice.Rnw:422-424
###################################################
data("Articles")
head(Articles,3)


###################################################
### code chunk number 3: Rchoice.Rnw:429-430 (eval = FALSE)
###################################################
## help(Articles)


###################################################
### code chunk number 4: Rchoice.Rnw:435-438
###################################################
poisson <- Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, 
                 family = poisson)
summary(poisson)


###################################################
### code chunk number 5: Rchoice.Rnw:447-449
###################################################
library(sandwich)
library(lmtest)


###################################################
### code chunk number 6: Rchoice.Rnw:452-453
###################################################
coeftest(poisson, vcov = sandwich)


###################################################
### code chunk number 7: Rchoice.Rnw:458-460
###################################################
vcov.stata <- vcovHC(poisson, type = "HC0") * nObs(poisson)/(nObs(poisson)-1)
coeftest(poisson, vcov = vcov.stata)


###################################################
### code chunk number 8: Rchoice.Rnw:467-468
###################################################
library(car)


###################################################
### code chunk number 9: Rchoice.Rnw:471-472
###################################################
deltaMethod(poisson, "phd/ment")


###################################################
### code chunk number 10: Rchoice.Rnw:489-493
###################################################
poisson.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
                       data = Articles, ranp = c(kid5 = "n", phd = "n", ment = "n"), 
                       family = poisson,
                       R = 10)


###################################################
### code chunk number 11: Rchoice.Rnw:499-500
###################################################
summary(poisson.ran)


###################################################
### code chunk number 12: Rchoice.Rnw:505-506
###################################################
pnorm(coef(poisson.ran)["mean.kid5"]/coef(poisson.ran)["sd.kid5"])


###################################################
### code chunk number 13: Rchoice.Rnw:514-516
###################################################
waldtest(poisson.ran, poisson)
lrtest(poisson.ran, poisson)


###################################################
### code chunk number 14: Rchoice.Rnw:521-523 (eval = FALSE)
###################################################
## poisson.ran2 <- update(poisson.ran, ranp = c(kid5 = "u", phd = "t" , ment = "cn"))
## summary(poisson.ran2)


###################################################
### code chunk number 15: Rchoice.Rnw:528-533
###################################################
poissonc.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, 
                        ranp = c(kid5 = "n", phd = "n", ment = "n"), 
                        family = poisson, correlation =  TRUE,
                        R = 10)
summary(poissonc.ran)


###################################################
### code chunk number 16: Rchoice.Rnw:538-540
###################################################
cov.Rchoice(poissonc.ran)
cor.Rchoice(poissonc.ran)


###################################################
### code chunk number 17: Rchoice.Rnw:557-564
###################################################
poissonH.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment | fem, 
                        data = Articles,
                        ranp = c(kid5="n", phd = "n", ment = "n"), 
                        family = poisson,
                        correlation =  TRUE,
                        R = 10)
summary(poissonH.ran)


###################################################
### code chunk number 18: Rchoice.Rnw:569-570
###################################################
lrtest(poissonH.ran, poissonc.ran)


###################################################
### code chunk number 19: Rchoice.Rnw:604-605
###################################################
plot(poissonH.ran, par = "ment", type = "histogram", bin=0.005)


###################################################
### code chunk number 20: Rchoice.Rnw:610-611
###################################################
plot(poissonH.ran, par = "ment")


###################################################
### code chunk number 21: Rchoice.Rnw:618-619
###################################################
plot(poissonH.ran, par = "ment", ind = TRUE, id = seq(1, 50, 1))


###################################################
### code chunk number 22: Rchoice.Rnw:631-638 (eval = FALSE)
###################################################
## data("Workmroz")
## probit.ran <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc, 
##                       ranp = c(k5 = "n", hc = "n"), 
##                       family = binomial('probit'), 
##                       data = Workmroz, 
##                       R = 100)
## summary(probit.ran)


###################################################
### code chunk number 23: Rchoice.Rnw:644-645
###################################################
data("Health")


###################################################
### code chunk number 24: Rchoice.Rnw:650-656 (eval = FALSE)
###################################################
## oprobit.ran <- Rchoice(newhsat ~ age + educ + hhinc + married + hhkids, 
##                       data = Health, family = ordinal('probit'), 
##                       subset = year == 1988, 
##                       ranp = c(age = "n", hhinc = "n"), 
##                       start = rep(0, 11))
## summary(oprobit.ran)


