### R code from vignette source 'Rchoice.Rnw'

###################################################
### code chunk number 1: Rchoice.Rnw:417-418
###################################################
library("Rchoice")


###################################################
### code chunk number 2: Rchoice.Rnw:423-425
###################################################
data("Articles")
head(Articles,3)


###################################################
### code chunk number 3: Rchoice.Rnw:430-431
###################################################
help(Articles)


###################################################
### code chunk number 4: Rchoice.Rnw:436-439
###################################################
poisson<-Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, 
                 link = "poisson")
summary(poisson)


###################################################
### code chunk number 5: Rchoice.Rnw:448-450
###################################################
require(sandwich)
require(lmtest)


###################################################
### code chunk number 6: Rchoice.Rnw:453-454
###################################################
coeftest(poisson, vcov = sandwich)


###################################################
### code chunk number 7: Rchoice.Rnw:459-461
###################################################
vcov.stata <- vcovHC(poisson, type = "HC0") * nObs(poisson)/(nObs(poisson)-1)
coeftest(poisson, vcov = vcov.stata)


###################################################
### code chunk number 8: Rchoice.Rnw:468-469
###################################################
require(car)


###################################################
### code chunk number 9: Rchoice.Rnw:472-473
###################################################
deltaMethod(poisson, "phd/ment")


###################################################
### code chunk number 10: Rchoice.Rnw:490-493
###################################################
poisson.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, 
                       data = Articles, ranp = c(kid5="n", phd = "n", ment = "n"), 
                       link = "poisson")


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
### code chunk number 14: Rchoice.Rnw:521-523
###################################################
poisson.ran2 <- update(poisson.ran, ranp = c(kid5 = "u", phd = "t" , ment = "cn"))
summary(poisson.ran2)


###################################################
### code chunk number 15: Rchoice.Rnw:528-532
###################################################
poissonc.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment, data = Articles, 
                        ranp = c(kid5="n", phd = "n", ment = "n"), 
                        link = "poisson", correlation =  TRUE)
summary(poissonc.ran)


###################################################
### code chunk number 16: Rchoice.Rnw:537-539
###################################################
cov.Rchoice(poissonc.ran)
cor.Rchoice(poissonc.ran)


###################################################
### code chunk number 17: Rchoice.Rnw:556-562
###################################################
poissonH.ran <- Rchoice(art ~ fem + mar + kid5 + phd + ment | fem, 
                        data = Articles,
                        ranp = c(kid5="n", phd = "n", ment = "n"), 
                        link = "poisson",
                        correlation =  TRUE)
summary(poissonH.ran)


###################################################
### code chunk number 18: Rchoice.Rnw:567-568
###################################################
lrtest(poissonH.ran, poissonc.ran)


###################################################
### code chunk number 19: Rchoice.Rnw:602-603
###################################################
plot(poissonH.ran, par = "ment", type = "histogram", bin=0.005)


###################################################
### code chunk number 20: Rchoice.Rnw:608-609
###################################################
plot(poissonH.ran, par = "ment")


###################################################
### code chunk number 21: Rchoice.Rnw:616-617
###################################################
plot(poissonH.ran, par = "ment", ind = TRUE, id = seq(1, 50, 1))


###################################################
### code chunk number 22: Rchoice.Rnw:630-631
###################################################
data("Workmroz")


###################################################
### code chunk number 23: Rchoice.Rnw:635-642
###################################################
probit.ran <- Rchoice(lfp ~ k5 + k618 + age + wc + hc + lwg + inc, 
                      ranp = c(k5 = "n", hc = "n"), 
                      link = "probit", 
                      print.level = 1, 
                      data = Workmroz, 
                      R = 100)
summary(probit.ran)


###################################################
### code chunk number 24: Rchoice.Rnw:647-649
###################################################
probit2.ran <- update(probit.ran, start = c(0,0,0,0,0,0,0,0,0,0))
summary(probit2.ran)


###################################################
### code chunk number 25: Rchoice.Rnw:655-656
###################################################
data("Health")


###################################################
### code chunk number 26: Rchoice.Rnw:661-667
###################################################
oprobit.ran<-Rchoice(newhsat ~ age + educ + hhinc + married + hhkids, 
                      data = Health, link = "ordered probit", 
                      subset = year == 1988, 
                      ranp = c(age = "n", hhinc = "n"), print.level=1, 
                      start = rep(0,11))
summary(oprobit.ran)


