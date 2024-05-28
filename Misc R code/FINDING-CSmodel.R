library(glmmTMB)
library(mgcv)
library(ggplot2)
library(dplyr)
data("chicago", package = "gamair")


sm_tmpd <- mgcv::smoothCon(s(tmpd),absorb.cons = TRUE,  data = chicago)[[1]]


re_tmpd <- mgcv::smooth2random(sm_tmpd, "", type = 2)


# Create a fake grouping variable
chicago$fake_group <- factor(rep(1, nrow(chicago)))



Xf_tmpd <- re_tmpd$Xf
Xr_tmpd <- re_tmpd$rand[[1]]


ftmb1 <- glmmTMB(formula = death ~Xf_tmpd + homdiag(0 + Xr_tmpd | fake_group), data = chicago, REML=TRUE)



summary(ftmb1)

gtmb1<-glmmTMB(death ~s(tmpd), data = chicago, REML=TRUE)




mgcv4<-gamm4(death~s(tmpd), data = chicago)


mgcv1<-gamm(death~s(tmpd), data = chicago, method="REML")


summary(mgcv1$lme)

stopifnot(all.equal(sigma(mgcv1$lme), sigma(ftmb1), tol = 1e-10))

## fitted values
stopifnot(all.equal(fitted(ftmb1), unname(fitted(mgcv1$lme)), tol = 1e-9))

c(mgcvlogLik = logLik(mgcv1$lme), gamm4logLik=logLik(mgcv4$mer) ,tprslogLik = logLik(gtmb1))

c(mgcvAIC = AIC(mgcv1$lme),gamm4AIC=AIC(mgcv4$mer), tprsAIC = AIC(gtmb1))

sigma(mgcv1$lme)
sigma(mgcv4$mer)
sigma(ftmb1)




## fitted values
stopifnot(all.equal(fitted(ftmb1), unname(fitted(mgcv4$mer)), tol = 2e-10))

c(mgcvlogLik = logLik(mgcv4$mer), tprslogLik = logLik(ftmb1))

c(mgcvAIC = AIC(mgcv1$lme), tprsAIC = AIC(ftmb1))






# cs -splines research





sm_tmpdcs <- mgcv::smoothCon(s(tmpd, bs="cs"),absorb.cons = TRUE,  data = chicago)[[1]]

re_tmpdcs <- mgcv::smooth2random(sm_tmpdcs, "", type = 2)



chicago$fake_group <- factor(rep(1, nrow(chicago)))



Xr_tmpdcs <- re_tmpdcs$rand[[1]]


ftmb1cs <- glmmTMB(formula = death ~ homdiag(0 + Xr_tmpdcs | fake_group), data = chicago, REML = TRUE)



mgcv1cs<-gamm(death~s(tmpd, bs="cs"), data = chicago,method="REML")



stopifnot(all.equal(sigma(mgcv1cs$lme), sigma(ftmb1cs), tol = 1e-10))

stopifnot(all.equal(fitted(ftmb1cs), unname(fitted(mgcv1cs$lme)), tol = 1e-12))

c(mgcvcslogLik = logLik(mgcv1cs$lme), cslogLik = logLik(ftmb1cs))

c(mgcvcsAIC = AIC(mgcv1cs$lme), csAIC = AIC(ftmb1cs))


sm_timecs <- mgcv::smoothCon(s(time, bs="cs"),absorb.cons = TRUE,  data = chicago)[[1]]

re_timecs <- mgcv::smooth2random(sm_timecs, "", type = 2)

Xr_timecs <- re_timecs$rand[[1]]

ftmb2cs <- glmmTMB(formula = death ~ homdiag(0 + Xr_tmpdcs | fake_group)+homdiag(0 + Xr_timecs | fake_group), data = chicago, REML = TRUE)



mgcv2cs<-gamm(death~s(tmpd, bs="cs")+s(time, bs="cs"), data = chicago,method="REML")

stopifnot(all.equal(sigma(mgcv2cs$lme), sigma(ftmb2cs), tol = 1e-10))

stopifnot(all.equal(fitted(ftmb2cs), unname(fitted(mgcv2cs$lme)), tol = 1e-12))

c(mgcvcslogLik = logLik(mgcv2cs$lme), cslogLik = logLik(ftmb2cs))

c(mgcvcsAIC = AIC(mgcv2cs$lme), csAIC = AIC(ftmb2cs))




# cc -splines research



sm_tmpdcc <- mgcv::smoothCon(s(tmpd, bs="cc"),absorb.cons = TRUE,  data = chicago)[[1]]

re_tmpdcc <- mgcv::smooth2random(sm_tmpdcc, "", type = 2)


chicago$fake_group <- factor(rep(1, nrow(chicago)))


Xr_tmpdcc <- re_tmpdcc$rand[[1]]


ftmb1cc <- glmmTMB(formula = death ~ homdiag(0 + Xr_tmpdcc | fake_group), data = chicago, REML = TRUE)

summary(ftmb1cc)






mgcv1cc<-gamm(death~s(tmpd, bs="cc"), data = chicago,method="REML")

summary(mgcv1cc$lme)


stopifnot(all.equal(sigma(mgcv1cc$lme), sigma(ftmb1cc), tol = 1e-10))


stopifnot(all.equal(fitted(ftmb1cc), unname(fitted(mgcv1cc$lme)), tol = 1e-10))


c(mgcvccLogLik = logLik(mgcv1cc$lme), ccLogLik = logLik(ftmb1cc))

c(mgcvccAIC = AIC(mgcv1cc$lme), ccAIC = AIC(ftmb1cc))


sm_timecc <- mgcv::smoothCon(s(time, bs="cc"),absorb.cons = TRUE,  data = chicago)[[1]]

re_timecc <- mgcv::smooth2random(sm_timecc, "", type = 2)

Xr_timecc <- re_timecc$rand[[1]]

ftmb2cc <- glmmTMB(formula = death ~ homdiag(0 + Xr_tmpdcc | fake_group)+homdiag(0 + Xr_timecc | fake_group), data = chicago, REML = TRUE)



mgcv2cc<-gamm(death~s(tmpd, bs="cc")+s(time, bs="cc"), data = chicago,method="REML")

stopifnot(all.equal(sigma(mgcv2cc$lme), sigma(ftmb2cc), tol = 1e-10))

stopifnot(all.equal(fitted(ftmb2cc), unname(fitted(mgcv2cc$lme)), tol = 1e-12))

c(mgcvcslogLik = logLik(mgcv2cc$lme), cslogLik = logLik(ftmb2cc))

c(mgcvcsAIC = AIC(mgcv2cc$lme), csAIC = AIC(ftmb2cc))


# re research
library(glmmTMB)
library(mgcv)

data("sleepstudy", package = "lme4")




sm_tmpdre <- mgcv::smoothCon(s(Subject, bs="re"),absorb.cons = TRUE,  data = sleepstudy)[[1]]


re_tmpdre <- mgcv::smooth2random(sm_tmpdre, "", type = 2)


# Create a fake grouping variable
sleepstudy$fake_group <- factor(rep(1, nrow(sleepstudy)))


Xr_tmpdre <- re_tmpdre$rand[[1]]


ftmb1re <- glmmTMB(formula = Reaction ~s(Days)+homdiag(0 + Xr_tmpdre | fake_group), data = sleepstudy, REML = TRUE)

summary(ftmb1re)


mgcv4re<-gamm4(Reaction~ s(Days)+s(Subject, bs="re"), data =sleepstudy)

summary(mgcv4re$mer)


stopifnot(all.equal(sigma(mgcv4re$mer), sigma(ftmb1re), tol = 1e-10))


stopifnot(all.equal(fitted(ftmb1re), unname(fitted(mgcv4re$mer)), tol = 2e-10))


c(mgcvreLogLik = logLik(mgcv4re$mer), reLogLik = logLik(ftmb1re))

c(mgcvreAIC = AIC(mgcv4re$mer), reAIC = AIC(ftmb1re))
sigma(ftmb1re)
sigma(mgcv4re$mer)

# Loglik and sigma not quite the same when using s() + homdiag(0 + Xr_tmpdre | fake_group)


ftmb2re <- glmmTMB(formula = Reaction ~homdiag(0 + Xr_tmpdre | fake_group), data = sleepstudy, REML = TRUE)

summary(ftmb2re)


mgcv2re<-gamm4(Reaction~ s(Subject, bs="re"), data =sleepstudy)

summary(mgcv2re$mer)

stopifnot(all.equal(sigma(mgcv2re$mer), sigma(ftmb2re), tol = 1e-10))


stopifnot(all.equal(fitted(ftmb2re), unname(fitted(mgcv2re$mer)), tol = 2e-10))


c(mgcvreLogLik = logLik(mgcv2re$mer), reLogLik = logLik(ftmb2re))

c(mgcvreAIC = AIC(mgcv2re$mer), reAIC = AIC(ftmb2re))
sigma(ftmb2re)
sigma(mgcv2re$mer)

# not the case when only using homdiag(0 + Xr_tmpdre | fake_group), no addition

