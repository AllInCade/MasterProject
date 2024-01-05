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


ftmb1 <- glmmTMB(formula = death ~Xf_tmpd + homdiag(0 + Xr_tmpd | fake_group), data = chicago, REML = TRUE)


mgcv1<-gam(death~s(tmpd), data = chicago)

summary(mgcv1)



stopifnot(all.equal(sigma(mgcv1), sigma(ftmb1), tol = 1e-5))

## fitted values
stopifnot(all.equal(fitted(ftmb1), unname(fitted(mgcv1)), tol = 2e-5))

c(mgcvlogLik = logLik(mgcv1), tprslogLik = logLik(ftmb1))

c(mgcvAIC = AIC(mgcv1), tprsAIC = AIC(ftmb1))






# cs -splines reasearch





sm_tmpdcs <- mgcv::smoothCon(s(tmpd, bs="cs"),absorb.cons = TRUE,  data = chicago)[[1]]

re_tmpdcs <- mgcv::smooth2random(sm_tmpdcs, "", type = 2)



chicago$fake_group <- factor(rep(1, nrow(chicago)))



Xr_tmpdcs <- re_tmpdcs$rand[[1]]


ftmb1cs <- glmmTMB(formula = death ~ homdiag(0 + Xr_tmpdcs | fake_group), data = chicago, REML = TRUE)



mgcv1cs<-gam(death~s(tmpd, bs="cs"), data = chicago)



stopifnot(all.equal(sigma(mgcv1cs), sigma(ftmb1cs), tol = 1e-8))


stopifnot(all.equal(fitted(ftmb1cs), unname(fitted(mgcv1cs)), tol = 2e-5))

c(mgcvcslogLik = logLik(mgcv1cs), cslogLik = logLik(ftmb1cs))

c(mgcvcsAIC = AIC(mgcv1cs), csAIC = AIC(ftmb1cs))



# cc -splines reasearch



sm_tmpdcc <- mgcv::smoothCon(s(tmpd, bs="cc"),absorb.cons = TRUE,  data = chicago)[[1]]

re_tmpdcc <- mgcv::smooth2random(sm_tmpdcc, "", type = 2)


chicago$fake_group <- factor(rep(1, nrow(chicago)))


Xr_tmpdcc <- re_tmpdcc$rand[[1]]


ftmb1cc <- glmmTMB(formula = death ~ homdiag(0 + Xr_tmpdcc | fake_group), data = chicago, REML = TRUE)




mgcv1cc<-gam(death~s(tmpd, bs="cc"), data = chicago)

summary(mgcv1cc)


stopifnot(all.equal(sigma(mgcv1cc), sigma(ftmb1cc), tol = 1e-8))


stopifnot(all.equal(fitted(ftmb1cc), unname(fitted(mgcv1cc)), tol = 2e-5))


c(mgcvLogLik = logLik(mgcv1cc), ccLogLik = logLik(ftmb1cc))

c(mgcvAIC = AIC(mgcv1cc), ccAIC = AIC(ftmb1cc))




