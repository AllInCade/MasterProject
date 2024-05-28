data("chicago", package = "gamair")

library(mgcv)

library(glmmTMB)


# Convert 'tmpd' to random effect representation

# absorb.cons	:
# Set to TRUE in order to have identifiability constraints absorbed into the basis.
sm_tmpd <- mgcv::smoothCon(s(tmpd),absorb.cons = TRUE, data = chicago)[[1]]
re_tmpd <- mgcv::smooth2random(sm_tmpd, "", type = 2)


# Extract Xf and Xr matrices
Xf_tmpd <- re_tmpd$Xf
Xr_tmpd <- re_tmpd$rand[[1]]

# Create a fake grouping variable
chicago$ID <- factor(rep(1, nrow(chicago)))


ftmb1 <- glmmTMB(formula = death ~ Xf_tmpd +homdiag(0 +Xr_tmpd | ID), data = chicago, REML=TRUE)

summary(ftmb1)


gtmb0S <- glmmTMB(formula = death ~ s(tmpd), data = chicago,REML=TRUE)

summary(gtmb0S)


