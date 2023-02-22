### "Sobel" test
# linear model
diffLR.test <- function(fit.l, fit.r, group) {
  coef.l <- summary(fit.l)$coefficients
  vcov.l <- vcov(fit.l)
  
  beta.l0 <- coef.l[rownames(coef.l)=="(Intercept)",1]
  beta.l1 <- coef.l[rownames(coef.l)==group,1]
  vcov.l <- vcov.l[rownames(vcov.l) %in% c("(Intercept)",group),colnames(vcov.l) %in% c("(Intercept)",group)]
  
  coef.r <- summary(fit.r)$coefficients
  vcov.r <- vcov(fit.r)
  
  beta.r0 <- coef.r[rownames(coef.r)=="(Intercept)",1]
  beta.r1 <- coef.r[rownames(coef.r)==group,1]
  vcov.r <- vcov.r[rownames(vcov.r) %in% c("(Intercept)",group),colnames(vcov.r) %in% c("(Intercept)",group)]
  
  # gratitude
  g <- c(-beta.r1, -beta.r0-beta.r1, -beta.l1, -beta.l0-beta.l1)
  
  # numerator of Z statistic
  delta <- beta.l0*beta.r0 - (beta.l0+beta.l1)*(beta.r0+beta.r1)
  # variance
  V2 <- t(g)%*%bdiag(vcov.l, vcov.r)%*%g
  
  chisq.stat <- as.numeric(delta^2/V2)
  pval <- 1 - pchisq(chisq.stat, df = 1)
  
  list(pval.cont=pval, delta.cont=-delta, chisq.cont=chisq.stat)
}

# logistic model
diffLR.test.p <- function(fit.l, fit.r, group) {
  coef.l <- summary(fit.l)$coefficients
  vcov.l <- vcov(fit.l)
  
  beta.l0 <- coef.l[rownames(coef.l)=="(Intercept)",1]
  beta.l1 <- coef.l[rownames(coef.l)==group,1]
  vcov.l <- vcov.l[rownames(vcov.l) %in% c("(Intercept)",group),colnames(vcov.l) %in% c("(Intercept)",group)]
  
  coef.r <- summary(fit.r)$coefficients
  vcov.r <- vcov(fit.r)
  
  beta.r0 <- coef.r[rownames(coef.r)=="(Intercept)",1]
  beta.r1 <- coef.r[rownames(coef.r)==group,1]
  vcov.r <- vcov.r[rownames(vcov.r) %in% c("(Intercept)",group),colnames(vcov.r) %in% c("(Intercept)",group)]
  
  g <- c(expit(beta.r0)*expit.d(beta.l0)-expit(beta.r0+beta.r1)*expit.d(beta.l0+beta.l1),
         -expit(beta.r0+beta.r1)*expit.d(beta.l0+beta.l1),
         expit(beta.l0)*expit.d(beta.r0)-expit(beta.l0+beta.l1)*expit.d(beta.r0+beta.r1),
         -expit(beta.l0+beta.l1)*expit.d(beta.r0+beta.r1)
  )
  
  delta <- expit(beta.l0)*expit(beta.r0) - expit(beta.l0+beta.l1)*expit(beta.r0+beta.r1)
  
  V2 <- t(g)%*%bdiag(vcov.l, vcov.r)%*%g
  
  chisq.stat <- as.numeric(delta^2/V2)
  pval <- 1 - pchisq(chisq.stat, df = 1)
  
  list(pval.disc=pval, delta.disc=-delta, chisq.disc=chisq.stat)
}

expit <- function(x) {
  p <- 1/(1 + exp(-x))
  p
}

# 1st derivative of expit function
expit.d <- function(x) {
  p <- expit(x)
  d <- p*(1-p)
  d
}
