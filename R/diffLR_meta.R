

### fit the model ###

diffLR.meta <- function(lig, rec, geData.L, geData.R, design_matrix,
                        re.cont, re.disc, 
                        add0.l = F,
                        add1.l = F,
                        add0.r = F,
                        add1.r = F
                        ) {
  
  ## data preparation
  samples <- 1:nrow(design_matrix)
  effect <- colnames(design_matrix)
  
  data4l_list <- vector("list", length(samples))
  data4r_list <- vector("list", length(samples))
  
  for (i in samples) {
    expr.l <- t(geData.L[[i]][geData.L[[i]]$gene==lig, -1])
    expr.r <- t(geData.R[[i]][geData.R[[i]]$gene==rec, -1])
    
    # just in case the L/R gene doesn't exist in some samples
    if (length(expr.l)==0) {
      expr.l <- data.frame(expr=rep(0, nrow(expr.l)))
    } else {
      expr.l <- data.frame(expr=expr.l[,1])
      rownames(expr.l) <- NULL
    }
    if (length(expr.r)==0) {
      expr.r <- data.frame(expr=rep(0, nrow(expr.r)))
    } else {
      expr.r <- data.frame(expr=expr.r[,1])
      rownames(expr.r) <- NULL
    }
  
    expr.l <- bind_cols(expr.l, design_matrix[i, , drop=FALSE])
    expr.r <- bind_cols(expr.r, design_matrix[i, , drop=FALSE])
    
    expr.l$sample.id <- i
    expr.r$sample.id <- i
    
    data4l_list[[i]] <- expr.l
    data4r_list[[i]] <- expr.r
    rm(expr.l, expr.r)
  }
  
  data4l <- bind_rows(data4l_list)
  data4r <- bind_rows(data4r_list)
  rm(data4l_list, data4r_list)  
  
  ## add pseudo-observations
  
  data4l <- data4l %>% 
    mutate(b = ifelse(expr > 0, 1, 0))
  data4r <- data4r %>% 
    mutate(b = ifelse(expr > 0, 1, 0))
  
  if (add1.l) {
    data4l.extra <- data.frame(expr = 0, sample.id = samples, b = 1) %>% 
      bind_cols(design_matrix)
    data4l <- bind_rows(data4l, data4l.extra)
  }
  if (add0.l) {
    data4l.extra <- data.frame(expr = 0, sample.id = samples, b = 0) %>% 
      bind_cols(design_matrix)
    data4l <- bind_rows(data4l, data4l.extra)
  }
  if (add1.r) {
    data4r.extra <- data.frame(expr = 0, sample.id = samples, b = 1) %>% 
      bind_cols(design_matrix)
    data4r <- bind_rows(data4r, data4r.extra)
  }
  if (add0.r) {
    data4r.extra <- data.frame(expr = 0, sample.id = samples, b = 0) %>% 
      bind_cols(design_matrix)
    data4r <- bind_rows(data4r, data4r.extra)
  }
  
  ## linear model fitting
  
  data4l.cont <- data4l %>% 
    filter(b == 1) %>% 
    select(-b)
  data4r.cont <- data4r %>% 
    filter(b == 1) %>% 
    select(-b)
  
  # just in case the gene isn't expressed in some sample
  data.sup <- data.frame(sample.id = samples) %>% 
    bind_cols(design_matrix)
  
  data4l.cont <- full_join(data4l.cont, data.sup, by = c(effect, "sample.id"))
  data4r.cont <- full_join(data4r.cont, data.sup, by = c(effect, "sample.id"))
  data4l.cont[is.na(data4l.cont)] <- 0
  data4r.cont[is.na(data4r.cont)] <- 0
  
  # just in case very few cells express the relevant gene (grouping factor must be < number of observations)
  if (nrow(data4l.cont) <= length(samples)) {
    form <- paste0("expr ~ 1 + ",effect)
    fit.cont.l <- lm(form, data = data4l.cont)
  } else if (re.cont) {
    form <- paste0("expr ~ 1 + (1|sample.id) + ",effect)
    fit.cont.l <- lmer(form, data = data4l.cont)
  } else {
    form <- paste0("expr ~ 1 + ",effect)
    fit.cont.l <- lm(form, data = data4l.cont)
  }
  
  if (nrow(data4r.cont) <= length(samples)) {
    form <- paste0("expr ~ 1 + ",effect)
    fit.cont.r <- lm(form, data = data4r.cont)
  } else if (re.cont) {
    form <- paste0("expr ~ 1 + (1|sample.id) + ",effect)
    fit.cont.r <- lmer(form, data = data4r.cont)
  } else {
    form <- paste0("expr ~ 1 + ",effect)
    fit.cont.r <- lm(form, data = data4r.cont)
  }
  
  # if (re.cont) {
  #   form <- paste0("expr ~ 1 + (1|sample.id) + ",effect)
  #   
  #   fit.cont.l <- lmer(form, data = data4l.cont)
  #   fit.cont.r <- lmer(form, data = data4r.cont)
  #   
  # } else {
  #   form <- paste0("expr ~ 1 + ",effect)
  #   
  #   fit.cont.l <- lm(form, data = data4l.cont)
  #   fit.cont.r <- lm(form, data = data4r.cont)
  # }
  
  ## linear model testing
  
  test.cont <- diffLR.test(fit.cont.l, fit.cont.r, effect = effect)
  
  ## logistic model fitting
  
  data4l.disc <- data4l %>% 
    select(-expr)
  data4r.disc <- data4r %>% 
    select(-expr)
    
  if (re.disc) {
    form <- paste0("b ~ 1 + (1|sample.id) + ",effect)
    
    fit.disc.l <- glmer(form, data = data4l.disc, family = binomial)
    fit.disc.r <- glmer(form, data = data4r.disc, family = binomial)
    
  } else {
    form <- paste0("b ~ 1 + ",effect)
    
    fit.disc.l <- glm(form, data = data4l.disc, family = binomial)
    fit.disc.r <- glm(form, data = data4r.disc, family = binomial)
  } 
  
  ## logistic model testing
  
  test.disc <- diffLR.test.p(fit.disc.l, fit.disc.r, effect = effect)
  
  # joint test
  test.both <- append(test.cont, test.disc)  
  
  chisq.stat <- test.both$chisq.cont + test.both$chisq.disc
  
  test.both[["pval.joint"]] <- 1 - pchisq(chisq.stat, df = 2)
  test.both <- within(test.both, rm(chisq.cont, chisq.disc))
  
  test.both
}


### "Sobel" test
# linear model
diffLR.test <- function(fit.l, fit.r, effect) {
  coef.l <- summary(fit.l)$coefficients
  vcov.l <- vcov(fit.l)
  
  beta.l0 <- coef.l[rownames(coef.l)=="(Intercept)",1]
  beta.l1 <- coef.l[rownames(coef.l)==effect,1]
  vcov.l <- vcov.l[rownames(vcov.l) %in% c("(Intercept)",effect),colnames(vcov.l) %in% c("(Intercept)",effect)]
  
  coef.r <- summary(fit.r)$coefficients
  vcov.r <- vcov(fit.r)
  
  beta.r0 <- coef.r[rownames(coef.r)=="(Intercept)",1]
  beta.r1 <- coef.r[rownames(coef.r)==effect,1]
  vcov.r <- vcov.r[rownames(vcov.r) %in% c("(Intercept)",effect),colnames(vcov.r) %in% c("(Intercept)",effect)]
  
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
diffLR.test.p <- function(fit.l, fit.r, effect) {
  coef.l <- summary(fit.l)$coefficients
  vcov.l <- vcov(fit.l)
  
  beta.l0 <- coef.l[rownames(coef.l)=="(Intercept)",1]
  beta.l1 <- coef.l[rownames(coef.l)==effect,1]
  vcov.l <- vcov.l[rownames(vcov.l) %in% c("(Intercept)",effect),colnames(vcov.l) %in% c("(Intercept)",effect)]
  
  coef.r <- summary(fit.r)$coefficients
  vcov.r <- vcov(fit.r)
  
  beta.r0 <- coef.r[rownames(coef.r)=="(Intercept)",1]
  beta.r1 <- coef.r[rownames(coef.r)==effect,1]
  vcov.r <- vcov.r[rownames(vcov.r) %in% c("(Intercept)",effect),colnames(vcov.r) %in% c("(Intercept)",effect)]
  
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
