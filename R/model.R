
# check the specification of arguments


# existing.lr is all of the ligand-receptor pairs that are expressed given the cell type pairs of interest
# geData.L is the gene expression matrix for ligand genes
# geData.R is the gene expression matrix for receptor genes
# the order of samples in geData.L and geData.R are matched

stat_analysis <- function(existing.lr, geData.L, geData.R, group.vec, ...) {

  message("Model fitting and statistical testing ...")

  result_list <- vector("list", length = nrow(existing.lr) )

  message("Fitting the models ... \n")

  for (j in 1:nrow(existing.lr)) {
    progress(j, max.value = nrow(existing.lr))

    lig <- existing.lr$ligand[j]
    rec <- existing.lr$receptor[j]

    # # create a df to store the statistics
    # data.sum <- data.frame(sample.id=samples,
    #                        count.cL=NA, count.cR=NA, # no of cells
    #                        prop.L=NA, prop.R=NA, # expression rates
    #                        mean.pos.L=NA, mean.pos.R=NA, # pos mean expression
    #                        mean.L=NA, mean.R=NA # overall mean expression
    # )
    #
    # for (i in samples) {
    #   expr.l <- t(geData.L[[i]][geData.L[[i]]$gene==lig, -1])
    #   expr.r <- t(geData.R[[i]][geData.R[[i]]$gene==rec, -1])
    #
    #   # no of cells
    #   data.sum$count.cL[i] <- ncol(geData.L[[i]])-1
    #   data.sum$count.cR[i] <- ncol(geData.R[[i]])-1
    #
    #   # just in case the L/R gene doesn't exist in some samples
    #   if (length(expr.l)==0) {
    #     expr.l <- rep(0,data.sum$count.cL[i])
    #   } else {
    #     expr.l <- c(expr.l)
    #   }
    #   if (length(expr.r)==0) {
    #     expr.r <- rep(0,data.sum$count.cR[i])
    #   } else {
    #     expr.r <- c(expr.r)
    #   }
    #
    #   # prop
    #   data.sum$prop.L[i] <- sum(expr.l>0)/length(expr.l)
    #   data.sum$prop.R[i] <- sum(expr.r>0)/length(expr.r)
    #
    #   # positive mean
    #   data.sum$mean.pos.L[i] <- mean(expr.l[expr.l>0])
    #   data.sum$mean.pos.R[i] <- mean(expr.r[expr.r>0])
    #
    #   # overall mean
    #   data.sum$mean.L[i] <- mean(expr.l)
    #   data.sum$mean.R[i] <- mean(expr.r)
    # }
    # # convert NaN to 0
    # data.sum[is.na(data.sum)] <- 0

    props.l <- numeric(length(samples))
    props.r <- numeric(length(samples))

    for (i in samples) {
      expr.l <- t(geData.L[[i]][geData.L[[i]]$gene==lig, -1])
      expr.r <- t(geData.R[[i]][geData.R[[i]]$gene==rec, -1])

      # just in case the L/R gene doesn't exist in some samples
      if (length(expr.l)==0) {
        expr.l <- rep(0,ncol(geData.L[[i]])-1)
      } else {
        expr.l <- c(expr.l)
      }
      if (length(expr.r)==0) {
        expr.r <- rep(0,ncol(geData.R[[i]])-1)
      } else {
        expr.r <- c(expr.r)
      }

      # prop
      props.l[i] <- sum(expr.l>0)/length(expr.l)
      props.r[i] <- sum(expr.r>0)/length(expr.r)

    }

    ### check if psuedo-observations are needed
    add0.l <- F
    add1.l <- F
    add0.r <- F
    add1.r <- F

    grps <- unique(group.vec)

    # ligand:
    # categorize the elements in props.l aka which group do they belong to?
    A <- all(props.l[group.vec == grps[1]] == 0)
    B <- all(props.l[group.vec == grps[2]] == 0)
    C <- all(props.l[group.vec == grps[1]] == 1)
    D <- all(props.l[group.vec == grps[2]] == 1)
    if (A | B) add1.l <- T # if one group only contains 0, then we should add 1 to this group
    if (C | D) add0.l <- T

    # receptor:

    A <- all(props.r[group.vec == grps[1]] == 0)
    B <- all(props.r[group.vec == grps[2]] == 0)
    C <- all(props.r[group.vec == grps[1]] == 1)
    D <- all(props.r[group.vec == grps[2]] == 1)
    if (A | B) add1.r <- T
    if (C | D) add0.r <- T

    ### fit the model and conduct tests

    #########################################################################

    ## data preparation

    data4l_list <- vector("list", length(samples))
    data4r_list <- vector("list", length(samples))

    for (i in samples) {
      expr.l <- t(geData.L[[i]][geData.L[[i]]$gene==lig, -1])
      expr.r <- t(geData.R[[i]][geData.R[[i]]$gene==rec, -1])

      # just in case the L/R gene doesn't exist in some samples
      # expr.l is a n*1 dataframe. Its nrow() always exists.
      # Its length() doesn't exist if the L/R gene doesn't exist.
      # length() = the number of elements in the matrix.
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

      expr.l <- bind_cols(expr.l, effect=group.vec[i])
      expr.r <- bind_cols(expr.r, effect=group.vec[i])

      expr.l$sample.id <- i
      expr.r$sample.id <- i

      data4l_list[[i]] <- expr.l
      data4r_list[[i]] <- expr.r
      rm(expr.l, expr.r)
    }

    data4l <- bind_rows(data4l_list)
    data4r <- bind_rows(data4r_list)
    #rm(data4l_list, data4r_list)

    ## add pseudo-observations

    data4l <- data4l %>%
      mutate(b = ifelse(expr > 0, 1, 0))
    data4r <- data4r %>%
      mutate(b = ifelse(expr > 0, 1, 0))

    if (add1.l) {
      data4l.extra <- data.frame(expr = 0, sample.id = samples, b = 1) %>%
        bind_cols(group.vec)
      data4l <- bind_rows(data4l, data4l.extra)
    }
    if (add0.l) {
      data4l.extra <- data.frame(expr = 0, sample.id = samples, b = 0) %>%
        bind_cols(group.vec)
      data4l <- bind_rows(data4l, data4l.extra)
    }
    if (add1.r) {
      data4r.extra <- data.frame(expr = 0, sample.id = samples, b = 1) %>%
        bind_cols(group.vec)
      data4r <- bind_rows(data4r, data4r.extra)
    }
    if (add0.r) {
      data4r.extra <- data.frame(expr = 0, sample.id = samples, b = 0) %>%
        bind_cols(group.vec)
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
      bind_cols(effect=group.vec)

    data4l.cont <- full_join(data4l.cont, data.sup, by = c("effect", "sample.id"))
    data4r.cont <- full_join(data4r.cont, data.sup, by = c("effect", "sample.id"))
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



    #########################################################################

    result_sublist[['interaction']] <- paste(lig, rec, sep = "_")

    result_list[[j]] <- result_sublist

  }
  message("\n Done. \n")

  ## includes pvals, deltas, interactions
  # list to df
  data.result <- as.data.frame(do.call(rbind, result_list))
  # convert each element of the df from list to character
  data.result <- as.data.frame(sapply(data.result, unlist))
  # reorder the cols
  data.result <- data.result[,c(6,1:5)]
  # character to numeric
  data.result[,2:6] <- sapply(data.result[,2:6], as.numeric)

  #
  message(paste0("Time collapsed: ",round(difftime(Sys.time(), start.time, units = "mins"), digits = 1)," mins \n"))

  ### return a data frame
  return(data.result)

}



