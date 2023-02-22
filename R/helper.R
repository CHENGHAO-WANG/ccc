
### check the dimensions and classes input

precheck <- function(geData.list, ctData.list, group.vec, ref, env) {

  message("Checking data dimensions ...")

  if (!is.vector(group.vec) & !is.factor(group.vec)) stop("group.vec must be a vector or a factor")

  n.grp <- length(unique(group.vec))

  if (n.grp > 6) warning(paste0("Too many group levels. <= 6 is recommended. Currently ",n.grp))
  if (n.grp < 2) stop("Less than 2 group levels. Check your group.vec input or try diffLR2")

  if (!(ref %in% group.vec)) stop("group.vec doesn't contain ref")

  if (!is.list(geData.list) | !is.list(ctData.list)) stop("geData.list and ctData.list must be lists of matrices or data.frames")

  if (!all(sapply(list(length(geData.list),length(ctData.list)), function(x) x == length(group.vec) ))) {
    stop("geData.list, ctData.list, group.vec indicate different number of samples")
  }

  samples <- 1:length(group.vec)

  for (i in samples) {
    if (ncol(geData.list[[i]]) != nrow(ctData.list[[i]])) stop(paste0("geData.list and ctData.list indicate different cell number in sample ",i))
  }

  assign("samples", samples, envir = env)
  assign("n.grp", n.grp, envir = env)

  # return(NULL)
}

### rename the cols of input and convert them to data.frames

REname <- function(geData.list, ctData.list, lrdb, env) {

  samples <- 1:length(geData.list)

  for (i in samples) {
    geData.list[[i]] <- as.data.frame(geData.list[[i]])
    geData.list[[i]] <- cbind(gene = rownames(geData.list[[i]]), data.frame(geData.list[[i]], row.names = NULL))

    ctData.list[[i]] <- as.data.frame(ctData.list[[i]])
    colnames(ctData.list[[i]]) <- c("cell", "celltype")
  }

  colnames(lrdb) <- c("ligand", "receptor")

  assign("geData.list", geData.list, envir = env)
  assign("ctData.list", ctData.list, envir = env)
  assign("lrdb", lrdb, envir = env)

  # return(NULL)
}

### separate the gene expression data of the cell types of interest
### only keep the ligand/receptor genes

preprocess <- function(geData.list, ctData.list, lrdb, celltypes.L, celltypes.R, env) {

  samples <- 1:length(geData.list)

  geData.L <- vector("list", length(samples))
  geData.R <- vector("list", length(samples))

  for (i in samples) {

    if ( length(intersect(ctData.list[[i]]$celltype, celltypes.L)) == 0 ) stop("celltypes.L don't exist in sample ",i)
    if ( length(intersect(ctData.list[[i]]$celltype, celltypes.R)) == 0 ) stop("celltypes.R don't exist in sample ",i)

    if (!setequal(intersect(ctData.list[[i]]$celltype, celltypes.L), celltypes.L)) warning(paste0("In sample ",i,", some cell type(s) in celltypes.L don't exist"))
    if (!setequal(intersect(ctData.list[[i]]$celltype, celltypes.R), celltypes.R)) warning(paste0("In sample ",i,", some cell type(s) in celltypes.R don't exist"))

    cell.L <- ctData.list[[i]][ctData.list[[i]]$celltype %in% celltypes.L, "cell"]

    cols <- c("gene", cell.L)
    geData.L[[i]] <- geData.list[[i]][,cols] %>%
      filter(gene %in% lrdb$ligand)

    cell.R <- ctData.list[[i]][ctData.list[[i]]$celltype %in% celltypes.R, "cell"]

    cols <- c("gene", cell.R)
    geData.R[[i]] <- geData.list[[i]][,cols] %>%
      filter(gene %in% lrdb$receptor)
  }

  assign("geData.L", geData.L, envir = env)
  assign("geData.R", geData.R, envir = env)

  rm(geData.list, ctData.list, envir = env)

  # return(NULL)
}

### filter ligand-receptor pairs that are expressed

prefilter <- function(lrdb, min.prop, filter.all, geData.L, geData.R) {

  lrdb <- lrdb %>% select(ligand, receptor)

  if (filter.all) {
    message("Finding L-R existed in all samples ...")
    # "exist" means both ligand and receptor are expressed > min.prop in all samples
    existing.lr <- t(apply(lrdb, 1, find_lr, geData.l=geData.L, geData.r=geData.R, min.prop = min.prop))

  } else {
    message("Finding L-R existed in at least one sample ...")
    # "exist" means both ligand and receptor are expressed > min.prop in at least one sample
    existing.lr <- t(apply(lrdb, 1, find_lr0, geData.l=geData.L, geData.r=geData.R, min.prop = min.prop))
  }

  existing.lr <- as.data.frame(existing.lr)
  colnames(existing.lr) <- c("ligand","receptor")
  existing.lr <- na.omit(existing.lr)

  return(existing.lr)
}

# find LR expressed in all samples (> min.prop)

find_lr <- function(row_lr, geData.l, geData.r, min.prop) {
  if ( all(sapply(geData.l, g_sample, g=row_lr[1], threshold = min.prop)) & all(sapply(geData.r, g_sample, g=row_lr[2], threshold = min.prop)) ) {
    return(row_lr)
  } else {
    rep(NA,2)
  }
}

# find LR expressed in at least one sample (> min.prop)

find_lr0 <- function(row_lr, geData.l, geData.r, min.prop) {
  if ( any(sapply(geData.l, g_sample, g=row_lr[1], threshold = min.prop) * sapply(geData.r, g_sample, g=row_lr[2], threshold = min.prop)) ) {
    return(row_lr)
  } else {
    rep(NA,2)
  }
}

# whether a gene g is expressed in the given sample

g_sample <- function(g, gedata.i, threshold) {
  if (g %in% gedata.i$gene) {
    if ( sum(gedata.i[gedata.i$gene==g,-1] > 0)/(ncol(gedata.i)-1) > threshold) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}


### obtain the list of all cell type pairs across the samples

get_celltype_pairs <- function(ctData.list) {
  # obtain the shared cell types across the samples
  celltype.list <- Reduce( intersect, get_celltype_list(ctData.list) )

  celltype.pairs <- combn(celltype.list, 2)

  celltype.pairs2 <- celltype.pairs[,2:1]

  celltype.pairs3 <- cbind(celltype.list, celltype.list)

  celltype.pairs4 <- rbind(celltype.pairs, celltype.pairs2, celltype.pairs3)

  return(celltype.pairs4)
}

# obtain the list of cell types in each sample

get_celltype_list <- function(ctData.list) {
    lapply(ctData.list, function(x) unique(x[,2,drop=TRUE]) )
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
