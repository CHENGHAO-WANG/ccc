
# geData: a list of dataframes, where each df storing the normalized gene expression matrix in each sample.
# For each df, 1st col for gene names, then each col for a cell (col name = cell id)
# ctData: a list of dataframes, where each df storing the cell type for each cell in each sample.
# For each df, col names = c("cell", "celltype")
# lr: a df storing the ligand-receptor pairs, with col names = c("ligand", "receptor")
# design_matrix: a df of sample level design matrix
# cdr_adj.cont: whether adjust for CDR in linear model (default is TRUE)
# cdr_adj.disc: whether adjust for CDR in logistic model (default is FALSE)
# re.cont: Should random intercept terms be included in the linear model?
# re.disc: Should random intercept terms be included in the logistic model?
# min.prop: threshold for highly expressed l-r pairs.

# return a data.frame



diffLR <- function(geData, ctData, lr, celltype.L, celltype.R, design_matrix,
                   re.cont=TRUE, re.disc=TRUE,
                   min.prop=0.10) {
  ######################
  data("ref", envir=environment())
  ######################
  start.time <- Sys.time()

  cat("Checking data structure ... \n")

  ### check packages loading
  if (!require(dplyr)) install.packages('dplyr')
  library(dplyr)
  if (!require(lme4)) install.packages('lme4')
  library(lme4)
  if (!require(svMisc)) install.packages('svMisc')
  library(svMisc)

  ### check the design matrix; there should be only 1 column
  if (ncol(design_matrix) != 1) {
    stop("Can only compare two groups at a time")
  }

  # ... and two groups
  grp.id <- design_matrix[[1]]
  grps <- unique(grp.id)
  if (length(grps) != 2) {
    stop("Can only compare two groups at a time")
  }

  ### check whether the number of samples are equal to the number of rows in design matrix

  if (!all(sapply(list(length(geData),length(ctData)), function(x) x == nrow(design_matrix) ))) {
    stop("geData, ctData, design_matrix indicate different number of samples")
  }

  samples <- 1:length(geData)

  ### check whether the number of cells are equal in each sample

  for (i in samples) {
    if (ncol(geData[[i]])-1 != nrow(ctData[[i]])) {
      stop(paste0("geData and ctData indicate different cell number in sample ",i))
    }
  }

  # check col names for ctData and lr
  for (i in samples) {
    if (!setequal(colnames(ctData[[i]]),c("cell", "celltype"))) {
      stop("colnames for ctData should be c(\"cell\", \"celltype\")")
    }
  }

  if (!setequal(colnames(lr),c("ligand", "receptor"))) {
    stop("colnames for lr should be c(\"ligand\", \"receptor\")")
  }

  ### rename the first column of geData
  for (i in samples) {
    colnames(geData[[i]])[1] <- "gene"
  }

  cat("Done. \n")

  ### compute cellular detection rate ###

  # cdr_list <- vector("list", length = length(samples))
  # for (i in samples) {
  #   cell.detected.counts <- summarise_all(geData[[i]], ~(if(is.numeric(.)) sum(.!=0) else "count"))
  #
  #   cdr <- as.numeric(cell.detected.counts[,-1])/nrow(geData[[i]])
  #   cell.ids <- names(cell.detected.counts)[-1]
  #
  #   cdr_list[[i]] <- data.frame(cell.id=cell.ids, cdr=cdr)
  #
  # }

  ### separate the gene expression data of the cell types of interest
  ### only keep the ligand/receptor genes

  geData.L <- vector("list", length(samples))
  geData.R <- vector("list", length(samples))
  for (i in samples) {
    cell.L <- ctData[[i]][ctData[[i]]$celltype==celltype.L, "cell"]

    cols <- c("gene", cell.L)
    geData.L[[i]] <- geData[[i]][,cols] %>%
      filter(gene %in% lr$ligand)

    cell.R <- ctData[[i]][ctData[[i]]$celltype==celltype.R, "cell"]

    cols <- c("gene", cell.R)
    geData.R[[i]] <- geData[[i]][,cols] %>%
      filter(gene %in% lr$receptor)
  }

  rm(geData, ctData)

  ### compute expression rate ###

  # prop_list.L <- vector("list", length = length(samples))
  # prop_list.R <- vector("list", length = length(samples))
  # for (i in samples) {
  #   gene.detected.counts <- summarise_all(as.data.frame(t(geData.L[[i]][,-1])), ~(if(is.numeric(.)) sum(.!=0) else "count"))
  #
  #   prop.g <- as.numeric(gene.detected.counts[,-1])/(ncol(geData.L[[i]])-1)
  #   gene.names <- names(gene.detected.counts)[-1]
  #
  #   prop_list.L[[i]] <- data.frame(gene=gene.names, prop=prop.g)
  #
  #   gene.detected.counts <- summarise_all(as.data.frame(t(geData.R[[i]][,-1])), ~(if(is.numeric(.)) sum(.!=0) else "count"))
  #
  #   prop.g <- as.numeric(gene.detected.counts[,-1])/(ncol(geData.R[[i]])-1)
  #   gene.names <- names(gene.detected.counts)[-1]
  #
  #   prop_list.R[[i]] <- data.frame(gene=gene.names, prop=prop.g)
  #
  # }

  # make sure the order of the columns is ligand and then receptor.
  lr <- lr %>% select(ligand, receptor)
  ### find all ligand-receptor pairs that exist in all samples
  # cat("Finding L-R existed in all samples ... \n")

  # high.lr <- t(apply(lr, 1, find_lr, geData.l=geData.L, geData.r=geData.R))
  #
  # high.lr <- as.data.frame(high.lr)
  # colnames(high.lr) <- c("ligand","receptor")
  #
  # # drop NAs
  # high.lr <- na.omit(high.lr)

  ### find all ligand-receptor pairs that exist in at least one sample
  # "exist" means both ligand and receptor are expressed > min.prop
  cat("Finding L-R existed in at least one sample ... \n")
  existing.lr <- t(apply(lr, 1, find_lr0, geData.l=geData.L, geData.r=geData.R, min.prop = min.prop))

  existing.lr <- as.data.frame(existing.lr)
  colnames(existing.lr) <- c("ligand","receptor")

  # drop NAs
  existing.lr <- na.omit(existing.lr)
  ########### we are here
  cat("Done. \n")

  result_list <- vector("list", length = nrow(existing.lr))

  ### fit the model ###

  cat("Fitting the models ... \n")
  for (j in 1:nrow(existing.lr)) {
    progress(j, max.value = nrow(existing.lr))

    lig <- existing.lr$ligand[j]
    rec <- existing.lr$receptor[j]

    # create a df to store the statistics
    data.sum <- data.frame(sample.id=samples,
                           count.cL=NA, count.cR=NA, # no of cells
                           prop.L=NA, prop.R=NA, # expression rates
                           mean.pos.L=NA, mean.pos.R=NA, # pos mean expression
                           mean.L=NA, mean.R=NA # overall mean expression
                           )

    for (i in samples) {
      expr.l <- t(geData.L[[i]][geData.L[[i]]$gene==lig, -1])
      expr.r <- t(geData.R[[i]][geData.R[[i]]$gene==rec, -1])

      # no of cells
      data.sum$count.cL[i] <- ncol(geData.L[[i]])-1
      data.sum$count.cR[i] <- ncol(geData.R[[i]])-1

      # just in case the L/R gene doesn't exist in some samples
      if (length(expr.l)==0) {
        expr.l <- rep(0,data.sum$count.cL[i])
      } else {
        expr.l <- c(expr.l)
      }
      if (length(expr.r)==0) {
        expr.r <- rep(0,data.sum$count.cR[i])
      } else {
        expr.r <- c(expr.r)
      }

      # prop
      data.sum$prop.L[i] <- sum(expr.l>0)/length(expr.l)
      data.sum$prop.R[i] <- sum(expr.r>0)/length(expr.r)

      # positive mean
      data.sum$mean.pos.L[i] <- mean(expr.l[expr.l>0])
      data.sum$mean.pos.R[i] <- mean(expr.r[expr.r>0])

      # overall mean
      data.sum$mean.L[i] <- mean(expr.l)
      data.sum$mean.R[i] <- mean(expr.r)
    }
    # convert NaN to 0
    data.sum[is.na(data.sum)] <- 0

    ### check if psuedo-observations are needed
    add0.l <- F
    add1.l <- F
    add0.r <- F
    add1.r <- F

    # ligand:
    props.l <- data.sum$prop.L
    # categorize the elements in props.l aka which group do they belong to?
    A <- all(props.l[grp.id == grps[1]] == 0)
    B <- all(props.l[grp.id == grps[2]] == 0)
    C <- all(props.l[grp.id == grps[1]] == 1)
    D <- all(props.l[grp.id == grps[2]] == 1)
    if (A | B) add1.l <- T # if one group only contains 0, then we should add 1 to this group
    if (C | D) add0.l <- T

    # receptor:
    props.r <- data.sum$prop.R

    A <- all(props.r[grp.id == grps[1]] == 0)
    B <- all(props.r[grp.id == grps[2]] == 0)
    C <- all(props.r[grp.id == grps[1]] == 1)
    D <- all(props.r[grp.id == grps[2]] == 1)
    if (A | B) add1.r <- T
    if (C | D) add0.r <- T

    ### fit the model and conduct tests

    result_sublist <- diffLR.meta(lig = lig, rec = rec, geData.L = geData.L, geData.R = geData.R, design_matrix = design_matrix,
                re.cont = re.cont, re.disc = re.disc,
                add0.l = add0.l, add1.l = add1.l, add0.r = add0.r, add1.r = add1.r)

    result_sublist[['interaction']] <- paste(lig, rec, sep = "_")

    result_list[[j]] <- result_sublist

  }
  cat("\n Done. \n")

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
  cat(paste0("Time collapsed: ",round(difftime(Sys.time(), start.time, units = "mins"), digits = 1)," mins \n"))

  ### return a data frame
  return(data.result)

}






### helper functions ###

# find LR expressed in all samples

find_lr <- function(row_lr, geData.l, geData.r) {
  if ( all(sapply(geData.l, g_sample, g=row_lr[1])) & all(sapply(geData.r, g_sample, g=row_lr[2])) ) {
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

#

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




