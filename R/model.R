
# check the specification of arguments

stat_analysis <- function(existing.lr, geData.L, geData.R) {
  
  message("Model fitting and statistical testing ...")
  
  result_list <- vector("list", length = nrow(existing.lr) )
  
  message("Fitting the models ... \n")
  
  if (n.grp == 2) {
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
  }
}