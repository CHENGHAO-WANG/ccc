#' 
#' 
#' 
#' 
#' @param multi_sub it 
#'  \itemize{
#'    \item
#'    \item
#'  }
#'  
#'  @param cell_type_padj adjust p-values for each sender-receiver pair or not

library(data.table)
library(lme4)
library(sandwich)
library(clubSandwich)


ccc_analysis <- function(expression_matrix, metadata,
                         cell_id_col = "cell_id", cell_type_col = "cell_type", group_col = "group", covar_col = NULL, cdr = TRUE,
                         id_col = NULL, lmm_re = TRUE, logmm_re = TRUE,
                         sender = NULL, receiver = NULL,
                         lr = c("omnipathr","ramilowski"),
                         multi_sub = c("minimum","arithmetic_mean","geometric_mean","min_avg_gene","min_rate_gene"),
                         contrast, sandwich = FALSE,
                         verbose = TRUE,
                         min_pct = 0.01, large_n = 2, min_avg_pct = 0,
                         min_cell = 10,
                         threshold = 0, sep_prop = 0, sep_n = 0, sep_detection = TRUE,
                         padj_method = "BH", cell_type_padj = TRUE,
                         control_logm = list(),
                         control_lmm = lme4::lmerControl() , control_logmm = list(),
                         chunk_size = 10,
                         nthreads = 1
) {
  old_nthreads <- getDTthreads()
  setDTthreads(nthreads)
  
  if (is.vector(contrast)) {
    contrast <- matrix(contrast, nrow = 1L, dimnames = list(NULL, names(v)))
  } else if (!is.matrix(contrast)) {
    stop("'contrast' must be either a named vector or a matrix with column names.")
  }
  if (is.null(colnames(contrast))) {
    stop("'contrast' must be either a named vector or a matrix with column names.")
  }
  if (qr(contrast)$rank < nrow(contrast)) {
    stop("'contrast' must be full row rank")
  }
  
  contrast <- contrast[, colSums(abs(contrast) > 0)]
  
  if (isFALSE(all(colnames(contrast) %in% unique(metadata[[group_col]])))) {
    stop(paste0("'contrast' contains group levels that are not present in the \"",group_col,"\" column of 'metadata'."))
  }
  
  multi_sub <- match.arg(multi_sub, choices = c("minimum","arithmetic_mean","geometric_mean","min_avg_gene","min_rate_gene"))
  
  err_msg <- " in 'covar_col'. Please use a different argument or rename it."
  if ("id" %in% covar_col) {
    stop(paste0("\"id\"",err_msg))
  }
  if ("group" %in% covar_col) {
    stop(paste0("\"group\"",err_msg))
  }
  if ("cdr" %in% covar_col) {
    stop(paste0("\"cdr\"",err_msg))
  }
  
  if(any(is.na(expression_matrix))) {
    stop("Missing values not allowed in expression matrix, please remove NA values.")
  }
  if(any(is.na(metadata))) {
    stop("Missing values not allowed in metadata, please remove NA values.")
  }
  if(any(is.na(covar_col))) {
    stop("Missing values not allowed in covar_col, please remove NA values.")
  }
  
  if ("y" %in% colnames(metadata)) {
    stop("\"y\" in column names of 'metadata'. Please rename it.")
  }
  if ("z" %in% colnames(metadata)) {
    stop("\"z\" in column names of 'metadata'. Please rename it.")
  }
  
  covar_exists <- covar_col %in% colnames(metadata)
  if(!is.null(covar_col) && all(covar_exists)){
    stop(paste0("The following columns specified in 'covar_col' do not exist in 'metadata': ", paste(covar_col[!covar_exists], collapse = ", ")))
  }
  if (isFALSE(group_col %in% colnames(metadata))) {
    stop(paste0("'group_col' \"",group_col, "\" does not exist in 'metadata'."))
  }
  if (isFALSE(cell_id_col %in% colnames(metadata))) {
    stop(paste0("'cell_id_col' \"",cell_id_col, "\" does not exist in 'metadata'."))
  }
  if (isFALSE(cell_type_col %in% colnames(metadata))) {
    stop(paste0("'cell_type_col' \"",cell_type_col, "\" does not exist in 'metadata'."))
  }
  
  if (isTRUE(lmm_re) || isTRUE(logmm_re)) {
    if (is.null(id_col)) {
      stop("'id_col' is not specified")
    }
    if (isFALSE(id_col %in% colnames(metadata))) {
      stop(paste0("'id_col' \"",id_col, "\" does not exist in 'metadata'."))
    }
  }
  if (isFALSE(lmm_re) && isFALSE(logmm_re) && !is.null(id_col)) {
    warning("'id_col' is not NULL. This input will be ignored, because 'lmm_re' and 'logmm_re' are FALSE")
    id_col <- NULL
  }
  
  padj_method <- match.arg(padj_method, p.adjust.methods)
  
  
  metadata <- rename_metadata(metadata = metadata,
                              cell_id_col = cell_id_col, cell_type_col = cell_type_col,
                              id_col = id_col, group_col = group_col)
  num_ids <- metadata[, uniqueN("id")]
  
  if (isTRUE(sep_detection)) {
    if (sep_prop < 0 || sep_prop > 1) {
      stop("'sep_prop' must be between 0 and 1 (inclusive).")
    }
    if (sep_n < 0 || sep_n > num_ids) {
      stop("'sep_n' must be between 0 and number of samples (inclusive).")
    }
  }
  
  
  
  
}

#'
#' @import data.table



rename_metadata <- function(metadata, cell_id_col, id_col, group_col, cell_type_col) {
  metadata <- as.data.table(metadata)
  if (is.null(id_col)) {
    setnames(metadata, old = c(cell_id_col, group_col, cell_type_col), new = c("cell_id", "group", "cell_type"))
    metadata[, id := group]
  } else {
    setnames(metadata, old = c(cell_id_col, id_col, group_col, cell_type_col), new = c("cell_id", "id", "group", "cell_type"))
  }
  return(metadata)
}



#' 
#' 
prep_lr <- function(lr) {
  if (is.data.frame(lr)) {
    stopifnot(colnames(lr) == c("ligand","receptor"))
    lr_name <- "user"
    lr_user <- lr
  } 
  
  lr_table <- switch(EXPR = lr_name,
                     "omnipathr" = {
                       data("ramilowski", envir = environment())
                       ramilowski
                       },
                     "ramilowski" = {
                       data("ramilowski", envir = environment())
                       ramilowski
                     },
                     "user" = lr_user,
                     stop("'lr' should be \"omnipathr\" or \"ramilowski\" or a data.frame of ligand-receptor pairs")
  )
  
  return(lr_table)
}



# intercell_network <- OmnipathR::intercell_network(ligand_receptor = TRUE, high_confidence = TRUE, simplify = TRUE)
# omnipathr <- data.frame(ligand=intercell_network$source_genesymbol,receptor=intercell_network$target_genesymbol)
# save(omnipathr, file = "omnipathr.rda")


# 
# oplan <- plan(multisession, workers = 2)
# on.exit(plan(oplan), add = TRUE)

#' 
#' @import data.table
#' 
filter_cell_type <- function(metadata, sender, receiver, min_cell, contrast) {
  
  if (is.null(sender)) {
    sender <- unique(metadata$cell_type)
  }
  if (is.null(receiver)) {
    receiver <- unique(metadata$cell_type)
  }
  
  metadata <- metadata[group %in% colnames(contrast)]
  
  # # Check if any rows remain after contrast filtering
  # if (nrow(metadata) == 0) {
  #   stop("No rows remain after applying contrast filtering.")
  # }
  
  # Find rows where cell_type appears in either sender OR receiver
  valid_cell_types <- union(sender, receiver)
  metadata_subset <- metadata[cell_type %in% valid_cell_types]
  
  # # Check if any rows remain
  # if (nrow(metadata_subset) == 0) {
  #   stop("No rows remain after filtering by sender or receiver cell types.")
  # }
  
  # Count occurrences of each id-cell_type combination
  counts <- metadata_subset[, .N, by = .(id, cell_type)]
  
  # Keep only cell_type where count >= min_cell for each id
  metadata_subset <- metadata_subset[counts[N >= min_cell], on = .(id, cell_type)]
  
  # # Check again if any rows remain after applying min_cell filter
  # if (nrow(metadata_subset) == 0) {
  #   stop("No rows remain after applying min_cell filter.")
  # }
  
  # Ensure each id has the same unique set of cell_type values
  valid_cell_types_by_id <- metadata_subset[, .(unique_cell_types = list(unique(cell_type))), by = id]
  
  # Find the intersection of cell_types across all ids
  common_cell_types <- Reduce(intersect, valid_cell_types_by_id$unique_cell_types)
  
  # Filter metadata_subset so that only these common cell_types remain for each id
  metadata_subset <- metadata_subset[cell_type %in% common_cell_types]
  
  # Check if any rows remain after ensuring consistent cell_types across ids
  if (nrow(metadata_subset) == 0) {
    stop("No cell types remain after 'min_cell' filtering.")
  }
  
  # Check if sender or receiver contain elements not in the final subset
  remaining_cell_types <- unique(metadata_subset$cell_type)
  
  missing_sender <- setdiff(sender, remaining_cell_types)
  missing_receiver <- setdiff(receiver, remaining_cell_types)
  
  sender <- intersect(sender, remaining_cell_types)
  receiver <- intersect(receiver, remaining_cell_types)
  
  if (length(sender) == 0) {
    stop("No sender cell types remain after 'min_cell' filtering.")
  }
  if (length(receiver) == 0) {
    stop("No receiver cell types remain after 'min_cell' filtering.")
  }
  
  if (length(missing_sender) > 0 || length(missing_receiver) > 0) {
    warning_msg <- "Some cell types in 'sender' or 'receiver' do not appear in the final subset."
    if (length(missing_sender) > 0) {
      warning_msg <- paste0(warning_msg, "\nMissing in sender: ", paste(missing_sender, collapse = ", "))
    }
    if (length(missing_receiver) > 0) {
      warning_msg <- paste0(warning_msg, "\nMissing in receiver: ", paste(missing_receiver, collapse = ", "))
    }
    warning(warning_msg)
  }
  
  return(list(metadata_subset = metadata_subset, sender = sender, receiver = receiver))
}

#' 
#' @import data.table
#' 

compute_cdr <- function(expression_matrix, metadata_subset, threshold) {
  # Ensure rownames and colnames are correctly set
  if (is.null(rownames(expression_matrix)) || is.null(colnames(expression_matrix))) {
    stop("expression_matrix must have rownames (genes) and colnames (cell ids).")
  }
  
  # Subset the expression matrix to only include the cell ids in metadata_subset
  common_cells <- intersect(colnames(expression_matrix), metadata_subset$cell_id)
  
  if (length(common_cells) == 0) {
    stop("No overlapping cell ids found between 'expression_matrix' and 'metadata' given 'sender' and 'receiver'.")
  }
  
  expression_subset <- expression_matrix[, common_cells, drop = FALSE]
  
  # Calculate cdr: fraction of genes with expression > threshold for each cell
  cdr_values <- colMeans(expression_subset > threshold)
  
  # Merge results with metadata_subset
  metadata_subset[, cdr := cdr_values[match(cell_id, names(cdr_values))]]
  
  return(metadata_subset)
}




run_analysis2 <- function(sender, receiver, lr_table, ){
  # Create all possible combinations of sender and receiver
  sender_receiver_combinations <- expand.grid(sender = sender, receiver = receiver)
  
  # Merge with lr_table to get all possible combinations
  pairs4analysis <- merge(sender_receiver_combinations, lr_table, by = NULL)
  
  # Convert to data.table
  setDT(pairs4analysis)
  npairs <- nrow(pairs4analysis)
  
  future_lapply(seq(1L, nrow(pairs4analysis), by = chunk_size), FUN = run_analysis(i))
  
  unique_ids <- unique(metadata_subset[,id])
  
  run_analysis <- function(i) {
    chunk <- pairs4analysis[i:min(i + chunk_size - 1L, npairs), ]
    results.i <- list()
    for (j in 1L:nrow(chunk)) {
      sender1 <- chunk$sender[i]
      ligand <- chunk$ligand[i]
      receiver1 <- chunk$receiver[i]
      receptor <- chunk$receptor[i]
      
      # Create copies of metadata_subset for sender and receiver
      data_sender_ligand <- metadata_subset[cell_type == sender]
      data_receiver_receptor <- metadata_subset[cell_type == receiver]
      
      # Handle single gene or multi-gene ligands and receptors
      ligand_genes <- unlist(strsplit(ligand, "_"))
      receptor_genes <- unlist(strsplit(receptor, "_"))
      
      # Subset expression matrix for ligands and receptors
      ligand_expr_values <- expression_matrix[ligand_genes, data_sender_ligand$cell_id, drop = TRUE]
      receptor_expr_values <- expression_matrix[receptor_genes, data_receiver_receptor$cell_id, drop = TRUE]
      
      # Function to compute expression value based on method
      compute_expression_value <- function(expr_values, multi_sub) {
        if (is.null(dim(expr_values))) return(expr_values)
        
        switch(
          multi_sub,
          "minimum" = apply(expr_values, 2L, min),
          "arithmetic_mean" = colMeans(expr_values),
          "geometric_mean" = apply(expr_values, 2L, FUN = function(x) prod(x)^(1 / length(x))),
          "min_avg_gene" = {
            gene_means <- rowMeans(expr_values)
            min_gene <- names(which.min(gene_means))
            expr_values[min_gene, , drop = TRUE]
          }
          "min_rate_gene" = {
            gene_rates <- rowMeans(expr_values > threshold)
            min_gene <- names(which.min(gene_rates))
            expr_values[min_gene, , drop = TRUE]
          }
        )
      }
      
      # Add expression column to metadata copies
      data_sender_ligand[, y := compute_expression_value(ligand_expr_values, multi_sub)]
      data_receiver_receptor[, y := compute_expression_value(receptor_expr_values, multi_sub)]
      
      # Compute expression rates
      compute_expression_rate <- function(data_subset) {
        sapply(metadata_subset$id, function(uid) {
          mean(data_subset[id == uid, y] > threshold)
        })
      }
      
      ligand_expression_rates <- compute_expression_rate(data_sender_ligand)
      receptor_expression_rates <- compute_expression_rate(data_receiver_receptor)
      
      # Check if the number of ids with both ligand and receptor expression rate >= min_pct is >= large_n
      valid_ids <- sum((ligand_expression_rates >= min_pct) & (receptor_expression_rates >= min_pct))
      if (valid_ids < large_n) {
        next
      } 

      ##
      # Add indicator column
      data_sender_ligand[, z := ifelse(y > threshold, 1, 0)]
      data_receiver_receptor[, z := ifelse(y > threshold, 1, 0)]
      
      # Subset data
      data_sender_ligand_1 <- data_sender_ligand[z == 1]
      data_receiver_receptor_1 <- data_receiver_receptor[z == 1]
      
      # Define covariates
      covar <- c(covar, if(isTRUE(cdr)) "cdr")

      center_covar(dt = data_sender_ligand_1, covar = covar) -> covariates # The 4 data sets are different. But the column names are the same.
      center_covar(dt = data_receiver_receptor_1, covar = covar)
      center_covar(dt = data_sender_ligand, covar = covar)
      center_covar(dt = data_receiver_receptor, covar = covar)
      
      covariates <- c("group", covar)
      
      # Define model formulas
      fixed_effects <- paste(covariates, collapse = " + ")
      formula_linear <- as.formula(paste("y ~ 0 +", fixed_effects))
      formula_logistic <- as.formula(paste("z ~ 0 +", fixed_effects))
      if (lmm_re) {
        formula_linear <- as.formula(paste("y ~ 0 +", fixed_effects, "+ (1|id)"))
      }
      if (logmm_re) {
        fixed_formula <- as.formula(paste("z ~ 0 +", fixed_effects))
        random_formula <- as.formula("~ 1 | id")
        formula_logistic <- list("fixed" = fixed_formula, "random" = random_formula)
      }
      
      
        
      # Fit models
      fit_linear <- function(data, formula) {
        tryCatch({
          cond <- detect_all_zeros(dt = data, id_col = "id", id = unique_ids)
          if (cond) {
            stop("Too few cells expressing the ligand/receptor gene for fitting a linear model.")
          } else {
            if (isTRUE(lmm_re)) {
              lmer(formula)
            } else {
              lm(formula)
            }
          }
        }, error = function(e) {
          return(error = e$message)
        })
      
      }
      fit_logistic <- function(data, formula) {
        tryCatch({
          cond <- isTRUE(sep_detection) && detect_re_separation(dt = data, z_col = "z", id_col = "id", num_ids = num_ids, sep_prop = sep_prop, sep_n = sep_n)
          if (cond) {
            stop("Complete or Quasi-complete separation detected.")
          } else {
            if (isTRUE(logmm_re)) {
              mixed_model(fixed = formula$fixed, random = formula$random, family = binomial())
            } else {
              glm(formula, family = binomial())
            }
          }
        }, error = function(e) {
          return(error = e$message)
        })
      }

      ##
      fit.l.linear <- fit_linear(data = data_sender_ligand_1, formula = formula_linear)
      fit.r.linear <- fit_linear(data = data_receiver_receptor_1, formula = formula_linear)
      fit.l.logistic <- fit_logistic(data = data_sender_ligand, formula = formula_logistic)
      fit.r.logistic <- fit_logistic(data = data_receiver_receptor, formula = formula_logistic)
      
      results.i[[length(results.i) + 1L]] <- ccc_test(fit.l.linear = fit.l.linear, fit.l.logistic = fit.l.logistic,
                                                      fit.r.linear = fit.r.linear, fit.r.logistic = fit.r.logistic,
                                                      contrast = contrast, re_lmm = re_lmm, re_logmm = re_logmm,
                                                      sandwich = sandwich)
      
    }
    rbindlist(results.i, fill = TRUE)
  }
 
  
  
}

ccc_test <- function(fit.l.linear, fit.l.logistic, fit.r.linear, fit.r.logistic,
                                        contrast, re_lmm, re_logmm, sandwich) {
  
  group_names <- colnames(contrast)
  
  if (!is.character(fit.l.linear) && !is.character(fit.r.linear)) {
    if (isTRUE(re_lmm)) {
      coef_l_lm <- fixef(fit.l.linear)
      coef_r_lm <- fixef(fit.r.linear)
      if (isTRUE(sandwich)) {
        vcov_l_lm_group <- vcovCR(fit.l.linear, type = "CR2")[group_names, group_names]
        vcov_r_lm_group <- vcovCR(fit.r.linear, type = "CR2")[group_names, group_names]
      }
    } else {
      coef_l_lm <- stats::coef(fit.l.linear)
      coef_r_lm <- stats::coef(fit.r.linear)
      if (isTRUE(sandwich)) {
        vcov_l_lm_group <- vcovHC(fit.l.linear, type = "HC3")[group_names, group_names]
        vcov_r_lm_group <- vcovHC(fit.r.linear, type = "HC3")[group_names, group_names]
      }
    }
    vcov_l_lm_group <- vcov(fit.l.linear)[group_names, group_names]
    vcov_r_lm_group <- vcov(fit.r.linear)[group_names, group_names]
    test.linear <- TRUE
  } else {
    test.linear <- FALSE
  }
  if (!is.character(fit.l.logistic) && !is.character(fit.r.logistic)) {
    if (isTRUE(re_logmm)) {
      m_l <- marginal_coefs(fit.l.logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      m_r <- marginal_coefs(fit.r.logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      coef_l_logm <- m_l$betas
      coef_r_logm <- m_l$betas
      vcov_l_logm_group <- m_l$var_betas[group_names, group_names] 
      vcov_r_logm_group <- m_r$var_betas[group_names, group_names]
    } else {
      coef_l_logm <- stats::coef(fit.l.logistic)
      coef_r_logm <- stats::coef(fit.r.logistic)
      if (isTRUE(sandwich)) {
        vcov_l_logm_group <- vcovHC(fit.l.logistic, type = "HC3")[group_names, group_names]
        vcov_r_logm_group <- vcovHC(fit.r.logistic, type = "HC3")[group_names, group_names]
      } else {
        vcov_l_logm_group <- vcov(fit.l.logistic)[group_names, group_names]
        vcov_r_logm_group <- vcov(fit.r.logistic)[group_names, group_names]
      }
    }
    test.logistic <- TRUE
  } else {
    test.logistic <- FALSE
  }

  # coef_l_lmm <- fixef(fit.l.lmm)
  # coef_l_logmm <- fixef(fit.l.logmm)
  # coef_r_lmm <- fixef(fit.r.lmm)
  # coef_r_logmm <- fixef(fit.r.logmm)
  
  if (isTRUE(test.linear)) {
    coef_l_lm <- coef_l_lm[names(coef_l_lm) %in% group_names]
    coef_r_lm <- coef_r_lm[names(coef_r_lm) %in% group_names]
    coef_l_lm <- coef_l_lm[group_names]
    coef_r_lm <- coef_r_lm[group_names]
  }
  if (isTRUE(test.logistic)) {
    coef_l_logm <- coef_l_logm[names(coef_l_logm) %in% group_names]
    coef_r_logm <- coef_r_logm[names(coef_r_logm) %in% group_names]
    coef_l_logm <- coef_l_logm[group_names]
    coef_r_logm <- coef_r_logm[group_names]
    
  }
  
  # Filter coefficients to only include group names
  # coef_l_lm <- coef_l_lm[names(coef_l_lm) %in% group_names]
  # coef_l_logm <- coef_l_logm[names(coef_l_logm) %in% group_names]
  # coef_r_lm <- coef_r_lm[names(coef_r_lm) %in% group_names]
  # coef_r_logm <- coef_r_logm[names(coef_r_logm) %in% group_names]
  
  # Ensure coefficients are named vectors
  # if (is.null(names(coef_l_lmm))) names(coef_l_lmm) <- names(fixef(fit.l.lmm)[names(fixef(fit.l.lmm)) %in% group_names])
  # if (is.null(names(coef_l_logmm))) names(coef_l_logmm) <- names(fixef(fit.l.logmm)[names(fixef(fit.l.logmm)) %in% group_names])
  # if (is.null(names(coef_r_lmm))) names(coef_r_lmm) <- names(fixef(fit.r.lmm)[names(fixef(fit.r.lmm)) %in% group_names])
  # if (is.null(names(coef_r_logmm))) names(coef_r_logmm) <- names(fixef(fit.r.logmm)[names(fixef(fit.r.logmm)) %in% group_names])
  
  # Reorder coefficients to match group_names order
  # coef_l_lm <- coef_l_lm[group_names]
  # coef_l_logm <- coef_l_logm[group_names]
  # coef_r_lm <- coef_r_lm[group_names]
  # coef_r_logm <- coef_r_logm[group_names]
  
  if (isTRUE(test.linear)) {
    product_vector_linear <- numeric(length(group_names))
    names(product_vector_linear) <- group_names
    product_vector_linear <- coef_l_lm * coef_r_lm
  }
  
  product_vector_hurdle <- product_vector_linear <- product_vector_logisitc <- numeric(length(group_names))
  names(product_vector_hurdle) <- names(product_vector_linear) <- names(product_vector_logisitc) <- group_names
  
  for (group in group_names) {
    product_vector_hurdle[group] <- coef_l_lm[group] * coef_r_lm[group] * plogis(coef_l_logm[group]) * plogis(coef_r_logm[group])
    product_vector_linear[group] <- coef_l_lm[group] * coef_r_lm[group]
    product_vector_logisitc[group] <- plogis(coef_l_logm[group]) * plogis(coef_r_logm[group])
  }
  
  effect_size_hurdle <- contrast %*% product_vector_hurdle
  effect_size_linear <- contrast %*% product_vector_linear
  effect_size_logistic <- contrast %*% product_vector_logisitc
  
  # Delta method for mean and covariance
  gradient_matrix_hurdle <- matrix(0, nrow = nrow(contrast), ncol = length(group_names) * 4L)
  
  for (i in seq_len(nrow(contrast))) {
    for (j in seq_along(group_names)) {
      group <- group_names[j]
      
      ind_l_lm <- which(names(coef_l_lm) == group)
      ind_l_logm <- which(names(coef_l_logm) == group) + length(group_names)
      ind_r_lm <- which(names(coef_r_lm) == group) + length(group_names) * 2L
      ind_r_logm <- which(names(coef_r_logm) == group) + length(group_names) * 3L
      
      gradient_matrix_hurdle[i, ind_l_lm] <- contrast[i, j] * coef_r_lm[group] * plogis(coef_l_logm[group]) * plogis(coef_r_logm[group])
      gradient_matrix_hurdle[i, ind_l_logm] <- contrast[i, j] * coef_l_lm[group] * coef_r_lm[group] * plogis(coef_r_logm[group]) * dlogis(coef_l_logm[group])
      gradient_matrix_hurdle[i, ind_r_lm] <- contrast[i, j] * coef_l_lm[group] * plogis(coef_l_logm[group]) * plogis(coef_r_logm[group])
      gradient_matrix_hurdle[i, ind_r_logm] <- contrast[i, j] * coef_l_lm[group] * coef_r_lm[group] * plogis(coef_l_logm[group]) * dlogis(coef_r_logm[group])
    }
  }
  
  # Extract relevant parts of vcov matrices
  # vcov_l_lmm_group <- vcov(fit.l.lmm)[group_names, group_names]
  # vcov_l_logmm_group <- vcov(fit.l.logmm)[group_names, group_names]
  # vcov_r_lmm_group <- vcov(fit.r.lmm)[group_names, group_names]
  # vcov_r_logmm_group <- vcov(fit.r.logmm)[group_names, group_names]
  
  vcov_hurdle <- bdiag(vcov_l_lm_group, vcov_l_logm_group, vcov_r_lm_group, vcov_r_logm_group)
  
  cov_effect_size_hurdle <- gradient_matrix_hurdle %*% as.matrix(vcov_hurdle) %*% t(gradient_matrix_hurdle)
  
  # Delta method | linear model
  gradient_matrix_linear <- matrix(0, nrow = nrow(contrast), ncol = length(group_names) * 2L)
  
  for (i in seq_len(nrow(contrast))) {
    for (j in seq_along(group_names)) {
      group <- group_names[j]
      
      ind_l_lm <- which(names(coef_l_lm) == group)
      ind_r_lm <- which(names(coef_r_lm) == group) + length(group_names)
      
      gradient_matrix_linear[i, ind_l_lm] <- contrast[i, j] * coef_r_lm[group] 
      gradient_matrix_linear[i, ind_r_lm] <- contrast[i, j] * coef_l_lm[group] 
    }
  }
  
  vcov_linear <- bdiag(vcov_l_lm_group, vcov_r_lm_group)
  
  cov_effect_size_linear <- gradient_matrix_linear %*% as.matrix(vcov_linear) %*% t(gradient_matrix_linear)
  
  # Delta method | logistic model
  gradient_matrix_logistic <- matrix(0, nrow = nrow(contrast), ncol = length(group_names) * 2L)
  
  for (i in seq_len(nrow(contrast))) {
    for (j in seq_along(group_names)) {
      group <- group_names[j]
      
      ind_l_logm <- which(names(coef_l_logm) == group)
      ind_r_logm <- which(names(coef_r_logm) == group) + length(group_names)
      
      gradient_matrix_logistic[i, ind_l_logm] <- contrast[i, j] * plogis(coef_r_logm[group]) * dlogis(coef_l_logm[group])
      gradient_matrix_logistic[i, ind_r_logm] <- contrast[i, j] * plogis(coef_l_logm[group]) * dlogis(coef_r_logm[group])
    }
  }
  
  vcov_logistic <- bdiag(vcov_l_logm_group, vcov_r_logm_group)
  
  cov_effect_size_logistic <- gradient_matrix_logistic %*% as.matrix(vcov_hurdle) %*% t(gradient_matrix_logistic)
  
  # test_statistic <- tryCatch({
  #   solve(cov_weighted_sum) %*% mean_weighted_sum
  # }, error = function(e) {
  #   MASS::ginv(cov_weighted_sum) %*% mean_weighted_sum
  # })
  
  ### testing
  
  test_stat_hurdle <- t(effect_size_hurdle) %*% solve(cov_effect_size_hurdle) %*% effect_size_hurdle
  p_value_hurdle <- pchisq(test_stat_hurdle, df = nrow(contrast), lower.tail = FALSE)
  test_stat_linear <- t(effect_size_linear) %*% solve(cov_effect_size_linear) %*% effect_size_linear
  p_value_linear <- pchisq(test_stat_linear, df = nrow(contrast), lower.tail = FALSE)
  test_stat_logistic <- t(effect_size_logistic) %*% solve(cov_effect_size_logistic) %*% effect_size_logistic
  p_value_logistic <- pchisq(test_stat_logistic, df = nrow(contrast), lower.tail = FALSE)
  test_stat_2part <- test_stat_linear + test_stat_logistic
  p_value_2part <- pchisq(test_stat_2part, df = nrow(contrast) * 2L, lower.tail = FALSE)
  
  # return(list(mean_weighted_sum = mean_weighted_sum, 
  #             test_statistic = test_statistic, 
  #             p_value = p_value))
  
  dt.result
  
  
}


library(future.apply)
plan(multisession)

library(progressr)
handlers(global = TRUE)
handlers("progress", "beepr")

my_fcn <- function(xs) {
  p <- progressor(along = xs)
  future_lapply(xs, function(x, ...) {
    Sys.sleep(6.0-x)
    p(sprintf("x=%g", x))
    sqrt(x)
  })
}

my_fcn(1:5)

library(data.table)

# Create a data.table with double vectors
dt <- data.table(
  id = 1:3, 
  values = list(c(1/3, 2/3), c(pi, exp(1)), c(sqrt(2), log(2)))
)

# Save as CSV
fwrite(dt, "output.csv")

# Read the CSV file
dt_read <- fread("output.csv")

# Check structure
str(dt_read)

dt_read[, values := lapply(values, function(x) as.numeric(strsplit(x, "\\|")[[1]]))]

# Check the updated structure
str(dt_read)

dt_read

dt <- data.table(
  id = 1:3, 
  values = list(c(1/3, 2/3), c(pi, exp(1)), c(sqrt(2), log(2), 5))
)
# Find the maximum length of vectors
max_len <- max(lengths(dt$values))

# Expand list-column into separate columns
dt_wide <- dt[, paste0("value_", seq_len(max_len)) := transpose(values)]

# Remove the original 'values' column
dt_wide[, values := NULL]

# Check the result
print(dt_wide)

library(lme4)
library(detectseparation)
# Generate example data
set.seed(123)
data <- data.frame(
  y = rbinom(100, 1, 0.5),    # Binary outcome (0 or 1)
  x = factor(rep(1:2, each = 50)),  # Categorical predictor
  group = factor(rep(1:4, each = 25))  # Random effect grouping (10 groups)
)

# Fit logistic mixed-effects model
model <- glmer(y ~ x + (1 | group), data = data, family = binomial)

# Print model summary
summary(model)

glm(y ~ x, data = data, family = binomial, method = 'detect_separation')



set.seed(123)

# Generate example data
data <- data.frame(
  y = rbinom(100, 1, 0.5),  # Binary outcome
  x = factor(rep(1:2, each = 50)),  # Categorical predictor
  group = factor(rep(1:10, each = 10))  # Random effect grouping (10 groups)
)

# Force three clusters (e.g., groups 1, 3, and 5) to have y = 0
data$y[data$group %in% c(1, 3, 5)] <- 0

# Fit logistic mixed-effects model
model <- glmer(y ~ x + (1 | group), data = data, family = binomial)

# Print model summary
summary(model)

library(data.table)



