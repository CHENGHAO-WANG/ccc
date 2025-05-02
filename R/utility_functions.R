utils::globalVariables(c(
  ".", "N", "V1", "cdr", "cell_id", "cell_type", "group", "id", "y", "z",
  "omnipathr", "ramilowski", "dispersion", "gene_id", "within_sample_correlation",
  "cell_continuous", "cell_discrete", "sample_continuous", "sample_discrete"
))

center_covar <- function(dt, covar) {
  # if (!is.data.table(dt)) stop("Input 'dt' must be a data.table.")

  # missing_covar <- covar[!covar %in% names(dt)]
  # if (length(missing_covar) > 0) {
  #   stop("The following columns do not exist in 'dt': ", paste(missing_covar, collapse = ", "))
  # }

  centered_column_names <- character(0)

  for (col_name in covar) {
    col <- dt[[col_name]]

    if (is.factor(col) || is.character(col)) {
      levels_col <- sort(unique(col)) # you could use levels(col) if it's a factor and you want to preserve order
      ref_level <- levels_col[1] # drop the first level
      dummy_levels <- levels_col[levels_col != ref_level]

      for (level in dummy_levels) {
        new_col_name <- paste0(col_name, "_", level)
        dt[, (new_col_name) := as.integer(col == level)]
        dt[, (new_col_name) := scale(.SD[[1]], center = TRUE, scale = FALSE), .SDcols = new_col_name]
        centered_column_names <- c(centered_column_names, new_col_name)
      }

      dt[, (col_name) := NULL] # remove original
    } else {
      dt[, (col_name) := scale(col, center = TRUE, scale = FALSE)]
      centered_column_names <- c(centered_column_names, col_name)
    }
  }

  return(centered_column_names)
}


detect_re_separation <- function(dt, z_col, id_col, num_ids = NULL, sep_prop = 0, sep_n = 0) {
  # Ensure dt is a data.table
  # setDT(dt)

  # Count the number of distinct id's
  if (is.null(num_ids)) {
    num_ids <- dt[, uniqueN(get(id_col))]
  }

  # Count the number of ids with all z values as either 0 or 1
  count_sep_ids <- dt[, all(get(z_col) == 0) || all(get(z_col) == 1), by = get(id_col)][, sum(V1)]

  # Check if the proportion exceeds the threshold
  detection <- count_sep_ids > sep_prop * num_ids && count_sep_ids > sep_n

  return(detection)
}


detect_all_zeros <- function(dt, id_col, id) {
  # setDT(dt)
  unique_id_col <- unique(dt[[id_col]]) # Get unique values from id_col
  !all(id %in% unique_id_col)
}

compute_group_stats <- function(dt, y_col = "y", group_col = "group", z_col = "z", prefix) {
  # Compute stats for groups that are present in the data
  stats_dt <- dt[, .(
    mean = mean(get(y_col)),
    positive_mean = if (any(get(z_col) == 1)) mean(get(y_col)[get(z_col) == 1]) else 0,
    expression_rate = mean(get(z_col) == 1)
  ), by = get(group_col)]

  setnames(stats_dt, "get", group_col)

  # # Ensure all groups from group_names are included
  # all_groups_dt <- data.table(group = group_names)
  # result <- merge(all_groups_dt, stats_dt, by = group_col, all.x = TRUE)
  #
  # # Replace NAs with 0 for missing groups
  # for (col in c("mean", "positive_mean", "expression_rate")) {
  #   result[is.na(get(col)), (col) := 0]
  # }
  #
  # return(result)
  stats_cols <- setdiff(names(stats_dt), group_col)
  setnames(stats_dt, stats_cols, paste0(prefix, stats_cols))
  stats_dt
}



#'
#'
prep_lr <- function(lr) {
  if (is.data.frame(lr)) {
    stopifnot(colnames(lr) == c("ligand", "receptor"))
    lr_name <- "user"
    lr_user <- lr
  }

  lr_table <- switch(EXPR = lr_name,
    "omnipathr" = {
      data("omnipathr", envir = environment())
      omnipathr
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



rename_metadata <- function(metadata, cell_id_col, id_col, group_col, cell_type_col) {
  # metadata0 <- as.data.table(metadata)
  # message("Current column names: ", paste(names(metadata), collapse = ", "))
  if (is.null(id_col)) {
    setnames(metadata, old = c(cell_id_col, group_col, cell_type_col), new = c("cell_id", "group", "cell_type"))
    metadata[, id := group]
  } else {
    data.table::setnames(metadata, old = c(cell_id_col, id_col, group_col, cell_type_col), new = c("cell_id", "id", "group", "cell_type")) # , skip_absent = T)
  }
  return(metadata)
}




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
      warning_msg <- paste0(warning_msg, "Missing in sender: ", paste(missing_sender, collapse = ", "))
    }
    if (length(missing_receiver) > 0) {
      warning_msg <- paste0(warning_msg, "Missing in receiver: ", paste(missing_receiver, collapse = ", "))
    }
    warning(warning_msg)
  }

  return(list(metadata_subset = metadata_subset, sender = sender, receiver = receiver))
}


compute_cdr <- function(expression_matrix, metadata_subset, threshold) {
  # # Ensure rownames and colnames are correctly set
  # if (is.null(rownames(expression_matrix)) || is.null(colnames(expression_matrix))) {
  #   stop("expression_matrix must have rownames (genes) and colnames (cell ids).")
  # }

  # Subset the expression matrix to only include the cell ids in metadata_subset
  common_cells <- intersect(colnames(expression_matrix), metadata_subset$cell_id)

  # if (length(common_cells) == 0) {
  #   stop("No overlapping cell ids found between 'expression_matrix' and 'metadata' given 'sender' and 'receiver'.")
  # }

  expression_subset <- expression_matrix[, common_cells, drop = FALSE]

  # Calculate cdr: fraction of genes with expression > threshold for each cell
  cdr_values <- Matrix::colMeans(expression_subset > threshold)

  # Merge results with metadata_subset
  metadata_subset[, cdr := cdr_values[match(cell_id, names(cdr_values))]]

  return(metadata_subset)
}

# Function to compute expression value based on method
compute_expression_value <- function(expr_values, multi_sub, threshold) {
  if (is.null(dim(expr_values))) {
    return(expr_values)
  }

  switch(multi_sub,
    "minimum" = apply(expr_values, 2L, min),
    "arithmetic_mean" = Matrix::colMeans(expr_values),
    "geometric_mean" = apply(expr_values, 2L, FUN = function(x) prod(x)^(1 / length(x))),
    "min_avg_gene" = {
      gene_means <- rowMeans(expr_values)
      min_gene <- names(which.min(gene_means))
      expr_values[min_gene, , drop = TRUE]
    },
    "min_rate_gene" = {
      gene_rates <- rowMeans(expr_values > threshold)
      min_gene <- names(which.min(gene_rates))
      expr_values[min_gene, , drop = TRUE]
    }
  )
}


ccc_test <- function(fit.l.linear, fit.l.logistic, fit.r.linear, fit.r.logistic,
                     contrast, lmm_re, logmm_re, sandwich,
                     sender, receiver, ligand, receptor) {
  dt.test <- data.table()
  group_names <- paste0("group", colnames(contrast))
  if (!is.null(fit.l.linear) && !is.null(fit.r.linear)) {
    if (isTRUE(lmm_re)) {
      coef_l_lm <- lme4::fixef(fit.l.linear)
      coef_r_lm <- lme4::fixef(fit.r.linear)
      if (isTRUE(sandwich)) {
        vcov_l_lm_group <- clubSandwich::vcovCR(fit.l.linear, type = "CR2")[group_names, group_names]
        vcov_r_lm_group <- clubSandwich::vcovCR(fit.r.linear, type = "CR2")[group_names, group_names]
      }
    } else {
      coef_l_lm <- stats::coef(fit.l.linear)
      coef_r_lm <- stats::coef(fit.r.linear)
      if (isTRUE(sandwich)) {
        vcov_l_lm_group <- sandwich::vcovHC(fit.l.linear, type = "HC3")[group_names, group_names]
        vcov_r_lm_group <- sandwich::vcovHC(fit.r.linear, type = "HC3")[group_names, group_names]
      }
    }
    if (isFALSE(sandwich)) {
      vcov_l_lm_group <- vcov(fit.l.linear)[group_names, group_names]
      vcov_r_lm_group <- vcov(fit.r.linear)[group_names, group_names]
    }
    test.linear <- TRUE
  } else {
    test.linear <- FALSE
  }

  if (!is.null(fit.l.logistic) && !is.null(fit.r.logistic)) {
    if (isTRUE(logmm_re)) {
      m_l <- GLMMadaptive::marginal_coefs(fit.l.logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      m_r <- GLMMadaptive::marginal_coefs(fit.r.logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      coef_l_logm <- m_l$betas
      coef_r_logm <- m_l$betas
      vcov_l_logm_group <- m_l$var_betas[group_names, group_names]
      vcov_r_logm_group <- m_r$var_betas[group_names, group_names]
    } else {
      coef_l_logm <- stats::coef(fit.l.logistic)
      coef_r_logm <- stats::coef(fit.r.logistic)
      if (isTRUE(sandwich)) {
        vcov_l_logm_group <- sandwich::vcovHC(fit.l.logistic, type = "HC3")[group_names, group_names]
        vcov_r_logm_group <- sandwich::vcovHC(fit.r.logistic, type = "HC3")[group_names, group_names]
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
  if (test.linear || test.logistic) {
    dt.test[, c("sender", "receiver", "ligand", "receptor") := list(sender, receiver, ligand, receptor)]
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
    effect_size_linear <- contrast %*% (coef_l_lm * coef_r_lm)
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

    vcov_linear <- Matrix::bdiag(vcov_l_lm_group, vcov_r_lm_group)
    cov_effect_size_linear <- gradient_matrix_linear %*% as.matrix(vcov_linear) %*% t(gradient_matrix_linear)

    test_stat_linear <- t(effect_size_linear) %*% chol2inv(chol(cov_effect_size_linear)) %*% effect_size_linear
    pvalue_linear <- pchisq(test_stat_linear, df = nrow(contrast), lower.tail = FALSE)

    dt.test[, c("effect_size_linear", "pvalue_linear") := list(list(effect_size_linear), pvalue_linear)]
    if (nrow(contrast) > 1L) {
      dt.test[, paste0("effect_size_linear_", seq_len(nrow(contrast))) := transpose(effect_size_linear)]
      dt.test[, effect_size_linear := NULL]
    }
  }

  if (isTRUE(test.logistic)) {
    effect_size_logistic <- contrast %*% (plogis(coef_l_logm) * plogis(coef_r_logm))
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

    vcov_logistic <- Matrix::bdiag(vcov_l_logm_group, vcov_r_logm_group)
    cov_effect_size_logistic <- gradient_matrix_logistic %*% as.matrix(vcov_logistic) %*% t(gradient_matrix_logistic)

    test_stat_logistic <- t(effect_size_logistic) %*% chol2inv(chol(cov_effect_size_logistic)) %*% effect_size_logistic
    pvalue_logistic <- pchisq(test_stat_logistic, df = nrow(contrast), lower.tail = FALSE)

    dt.test[, c("effect_size_logistic", "pvalue_logistic") := list(list(effect_size_logistic), pvalue_logistic)]
    if (nrow(contrast) > 1L) {
      dt.test[, paste0("effect_size_logistic_", seq_len(nrow(contrast))) := transpose(effect_size_logistic)]
      dt.test[, effect_size_logistic := NULL]
    }
  }

  if (isTRUE(test.linear) && isTRUE(test.logistic)) {
    effect_size_hurdle <- contrast %*% (coef_l_lm * coef_r_lm * plogis(coef_l_logm) * plogis(coef_r_logm))
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
    vcov_hurdle <- Matrix::bdiag(vcov_l_lm_group, vcov_l_logm_group, vcov_r_lm_group, vcov_r_logm_group)
    cov_effect_size_hurdle <- gradient_matrix_hurdle %*% as.matrix(vcov_hurdle) %*% t(gradient_matrix_hurdle)

    test_stat_hurdle <- t(effect_size_hurdle) %*% chol2inv(chol(cov_effect_size_hurdle)) %*% effect_size_hurdle
    pvalue_hurdle <- pchisq(test_stat_hurdle, df = nrow(contrast), lower.tail = FALSE)
    test_stat_2part <- test_stat_linear + test_stat_logistic
    pvalue_2part <- pchisq(test_stat_2part, df = nrow(contrast) * 2L, lower.tail = FALSE)

    dt.test[, c("effect_size_hurdle", "pvalue_hurdle", "pvalue_2part") := list(list(effect_size_hurdle), pvalue_hurdle, pvalue_2part)]
    if (nrow(contrast) > 1L) {
      dt.test[, paste0("effect_size_hurdle_", seq_len(nrow(contrast))) := transpose(effect_size_hurdle)]
      dt.test[, effect_size_hurdle := NULL]
    }
  }




  # test_statistic <- tryCatch({
  #   solve(cov_weighted_sum) %*% mean_weighted_sum
  # }, error = function(e) {
  #   MASS::ginv(cov_weighted_sum) %*% mean_weighted_sum
  # })



  # return(list(mean_weighted_sum = mean_weighted_sum,
  #             test_statistic = test_statistic,
  #             p_value = p_value))

  dt.test
}


## truncated normal
rtrunc_norm <- function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  lower_p <- pnorm(lower, mean = mean, sd = sd)
  upper_p <- pnorm(upper, mean = mean, sd = sd)

  u_samples <- runif(n, min = lower_p, max = upper_p)

  normal_samples <- qnorm(u_samples, mean = mean, sd = sd)

  normal_samples
}

## inverse_gamma
rinvgamma <- function(n, shape, rate = 1, scale = 1 / rate) {
  1 / rgamma(n, shape, scale)
}
