utils::globalVariables(c(
  ".", "N", "V1", "cdr", "cell_id", "cell_type", "group", "id", "y", "z",
  "omnipathr", "ramilowski", "dispersion", "gene_id", "within_sample_correlation",
  "cell_continuous", "cell_discrete", "sample_continuous", "sample_discrete",
  "sender", "receiver", "background_sep", "target_sep", "n_classes" 
))

center_covar <- function(dt, covar) {
  centered_column_names <- character(0)

  for (col_name in covar) {
    col <- dt[[col_name]]

    if (is.factor(col) || is.character(col)) {
      levels_col <- sort(unique(col))
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

detect_re_separation2 <- function(dt, z_col, id_col, num_ids = NULL, sep_sample_prop = 0, sep_sample_n = 0) {
  # Count the number of distinct id's
  if (is.null(num_ids)) {
    num_ids <- dt[, uniqueN(get(id_col))]
  }
  
  # Function to check if all z values are 0 or all are 1
  is_separated <- function(z) all(z == 0) || all(z == 1)
  
  # Check per id whether it's separated within target or background
  count_sep_ids <- dt[
    , .(target_sep = is_separated(get(z_col)[class == "target"]),
        background_sep = is_separated(get(z_col)[class == "background"])),
    by = get(id_col)
  ][, sum(target_sep | background_sep)]
  
  # Check if the proportion exceeds the threshold
  detection <- count_sep_ids > sep_sample_prop * num_ids && count_sep_ids > sep_sample_n
  
  return(detection)
}

detect_all_zeros <- function(dt, id_col, id) {
  unique_id_col <- unique(dt[[id_col]]) # Get unique values from id_col
  !all(id %in% unique_id_col)
}

detect_all_zeros2 <- function(dt, id_col, id) {
  # Get ids that have both 'target' and 'background'
  id_with_both_classes <- dt[, .(n_classes = uniqueN(class)), by = get(id_col)][n_classes == 2, get]
  
  # Check if all input ids are in that set
  !all(id %in% id_with_both_classes)
}


compute_group_stats <- function(dt, y_col = "y", group_col = "group", z_col = "z", prefix) {
  # Compute stats for groups that are present in the data
  stats_dt <- dt[, .(
    mean = mean(get(y_col)),
    mean_se = sd(get(y_col))/sqrt(.N),
    positive_mean = if (any(get(z_col) == 1)) mean(get(y_col)[get(z_col) == 1]) else 0,
    positive_mean_se = if (any(get(z_col) == 1)) sd(get(y_col)[get(z_col) == 1])/sqrt(sum(get(z_col) == 1)) else 0,
    expression_rate = mean(get(z_col) == 1),
    expression_rate_se = sqrt(mean(get(z_col) == 1) * (1 - mean(get(z_col) == 1))/sqrt(.N))
  ), by = get(group_col)]

  setnames(stats_dt, "get", group_col)

  stats_cols <- setdiff(names(stats_dt), group_col)
  setnames(stats_dt, stats_cols, paste0(prefix, stats_cols))
  stats_dt
}


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
  if (is.null(id_col)) {
    setnames(metadata, old = c(cell_id_col, group_col, cell_type_col), new = c("cell_id", "group", "cell_type"))
    metadata[, id := group]
  } else {
    data.table::setnames(metadata, old = c(cell_id_col, id_col, group_col, cell_type_col), new = c("cell_id", "id", "group", "cell_type")) # , skip_absent = T)
  }
  return(metadata)
}

rename_metadata2 <- function(metadata, cell_id_col, id_col, cell_type_col) {
  if (is.null(id_col)) {
    setnames(metadata, old = c(cell_id_col, cell_type_col), new = c("cell_id", "cell_type"))
    metadata[, id := "constant"]
  } else {
    data.table::setnames(metadata, old = c(cell_id_col, id_col, cell_type_col), new = c("cell_id", "id", "cell_type")) # , skip_absent = T)
  }
  return(metadata)
}

filter_cell_type <- function(metadata, sender, receiver, min_cell) {
  if (is.null(sender)) {
    sender <- unique(metadata$cell_type)
  }
  if (is.null(receiver)) {
    receiver <- unique(metadata$cell_type)
  }
  
  # To keep all group levels in model fittings
  # metadata <- metadata[group %in% colnames(contrast)]

  # Find rows where cell_type appears in either sender OR receiver
  valid_cell_types <- union(sender, receiver)
  metadata_subset <- metadata[cell_type %in% valid_cell_types]

  # Count occurrences of each id-cell_type combination
  counts <- metadata_subset[, .N, by = .(id, cell_type)]

  # Keep only cell_type where count >= min_cell for each id
  metadata_subset <- metadata_subset[counts[N >= min_cell], on = .(id, cell_type)]

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
  # Subset the expression matrix to only include the cell ids in metadata_subset
  common_cells <- intersect(colnames(expression_matrix), metadata_subset$cell_id)

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

stouffer_combine_pvalues <- function(pvalues) {
  # Convert p-values to z-scores
  z_scores <- qnorm(pvalues, lower.tail = FALSE)
  # Combine z-scores using Stouffer's method
  combined_z <- sum(z_scores) / sqrt(length(z_scores))
  # Convert back to p-value
  pnorm(combined_z, lower.tail = FALSE)
}

fisher_combine_pvalues <- function(pvalues) {
  # Fisher's method: -2 * sum(log(p))
  chi_sq <- -2 * sum(log(pvalues))
  # Degrees of freedom is 2 * number of p-values
  df <- 2 * length(pvalues)
  # Convert to p-value
  pchisq(chi_sq, df = df, lower.tail = FALSE)
}

ccc_estimate <- function(fit.l.linear, fit.l.logistic, fit.r.linear, fit.r.logistic,
                     unique_levels, lmm_re, logmm_re, sandwich,
                     sender, receiver, ligand, receptor, var_to_test = c("group", "class"),
                     marginal_cores = 1L, marginal = FALSE, approx = TRUE, num_ids = NULL) {
  dt.est <- data.table()
  var_to_test <- var_to_test[1L]
  var_names <- paste0(var_to_test, unique_levels)
  
  test.linear <- FALSE
  test.logistic <- FALSE
  if (!is.null(fit.l.linear) && !is.null(fit.r.linear)) {
    if (isTRUE(lmm_re)) {
      coef_l_lm <- lme4::fixef(fit.l.linear)
      coef_r_lm <- lme4::fixef(fit.r.linear)
      if (isTRUE(sandwich)) {
        vcov_l_lm <- clubSandwich::vcovCR(fit.l.linear, type = "CR1p")[var_names, var_names]
        vcov_r_lm <- clubSandwich::vcovCR(fit.r.linear, type = "CR1p")[var_names, var_names]
      }
    } else {
      coef_l_lm <- stats::coef(fit.l.linear)
      coef_r_lm <- stats::coef(fit.r.linear)
      if (isTRUE(sandwich)) {
        vcov_l_lm <- sandwich::vcovHC(fit.l.linear, type = "HC3")[var_names, var_names]
        vcov_r_lm <- sandwich::vcovHC(fit.r.linear, type = "HC3")[var_names, var_names]
      }
    }
    if (isFALSE(sandwich)) {
      vcov_l_lm <- vcov(fit.l.linear)[var_names, var_names]
      vcov_r_lm <- vcov(fit.r.linear)[var_names, var_names]
    }
    coef_l_lm <- coef_l_lm[var_names]
    coef_r_lm <- coef_r_lm[var_names]
    test.linear <- TRUE
  }
  
  if (!is.null(fit.l.logistic) && !is.null(fit.r.logistic)) {
    if (isTRUE(logmm_re)) {
      if (isTRUE(marginal)) {
        if (isTRUE(approx)) {
          # Use attenuation approximation for marginal coefficients
          # First get conditional coefficients
          coef_l_logm <- lme4::fixef(fit.l.logistic)
          coef_r_logm <- lme4::fixef(fit.r.logistic)
          vcov_l_logm <- vcov(fit.l.logistic, sandwich = sandwich)[var_names, var_names]
          vcov_r_logm <- vcov(fit.r.logistic, sandwich = sandwich)[var_names, var_names]
          
          # Get random intercept variances
          var_l_random <- fit.l.logistic$D[1, 1]
          var_r_random <- fit.r.logistic$D[1, 1]
          
          # Calculate attenuation factors
          # Factor = 1 / sqrt(1 + (16/15) * (sqrt(3)/pi) * sigma^2)
          factor_l <- 1 / sqrt(1 + (16/15) * (sqrt(3)/pi) * var_l_random)
          factor_r <- 1 / sqrt(1 + (16/15) * (sqrt(3)/pi) * var_r_random)
          
          # Apply attenuation correction to coefficients
          coef_l_logm <- coef_l_logm * factor_l
          coef_r_logm <- coef_r_logm * factor_r
          
          # Apply attenuation correction to covariance matrices
          vcov_l_logm <- vcov_l_logm * (factor_l^2)
          vcov_r_logm <- vcov_r_logm * (factor_r^2)
        } else {
          # Use exact marginal coefficients
          m_l <- GLMMadaptive::marginal_coefs(fit.l.logistic, std_errors = TRUE, cores = marginal_cores, sandwich = sandwich)
          m_r <- GLMMadaptive::marginal_coefs(fit.r.logistic, std_errors = TRUE, cores = marginal_cores, sandwich = sandwich)
          coef_l_logm <- m_l$betas
          coef_r_logm <- m_r$betas
          vcov_l_logm <- m_l$var_betas[var_names, var_names]
          vcov_r_logm <- m_r$var_betas[var_names, var_names]
        }
      } else {
        # Use conditional coefficients (default behavior)
        coef_l_logm <- lme4::fixef(fit.l.logistic)
        coef_r_logm <- lme4::fixef(fit.r.logistic)
        vcov_l_logm <- vcov(fit.l.logistic, sandwich = sandwich)[var_names, var_names]
        vcov_r_logm <- vcov(fit.r.logistic, sandwich = sandwich)[var_names, var_names]
      }
      
      if (isTRUE(sandwich) && isTRUE(marginal)) {
        # small cluster number correction for marginal coefficients
        vcov_l_logm <- vcov_l_logm * (num_ids / (num_ids - 1))
        vcov_r_logm <- vcov_r_logm * (num_ids / (num_ids - 1))
      }
    } else {
      # For non-mixed models, marginal and approx arguments are ignored
      coef_l_logm <- stats::coef(fit.l.logistic)
      coef_r_logm <- stats::coef(fit.r.logistic)
      if (isTRUE(sandwich)) {
        vcov_l_logm <- sandwich::vcovHC(fit.l.logistic, type = "HC3")[var_names, var_names]
        vcov_r_logm <- sandwich::vcovHC(fit.r.logistic, type = "HC3")[var_names, var_names]
      } else {
        vcov_l_logm <- vcov(fit.l.logistic)[var_names, var_names]
        vcov_r_logm <- vcov(fit.r.logistic)[var_names, var_names]
      }
    }
    coef_l_logm <- coef_l_logm[var_names]
    coef_r_logm <- coef_r_logm[var_names]
    test.logistic <- TRUE
  }

  # if (isTRUE(test.linear)) {
  #   coef_l_lm <- coef_l_lm[names(coef_l_lm) %in% var_names]
  #   coef_r_lm <- coef_r_lm[names(coef_r_lm) %in% var_names]
  #   coef_l_lm <- coef_l_lm[var_names]
  #   coef_r_lm <- coef_r_lm[var_names]
  # }
  # if (isTRUE(test.logistic)) {
  #   coef_l_logm <- coef_l_logm[names(coef_l_logm) %in% var_names]
  #   coef_r_logm <- coef_r_logm[names(coef_r_logm) %in% var_names]
  #   coef_l_logm <- coef_l_logm[var_names]
  #   coef_r_logm <- coef_r_logm[var_names]
  # }
  if (test.linear || test.logistic) {
    dt.est[, c("sender", "receiver", "ligand", "receptor") := list(sender, receiver, ligand, receptor)]
  }
  
  if(test.linear) {
    dt.est[, c("coef_l_lm", "coef_r_lm") := list(list(coef_l_lm), list(coef_r_lm))]
    dt.est[, c("vcov_l_lm", "vcov_r_lm") := list(list(vcov_l_lm), list(vcov_r_lm))]
  }
  if (test.logistic) {
    dt.est[, c("coef_l_logm", "coef_r_logm") := list(list(coef_l_logm), list(coef_r_logm))]
    dt.est[, c("vcov_l_logm", "vcov_r_logm") := list(list(vcov_l_logm), list(vcov_r_logm))]
  }

  dt.est
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
rinvgamma <- function(n, shape, rate, scale = 1 / rate) {
  1 / rgamma(n, shape = shape, rate = scale)
}

## arguments that are different from defaults
args_differ_from_defaults <- function(...) {
  arg_names <- c(...)
  f <- sys.function(sys.parent())
  env <- parent.frame()
  defaults <- formals(f)
  differing <- vapply(arg_names, function(arg) {
    current_val <- get(arg, envir = env)
    default_val <- eval(defaults[[arg]], envir = env)
    !identical(current_val, default_val)
  }, logical(1))
  
  names(differing)[differing]
}
