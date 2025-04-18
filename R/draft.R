center_covar <- function(dt, covar) {
  dt_copy <- copy(dt) # Avoid modifying original data.table
  
  categorical_covar <- covar[sapply(dt_copy[, covar, with = FALSE], function(col) is.factor(col) || is.character(col))]
  
  # Convert categorical covar columns to numeric (one-hot encoding)
  for (col_name in categorical_covar) {
    col <- dt_copy[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt_copy[, (new_col_name) := as.integer(col == level)]
    }
    dt_copy[, (col_name) := NULL] # Remove original categorical column
    #Add new one hot encoded columns to the covar list, and remove the original column.
    covar <- c(covar[covar != col_name], names(dt_copy)[startsWith(names(dt_copy), paste0(col_name, "_"))])
  }
  
  # Center specified columns
  for (col_name in covar) {
    dt_copy[, (col_name) := scale(dt_copy[[col_name]], center = TRUE, scale = FALSE)]
  }
  
  return(dt_copy)
}

preprocess_dt_center_covar <- function(dt, covar, make_copy = TRUE) {
  if (make_copy) {
    dt_copy <- copy(dt)
  } else {
    dt_copy <- dt
  }
  
  covar_exists <- covar %in% names(dt_copy)
  if (!all(covar_exists)) {
    stop(paste0("The following columns specified in covar do not exist in the data.table: ", paste(covar[!covar_exists], collapse = ", ")))
  }
  
  categorical_covar <- covar[sapply(dt_copy[, covar, with = FALSE], function(col) is.factor(col) || is.character(col))]
  
  for (col_name in categorical_covar) {
    col <- dt_copy[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt_copy[, (new_col_name) := as.integer(col == level)]
    }
    dt_copy[, (col_name) := NULL]
    covar <- c(covar[covar != col_name], names(dt_copy)[startsWith(names(dt_copy), paste0(col_name, "_"))])
  }
  
  for (col_name in covar) {
    dt_copy[, (col_name) := scale(dt_copy[[col_name]], center = TRUE, scale = FALSE)]
  }
  
  return(dt_copy)
}

# Example Usage:
dt <- data.table(
  A = factor(c("a", "b", "a")),
  B = c(1, 2, 3),
  C = c("x", "y", "z"),
  D = 4:7,
  E = 8:11
)

covar_cols <- c("A", "B", "E")

# Default: make_copy = TRUE
preprocessed_dt_copy <- preprocess_dt_center_covar(dt, covar_cols)
print(preprocessed_dt_copy)
print(dt) # Original 'dt' is unchanged

# make_copy = FALSE
preprocessed_dt_no_copy <- preprocess_dt_center_covar(dt, covar_cols, make_copy = FALSE)
print(preprocessed_dt_no_copy)
print(dt) # Original 'dt' is now modified

vcov_lm <- function(fit, group_names) {
  
}


vcov_logm <- function(fit, group_names) {
  
}








#' 
#' @import data.table
#' 
filter_lr <- function(expression_matrix, metadata_subset,
                      cdr,
                      sender, receiver,
                      lr_table,
                      multi_sub = c("minimum","arithmetic_mean","geometric_mean","min_avg_gene","min_rate_gene"),
                      contrast,
                      verbose = TRUE,
                      min_pct = 0.01, large_n = 2, min_avg_pct = 0,
                      min_cell = 10,
                      threshold = 0,
                      ...) {
  if (is.null(id)) id <- group
  
}

# compute_expression_rates <- function(chunk, metadata_subset, expression_matrix, threshold, min_pct, large_n, method = "min") {
#   library(data.table)
#   
#   for (i in 1:nrow(chunk)) {
#     sender <- chunk$sender[i]
#     ligand <- chunk$ligand[i]
#     receiver <- chunk$receiver[i]
#     receptor <- chunk$receptor[i]
#     
#     # Create copies of metadata_subset for sender and receiver
#     data_sender_ligand <- copy(metadata_subset)[cell_type == sender]
#     data_receiver_receptor <- copy(metadata_subset)[cell_type == receiver]
#     
#     # if (nrow(data_sender_ligand) == 0 || nrow(data_receiver_receptor) == 0) next  # Skip if no matching cells
#     
#     # Handle single gene or multi-gene ligands and receptors
#     ligand_genes <- unlist(strsplit(ligand, "_"))
#     receptor_genes <- unlist(strsplit(receptor, "_"))
#     
#     # Subset expression matrix for ligands and receptors
#     ligand_expr_values <- expression_matrix[ligand_genes, data_sender_ligand$cell_id, drop = FALSE]
#     receptor_expr_values <- expression_matrix[receptor_genes, data_receiver_receptor$cell_id, drop = FALSE]
#     
#     # Function to compute expression value based on method
#     compute_expression_value <- function(expr_values, method) {
#       if (length(dim(expr_values)) == 0) return(expr_values)
#       
#       switch(
#         method,
#         "min" = apply(expr_values, 2, min, na.rm = TRUE),
#         "mean" = colMeans(expr_values, na.rm = TRUE),
#         "lowest" = {
#           gene_rates <- colSums(expr_values > threshold, na.rm = TRUE) / ncol(expr_values)
#           min_gene <- names(which.min(gene_rates))
#           expr_values[min_gene, ]
#         },
#         stop("Invalid 'method'. Choose from 'min', 'mean', or 'lowest'.")
#       )
#     }
#     
#     # Add expression column to metadata copies
#     data_sender_ligand[, expression := compute_expression_value(ligand_expr_values, method)]
#     data_receiver_receptor[, expression := compute_expression_value(receptor_expr_values, method)]
#     
#     # Compute expression rates
#     unique_ids <- unique(c(data_sender_ligand$id, data_receiver_receptor$id))
#     
#     compute_expression_rate <- function(data_subset) {
#       sapply(unique_ids, function(uid) {
#         sum(data_subset[id == uid, expression] > threshold, na.rm = TRUE) / nrow(data_subset[id == uid])
#       })
#     }
#     
#     ligand_expression_rates <- compute_expression_rate(data_sender_ligand)
#     receptor_expression_rates <- compute_expression_rate(data_receiver_receptor)
#     
#     # Check if the number of ids with both ligand and receptor expression rate >= min_pct is >= large_n
#     valid_ids <- sum((ligand_expression_rates >= min_pct) & (receptor_expression_rates >= min_pct), na.rm = TRUE)
#     if (valid_ids < large_n) next
#     
#     # Execute some code (Replace this with your actual code)
#     print(paste("Processing ligand:", ligand, "and receptor:", receptor))
#   }
# }

# library(data.table)
# library(lme4)
# 
# fit_models <- function(data_sender_ligand, data_receiver_receptor, threshold, id = NULL, 
#                        lmm_re = FALSE, logmm_re = FALSE, covar = NULL) {
#   
#   # Add indicator column
#   data_sender_ligand[, indicator := ifelse(expression > threshold, 1, 0)]
#   data_receiver_receptor[, indicator := ifelse(expression > threshold, 1, 0)]
#   
#   # Subset data
#   data_sender_ligand_1 <- data_sender_ligand[indicator == 1]
#   data_receiver_receptor_1 <- data_receiver_receptor[indicator == 1]
#   
#   # Determine model types
#   if (is.null(id_col)) {
#     lmm_re <- FALSE
#     logmm_re <- FALSE
#   }
#   
#   # Define covariates
#   covariates <- c("group", covar)
#   covariates <- covariates[!is.na(covariates)]
#   
#   # Center covariates
#   center_covar <- function(dt, covars) {
#     for (cov in covars) {
#       if (cov %in% names(dt) && is.numeric(dt[[cov]])) {
#         dt[, (cov) := scale(get(cov), center = TRUE, scale = FALSE)]
#       }
#     }
#   }
#   center_covariates(data_sender_ligand_1, covariates)
#   center_covariates(data_receiver_receptor_1, covariates)
#   center_covariates(data_sender_ligand, covariates)
#   center_covariates(data_receiver_receptor, covariates)
#   
#   # Define model formulas
#   fixed_effects <- paste(covariates, collapse = " + ")
#   formula_linear <- as.formula(paste("expression ~ 0 +", fixed_effects))
#   formula_logistic <- as.formula(paste("indicator ~ 0 +", fixed_effects))
#   
#   if (lmm_re) {
#     formula_linear <- as.formula(paste("expression ~ 0 +", fixed_effects, "+ (1|id)"))
#   }
#   if (logmm_re) {
#     formula_logistic <- as.formula(paste("indicator ~ 0 +", fixed_effects, "+ (1|id)"))
#   }
#   
#   # Fit models
#   fit_model <- function(data, formula, family = NULL) {
#     tryCatch({
#       if (is.null(family)) {
#         if ("id" %in% names(data) && lmm_re) {
#           lmer(formula, data = data)
#         } else {
#           lm(formula, data = data)
#         }
#       } else {
#         if ("id" %in% names(data) && logmm_re) {
#           glmer(formula, data = data, family = family)
#         } else {
#           glm(formula, data = data, family = family)
#         }
#       }
#     }, error = function(e) {
#       return(list(error = e$message))
#     })
#   }
#   
#   # Run models
#   results <- list(
#     linear_sender = fit_model(data_sender_ligand_1, formula_linear),
#     linear_receiver = fit_model(data_receiver_receptor_1, formula_linear),
#     logistic_sender = fit_model(data_sender_ligand, formula_logistic, binomial),
#     logistic_receiver = fit_model(data_receiver_receptor, formula_logistic, binomial)
#   )
#   
#   return(results)
# }

library(numDeriv)
func2 <- function(x) c(sin(x), cos(x))
x <- (0:1)*2*pi
jacobian(func2, x)
jacobian(func2, x, "complex")

func3 <- function(x) {
  x[1]*x[2]-x[3]*x[4]
}

jacobian(func3, 1:4)


library(glmmTMB)  # Ensure necessary libraries are loaded
library(Matrix)
library(numDeriv)  # For delta method calculations

compute_test_stat <- function(fit.l.lmm, fit.l.logmm, fit.r.lmm, fit.r.logmm, contrast) {
  # Extract fixed effect estimates for 'group'
  beta_l_lmm <- fixef(fit.l.lmm)["group"]
  beta_r_lmm <- fixef(fit.r.lmm)["group"]
  beta_l_logmm <- fixef(fit.l.logmm)["group"]
  beta_r_logmm <- fixef(fit.r.logmm)["group"]
  
  # Ensure matching order with contrast columns
  group_levels <- colnames(contrast)
  beta_l_lmm <- beta_l_lmm[group_levels]
  beta_r_lmm <- beta_r_lmm[group_levels]
  beta_l_logmm <- beta_l_logmm[group_levels]
  beta_r_logmm <- beta_r_logmm[group_levels]
  
  # Compute transformed parameter values
  inv_logit <- function(x) exp(x) / (1 + exp(x))
  param_vector <- beta_l_lmm * beta_r_lmm * inv_logit(beta_l_logmm) * inv_logit(beta_r_logmm)
  
  # Compute weighted sum
  weighted_sum <- contrast %*% param_vector
  
  # Compute mean of weighted sum
  mean_vector <- colSums(contrast * param_vector)
  
  # Compute covariance using the delta method
  cov_matrix <- matrix(0, nrow = nrow(contrast), ncol = nrow(contrast))
  jacobian_fun <- function(params) contrast %*% params
  
  param_estimates <- c(beta_l_lmm, beta_r_lmm, beta_l_logmm, beta_r_logmm)
  vcov_combined <- bdiag(vcov(fit.l.lmm), vcov(fit.r.lmm), vcov(fit.l.logmm), vcov(fit.r.logmm))
  jacobian <- jacobian(jacobian_fun, param_estimates)
  cov_weighted_sum <- jacobian %*% vcov_combined %*% t(jacobian)
  
  # Compute test statistic and p-value
  test_stat <- t(mean_vector) %*% solve(cov_weighted_sum) %*% mean_vector
  p_value <- 1 - pchisq(test_stat, df = length(mean_vector))
  
  list(mean = mean_vector, test_stat = test_stat, p_value = p_value)
}

library(lme4)
library(MASS)

calculate_weighted_sum_test <- function(fit.l.lmm, fit.l.logmm, fit.r.lmm, fit.r.logmm, contrast) {
  
  group_names <- colnames(contrast)
  
  coef_l_lmm <- fixef(fit.l.lmm)
  coef_l_logmm <- fixef(fit.l.logmm)
  coef_r_lmm <- fixef(fit.r.lmm)
  coef_r_logmm <- fixef(fit.r.logmm)
  
  # Filter coefficients to only include group names
  coef_l_lmm <- coef_l_lmm[names(coef_l_lmm) %in% group_names]
  coef_l_logmm <- coef_l_logmm[names(coef_l_logmm) %in% group_names]
  coef_r_lmm <- coef_r_lmm[names(coef_r_lmm) %in% group_names]
  coef_r_logmm <- coef_r_logmm[names(coef_r_logmm) %in% group_names]
  
  # Ensure coefficients are named vectors
  if (is.null(names(coef_l_lmm))) names(coef_l_lmm) <- names(fixef(fit.l.lmm)[names(fixef(fit.l.lmm)) %in% group_names])
  if (is.null(names(coef_l_logmm))) names(coef_l_logmm) <- names(fixef(fit.l.logmm)[names(fixef(fit.l.logmm)) %in% group_names])
  if (is.null(names(coef_r_lmm))) names(coef_r_lmm) <- names(fixef(fit.r.lmm)[names(fixef(fit.r.lmm)) %in% group_names])
  if (is.null(names(coef_r_logmm))) names(coef_r_logmm) <- names(fixef(fit.r.logmm)[names(fixef(fit.r.logmm)) %in% group_names])
  
  # Reorder coefficients to match group_names order
  coef_l_lmm <- coef_l_lmm[group_names]
  coef_l_logmm <- coef_l_logmm[group_names]
  coef_r_lmm <- coef_r_lmm[group_names]
  coef_r_logmm <- coef_r_logmm[group_names]
  
  product_vector <- numeric(length(group_names))
  names(product_vector) <- group_names
  
  for (group in group_names) {
    product_vector[group] <- coef_l_lmm[group] * coef_r_lmm[group] * plogis(coef_l_logmm[group]) * plogis(coef_r_logmm[group])
  }
  
  weighted_sum <- contrast %*% product_vector
  
  # Delta method for mean and covariance
  gradient_matrix <- matrix(0, nrow = nrow(contrast), ncol = length(group_names) * 4)
  
  for (i in seq_len(nrow(contrast))) {
    for (j in seq_along(group_names)) {
      group <- group_names[j]
      
      ind_l_lmm <- which(names(coef_l_lmm) == group)
      ind_l_logmm <- which(names(coef_l_logmm) == group) + length(group_names)
      ind_r_lmm <- which(names(coef_r_lmm) == group) + length(group_names) * 2
      ind_r_logmm <- which(names(coef_r_logmm) == group) + length(group_names) * 3
      
      gradient_matrix[i, ind_l_lmm] <- contrast[i, j] * coef_r_lmm[group] * plogis(coef_l_logmm[group]) * plogis(coef_r_logmm[group])
      gradient_matrix[i, ind_l_logmm] <- contrast[i, j] * coef_l_lmm[group] * coef_r_lmm[group] * plogis(coef_r_logmm[group]) * dlogis(coef_l_logmm[group])
      gradient_matrix[i, ind_r_lmm] <- contrast[i, j] * coef_l_lmm[group] * plogis(coef_l_logmm[group]) * plogis(coef_r_logmm[group])
      gradient_matrix[i, ind_r_logmm] <- contrast[i, j] * coef_l_lmm[group] * coef_r_lmm[group] * plogis(coef_l_logmm[group]) * dlogis(coef_r_logmm[group])
    }
  }
  
  # Extract relevant parts of vcov matrices
  vcov_l_lmm_group <- vcov(fit.l.lmm)[group_names, group_names]
  vcov_l_logmm_group <- vcov(fit.l.logmm)[group_names, group_names]
  vcov_r_lmm_group <- vcov(fit.r.lmm)[group_names, group_names]
  vcov_r_logmm_group <- vcov(fit.r.logmm)[group_names, group_names]
  
  vcov_combined <- bdiag(vcov_l_lmm_group, vcov_l_logmm_group, vcov_r_lmm_group, vcov_r_logmm_group)
  
  cov_weighted_sum <- gradient_matrix %*% as.matrix(vcov_combined) %*% t(gradient_matrix)
  
  mean_weighted_sum <- weighted_sum
  
  test_statistic <- tryCatch({
    solve(cov_weighted_sum) %*% mean_weighted_sum
  }, error = function(e) {
    MASS::ginv(cov_weighted_sum) %*% mean_weighted_sum
  })
  
  p_value <- pchisq(test_statistic, df = nrow(contrast), lower.tail = FALSE)
  
  return(list(mean_weighted_sum = mean_weighted_sum, 
              test_statistic = test_statistic, 
              p_value = p_value))
}

### check required arguments
check_required_args <- function(f=sys.function(), mc=match.call()) {
  fmls <- formals(f)
  required_args <- names(fmls)[sapply(fmls, is.symbol)]
  user_args <- names(mc)[-1]
  missing_required <- setdiff(required_args, user_args)
  
  missing_required
}

h <- function(x,y) {
  missing_args <- check_required_args(f=sys.function(),mc=match.call())
  if (length(missing_args)>0) {
    stop("no")
  }
  #x+y
}
h(1)

hh <- function(x,y) {
  x+y
}

# intercell_network <- OmnipathR::intercell_network(ligand_receptor = TRUE, high_confidence = TRUE, simplify = TRUE)
# omnipathr <- data.frame(ligand=intercell_network$source_genesymbol,receptor=intercell_network$target_genesymbol)
# save(omnipathr, file = "omnipathr.rda")


# 
# oplan <- plan(multisession, workers = 2)
# on.exit(plan(oplan), add = TRUE)