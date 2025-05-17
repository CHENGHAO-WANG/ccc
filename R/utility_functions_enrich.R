#' One-sided Test for Enrichment Analysis
#'
#' This function performs one-sided tests for enrichment analysis, comparing ligand vs receptor
#' expression in target vs background cell type pairs. It calculates effect sizes and p-values 
#' for linear and logistic models, and combines them using various methods.
#'
#' @param fit_ligand_linear Linear model fit for ligand expression
#' @param fit_ligand_logistic Logistic model fit for ligand expression
#' @param fit_receptor_linear Linear model fit for receptor expression
#' @param fit_receptor_logistic Logistic model fit for receptor expression
#' @param lmm_re Logical, whether linear models include random effects
#' @param logmm_re Logical, whether logistic models include random effects
#' @param sandwich Logical, whether to use sandwich variance estimators
#' @param sender Sender cell type
#' @param receiver Receiver cell type
#' @param ligand Ligand gene
#' @param receptor Receptor gene
#'
#' @return A data.table with test results
#' @importFrom stats pnorm
#' @importFrom stats pt
#' @importFrom stats vcov
#' @importFrom stats coef
#' @importFrom stats plogis
#' @importFrom stats dlogis
#' @importFrom data.table data.table

ccc_test_enrich <- function(fit_ligand_linear, fit_ligand_logistic, 
                           fit_receptor_linear, fit_receptor_logistic,
                           lmm_re, logmm_re, sandwich,
                           sender, receiver, ligand, receptor) {
  dt.test <- data.table()
  
  # Initialize test flags
  test.linear <- FALSE
  test.logistic <- FALSE
  
  class_names <- c("classtarget", "classbackground")
  contrast <- matrix(c(1, -1), nrow = 1L)
  colnames(contrast) <- class_names
  
  # Check if we can perform linear model tests
  if (!is.null(fit_ligand_linear) && !is.null(fit_receptor_linear)) {
    # Extract coefficients for ligand and receptor
    if (isTRUE(lmm_re)) {
      coef_ligand_lm <- lme4::fixef(fit_ligand_linear)
      coef_receptor_lm <- lme4::fixef(fit_receptor_linear)
      
      if (isTRUE(sandwich)) {
        vcov_ligand_lm <- clubSandwich::vcovCR(fit_ligand_linear, type = "CR2")[class_names, class_names]
        vcov_receptor_lm <- clubSandwich::vcovCR(fit_receptor_linear, type = "CR2")[class_names, class_names]
      } else {
        vcov_ligand_lm <- vcov(fit_ligand_linear)[class_names, class_names]
        vcov_receptor_lm <- vcov(fit_receptor_linear)[class_names, class_names]
      }
    } else {
      coef_ligand_lm <- stats::coef(fit_ligand_linear)
      coef_receptor_lm <- stats::coef(fit_receptor_linear)
      
      if (isTRUE(sandwich)) {
        vcov_ligand_lm <- sandwich::vcovHC(fit_ligand_linear, type = "HC3")[class_names, class_names]
        vcov_receptor_lm <- sandwich::vcovHC(fit_receptor_linear, type = "HC3")[class_names, class_names]
      } else {
        vcov_ligand_lm <- vcov(fit_ligand_linear)[class_names, class_names]
        vcov_receptor_lm <- vcov(fit_receptor_linear)[class_names, class_names]
      }
    }
    
    # # Calculate effect size (difference between target and background for both ligand and receptor)
    # effect_size_linear <- coef_ligand_lm - coef_receptor_lm
    # 
    # # Calculate standard error of the difference
    # se_linear <- sqrt(vcov_ligand_lm + vcov_receptor_lm)
    # 
    # # Calculate test statistic (t-statistic for one-sided test)
    # test_stat_linear <- effect_size_linear / se_linear
    # 
    # # Calculate one-sided p-value (testing if ligand > receptor in target vs background)
    # pvalue_linear <- pnorm(test_stat_linear, lower.tail = FALSE)
    
    test.linear <- TRUE
  }
  
  # Check if we can perform logistic model tests
  if (!is.null(fit_ligand_logistic) && !is.null(fit_receptor_logistic)) {
    # Extract coefficients for ligand and receptor
    if (isTRUE(logmm_re)) {
      m_ligand <- GLMMadaptive::marginal_coefs(fit_ligand_logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      m_receptor <- GLMMadaptive::marginal_coefs(fit_receptor_logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      
      coef_ligand_logm <- m_ligand$betas
      coef_receptor_logm <- m_receptor$betas
      
      vcov_ligand_logm <- m_ligand$var_betas[class_names, class_names]
      vcov_receptor_logm <- m_receptor$var_betas[class_names, class_names]
    } else {
      coef_ligand_logm <- stats::coef(fit_ligand_logistic)
      coef_receptor_logm <- stats::coef(fit_receptor_logistic)
      
      if (isTRUE(sandwich)) {
        vcov_ligand_logm <- sandwich::vcovHC(fit_ligand_logistic, type = "HC3")[class_names, class_names]
        vcov_receptor_logm <- sandwich::vcovHC(fit_receptor_logistic, type = "HC3")[class_names, class_names]
      } else {
        vcov_ligand_logm <- vcov(fit_ligand_logistic)[class_names, class_names]
        vcov_receptor_logm <- vcov(fit_receptor_logistic)[class_names, class_names]
      }
    }
    
    # # Calculate effect size (difference in probabilities)
    # effect_size_logistic <- plogis(coef_ligand_logm) - plogis(coef_receptor_logm)
    # 
    # # Calculate gradient for delta method
    # gradient_ligand <- dlogis(coef_ligand_logm)
    # gradient_receptor <- -dlogis(coef_receptor_logm)
    # 
    # # Calculate variance using delta method
    # var_logistic <- (gradient_ligand^2 * vcov_ligand_logm) + (gradient_receptor^2 * vcov_receptor_logm)
    # se_logistic <- sqrt(var_logistic)
    # 
    # # Calculate test statistic
    # test_stat_logistic <- effect_size_logistic / se_logistic
    # 
    # # Calculate one-sided p-value (testing if ligand > receptor in target vs background)
    # pvalue_logistic <- pnorm(test_stat_logistic, lower.tail = FALSE)
    
    test.logistic <- TRUE
  }
  
  if (isTRUE(test.linear)) {
    coef_ligand_lm <- coef_ligand_lm[names(coef_ligand_lm) %in% class_names]
    coef_receptor_lm <- coef_receptor_lm[names(coef_receptor_lm) %in% class_names]
    coef_ligand_lm <- coef_ligand_lm[class_names]
    coef_receptor_lm <- coef_receptor_lm[class_names]
  }
  if (isTRUE(test.logistic)) {
    coef_ligand_logm <- coef_ligand_logm[names(coef_ligand_logm) %in% class_names]
    coef_receptor_logm <- coef_receptor_logm[names(coef_receptor_logm) %in% class_names]
    coef_ligand_logm <- coef_ligand_logm[class_names]
    coef_receptor_logm <- coef_receptor_logm[class_names]
  }
  
  # Add basic information to results
  if (test.linear || test.logistic) {
    dt.test[, c("sender", "receiver", "ligand", "receptor") := list(sender, receiver, ligand, receptor)]
  }
  
  # Add linear model results if available
  if (isTRUE(test.linear)) {
    effect_size_linear <- contrast %*% (coef_ligand_lm * coef_receptor_lm)
    gradient_matrix_linear <- matrix(0, nrow = 1L, ncol = 4L)
    
    for (j in seq_along(class_names)) {
      group <- group_names[j]
      
      ind_l_lm <- which(names(coef_l_lm) == group)
      ind_r_lm <- which(names(coef_r_lm) == group) + length(group_names)
      
      gradient_matrix_linear[i, ind_l_lm] <- contrast[i, j] * coef_r_lm[group]
      gradient_matrix_linear[i, ind_r_lm] <- contrast[i, j] * coef_l_lm[group]
    }
    
    dt.test[, c("effect_size_linear", "se_linear", "test_stat_linear", "pvalue_linear") := 
            list(effect_size_linear, se_linear, test_stat_linear, pvalue_linear)]
  }
  
  # Add logistic model results if available
  if (isTRUE(test.logistic)) {
    dt.test[, c("effect_size_logistic", "se_logistic", "test_stat_logistic", "pvalue_logistic") := 
            list(effect_size_logistic, se_logistic, test_stat_logistic, pvalue_logistic)]
  }
  
  # Combine p-values if both tests are available
  if (isTRUE(test.linear) && isTRUE(test.logistic)) {
    # Calculate hurdle effect size (product of linear and logistic effects)
    effect_size_hurdle <- coef_ligand_lm * plogis(coef_ligand_logm) - coef_receptor_lm * plogis(coef_receptor_logm)
    
    # Combine p-values using Stouffer's method
    pvalue_stouffer <- stouffer_combine_pvalues(c(pvalue_linear, pvalue_logistic))
    
    # Combine p-values using Fisher's method
    pvalue_fisher <- fisher_combine_pvalues(c(pvalue_linear, pvalue_logistic))
    
    dt.test[, c("effect_size_hurdle", "pvalue_stouffer", "pvalue_fisher") := 
            list(effect_size_hurdle, pvalue_stouffer, pvalue_fisher)]
  }
  
  return(dt.test)
}


has_separation <- function(subdt, z_col, class_col) {
  any(sapply(unique(subdt[[class_col]]), function(val) {
    z_vals <- subdt[[z_col]][subdt[[class_col]] == val]
    all(z_vals == 0) || all(z_vals == 1)
  }))
}

# for ccc_enrich
detect_separation <- function(dt, id_col, z_col, class_col, num_ids = NULL, sep_prop = 0, sep_n = 0) {
  if (is.null(num_ids)) {
    num_ids <- dt[, uniqueN(get(id_col))]
  }
  
  count_sep_ids <- dt[, .(separation = has_separation(.SD, z_col = z_col, class_col = class_col)), by = get(id_col)][separation == TRUE, .N]
  detection <- count_sep_ids > sep_prop * num_ids && count_sep_ids > sep_n
  
  detection
}

# dt <- data.table(sample_id = rep(1:3, each = 6),
#                  response = c(0,0,1,1,0,1,
#                               1,1,1,1,1,1,
#                               1,1,0,0,1,0),
#                  cell_group = rep(c("target", "background"), each = 3, times = 3))
# dt
# detect_separation(dt, "sample_id", "response", "cell_group")

# for ccc_enrich
# detect_zeros <- function(dt, id_col, z_col, class_col) {
#   dt[, has_nonzero := {
#     all(sapply(unique(get(class_col)), function(val) {
#       any(get(z_col)[get(class_col) == val] == 1)
#     }))
#   }, by = get(id_col)][, !all(has_nonzero)]
# }

# detect_zeros(dt, "sample_id", "response", "cell_group")

extract_id_class_reference <- function(dt, id_col, class_col) {
  unique(dt[, .(id = get(id_col), class = get(class_col))])
}

# for ccc_enrich
detect_zeros <- function(reference, dt1, id_col, class_col) {
  present <- unique(dt1[, .(id = get(id_col), class = get(class_col))])
  present[, is_present := TRUE]
  merged <- merge(reference, present, by = c("id", "class"), all.x = TRUE)
  any(is.na(merged$is_present))
}

# ref <- extract_id_class_reference(dt, "sample_id", "cell_group")
# dt1 <- dt[response == 1, ]
# dt1
# 
# detect_zeros(reference = ref, dt1 = dt1, "sample_id", "cell_group")
