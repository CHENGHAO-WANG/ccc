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
  
  # Check if we can perform linear model tests
  if (!is.null(fit_ligand_linear) && !is.null(fit_receptor_linear)) {
    # Extract coefficients for ligand and receptor
    if (isTRUE(lmm_re)) {
      coef_ligand_lm <- lme4::fixef(fit_ligand_linear)["grouptarget"]
      coef_receptor_lm <- lme4::fixef(fit_receptor_linear)["grouptarget"]
      
      if (isTRUE(sandwich)) {
        vcov_ligand_lm <- clubSandwich::vcovCR(fit_ligand_linear, type = "CR2")["grouptarget", "grouptarget"]
        vcov_receptor_lm <- clubSandwich::vcovCR(fit_receptor_linear, type = "CR2")["grouptarget", "grouptarget"]
      } else {
        vcov_ligand_lm <- vcov(fit_ligand_linear)["grouptarget", "grouptarget"]
        vcov_receptor_lm <- vcov(fit_receptor_linear)["grouptarget", "grouptarget"]
      }
    } else {
      coef_ligand_lm <- stats::coef(fit_ligand_linear)["grouptarget"]
      coef_receptor_lm <- stats::coef(fit_receptor_linear)["grouptarget"]
      
      if (isTRUE(sandwich)) {
        vcov_ligand_lm <- sandwich::vcovHC(fit_ligand_linear, type = "HC3")["grouptarget", "grouptarget"]
        vcov_receptor_lm <- sandwich::vcovHC(fit_receptor_linear, type = "HC3")["grouptarget", "grouptarget"]
      } else {
        vcov_ligand_lm <- vcov(fit_ligand_linear)["grouptarget", "grouptarget"]
        vcov_receptor_lm <- vcov(fit_receptor_linear)["grouptarget", "grouptarget"]
      }
    }
    
    # Calculate effect size (difference between target and background for both ligand and receptor)
    effect_size_linear <- coef_ligand_lm - coef_receptor_lm
    
    # Calculate standard error of the difference
    se_linear <- sqrt(vcov_ligand_lm + vcov_receptor_lm)
    
    # Calculate test statistic (t-statistic for one-sided test)
    test_stat_linear <- effect_size_linear / se_linear
    
    # Calculate one-sided p-value (testing if ligand > receptor in target vs background)
    pvalue_linear <- pnorm(test_stat_linear, lower.tail = FALSE)
    
    test.linear <- TRUE
  }
  
  # Check if we can perform logistic model tests
  if (!is.null(fit_ligand_logistic) && !is.null(fit_receptor_logistic)) {
    # Extract coefficients for ligand and receptor
    if (isTRUE(logmm_re)) {
      m_ligand <- GLMMadaptive::marginal_coefs(fit_ligand_logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      m_receptor <- GLMMadaptive::marginal_coefs(fit_receptor_logistic, std_errors = TRUE, cores = 1L, sandwich = sandwich)
      
      coef_ligand_logm <- m_ligand$betas["grouptarget"]
      coef_receptor_logm <- m_receptor$betas["grouptarget"]
      
      vcov_ligand_logm <- m_ligand$var_betas["grouptarget", "grouptarget"]
      vcov_receptor_logm <- m_receptor$var_betas["grouptarget", "grouptarget"]
    } else {
      coef_ligand_logm <- stats::coef(fit_ligand_logistic)["grouptarget"]
      coef_receptor_logm <- stats::coef(fit_receptor_logistic)["grouptarget"]
      
      if (isTRUE(sandwich)) {
        vcov_ligand_logm <- sandwich::vcovHC(fit_ligand_logistic, type = "HC3")["grouptarget", "grouptarget"]
        vcov_receptor_logm <- sandwich::vcovHC(fit_receptor_logistic, type = "HC3")["grouptarget", "grouptarget"]
      } else {
        vcov_ligand_logm <- vcov(fit_ligand_logistic)["grouptarget", "grouptarget"]
        vcov_receptor_logm <- vcov(fit_receptor_logistic)["grouptarget", "grouptarget"]
      }
    }
    
    # Calculate effect size (difference in probabilities)
    effect_size_logistic <- plogis(coef_ligand_logm) - plogis(coef_receptor_logm)
    
    # Calculate gradient for delta method
    gradient_ligand <- dlogis(coef_ligand_logm)
    gradient_receptor <- -dlogis(coef_receptor_logm)
    
    # Calculate variance using delta method
    var_logistic <- (gradient_ligand^2 * vcov_ligand_logm) + (gradient_receptor^2 * vcov_receptor_logm)
    se_logistic <- sqrt(var_logistic)
    
    # Calculate test statistic
    test_stat_logistic <- effect_size_logistic / se_logistic
    
    # Calculate one-sided p-value (testing if ligand > receptor in target vs background)
    pvalue_logistic <- pnorm(test_stat_logistic, lower.tail = FALSE)
    
    test.logistic <- TRUE
  }
  
  # Add basic information to results
  if (test.linear || test.logistic) {
    dt.test[, c("sender", "receiver", "ligand", "receptor") := list(sender, receiver, ligand, receptor)]
  }
  
  # Add linear model results if available
  if (isTRUE(test.linear)) {
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