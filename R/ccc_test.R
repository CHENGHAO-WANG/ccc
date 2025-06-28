#' Identify Differential or Enriched Cell-Cell Communication
#' 
#' Perform Wald tests on the product of ligand expression levels in sender and receptor expression levels in receiver for differential or enriched cell-cell communication.
#' 
#' @param ccc_obj output of [ccc::ccc_diff()] or [ccc::ccc_enrich()].
#' @param contrast a named numeric vector or a numeric matrix with column names. The names should match the levels of the variable being tested (specified by `group_col` argument of [ccc::ccc_diff()]; the testings are performed based on this). Only used if `ccc_obj` is the output of [ccc::ccc_diff()].
#' @param test_type a character string for the type of test. Either "chisq" or "z". If "chisq", perform Wald Chisq tests. If "z", perform Wald Z tests. Defaults to "chisq" for `ccc_diff` and "z" for `ccc_enrich`. Must be "chisq" if `contrast` has more than one row.
#' @param ha a character string specifying the alternative hypothesis. Let \eqn{\theta} denote the difference in the product of ligand and receptor expression levels comparing target to background cell types (`ccc_enrich`) or specified by `contrast` (`ccc_diff`).
#'  \itemize{
#'   \item \dQuote{\code{greater}}: \eqn{\theta > c}
#'   \item \dQuote{\code{less}}: \eqn{\theta < c}
#'   \item \dQuote{\code{greater.abs}}: \eqn{|\theta| > c}
#'   \item \dQuote{\code{less.abs}}: \eqn{|\theta| < c}
#'  }
#' Only used if `test_type = "z"`. Defaults to "greater" for `ccc_enrich`. \eqn{c} is specified using `c_linear`, `c_logistic`, and `c_hurdle` arguments.
#' @param c_linear,c_logistic,c_hurdle numeric scalar. `c_linear` is for the tests on the linear component, `c_logistic` is for the tests on the logistic component, and `c_hurdle` is for the tests on the hurdle model. Default to 0. Must be non-negative for `ha = "greater.abs"` and `ha = "less.abs"`.
#' @param score a character string specifying the method for computing effect sizes.
#'  \itemize{
#'   \item \dQuote{\code{product}}: (the default) effect sizes are computed as the product of ligand and receptor expressions. For linear component: product of conditional mean expressions. For logistic component: product of expression rates (probabilities). For hurdle: product of all four components.
#'   \item \dQuote{\code{sum}}: effect sizes are computed as the sum of ligand and receptor expressions. For linear component: sum of conditional mean expressions. For logistic component: sum of log-odds (coefficients). For hurdle: sum of (probability × conditional mean) for ligand and receptor separately.
#'  }
#' @param verbose logical scalar. If `TRUE` (the default), display a progress bar. The default handler is "progress". This package uses the \pkg{progressr} framework for progress reporting, so users can customize the progress bar. See [progressr::handlers()] for customizing progress bar behavior. Note that progress bars are not displayed in non-interactive mode or R markdown.
#' @param padj_method a character string for multiple testing correction method. This is passed to [stats::p.adjust()]. Defaults to "BH".
#' @param cell_type_padj logical scalar. If `TRUE` (the default), adjust p-values for each sender-receiver pair.
#' @param chunk_size integer scalar. The number of communication pairs (each defined by a distinct combination of sender, receiver, ligand, and receptor) per chunk. Passed to the `future.chunk.size` argument of [future.apply::future_mapply()]. Defaults to 10. To enable parallelization, users should use the \pkg{future} package.
#' 
#' @details
#' This function performs Wald tests based on the output of [ccc::ccc_diff()] or [ccc::ccc_enrich()].
#' By default, Wald Chisq tests are performed for differential ccc analysis and one-sided Wald Z tests are performed for enriched ccc analysis (alternative hypothesis: communication event in target cell types > background cell types).
#' 
#' If the linear component is fitted successfully, this function performs a Wald test on the product (or sum, if `score = "sum"`) of ligand and receptor conditional mean expressions (condition on expression levels > `threshold`).
#' If the logistic componenet is fitted successfully, this function performs a Wald test on the product of ligand and receptor expression rates (`score = "product"`) or the sum of ligand and receptor log-odds (`score = "sum"`).
#' If both components are fitted successfully, this function performs a Wald test on the combined effect. For `score = "product"`: conditional mean of ligand \verb{*} expression rate of ligand \verb{*} conditional mean of receptor \verb{*} expression rate of receptor. For `score = "sum"`: (expression rate of ligand \verb{*} conditional mean of ligand) + (expression rate of receptor \verb{*} conditional mean of receptor).
#' The delta method is applied to obtain the relevant standard errors, using the estimates from the output of `ccc_*`.
#' The effect sizes are computed as the product or sum (depending on `score`) of ligand and receptor expressions. For `score = "product"`: conditional mean expressions, expression rates, and mean expression levels are multiplied. For `score = "sum"`: conditional mean expressions and log-odds are summed (with appropriate transformations for hurdle model). Each element corresponds to a level of `group_col`, multiplied by `contrast` (`ccc_diff`), or the difference comparing target against background cell types (`ccc_enrich`). 
#' 
#' If `ha = "greater.abs"`, the test statistic is \eqn{(|\hat{\theta}| - c)/se(\hat{\theta})}, and the p-value is two-sided.
#' If `ha = "less.abs"`, the test statistics are \eqn{(\hat{\theta} - c)/se(\hat{\theta})} and \eqn{(\hat{\theta} + c)/se(\hat{\theta})}. The p-value is the maximum of the two one-sided tests.
#' 
#' When both linear and logistic components are fitted successfully in `ccc_*`, three different methods are used to combine their p-values:
#' \itemize{
#'   \item{\code{"2-part"}}: a chi-square test statistic is computed as the sum of the linear and logistic test statistics, with the degrees of freedom added. This is only available if `test_type = "chisq"`.
#'   \item{\code{"Stouffer's method"}}: p-values from the linear and logistic components are converted to Z-scores, summed up, scaled by \eqn{\sqrt{2}} and converted back to p-values.
#'   \item{\code{"Fisher's method"}}: a chi-square test statistic is computed as -2 times the sum of the logarithms of the p-values from the linear and logistic components, with the degrees of freedom equal to 4.
#' }
#' 
#' Parallelization is supported via \pkg{future}. Progress bars are customized using \pkg{progressr}.
#' 
#' @returns a data frame contains effect sizes, p-values and adjusted p-values for each cell-cell communication event.
#' 
#' @examples
#' \dontrun{
#' ## Data simulation
#' dat.sim <- sim_count(seed = 123)
#' expression_matrix <- log_normalize(dat.sim$counts)
#' metadata <- dat.sim$metadata
#' lrdb.sim <- sim_lr(seed = 456, n_lr = 10)
#' 
#' ## Differential ccc analysis
#' a <- ccc_diff(expression_matrix = expression_matrix, metadata = metadata,
#'    id_col = "sample", lr = lrdb.sim, sender = "CT1", receiver = "CT3")
#' 
#' contrast1 <- matrix(c(-1, 1, 0,
#'                       -1, 0, 1), nrow = 2L, byrow = TRUE,
#'                       dimnames = list(NULL, c("grp1","grp2","grp3")))
#' a.result1 <- ccc_test(a, contrast = contrast1)
#' head(a.result1)
#' 
#' contrast2 <- c(grp1 = -1, grp2 = 1)
#' a.result2 <- ccc_test(a, contrast = contrast2, test_type = "z", ha = "greater")
#' head(a.result2)
#' 
#' ## Enriched ccc analysis
#' b <- ccc_enrich(expression_matrix = expression_matrix, metadata = metadata,
#'     id_col = "sample", lr = lrdb.sim, sender = "CT1", receiver = "CT3")
#' 
#' b.result1 <- ccc_test(b)
#' head(b.result1)
#' 
#' b.result2 <- ccc_test(b, test_type = "z", ha = "greater",
#'     c_linear = 0.1, c_logistic = 0.1, c_hurdle = 0.1)
#' head(b.result2)
#' 
#' ## Using sum scoring method
#' a.result.sum <- ccc_test(a, contrast = contrast2, score = "sum")
#' head(a.result.sum)
#' 
#' b.result.sum <- ccc_test(b, score = "sum")
#' head(b.result.sum)
#' }
#' 
#' @export

ccc_test <- function(ccc_obj, contrast = NULL, test_type = NULL, ha = NULL,
                     c_linear = 0, c_logistic = 0, c_hurdle = 0,
                     score = c("product", "sum"), verbose = TRUE, padj_method = "BH", 
                     cell_type_padj = TRUE, chunk_size = 10) {
  # verbose
  assertthat::assert_that(assertthat::is.flag(verbose))
  
  if (verbose) {
    if (interactive() && isFALSE(getOption("rstudio.notebook.executing"))) {
      if (!handlers(global = NA)) {
        handlers(global = TRUE)
        handlers("progress")
        message(
          "Info: No global progress bars were found; the 'progress' handler has been enabled. ",
          "See `vignette('progressr-intro')` for how to customize the progress bar settings."
        )
      }
    }
  }
  
  assertthat::assert_that(assertthat::is.flag(cell_type_padj))
  padj_method <- match.arg(padj_method, stats::p.adjust.methods)
  score <- match.arg(score)
  
  if (!is.null(test_type)) {
    test_type <- match.arg(test_type, choices = c("chisq", "z"))
  }
  if (!is.null(ha)) {
    ha <- match.arg(ha, choices = c("greater.abs", "less.abs", "greater", "less"))
  }
  if (ccc_obj$func == "ccc_diff") {
    var_name <- "group"
    if (is.null(contrast)) {
      stop("`contrast` is required for differential ccc analysis.")
    }
    if (is.vector(contrast)) {
      contrast <- matrix(contrast, nrow = 1L, dimnames = list(NULL, names(contrast)))
    } else if (!is.matrix(contrast)) {
      stop("`contrast` must be either a named vector or a matrix with column names.")
    }
    if (is.null(colnames(contrast))) {
      stop("`contrast` must be either a named vector or a matrix with column names.")
    }
    if (qr(contrast)$rank < nrow(contrast)) {
      stop("`contrast` must be full row rank")
    }
    if (is.null(test_type)) {
      test_type <- "chisq"
      changed_args <- args_differ_from_defaults('ha', 'c_linear', 'c_logistic', 'c_hurdle')
      if (length(changed_args) > 0L) {
        message("Perform Wald Chisq tests for differential ccc analysis by default.")
      }
    }
    if (nrow(contrast) > 1L) {
      if (test_type == "z") {
        test_type <- "chisq"
        message("`test_type` is set to 'chisq' for multiple-row contrasts.")
      }
    }
    if (test_type == "z") {
      if (is.null(ha)) {
        stop("`ha` is required for Wald Z tests.")
      }
    }
  } else if (ccc_obj$func == "ccc_enrich") {
    var_name <- "class"
    if (!is.null(contrast)) {
      message("`contrast` is ignored for enriched ccc analysis. Evaluating the difference of target cell type pairs relative to background by design.")
    }
    contrast <- matrix(c(1, -1), nrow = 1L)
    colnames(contrast) <- c("target", "background")
    if (is.null(test_type)) {
      test_type <- "z"
    }
    if (is.null(ha)) {
      ha <- "greater"
    }
  } else {
    stop("`ccc_obj` must be the output of `ccc_diff()` or `ccc_enrich()`.")
  }
  if (test_type == "chisq") {
    changed_args <- args_differ_from_defaults('ha', 'c_linear', 'c_logistic', 'c_hurdle')
    if (length(changed_args) > 0L) {
      message(paste0("`", paste(changed_args, collapse = "`"), "`", " ignored for Wald Chisq tests. Perform two-sided tests of simple null hypothesis."))
    }
  }
  if (test_type == "z") {
    if (ha == "greater.abs" || ha == "less.abs") {
      if (c_linear < 0 || c_logistic < 0 || c_hurdle < 0) {
        stop("`c_linear`, `c_logistic`, and `c_hurdle` must be non-negative for Z tests of composite null hypotheses.")
      }
    }
  }
  
  dt.est.all <- as.data.table(ccc_obj$estimate)
  unique_levels <- colnames(contrast)
  var_names <- paste0(var_name, unique_levels)
  if (verbose) {
    p <- progressr::progressor(along = seq_len(nrow(dt.est.all)))
    # message("Starting statistical analysis...")
  }
  
  wald.test <- function(coef_l_lm, coef_r_lm, coef_l_logm, coef_r_logm,
                   vcov_l_lm, vcov_r_lm, vcov_l_logm, vcov_r_logm) {
    test.linear <- FALSE
    test.logistic <- FALSE
    if (!is.null(coef_l_lm) && !is.null(coef_r_lm)) {
      coef_l_lm <- coef_l_lm[var_names]
      coef_r_lm <- coef_r_lm[var_names]
      vcov_l_lm <- vcov_l_lm[var_names, var_names]
      vcov_r_lm <- vcov_r_lm[var_names, var_names]
      test.linear <- TRUE
    }
    if (!is.null(coef_l_logm) && !is.null(coef_r_logm)) {
      coef_l_logm <- coef_l_logm[var_names]
      coef_r_logm <- coef_r_logm[var_names]
      vcov_l_logm <- vcov_l_logm[var_names, var_names]
      vcov_r_logm <- vcov_r_logm[var_names, var_names]
      test.logistic <- TRUE
    }
    dt.test <- data.table()
    
    if (isTRUE(test.linear)) {
      if (score == "sum") {
        # EXACT COMPUTATION: Linear combination
        effect_size_linear <- contrast %*% (coef_l_lm + coef_r_lm)
        # Exact variance for linear combination (no delta method needed)
        cov_effect_size_linear <- contrast %*% (vcov_l_lm + vcov_r_lm) %*% t(contrast)
      } else if (score == "product") {
        # EXACT COMPUTATION: Product of independent coefficients
        effect_size_linear <- contrast %*% (coef_l_lm * coef_r_lm)
        # Exact product covariance using outer products (including Hadamard term)
        product_cov_linear <- (coef_l_lm %o% coef_l_lm) * vcov_r_lm + 
                              (coef_r_lm %o% coef_r_lm) * vcov_l_lm + 
                              vcov_l_lm * vcov_r_lm
        # Apply contrast exactly
        cov_effect_size_linear <- contrast %*% product_cov_linear %*% t(contrast)
      }
      
      if (test_type == "chisq") {
        test_stat_linear <- t(effect_size_linear) %*% chol2inv(chol(cov_effect_size_linear)) %*% effect_size_linear
        pvalue_linear <- pchisq(test_stat_linear, df = nrow(contrast), lower.tail = FALSE)
      } else if (test_type == "z") {
        theta <- as.numeric(effect_size_linear)
        se <- as.numeric(sqrt(cov_effect_size_linear))
        c <- c_linear
        pvalue_linear <- switch(ha,
          "greater" = pnorm(as.numeric((theta - c)/se), lower.tail = FALSE),
          "less" = pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
          "greater.abs" = min(1, 2 * pnorm(as.numeric((abs(theta) - c)/se), lower.tail = FALSE)),
          "less.abs" = max(pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
          pnorm(as.numeric((theta + c)/se), lower.tail = FALSE))
        )
      }
      
      effect_size_linear <- as.numeric(effect_size_linear)
      dt.test[, c("effect_size_linear", "pvalue_linear") := list(list(effect_size_linear), pvalue_linear)]
      dt.test[, paste0("effect_size_linear_", seq_len(nrow(contrast))) := transpose(effect_size_linear)]
      
      # Calculate standard errors for linear effect size
      std_err_linear <- sqrt(diag(cov_effect_size_linear))
      dt.test[, paste0("std_err_linear_", seq_len(nrow(contrast))) := as.list(std_err_linear)]
      
      # Get column names
      cols_all <- names(dt.test)
      
      # Find the index of col0
      col0_index <- which(cols_all == "effect_size_linear")
      dt.test[, effect_size_linear := NULL]
      new_cols <- paste0("effect_size_linear_", seq_len(nrow(contrast)))
      std_err_cols <- paste0("std_err_linear_", seq_len(nrow(contrast)))
      cols_new <- names(dt.test)
      cols_reordered <- append(cols_new[!cols_new %in% c(new_cols, std_err_cols)], 
                             c(new_cols, std_err_cols), 
                             after = col0_index - 1)
      
      # Apply new column order
      setcolorder(dt.test, cols_reordered)  
    }
    
    if (isTRUE(test.logistic)) {
      if (score == "sum") {
        # EXACT COMPUTATION: Linear combination of log-odds
        effect_size_logistic <- contrast %*% (coef_l_logm + coef_r_logm)
        # Exact variance for linear combination (no delta method needed)
        cov_effect_size_logistic <- contrast %*% (vcov_l_logm + vcov_r_logm) %*% t(contrast)
      } else if (score == "product") {
        # HYBRID COMPUTATION: Delta method for probabilities, then exact for products
        # Step 1: Compute probability means and covariances using delta method
        p_l_mean <- plogis(coef_l_logm)
        p_r_mean <- plogis(coef_r_logm)
        
        # Delta method for probability covariances
        grad_l <- diag(dlogis(coef_l_logm))
        grad_r <- diag(dlogis(coef_r_logm))
        p_l_cov <- grad_l %*% vcov_l_logm %*% grad_l
        p_r_cov <- grad_r %*% vcov_r_logm %*% grad_r
        
        # Step 2: Exact computation for probability products using outer products
        effect_size_logistic <- contrast %*% (p_l_mean * p_r_mean)
        product_cov_logistic <- (p_l_mean %o% p_l_mean) * p_r_cov + 
                               (p_r_mean %o% p_r_mean) * p_l_cov + 
                               p_l_cov * p_r_cov
        # Apply contrast exactly
        cov_effect_size_logistic <- contrast %*% product_cov_logistic %*% t(contrast)
      }
      
      if (test_type == "chisq") {
        test_stat_logistic <- t(effect_size_logistic) %*% chol2inv(chol(cov_effect_size_logistic)) %*% effect_size_logistic
        pvalue_logistic <- pchisq(test_stat_logistic, df = nrow(contrast), lower.tail = FALSE)
      } else if (test_type == "z") {
        theta <- as.numeric(effect_size_logistic)
        se <- as.numeric(sqrt(cov_effect_size_logistic))
        c <- c_logistic
        pvalue_logistic <- switch(ha,
                                  "greater" = pnorm(as.numeric((theta - c)/se), lower.tail = FALSE),
                                  "less" = pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
                                  "greater.abs" = min(1, 2 * pnorm(as.numeric((abs(theta) - c)/se), lower.tail = FALSE)),
                                  "less.abs" = max(pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
                                                   pnorm(as.numeric((theta + c)/se), lower.tail = FALSE))
        )
      }
      
      effect_size_logistic <- as.numeric(effect_size_logistic)
      dt.test[, c("effect_size_logistic", "pvalue_logistic") := list(list(effect_size_logistic), pvalue_logistic)]
      dt.test[, paste0("effect_size_logistic_", seq_len(nrow(contrast))) := transpose(effect_size_logistic)]
      
      # Calculate standard errors for logistic effect size
      std_err_logistic <- sqrt(diag(cov_effect_size_logistic))
      dt.test[, paste0("std_err_logistic_", seq_len(nrow(contrast))) := as.list(std_err_logistic)]
      
      # Get column names
      cols_all <- names(dt.test)
      
      # Find the index of col0
      col0_index <- which(cols_all == "effect_size_logistic")
      dt.test[, effect_size_logistic := NULL]
      new_cols <- paste0("effect_size_logistic_", seq_len(nrow(contrast)))
      std_err_cols <- paste0("std_err_logistic_", seq_len(nrow(contrast)))
      cols_new <- names(dt.test)
      cols_reordered <- append(cols_new[!cols_new %in% c(new_cols, std_err_cols)], 
                             c(new_cols, std_err_cols), 
                             after = col0_index - 1)
      
      # Apply new column order
      setcolorder(dt.test, cols_reordered)  
    }
    
    if (isTRUE(test.linear) && isTRUE(test.logistic)) {
      if (score == "sum") {
        # HYBRID COMPUTATION: Exact for p*μ products, then exact sum
        # Get probability components (already computed above if logistic was processed)
        if (!exists("p_l_mean")) {
          p_l_mean <- plogis(coef_l_logm)
          p_r_mean <- plogis(coef_r_logm)
          grad_l <- diag(dlogis(coef_l_logm))
          grad_r <- diag(dlogis(coef_r_logm))
          p_l_cov <- grad_l %*% vcov_l_logm %*% grad_l
          p_r_cov <- grad_r %*% vcov_r_logm %*% grad_r
        }
        
        # Exact computation for p*μ products
        term_l_mean <- p_l_mean * coef_l_lm  # E[p_l * μ_l]
        term_r_mean <- p_r_mean * coef_r_lm  # E[p_r * μ_r]
        
        # Exact covariance for p*μ products using outer products (including Hadamard terms)
        term_l_cov <- (p_l_mean %o% p_l_mean) * vcov_l_lm + 
                      (coef_l_lm %o% coef_l_lm) * p_l_cov + 
                      vcov_l_lm * p_l_cov
        term_r_cov <- (p_r_mean %o% p_r_mean) * vcov_r_lm + 
                      (coef_r_lm %o% coef_r_lm) * p_r_cov + 
                      vcov_r_lm * p_r_cov
        
        # Exact computation for sum
        effect_size_hurdle <- contrast %*% (term_l_mean + term_r_mean)
        cov_effect_size_hurdle <- contrast %*% (term_l_cov + term_r_cov) %*% t(contrast)
      } else if (score == "product") {
        # HYBRID COMPUTATION: Build up from components
        # Linear product component (already computed if linear was processed)
        if (!exists("product_cov_linear")) {
          linear_product_mean <- coef_l_lm * coef_r_lm
          product_cov_linear <- (coef_l_lm %o% coef_l_lm) * vcov_r_lm + 
                               (coef_r_lm %o% coef_r_lm) * vcov_l_lm + 
                               vcov_l_lm * vcov_r_lm
        } else {
          linear_product_mean <- coef_l_lm * coef_r_lm
        }
        
        # Probability product component (already computed if logistic was processed)
        if (!exists("product_cov_logistic")) {
          p_l_mean <- plogis(coef_l_logm)
          p_r_mean <- plogis(coef_r_logm)
          grad_l <- diag(dlogis(coef_l_logm))
          grad_r <- diag(dlogis(coef_r_logm))
          p_l_cov <- grad_l %*% vcov_l_logm %*% grad_l
          p_r_cov <- grad_r %*% vcov_r_logm %*% grad_r
          prob_product_mean <- p_l_mean * p_r_mean
          product_cov_logistic <- (p_l_mean %o% p_l_mean) * p_r_cov + 
                                 (p_r_mean %o% p_r_mean) * p_l_cov + 
                                 p_l_cov * p_r_cov
        } else {
          prob_product_mean <- p_l_mean * p_r_mean
        }
        
        # Exact computation for hurdle = (linear product) * (probability product)
        effect_size_hurdle <- contrast %*% (linear_product_mean * prob_product_mean)
        hurdle_product_cov <- (linear_product_mean %o% linear_product_mean) * product_cov_logistic + 
                             (prob_product_mean %o% prob_product_mean) * product_cov_linear + 
                             product_cov_linear * product_cov_logistic
        cov_effect_size_hurdle <- contrast %*% hurdle_product_cov %*% t(contrast)
      }
      
      if (test_type == "chisq") {
        test_stat_hurdle <- t(effect_size_hurdle) %*% chol2inv(chol(cov_effect_size_hurdle)) %*% effect_size_hurdle
        pvalue_hurdle <- pchisq(test_stat_hurdle, df = nrow(contrast), lower.tail = FALSE)
        
        # 2-part p-values only available for chisq test
        if (exists("test_stat_linear") && exists("test_stat_logistic")) {
          test_stat_2part <- test_stat_linear + test_stat_logistic
          pvalue_2part <- pchisq(test_stat_2part, df = nrow(contrast) * 2L, lower.tail = FALSE)
        } else {
          pvalue_2part <- NA
        }
      } else if (test_type == "z") {
        theta <- as.numeric(effect_size_hurdle)
        se <- as.numeric(sqrt(cov_effect_size_hurdle))
        c <- c_hurdle
        pvalue_hurdle <- switch(ha,
                                "greater" = pnorm(as.numeric((theta - c)/se), lower.tail = FALSE),
                                "less" = pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
                                "greater.abs" = min(1, 2 * pnorm(as.numeric((abs(theta) - c)/se), lower.tail = FALSE)),
                                "less.abs" = max(pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
                                                 pnorm(as.numeric((theta + c)/se), lower.tail = FALSE))
        )
        pvalue_2part <- NA  # Not applicable for z-test
      }
      
      # Combine p-values using existing methods
      pvalues_to_combine <- c()
      if (exists("pvalue_linear")) pvalues_to_combine <- c(pvalues_to_combine, pvalue_linear)
      if (exists("pvalue_logistic")) pvalues_to_combine <- c(pvalues_to_combine, pvalue_logistic)
      
      if (length(pvalues_to_combine) == 2) {
        pvalue_stouffer <- stouffer_combine_pvalues(pvalues_to_combine)
        pvalue_fisher <- fisher_combine_pvalues(pvalues_to_combine)
      } else {
        pvalue_stouffer <- NA
        pvalue_fisher <- NA
      }
      
      effect_size_hurdle <- as.numeric(effect_size_hurdle)
      dt.test[, c("effect_size_hurdle", "pvalue_hurdle", "pvalue_2part", "pvalue_stouffer", "pvalue_fisher") := 
                list(list(effect_size_hurdle), pvalue_hurdle, pvalue_2part, pvalue_stouffer, pvalue_fisher)]
      dt.test[, paste0("effect_size_hurdle_", seq_len(nrow(contrast))) := transpose(effect_size_hurdle)]
      
      # Calculate standard errors for hurdle effect size
      std_err_hurdle <- sqrt(diag(cov_effect_size_hurdle))
      dt.test[, paste0("std_err_hurdle_", seq_len(nrow(contrast))) := as.list(std_err_hurdle)]
      
      # Get column names
      cols_all <- names(dt.test)
      
      # Find the index of col0
      col0_index <- which(cols_all == "effect_size_hurdle")
      dt.test[, effect_size_hurdle := NULL]
      new_cols <- paste0("effect_size_hurdle_", seq_len(nrow(contrast)))
      std_err_cols <- paste0("std_err_hurdle_", seq_len(nrow(contrast)))
      cols_new <- names(dt.test)
      cols_reordered <- append(cols_new[!cols_new %in% c(new_cols, std_err_cols)], 
                             c(new_cols, std_err_cols), 
                             after = col0_index - 1)
      
      # Apply new column order
      setcolorder(dt.test, cols_reordered)  
    }
    
    if (verbose) {
      p()
    }
    dt.test
  }
  
  future_mapply(FUN = wald.test,
         coef_l_lm = dt.est.all$coef_l_lm, coef_r_lm = dt.est.all$coef_r_lm,
         coef_l_logm = dt.est.all$coef_l_logm, coef_r_logm = dt.est.all$coef_r_logm,
         vcov_l_lm = dt.est.all$vcov_l_lm, vcov_r_lm = dt.est.all$vcov_r_lm,
         vcov_l_logm = dt.est.all$vcov_l_logm, vcov_r_logm = dt.est.all$vcov_r_logm,
         SIMPLIFY = FALSE, future.chunk.size = chunk_size
  ) -> test.list
  
  cols <- grep("^(coef_|vcov_)", colnames(dt.est.all), value = TRUE)
  dt.test.all <- dt.est.all[, (cols) := NULL]
  test.results <- rbindlist(test.list, fill = TRUE)
  dt.test.all <- dt.test.all[, names(test.results) := test.results]
  
  pval_cols <- grep("^pvalue", names(dt.test.all), value = TRUE)
  if (cell_type_padj) {
    dt.test.all[, (paste0("padj_", sub("pvalue\\_", "", pval_cols))) :=
                  lapply(.SD, function(p) p.adjust(p, method = padj_method)),
                by = .(sender, receiver), .SDcols = pval_cols
    ]
  } else {
    dt.test.all[, (paste0("padj_", sub("pvalue\\_", "", pval_cols))) :=
                  lapply(.SD, function(p) p.adjust(p, method = padj_method)),
                .SDcols = pval_cols
    ]
  }
  ordered_cols <- names(dt.test.all)
  for (pval in rev(pval_cols)) {
    padj <- sub("pvalue\\_", "padj_", pval)
    if (padj %in% names(dt.test.all)) {
      idx <- match(pval, ordered_cols)
      ordered_cols <- append(ordered_cols[ordered_cols != padj], padj, after = idx)
    }
  }
  setcolorder(dt.test.all, ordered_cols)
  
  as.data.frame(dt.test.all)
}


