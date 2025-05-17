#' Identify Differential or Enriched Cell-Cell Communication
#' 
#' Perform Wald tests on the product of ligand expression levels in sender and receptor expression levels in receiver for differential or enriched cell-cell communication.
#' 
#' 

ccc_test <- function(ccc_obj, contrast = NULL, test_type = NULL, ha = NULL,
                     c_linear = 0, c_logisitc = 0, c_hurdel = 0,
                     verbose = TRUE, padj_method = "BH", cell_type_padj = TRUE) {
  if (!is.null(test_type)) {
    test_type <- match.arg(test_type, choices = c("chisq", "z"))
  }
  if (is.null(ha)) {
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
      stop("'contrast' must be either a named vector or a matrix with column names.")
    }
    if (is.null(colnames(contrast))) {
      stop("'contrast' must be either a named vector or a matrix with column names.")
    }
    if (qr(contrast)$rank < nrow(contrast)) {
      stop("'contrast' must be full row rank")
    }
    if (is.null(test_type)) {
      test_type <- "chisq"
      changed_args <- args_differ_from_defaults('ha', 'c_linear', 'c_logisitc', 'c_hurdel')
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
    colnames(contrast) <- c("classtarget", "classbackground")
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
    changed_args <- args_differ_from_defaults('ha', 'c_linear', 'c_logisitc', 'c_hurdel')
    if (length(changed_args) > 0L) {
      message(paste0("`", paste(changed_args, collapse = "`"), "`", " ignored for Wald Chisq tests. Perform two-sided tests of simple null hypothesis."))
    }
  }
  if (test_type == "z") {
    if (ha == "greater.abs" || ha == "less.abs") {
      if (c_linear < 0 || c_logisitc < 0 || c_hurdel < 0) {
        stop("`c_linear`, `c_logisitc`, and `c_hurdel` must be non-negative for Z tests of composite null hypotheses.")
      }
    }
  }
  
  dt.test.all <- as.data.table(ccc_obj$estimate)
  unique_levels <- colnames(contrast)
  var_names <- paste0(var_name, unique_levels)
  
  wald.test <- function(coef_l_lm, coef_r_lm, coef_llogm, coef_r_logm,
                   vcov_l_lm, vcov_r_lm, vcov_l_logm, vcov_r_logm) {
    test.linear <- FALSE
    test.logistic <- FALSE
    if (!is.null(coef_l_lm) && !is.null(coef_r_lm)) {
      test.linear <- TRUE
    }
    if (!is.null(coef_l_logm) && !is.null(coef_r_logm)) {
      test.logistic <- TRUE
    }
    dt.test <- data.table()
    
    if (isTRUE(test.linear)) {
      effect_size_linear <- contrast %*% (coef_l_lm * coef_r_lm)
      gradient_matrix_linear <- matrix(0, nrow = nrow(contrast), ncol = length(var_names) * 2L)
      
      for (i in seq_len(nrow(contrast))) {
        for (j in seq_along(var_names)) {
          var_level <- var_names[j]
          level <- unique_levels[j]
          
          ind_l_lm <- which(names(coef_l_lm) == var_level)
          ind_r_lm <- which(names(coef_r_lm) == var_level) + length(var_names)
          
          gradient_matrix_linear[i, ind_l_lm] <- contrast[i, level] * coef_r_lm[var_level]
          gradient_matrix_linear[i, ind_r_lm] <- contrast[i, level] * coef_l_lm[var_level]
        }
      }
      
      vcov_linear <- Matrix::bdiag(vcov_l_lm, vcov_r_lm)
      cov_effect_size_linear <- gradient_matrix_linear %*% as.matrix(vcov_linear) %*% t(gradient_matrix_linear)
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
      # Get column names
      cols_all <- names(dt.test)
      
      # Find the index of col0
      col0_index <- which(cols_all == "effect_size_linear")
      dt.test[, effect_size_linear := NULL]
      new_cols <- paste0("effect_size_linear_", seq_len(nrow(contrast)))
      cols_new <- names(dt.test)
      cols_reordered <- append(cols_new[cols_new != new_cols], new_cols, after = col0_index - 1)
      
      # Apply new column order
      setcolorder(dt.test, cols_reordered)  
    }
    
    if (isTRUE(test.logistic)) {
      effect_size_logistic <- contrast %*% (plogis(coef_l_logm) * plogis(coef_r_logm))
      gradient_matrix_logistic <- matrix(0, nrow = nrow(contrast), ncol = length(var_names) * 2L)
      
      for (i in seq_len(nrow(contrast))) {
        for (j in seq_along(var_names)) {
          var_level <- var_names[j]
          level <- unique_levels[j]
          
          ind_l_logm <- which(names(coef_l_logm) == var_level)
          ind_r_logm <- which(names(coef_r_logm) == var_level) + length(var_names)
          
          gradient_matrix_logistic[i, ind_l_logm] <- contrast[i, level] * plogis(coef_r_logm[var_level]) * dlogis(coef_l_logm[var_level])
          gradient_matrix_logistic[i, ind_r_logm] <- contrast[i, level] * plogis(coef_l_logm[var_level]) * dlogis(coef_r_logm[var_level])
        }
      }
      
      vcov_logistic <- Matrix::bdiag(vcov_l_logm, vcov_r_logm)
      cov_effect_size_logistic <- gradient_matrix_logistic %*% as.matrix(vcov_logistic) %*% t(gradient_matrix_logistic)
      
      if (test_type == "chisq") {
        test_stat_logistic <- t(effect_size_logistic) %*% chol2inv(chol(cov_effect_size_logistic)) %*% effect_size_logistic
        pvalue_logistic <- pchisq(test_stat_logistic, df = nrow(contrast), lower.tail = FALSE)
      } else if (test_type == "z") {
        theta <- as.numeric(effect_size_logistic)
        se <- as.numeric(sqrt(cov_effect_size_logistic))
        c <- c_logisitc
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
      # Get column names
      cols_all <- names(dt.test)
      
      # Find the index of col0
      col0_index <- which(cols_all == "effect_size_logistic")
      dt.test[, effect_size_logistic := NULL]
      new_cols <- paste0("effect_size_logistic_", seq_len(nrow(contrast)))
      cols_new <- names(dt.test)
      cols_reordered <- append(cols_new[cols_new != new_cols], new_cols, after = col0_index - 1)
      
      # Apply new column order
      setcolorder(dt.test, cols_reordered)  
    }
    
    if (isTRUE(test.linear) && isTRUE(test.logistic)) {
      effect_size_hurdle <- contrast %*% (coef_l_lm * coef_r_lm * plogis(coef_l_logm) * plogis(coef_r_logm))
      gradient_matrix_hurdle <- matrix(0, nrow = nrow(contrast), ncol = length(var_names) * 4L)
      
      for (i in seq_len(nrow(contrast))) {
        for (j in seq_along(var_names)) {
          var_level <- var_names[j]
          level <- unique_levels[j]
          
          ind_l_lm <- which(names(coef_l_lm) == var_level)
          ind_l_logm <- which(names(coef_l_logm) == var_level) + length(var_names)
          ind_r_lm <- which(names(coef_r_lm) == var_level) + length(var_names) * 2L
          ind_r_logm <- which(names(coef_r_logm) == var_level) + length(var_names) * 3L
          
          gradient_matrix_hurdle[i, ind_l_lm] <- contrast[i, level] * coef_r_lm[var_level] * plogis(coef_l_logm[var_level]) * plogis(coef_r_logm[var_level])
          gradient_matrix_hurdle[i, ind_l_logm] <- contrast[i, level] * coef_l_lm[var_level] * coef_r_lm[var_level] * plogis(coef_r_logm[var_level]) * dlogis(coef_l_logm[var_level])
          gradient_matrix_hurdle[i, ind_r_lm] <- contrast[i, level] * coef_l_lm[var_level] * plogis(coef_l_logm[var_level]) * plogis(coef_r_logm[var_level])
          gradient_matrix_hurdle[i, ind_r_logm] <- contrast[i, level] * coef_l_lm[var_level] * coef_r_lm[var_level] * plogis(coef_l_logm[var_level]) * dlogis(coef_r_logm[var_level])
        }
      }
      vcov_hurdle <- Matrix::bdiag(vcov_l_lm, vcov_l_logm, vcov_r_lm, vcov_r_logm)
      cov_effect_size_hurdle <- gradient_matrix_hurdle %*% as.matrix(vcov_hurdle) %*% t(gradient_matrix_hurdle)
      
      if (test_type == "chisq") {
        test_stat_hurdle <- t(effect_size_hurdle) %*% chol2inv(chol(cov_effect_size_hurdle)) %*% effect_size_hurdle
        pvalue_hurdle <- pchisq(test_stat_hurdle, df = nrow(contrast), lower.tail = FALSE)
      } else if (test_type == "z") {
        theta <- as.numeric(effect_size_hurdle)
        se <- as.numeric(sqrt(cov_effect_size_hurdle))
        c <- c_hurdel
        pvalue_hurdle <- switch(ha,
                                "greater" = pnorm(as.numeric((theta - c)/se), lower.tail = FALSE),
                                "less" = pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
                                "greater.abs" = min(1, 2 * pnorm(as.numeric((abs(theta) - c)/se), lower.tail = FALSE)),
                                "less.abs" = max(pnorm(as.numeric((theta - c)/se), lower.tail = TRUE),
                                                 pnorm(as.numeric((theta + c)/se), lower.tail = FALSE))
        )
      }
      if (test_type == "chisq") {
        # 2-part p-values only available for chisq test
        test_stat_2part <- test_stat_linear + test_stat_logistic
        pvalue_2part <- pchisq(test_stat_2part, df = nrow(contrast) * 2L, lower.tail = FALSE)
        dt.test[, pvalue_2part := pvalue_2part]
      }
      pvalue_stouffer <- stouffer_combine_pvalues(c(pvalue_linear, pvalue_logistic))
      pvalue_fisher <- fisher_combine_pvalues(c(pvalue_linear, pvalue_logistic))
      
      effect_size_hurdle <- as.numeric(effect_size_hurdle)
      dt.test[, c("effect_size_hurdle", "pvalue_hurdle", "pvalue_stouffer", "pvalue_fisher") := 
                list(list(effect_size_hurdle), pvalue_hurdle, pvalue_stouffer, pvalue_fisher)]
      dt.test[, paste0("effect_size_hurdle_", seq_len(nrow(contrast))) := transpose(effect_size_hurdle)]
      # Get column names
      cols_all <- names(dt.test)
      
      # Find the index of col0
      col0_index <- which(cols_all == "effect_size_hurdle")
      dt.test[, effect_size_hurdle := NULL]
      new_cols <- paste0("effect_size_hurdle_", seq_len(nrow(contrast)))
      cols_new <- names(dt.test)
      cols_reordered <- append(cols_new[cols_new != new_cols], new_cols, after = col0_index - 1)
      
      # Apply new column order
      setcolorder(dt.test, cols_reordered)  
    }
    
    dt.test
  }
  
  mapply(FUN = wald.test,
         coef_l_lm = dt.test.all$coef_l_lm, coef_r_lm = dt.test.all$coef_r_lm,
         coef_l_logm = dt.test.all$coef_l_logm, coef_r_logm = dt.test.all$coef_r_logm,
         vcov_l_lm = dt.test.all$vcov_l_lm, vcov_r_lm = dt.test.all$vcov_r_lm,
         vcov_l_logm = dt.test.all$vcov_l_logm, vcov_r_logm = dt.test.all$vcov_r_logm,
         SIMPLIFY = FALSE
  ) -> test.list
  
  test.results <- rbindlist(test.list)
  dt.test.all <- dt.test.all[, names(test.results) := test.reslts]
  
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

