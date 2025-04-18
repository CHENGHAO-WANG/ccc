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
library(assertthat)

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
  on.exit(setDTthreads(old_nthreads), add = TRUE)
  
  assertthat::assert_that(assertthat::is.flag(verbose))
  
  if (verbose) {
    if (!progressr::handlers(global = NA) && interactive()) {
      # If no progressr bar settings are configured, then set cli as the default.
      progressr::handlers(global = TRUE)
      ohandlers <- progressr::handlers('cli')
      on.exit(progressr::handlers(ohandlers), add = TRUE)
      on.exit(progressr::handlers(global = FALSE), add = TRUE)
      message(
        "Info: No global progress bars were found; the cli handler has been enabled. ",
        "See `vignette('ccc-intro')` for how to customize the progress bar settings."
      )
    }
  }
  
  assertthat::assert_that(assertthat::is.flag(cdr))
  assertthat::assert_that(assertthat::is.flag(lmm_re))
  assertthat::assert_that(assertthat::is.flag(logmm_re))
  assertthat::assert_that(assertthat::is.flag(sandwich))
  assertthat::assert_that(assertthat::is.flag(sep_detection))
  assertthat::assert_that(assertthat::is.flag(cell_type_padj))
  
  required_args <- c('expression_matrix', 'metadata', 'contrast')
  passed_args <- names(match.call())[-1]
  missing_args <- setdiff(required_args, passed_args)
  if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
  }
  
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
  
  # contrast <- contrast[, colSums(abs(contrast) > 0)]
  if(any(is.na(expression_matrix))) {
    stop("Missing values not allowed in expression matrix, please remove NA values.")
  }
  if(any(is.na(metadata))) {
    stop("Missing values not allowed in metadata, please remove NA values.")
  }
  # if(any(is.na(covar_col))) {
  #   stop("Missing values not allowed in covar_col, please remove NA values.")
  # }
  
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
  
  if ("y" %in% colnames(metadata)) {
    stop("\"y\" in column names of 'metadata'. Please rename it.")
  }
  if ("z" %in% colnames(metadata)) {
    stop("\"z\" in column names of 'metadata'. Please rename it.")
  }
  
  covar_exists <- covar_col %in% colnames(metadata)
  if(!is.null(covar_col) && !all(covar_exists)){
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
  
  if (!all(colnames(contrast) %in% unique(metadata[[group_col]]))) {
    stop(paste0("'contrast' contains group levels that are not present in the \"",group_col,"\" column of 'metadata'."))
  }
  metadata <- metadata[metadata[[group_col]] %in% colnames(contrast), ]
  
  if (isTRUE(lmm_re) || isTRUE(logmm_re)) {
    if (is.null(id_col)) {
      stop("'id_col' is not specified, while 'lmm_re' or 'logmm_re' is TRUE")
    }
    if (isFALSE(id_col %in% colnames(metadata))) {
      stop(paste0("'id_col' \"",id_col, "\" does not exist in 'metadata'."))
    }
  }
  if (isFALSE(lmm_re) && isFALSE(logmm_re) && !is.null(id_col)) {
    warning("'id_col' is not NULL. This input will be ignored, because 'lmm_re' and 'logmm_re' are FALSE")
    id_col <- NULL
  }
  
  if (verbose) {
    if (is.null(sender)) {
      message("'sender' is not specified. All cell types will be considered as potential senders in the analysis.")
    }
    if (is.null(receiver)) {
      message("'receiver' is not specified. All cell types will be considered as potential receivers in the analysis.")
    }
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
  
  ### lr table
  lr_table <- prep_lr(lr = lr)
  
  ### rename metadata
  metadata <- rename_metadata(metadata = metadata, cell_id_col = cell_id_col,
                              id_col = id_col, group_col = group_col, cell_type_col = cell_type_col)
  
  ### filter cell type
  filtered_obj <- filter_cell_type(metadata = metadata, sender = sender, receiver = receiver,
                                   min_cell = min_cell, contrast = contrast)
  metadata_subset <- filtered_obj$metadata_subset
  sender <- filtered_obj$sender
  receiver <- filtered_obj$receiver
  
  rm(metadata, filtered_obj)
  gc()
  
  ### compute cdr if specified by the user
  if (isTRUE(cdr)) {
    metadata_subset <- compute_cdr(expression_matrix = expression_matrix,
                                   metadata_subset = metadata_subset, threshold = threshold)
  }
  
  ### filter lr; fit models; conduct tests
  
  # Create all possible combinations of sender and receiver
  sender_receiver_combinations <- expand.grid(sender = sender, receiver = receiver)
  # Merge with lr_table to get all possible combinations
  pairs4analysis <- base::merge(sender_receiver_combinations, lr_table, by = NULL)
  
  # Convert to data.table
  setDT(pairs4analysis)
  npairs <- nrow(pairs4analysis)
  
  unique_ids <- unique(metadata_subset[,id])
  
  i_s <- seq(1L, nrow(pairs4analysis), by = chunk_size)
  p <- progressr::progressor(along = i_s)
  
  run_analysis <- function(i) {
    chunk <- pairs4analysis[i:min(i + chunk_size - 1L, npairs), ]
    p()
    
    results.summary <- results.test <- results.error <- results.warning <- list()
    
    for (j in 1L:nrow(chunk)) {
      sender <- chunk$sender[j]
      ligand <- chunk$ligand[j]
      receiver <- chunk$receiver[j]
      receptor <- chunk$receptor[j]
      
      # Create copies of metadata_subset for sender and receiver
      data_sender_ligand <- metadata_subset[cell_type == sender]
      data_receiver_receptor <- metadata_subset[cell_type == receiver]
      
      # Handle single gene or multi-gene ligands and receptors
      ligand_genes <- unlist(strsplit(ligand, "_"))
      receptor_genes <- unlist(strsplit(receptor, "_"))
      
      # Subset expression matrix for ligands and receptors
      ligand_expr_values <- expression_matrix[ligand_genes, data_sender_ligand$cell_id, drop = TRUE]
      receptor_expr_values <- expression_matrix[receptor_genes, data_receiver_receptor$cell_id, drop = TRUE]
      
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
      
      # results.j[ligand.mean := mean(data_sender_ligand$y)]
      # positive_expression <- data_sender_ligand_1$y
      # results.j[ligand.positive_mean := ifelse(length(data_sender_ligand_1$y) > 0, mean(positive_expression))]
      # results.j[receptor.mean := mean(data_receiver_receptor$y)]
      
      # descriptive statistics summary
      dt.summary.ligand <- compute_group_stats(dt = data_sender_ligand, prefix = "ligand.")
      dt.summary.receptor <- compute_group_stats(dt = data_receiver_receptor, prefix = "receptor.")
      dt.summary <- data.table::merge(dt.summary.ligand, dt.summary.receptor, by = 'group')
      dt.summary[ , c("sender", "receiver", "ligand", "receptor") := list(sender, receiver, ligand, receptor)]
      results.summary[[length(results.summary) + 1L]] <- dt.summary
      
      # Define covariates
      if (is.null(covar_col) && isFALSE(cdr)) {
        covariates <- "group"
      } else {
        covar <- c(covar_col, if(isTRUE(cdr)) "cdr")
        
        center_covar(dt = data_sender_ligand_1, covar = covar) -> covariates # The 4 data sets are different. But the column names are the same.
        center_covar(dt = data_receiver_receptor_1, covar = covar)
        center_covar(dt = data_sender_ligand, covar = covar)
        center_covar(dt = data_receiver_receptor, covar = covar)
        
        covariates <- c("group", covariates)
      }
      
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
      fit_linear <- function(data, formula, name) {
        warnings. <- list()
        fit <- tryCatch(
          withCallingHandlers({
          cond <- detect_all_zeros(dt = data, id_col = "id", id = unique_ids)
          if (cond) {
            stop("Too few cells expressing the ligand/receptor gene for fitting a linear model.")
          } else {
            if (isTRUE(lmm_re)) {
              lmer(formula, data = data, control = control_lmm)
            } else {
              lm(formula, data = data)
            }
          }
        }, warning = function(w) {
          warnings.[[length(warnings.) + 1L]] <<- w$message
          invokeRestart("muffleWarning")
        }), error = function(e) {
          error_messages[[name]] <<- e$message
          return(NULL)
        })
        if (length(warnings.) > 0) {
          warning_messages[[name]] <<- warnings.
        }
        fit
      }
      
      fit_logistic <- function(data, formula, name) {
        warnings. <- list()
        fit <- tryCatch(
          withCallingHandlers({
          cond <- isTRUE(sep_detection) && detect_re_separation(dt = data, z_col = "z", id_col = "id", num_ids = num_ids, sep_prop = sep_prop, sep_n = sep_n)
          if (cond) {
            stop("Complete or Quasi-complete separation detected.")
          } else {
            if (isTRUE(logmm_re)) {
              mixed_model(fixed = formula$fixed, random = formula$random, family = binomial(), data = data, control = control_logmm)
            } else {
              glm(formula, family = binomial(), data = data, control = control_logm)
            }
          }
        }, warning = function(w) {
          warnings.[[length(warnings.) + 1L]] <<- w$message
          invokeRestart("muffleWarning")
        }), error = function(e) {
          error_messages[[name]] <<- e$message
          return(NULL)
        })
        if (length(warnings.) > 0) {
          warning_messages[[name]] <<- warnings.
        }
        fit
      }
      
      ##
      warning_messages <- error_messages <- list()
      fit.l.linear <- fit_linear(data = data_sender_ligand_1, formula = formula_linear, name = "ligand.linear")
      fit.r.linear <- fit_linear(data = data_receiver_receptor_1, formula = formula_linear, name = "ligand.logistic")
      fit.l.logistic <- fit_logistic(data = data_sender_ligand, formula = formula_logistic, name = "receptor.linear")
      fit.r.logistic <- fit_logistic(data = data_receiver_receptor, formula = formula_logistic, name = "receptor.logistic")
      if (length(warning_messages) > 0) {
        results.warning[[paste(sender, receiver, ligand, receptor, sep = "-")]] <- warning_messages
      }
      if (length(error_messages) > 0) {
        results.error[[paste(sender, receiver, ligand, receptor, sep = "-")]] <- error_messages
      }
      
      ## error message
      # model_objects <- list(ligand.linear = fit.l.linear,
      #                       ligand.logistic = fit.l.logistic,
      #                       receptor.linear = fit.r.linear,
      #                       receptor.logistic = fit.r.logistic)
      # is_errors <- sapply(model_objects, is.character)
      # error_messages <- model_objects[is_errors]
      # if (length(error_messages) > 0) {
      #   results.error[[paste(sender, receiver, ligand, receptor, sep = "-")]] <- error_messages
      # }
      
      results.test[[length(results.test) + 1L]] <- ccc_test(fit.l.linear = fit.l.linear, fit.l.logistic = fit.l.logistic,
                                                            fit.r.linear = fit.r.linear, fit.r.logistic = fit.r.logistic,
                                                            contrast = contrast, re_lmm = re_lmm, re_logmm = re_logmm,
                                                            sandwich = sandwich, results.j = results.j)
      
    }
    #rbindlist(results.i, fill = TRUE)
    list(descriptive_stats = rbindlist(results.summary),
         test_results = rbindlist(results.test, fill = TRUE),
         errors = results.error,
         warnings = results.warning)
  }
  
  results_obj <- future_lapply(i_s, FUN = run_analysis(i))
  
  list.descriptive_stats <- lapply(results_obj, \(x) x$descriptive_stats)
  list.test_results <- lapply(results_obj, \(x) x$test_results)
  list.errors <- lapply(results_obj, \(x) x$errors)
  list.warnings <- lapply(results_obj, \(x) x$warnings)
  
  list(summary = rbindlist(list.descriptive_stats),
       test = rbindlist(list.test_results, fill = TRUE),
       errors = do.call(c, list.errors),
       warnings = do.call(c, list.warnings))
  
}







