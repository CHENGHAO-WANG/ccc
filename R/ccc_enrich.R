#' Enrichment Analysis for Cell-Cell Communication
#'
#' This function performs enrichment analysis for cell-cell communication by comparing the expression of ligand-receptor pairs
#' in target cell type pairs against background cell types. For each target cell type pair, it calculates:
#' \itemize{
#'   \item Expression metrics (mean, standard deviation, percentage above threshold)
#'   \item Statistical tests comparing target vs background expression using linear and logistic models
#'   \item Fold change between target and background
#'   \item Multiple p-value combinations (Hurdle, Stouffer's method, Fisher's method)
#' }
#'
#' @param expression_matrix a numeric matrix of normalized counts, with rows corresponding to genes and columns corresponding to cells. Both row names (gene symbols) and column names (cell identifiers) must be provided.
#' @param metadata a data frame containing cell-level metadata (e.g., cell type, group, id, covariates).
#' @param cell_id_col a character string specifying the name of the column in `metadata` that contains cell identifiers. These ID's should match the column names of `expression_matrix`. Defaults to `"cell_id"`.
#' @param cell_type_col a character string specifying the name of the column in `metadata` that contains cell type annotations. Defaults to `"cell_type"`.
#' @param id_col a character string specifying the name of the column in `metadata` that contains individual-level (sample-level) ID's. Defaults to `"id"`.
#' @param sender a character string or a character vector specifying the cell types as sender (expressing ligands in cell-cell communication). Defaults to all cell types.
#' @param receiver a character string or a character vector specifying the cell types as receiver (expressing receptors in cell-cell communication). Defaults to all cell types.
#' @param lr specifies the ligand-receptor database to use. Can be `"omnipathr"` (the default), `"ramilowski"`, or a user-supplied data frame. If `"omnipathr"` or `"ramilowski"` is provided, the corresponding built-in dataset is used. If a data frame is provided, it must contain exactly two columns named `"ligand"` and `"receptor"`, with gene symbols as entries. For multi-subunit ligands/receptors, the gene symbols of all subunits must be joined by `_`. (e.g., `"CLCF1_CRLF1"` for a ligand composed of CLCF1 and CRLF1)
#' @param multi_sub a character string specifying how to handle multi-subunit ligands/receptors.
#'  \itemize{
#'    \item \dQuote{\code{minimum}}: (the default) the expression level for each cell is defined as the minimum expression across all subunit genes.
#'    \item \dQuote{\code{arithmetic_mean}}: the expression level for each cell is defined as the arithmetic mean expression across all subunit genes.
#'    \item \dQuote{\code{geometric_mean}}: the expression level for each cell is defined as the geometric mean expression across all subunit genes.
#'    \item \dQuote{\code{min_avg_gene}}: the subunit gene with the minimum average expression is selected.
#'    \item \dQuote{\code{min_rate_gene}}: the subunit gene with the minimum expression rate is selected. The expression rate for a gene is calculated as the number of cells with expression level above `threshold` divided by the total number of cells.
#'  }
#' @param verbose logical scalar. If `TRUE` (the default), display a progress bar. The default handler is "cli". This package uses the \pkg{progressr} framework for progress reporting, so users can customize the progress bar. See [progressr::handlers()] for customizing progress bar behavior.
#' @param min_cell integer scalar. Filter out cell types with fewer than `min_cell` cells. Defaults to 10.
#' @param min_pct numeric scalar. Only test ligand-receptor pairs that are expressed above `threshold` in a minimum fraction of `min_pct` cells for `large_n` individuals/samples in sender and receiver cell types respectively. Defaults to 0.01.
#' @param large_n integer scalar. Number of individuals/samples that are considered to be "large". Defaults to 2.
#' @param min_total_pct numeric scalar. Only test ligand-receptor pairs that are detected (expression level above `threshold`) in a minimum fraction of `min_total_pct` cells across all individuals/samples in sender and receiver cell types respectively. Defaults to 0.
#' @param threshold numeric scalar. A gene is considered expressed in a cell if its expression level is greater than `threshold`. Defaults to 0.
#' @param sep_detection logical scalar. If `TRUE` (the default), detect complete or quasi-complete separation in logistic models.
#' @param sep_prop numeric scalar. For each ligand/receptor gene, if it is expressed above/below `threshold` in all cells of more than a `sep_prop` fraction of individuals/samples, this is considered complete or quasi-complete separation and the logistic model for that gene is skipped.
#' @param sep_n numeric scalar. For each ligand/receptor gene, if it is expressed above/below `threshold` in all cells of more than `sep_n` individuals/samples, this is considered complete or quasi-complete separation and the logistic model for that gene is skipped.
#' @param padj_method a character string for multiple testing correction method. This is passed to [stats::p.adjust()]. Defaults to "BH".
#' @param cell_type_padj logical scalar. If `TRUE` (the default), adjust p-values for each sender-receiver pair.
#' @param chunk_size integer scalar. The number of communication pairs (each defined by a distinct combination of sender, receiver, ligand, and receptor) to be sent to each parallel environment. Defaults to 10. To enable parallelization, users should use the \pkg{future} package.
#' @param lmm_re logical scalar. If `TRUE`, use linear mixed models with random effects for the linear model component. Defaults to `FALSE`.
#' @param logmm_re logical scalar. If `TRUE`, use generalized linear mixed models with random effects for the logistic model component. Defaults to `FALSE`.
#' @param sandwich logical scalar. If `TRUE`, use sandwich variance estimators for robust inference. Defaults to `FALSE`.
#'
#' @details
#' This function performs enrichment analysis for cell-cell communication. For each target cell type pair (sender-receiver), 
#' it compares the expression of ligand-receptor pairs against background cell types. The background for a given target pair 
#' consists of all other cell types: for sender cell type X and receiver cell type Y, the background sender includes all cell types 
#' except X, and the background receiver includes all cell types except Y.
#'
#' The analysis is one-sided, testing whether the target cell type pair has higher expression of ligand-receptor pairs 
#' compared to the background. Multiple statistical methods are used to combine p-values from different tests.
#'
#' Parallelization is supported via \pkg{future}. Progress bars are customized using \pkg{progressr}.
#'
#' @returns A list with the following elements:
#'  \itemize{
#'    \item{\code{summary}}: a data frame of descriptive statistics for ligand/receptor gene expressions in target and background cell types.
#'    \item{\code{test}}: a data frame containing the results of enrichment analysis, including fold changes, z-scores, p-values, adjusted p-values, etc.
#'    \item(\code{errors}): a list of communication pairs for which analysis failed, along with corresponding error messages.
#'    \item{\code{warnings}}: a list of communication pairs for which warnings are issued during analysis, along with corresponding warning messages.
#'    \item{\code{messages}}: a list of communication pairs for which diagnostic messages are generated during analysis, along with corresponding messages.
#'  }
#'
#' @import data.table
#' @importFrom future.apply future_lapply
#' @importFrom progressr progressor
#' @importFrom utils data
#' @importFrom stats as.formula p.adjust pchisq pnorm qnorm sd
#'
#' @export

ccc_enrich <- function(expression_matrix, metadata,
                      cell_id_col = "cell_id", cell_type_col = "cell_type",
                      id_col = "id", sender = NULL, receiver = NULL,
                      lr = "omnipathr", multi_sub = "minimum",
                      verbose = TRUE, min_cell = 10,
                      min_pct = 0.01, large_n = 2, min_total_pct = 0,
                      threshold = 0, sep_detection = TRUE, sep_prop = 0, sep_n = 0,
                      padj_method = "BH", cell_type_padj = TRUE,
                      chunk_size = 10, lmm_re = FALSE, logmm_re = FALSE, sandwich = FALSE) {
  old_nthreads <- getDTthreads()
  on.exit(setDTthreads(old_nthreads), add = TRUE)

  assertthat::assert_that(assertthat::is.flag(verbose))

  if (verbose) {
    if (interactive()) {
      if (!handlers(global = NA)) {
        handlers(global = TRUE)
        handlers("cli")
        message(
          "Info: No global progress bars were found; the cli handler has been enabled. ",
          "See `vignette('progressr-intro')` for how to customize the progress bar settings."
        )
      }
    }
  }

  assertthat::assert_that(assertthat::is.flag(sep_detection))
  assertthat::assert_that(assertthat::is.flag(cell_type_padj))

  required_args <- c("expression_matrix", "metadata")
  passed_args <- names(match.call())[-1]
  missing_args <- setdiff(required_args, passed_args)
  if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
  }

  if (is.null(rownames(expression_matrix)) || is.null(colnames(expression_matrix))) {
    stop("'expression_matrix' must have rownames (genes) and colnames (cell ids).")
  }
  
  if (any(is.na(expression_matrix))) {
    stop("Missing values not allowed in expression matrix, please remove NA values.")
  }
  if (any(is.na(metadata))) {
    stop("Missing values not allowed in metadata, please remove NA values.")
  }

  multi_sub <- match.arg(multi_sub, choices = c("minimum", "arithmetic_mean", "geometric_mean", "min_avg_gene", "min_rate_gene"))

  if ("y" %in% colnames(metadata)) {
    stop("\"y\" in column names of 'metadata'. Please rename it.")
  }
  if ("z" %in% colnames(metadata)) {
    stop("\"z\" in column names of 'metadata'. Please rename it.")
  }

  if (isFALSE(cell_id_col %in% colnames(metadata))) {
    stop(paste0("'cell_id_col' \"", cell_id_col, "\" does not exist in 'metadata'."))
  }
  if (isFALSE(cell_type_col %in% colnames(metadata))) {
    stop(paste0("'cell_type_col' \"", cell_type_col, "\" does not exist in 'metadata'."))
  }
  if (isFALSE(id_col %in% colnames(metadata))) {
    stop(paste0("'id_col' \"", id_col, "\" does not exist in 'metadata'."))
  }

  if (is.null(sender)) {
    message("'sender' is not specified. All cell types will be considered as potential senders in the analysis.")
  }
  if (is.null(receiver)) {
    message("'receiver' is not specified. All cell types will be considered as potential receivers in the analysis.")
  }

  padj_method <- match.arg(padj_method, stats::p.adjust.methods)

  metadata <- as.data.table(metadata)
  metadata <- rename_metadata(
    metadata = metadata,
    cell_id_col = cell_id_col, cell_type_col = cell_type_col,
    id_col = id_col, group_col = id_col  # Use id_col as group_col for compatibility with utility functions
  )

  num_ids <- metadata[, uniqueN(id)]

  if (isTRUE(sep_detection)) {
    if (sep_prop < 0 || sep_prop > 1) {
      stop("'sep_prop' must be between 0 and 1 (inclusive).")
    }
    if (sep_n < 0 || sep_n > num_ids) {
      stop("'sep_n' must be between 0 and number of samples (inclusive).")
    }
  }

  #####################################################
  if (length(setdiff(metadata$cell_id, colnames(expression_matrix))) > 0) {
    stop("Cell ids from 'metadata' and 'expression_matrix' do not match, or the column names of 'expression_matrix' are not cell ids.")
  }
  if (!is.null(sender)) {
    missing_ct <- sender[!(sender %in% metadata$cell_type)]
    if (length(missing_ct) > 0) {
      stop(paste0("These sender cell types are missing in 'metadata': ", paste(missing_ct, collapse = ", ")))
    }
  }
  if (!is.null(receiver)) {
    missing_ct <- receiver[!(receiver %in% metadata$cell_type)]
    if (length(missing_ct) > 0) {
      stop(paste0("These receiver cell types are missing in 'metadata': ", paste(missing_ct, collapse = ", ")))
    }
  }

  ### lr table
  lr_table <- prep_lr(lr = lr)

  ### filter cell type
  # Create a dummy contrast matrix for compatibility with filter_cell_type
  dummy_contrast <- matrix(1, nrow = 1, ncol = length(unique(metadata$id)))
  colnames(dummy_contrast) <- unique(metadata$id)
  
  # Store original metadata for background cell types
  metadata_original <- copy(metadata)
  
  # Only filter cell types for sender and receiver, but keep all cells in the dataset
  filtered_obj <- filter_cell_type(
    metadata = metadata, sender = sender, receiver = receiver,
    min_cell = min_cell, contrast = dummy_contrast
  )
  
  # Update sender and receiver lists with filtered cell types
  sender <- filtered_obj$sender
  receiver <- filtered_obj$receiver
  
  # Keep the original metadata instead of using the filtered subset
  metadata_subset <- metadata_original
  
  rm(metadata, filtered_obj, metadata_original)
  gc()

  ### compute cdr if needed for future extensions
  metadata_subset <- compute_cdr(
    expression_matrix = expression_matrix,
    metadata_subset = metadata_subset, threshold = threshold
  )

  ### filter lr; prepare for analysis
  # Create all possible combinations of sender and receiver
  sender_receiver_combinations <- expand.grid(sender = sender, receiver = receiver)

  pairs4analysis <- base::merge(sender_receiver_combinations, lr_table, by = NULL)

  setDT(pairs4analysis)
  npairs <- nrow(pairs4analysis)
  i_s <- seq(1L, nrow(pairs4analysis), by = chunk_size)
  if (verbose) {
    p <- progressr::progressor(along = i_s)
    message("Starting enrichment analysis...")
  }

  # Function to analyze a chunk of pairs
  analyze_chunk <- function(i) {
    chunk <- pairs4analysis[i:min(i + chunk_size - 1L, npairs), ]

    results.summary <- results.test <- results.error <- results.warning <- results.message <- list()
    for (j in 1L:nrow(chunk)) {
      target_sender <- chunk$sender[j]
      ligand <- chunk$ligand[j]
      target_receiver <- chunk$receiver[j]
      receptor <- chunk$receptor[j]

      # Define background cell types (all except the target)
      background_sender <- setdiff(sender, target_sender)
      background_receiver <- setdiff(receiver, target_receiver)

      # Create copies of metadata_subset for target and background
      data_target_sender <- metadata_subset[cell_type == target_sender]
      data_target_receiver <- metadata_subset[cell_type == target_receiver]
      # For background, include ALL cells from non-target cell types, not just from filtered sender/receiver
      data_background_sender <- metadata_subset[cell_type != target_sender]
      data_background_receiver <- metadata_subset[cell_type != target_receiver]

      # Handle single gene or multi-gene ligands and receptors
      ligand_genes <- unlist(strsplit(ligand, "_"))
      receptor_genes <- unlist(strsplit(receptor, "_"))

      # Check if all genes exist in the expression matrix
      missing_genes <- setdiff(c(ligand_genes, receptor_genes), rownames(expression_matrix))
      if (length(missing_genes) > 0) {
        results.error[[paste(target_sender, target_receiver, ligand, receptor, sep = "_")]] <- 
          paste("Missing genes in expression matrix:", paste(missing_genes, collapse = ", "))
        next
      }

      # Try to process this pair with error handling
      tryCatch({
        # Extract expression for ligand genes in target sender
        ligand_expr_target <- if (length(ligand_genes) == 1) {
          expression_matrix[ligand_genes, data_target_sender$cell_id, drop = FALSE]
        } else {
          expression_matrix[ligand_genes, data_target_sender$cell_id, drop = FALSE]
        }
        
        # Extract expression for receptor genes in target receiver
        receptor_expr_target <- if (length(receptor_genes) == 1) {
          expression_matrix[receptor_genes, data_target_receiver$cell_id, drop = FALSE]
        } else {
          expression_matrix[receptor_genes, data_target_receiver$cell_id, drop = FALSE]
        }
        
        # Compute expression values for target cell types
        ligand_expr_target_value <- compute_expression_value(ligand_expr_target, multi_sub, threshold)
        receptor_expr_target_value <- compute_expression_value(receptor_expr_target, multi_sub, threshold)
        
        # Compute expression rates for filtering
        compute_expression_rate <- function(data_subset, expr_values) {
          sapply(unique(data_subset$id), function(uid) {
            cells_for_id <- data_subset[id == uid, cell_id]
            mean(expr_values[cells_for_id] > threshold)
          })
        }
        
        ligand_expression_rates <- compute_expression_rate(data_target_sender, ligand_expr_target_value)
        receptor_expression_rates <- compute_expression_rate(data_target_receiver, receptor_expr_target_value)
        
        # Check if the number of ids with both ligand and receptor expression rate >= min_pct is >= large_n
        valid_ids <- sum((ligand_expression_rates >= min_pct) & (receptor_expression_rates >= min_pct))
        if (valid_ids < large_n) {
          next
        }
        
        # Another filtering step based on total percentage
        total_pct_ligand <- mean(ligand_expr_target_value > threshold)
        total_pct_receptor <- mean(receptor_expr_target_value > threshold)
        if (min(total_pct_ligand, total_pct_receptor) < min_total_pct) {
          next
        }
        
        # Extract expression for ligand genes in background sender
        ligand_expr_background <- if (length(ligand_genes) == 1) {
          expression_matrix[ligand_genes, data_background_sender$cell_id, drop = FALSE]
        } else {
          expression_matrix[ligand_genes, data_background_sender$cell_id, drop = FALSE]
        }
        
        # Extract expression for receptor genes in background receiver
        receptor_expr_background <- if (length(receptor_genes) == 1) {
          expression_matrix[receptor_genes, data_background_receiver$cell_id, drop = FALSE]
        } else {
          expression_matrix[receptor_genes, data_background_receiver$cell_id, drop = FALSE]
        }

        # Compute expression values for background data
        ligand_expr_background_value <- compute_expression_value(ligand_expr_background, multi_sub, threshold)
        receptor_expr_background_value <- compute_expression_value(receptor_expr_background, multi_sub, threshold)

        # Prepare data for model fitting
        # Target data
        target_data <- data.table(
          ligand_expr = ligand_expr_target_value,
          receptor_expr = receptor_expr_target_value,
          product = ligand_expr_target_value * receptor_expr_target_value,
          group = "target",
          id = data_target_sender$id,  # Using sender's id for simplicity
          binary = (ligand_expr_target_value * receptor_expr_target_value) > threshold
        )
        
        # Background data
        background_data <- data.table(
          ligand_expr = ligand_expr_background_value,
          receptor_expr = receptor_expr_background_value,
          product = ligand_expr_background_value * receptor_expr_background_value,
          group = "background",
          id = data_background_sender$id,  # Using sender's id for simplicity
          binary = (ligand_expr_background_value * receptor_expr_background_value) > threshold
        )
        
        # Calculate expression metrics
        target_mean <- mean(target_data$product)
        target_sd <- sd(target_data$product)
        target_n <- nrow(target_data)
        target_pct_above <- mean(target_data$binary)

        background_mean <- mean(background_data$product)
        background_sd <- sd(background_data$product)
        background_n <- nrow(background_data)
        background_pct_above <- mean(background_data$binary)

        # Calculate fold change
        fold_change <- if (background_mean > 0) target_mean / background_mean else Inf
        
        # Prepare data for ligand and receptor models
        # Ligand data
        ligand_data <- data.table(
          y = c(ligand_expr_target_value, ligand_expr_background_value),
          z = c(ligand_expr_target_value > threshold, ligand_expr_background_value > threshold),
          group = c(rep("target", length(ligand_expr_target_value)), rep("background", length(ligand_expr_background_value))),
          id = c(data_target_sender$id, data_background_sender$id)
        )
        
        # Receptor data
        receptor_data <- data.table(
          y = c(receptor_expr_target_value, receptor_expr_background_value),
          z = c(receptor_expr_target_value > threshold, receptor_expr_background_value > threshold),
          group = c(rep("target", length(receptor_expr_target_value)), rep("background", length(receptor_expr_background_value))),
          id = c(data_target_receiver$id, data_background_receiver$id)
        )
        
        # Subset data for linear models (only values above threshold)
        ligand_data_1 <- ligand_data[z == TRUE]
        receptor_data_1 <- receptor_data[z == TRUE]
        
        # Initialize model fits as NULL
        fit_l_linear <- fit_l_logistic <- fit_r_linear <- fit_r_logistic <- NULL
        
        # Define model formulas
        formula_linear <- if (isTRUE(lmm_re)) {
          as.formula("y ~ 0 + group + (1|id)")
        } else {
          as.formula("y ~ 0 + group")
        }
        
        formula_logistic <- if (isTRUE(logmm_re)) {
          list(
            "fixed" = as.formula("z ~ 0 + group"),
            "random" = as.formula("~ 1 | id")
          )
        } else {
          as.formula("z ~ 0 + group")
        }
        
        # Function to fit models with error handling
        fit_model <- function(part, data, formula, name) {
          stopifnot(part == "linear" || part == "logistic")
          warnings. <- list()
          messages. <- list()
          fit <- tryCatch(
            withCallingHandlers(
              {
                if (part == "linear") {
                  if (nrow(data) == 0) {
                    stop("No data above threshold for fitting a linear model.")
                  } else {
                    if (isTRUE(lmm_re)) {
                      lme4::lmer(formula, data = data)
                    } else {
                      lm(formula, data = data)
                    }
                  }
                } else if (part == "logistic") {
                  if (var(data$z) == 0) {
                    stop("No variation in binary outcome for logistic model.")
                  } else if (isTRUE(sep_detection) && detect_re_separation(dt = data, z_col = "z", id_col = "id", num_ids = num_ids, sep_prop = sep_prop, sep_n = sep_n)) {
                    stop("Complete or quasi-complete separation detected in more than the threshold number of samples.")
                  } else {
                    if (isTRUE(logmm_re)) {
                      GLMMadaptive::mixed_model(
                        fixed = formula$fixed, random = formula$random,
                        family = binomial(), data = data
                      )
                    } else {
                      glm(formula, family = binomial(), data = data)
                    }
                  }
                }
              },
              warning = function(w) {
                warnings.[[length(warnings.) + 1L]] <<- w$message
                invokeRestart("muffleWarning")
              },
              message = function(m) {
                messages.[[length(messages.) + 1L]] <<- m$message
                invokeRestart("muffleMessage")
              }
            ),
            error = function(e) {
              results.warning[[paste(target_sender, target_receiver, ligand, receptor, sep = "_")]] <<- 
                paste("Error fitting", name, "model:", conditionMessage(e))
              return(NULL)
            }
          )
          fit
        }
        
        # Fit models for ligand and receptor
        fit_l_linear <- fit_model(part = "linear", data = ligand_data_1, formula = formula_linear, name = "ligand.linear")
        fit_r_linear <- fit_model(part = "linear", data = receptor_data_1, formula = formula_linear, name = "receptor.linear")
        fit_l_logistic <- fit_model(part = "logistic", data = ligand_data, formula = formula_logistic, name = "ligand.logistic")
        fit_r_logistic <- fit_model(part = "logistic", data = receptor_data, formula = formula_logistic, name = "receptor.logistic")
        
        # Perform statistical tests using ccc_test_enrich
        test_results <- ccc_test_enrich(
          fit_ligand_linear = fit_l_linear,
          fit_ligand_logistic = fit_l_logistic,
          fit_receptor_linear = fit_r_linear,
          fit_receptor_logistic = fit_r_logistic,
          lmm_re = lmm_re,
          logmm_re = logmm_re,
          sandwich = sandwich,
          sender = target_sender,
          receiver = target_receiver,
          ligand = ligand,
          receptor = receptor
        )
        
        # Extract p-values from test results
        p_value_linear <- if ("pvalue_linear" %in% names(test_results)) test_results$pvalue_linear else NA
        p_value_logistic <- if ("pvalue_logistic" %in% names(test_results)) test_results$pvalue_logistic else NA
        p_value_stouffer <- if ("pvalue_stouffer" %in% names(test_results)) test_results$pvalue_stouffer else NA
        p_value_fisher <- if ("pvalue_fisher" %in% names(test_results)) test_results$pvalue_fisher else NA
        
        # If we couldn't fit models, fall back to simple Z-test
        if (is.na(p_value_linear) && is.na(p_value_logistic)) {
          # Perform one-sided Z-test (target > background)
          z_score <- (target_mean - background_mean) / 
            sqrt((target_sd^2 / target_n) + (background_sd^2 / background_n))
          p_value_linear <- pnorm(z_score, lower.tail = FALSE)
          
          # Fisher's exact test for expression rates
          fisher_table <- matrix(
            c(
              sum(target_data$binary), sum(!target_data$binary),
              sum(background_data$binary), sum(!background_data$binary)
            ),
            nrow = 2
          )
          p_value_logistic <- fisher.test(fisher_table, alternative = "greater")$p.value
          
          # Combine p-values using different methods
          # Stouffer's method
          z_values <- c(qnorm(p_value_linear, lower.tail = FALSE), qnorm(p_value_logistic, lower.tail = FALSE))
          p_value_stouffer <- pnorm(sum(z_values) / sqrt(length(z_values)), lower.tail = FALSE)
          
          # Fisher's method
          p_value_fisher <- pchisq(-2 * sum(log(c(p_value_linear, p_value_logistic))), df = 4, lower.tail = FALSE)
        }

        # Store summary statistics
        results.summary[[paste(target_sender, target_receiver, ligand, receptor, sep = "_")]] <- data.table(
          sender = target_sender,
          receiver = target_receiver,
          ligand = ligand,
          receptor = receptor,
          target_mean = target_mean,
          target_sd = target_sd,
          target_n = target_n,
          target_pct_above = target_pct_above,
          background_mean = background_mean,
          background_sd = background_sd,
          background_n = background_n,
          background_pct_above = background_pct_above
        )

        # Store test results
        results.test[[paste(target_sender, target_receiver, ligand, receptor, sep = "_")]] <- data.table(
          sender = target_sender,
          receiver = target_receiver,
          ligand = ligand,
          receptor = receptor,
          fold_change = fold_change,
          p_value_linear = p_value_linear,
          p_value_logistic = p_value_logistic,
          p_value_stouffer = p_value_stouffer,
          p_value_fisher = p_value_fisher,
          effect_size_linear = if ("effect_size_linear" %in% names(test_results)) test_results$effect_size_linear else NA,
          effect_size_logistic = if ("effect_size_logistic" %in% names(test_results)) test_results$effect_size_logistic else NA,
          effect_size_hurdle = if ("effect_size_hurdle" %in% names(test_results)) test_results$effect_size_hurdle else NA
        )
      },
      error = function(e) {
        results.error[[paste(target_sender, target_receiver, ligand, receptor, sep = "_")]] <<- conditionMessage(e)
      },
      warning = function(w) {
        results.warning[[paste(target_sender, target_receiver, ligand, receptor, sep = "_")]] <<- conditionMessage(w)
      },
      message = function(m) {
        results.message[[paste(target_sender, target_receiver, ligand, receptor, sep = "_")]] <<- conditionMessage(m)
      })
    }

    if (verbose) p()

    list(
      summary = results.summary,
      test = results.test,
      error = results.error,
      warning = results.warning,
      message = results.message
    )
  }

  # Run analysis in parallel or sequentially
  if (requireNamespace("future.apply", quietly = TRUE)) {
    results_list <- future.apply::future_lapply(i_s, analyze_chunk)
  } else {
    results_list <- lapply(i_s, analyze_chunk)
  }

  # Combine results
  summary_list <- lapply(results_list, function(x) x$summary)
  test_list <- lapply(results_list, function(x) x$test)
  error_list <- lapply(results_list, function(x) x$error)
  warning_list <- lapply(results_list, function(x) x$warning)
  message_list <- lapply(results_list, function(x) x$message)

  summary_dt <- rbindlist(unlist(summary_list, recursive = FALSE), fill = TRUE)
  test_dt <- rbindlist(unlist(test_list, recursive = FALSE), fill = TRUE)
  error_list <- unlist(error_list, recursive = FALSE)
  warning_list <- unlist(warning_list, recursive = FALSE)
  message_list <- unlist(message_list, recursive = FALSE)

  # Adjust p-values
  if (cell_type_padj) {
    # Adjust p-values for each sender-receiver pair
    test_dt[, c("padj_linear", "padj_logistic", "padj_stouffer", "padj_fisher") := 
              .(p.adjust(p_value_linear, method = padj_method),
                p.adjust(p_value_logistic, method = padj_method),
                p.adjust(p_value_stouffer, method = padj_method),
                p.adjust(p_value_fisher, method = padj_method)),
            by = .(sender, receiver)]
  } else {
    # Adjust p-values globally
    test_dt[, c("padj_linear", "padj_logistic", "padj_stouffer", "padj_fisher") := 
              .(p.adjust(p_value_linear, method = padj_method),
                p.adjust(p_value_logistic, method = padj_method),
                p.adjust(p_value_stouffer, method = padj_method),
                p.adjust(p_value_fisher, method = padj_method))]
  }

  # Return results
  list(
    summary = summary_dt,
    test = test_dt,
    errors = error_list,
    warnings = warning_list,
    messages = message_list
  )
}