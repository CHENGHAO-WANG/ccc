#' Identify Differential Cell-Cell Communication
#'
#' For each communication pair (defined by a distinct combination of sender, receiver, ligand, and receptor), fit two gene-wise hurdle models (linear for expression levels > `threshold`; logistic for expression levels > `threshold` vs. expression levels <= `threshold`).
#' Then perform Wald tests on the product of ligand expression levels in sender and receptor expression levels in receiver for differential cell-cell communication.
#'
#' @param expression_matrix a numeric matrix of normalized counts, with rows corresponding to genes and columns corresponding to cells. Both row names (gene symbols) and column names (cell identifiers) must be provided.
#' @param metadata a data frame containing cell-level metadata (e.g., cell type, group, id, covariates).
#' @param contrast a named numeric vector or a numeric matrix with column names. The names should match the levels of the variable being tested (specified by `group_col`). The testings are performed based on this.
#' @param cell_id_col a character string specifying the name of the column in `metadata` that contains cell identifiers. These ID's should match the column names of `expression_matrix`. Defaults to `"cell_id"`.
#' @param cell_type_col a character string specifying the name of the column in `metadata` that contains cell type annotations. Defaults to `"cell_type"`.
#' @param group_col a character string specifying the name of the column in `metadata` that represents the variable to be tested. The values in this column should match the names specified in `contrast`. Defaults to `group`.
#' @param covar_col a character string or a character vector specifying the column(s) in `metadata` that represent covariates to include in the model. Defaults to `NULL`, meaning no covariates are adjusted for.
#' @param cdr logical scalar. If `TRUE` (the default), calculate and adjust for cellular detection rates (CDR). The CDR of a cell is defined as the number of genes with expression level above `threshold` divided by the total number of genes.
#' @param id_col a character string specifying the name of the column in `metadata` that contains individual-level (sample-level) ID's. Used for random effect modeling. Must be provided if `lmm_re == TRUE` or `logmm_re == TRUE`; otherwise, this can be omitted. Defaults to `"id"`.
#' @param lmm_re logical scalar. Should a random effect be included in the linear component of the hurdle model? If `TRUE` (the default), fit a linear mixed-effects model with a random intercept based on `id_col` for each ligand/receptor gene; if `FALSE`, fit a linear model without random effects.
#' @param logmm_re logical scalar. Should a random effect be included in the logistic component of the hurdle model? If `TRUE` (the default), fit a logistic mixed-effects model with a random intercept based on `id_col` for each ligand/receptor gene; if `FALSE`, fit a logistic model without random effects.
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
#' @param sandwich logical scalar. If `TRUE`, sandwich standard errors are used in the calculations. Defaults to `FALSE`.
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
#' @param control_logm control parameters for optimization in [stats::glm()].
#' @param control_lmm control parameters for optimization in [lme4::lmer()].
#' @param control_logmm control parameters for optimization in [GLMMadaptive::mixed_model].
#' @param chunk_size integer scalar. The number of communication pairs (each defined by a distinct combination of sender, receiver, ligand, and receptor) to be sent to each parallel environment. Defaults to 10. To enable parallelization, users should use the \pkg{future} package.
#' @details
#' This function perform differential cell-cell communication analysis. For each communication pair, a hurdle model is fitted to ligand expression data in sender and another hurdle model is fitted to receptor expression data in receiver.
#' The delta method is applied to obtain appropriate standard errors for the product of ligand and receptor expression levels, using estimates from the fitted hurdle models. Then Wald tests are performed.
#'
#' Users can either specify the relevant column names of `metadata` using the arguments `cell_id_col`, `cell_type_col`, `group_col`, `covar_col` (can be omitted if no covariates are to be adjusted for), and `id_col`.
#' To adjust for CDR, simply set `cdr = TRUE` (This is also the default setting); do not calculate and add CDR manually to `metadata`.
#'
#' For each ligand/receptor gene, if the model fitting fails or produce warnings or diagnostic messages, this function will return the relevant information in the output. If some p-values are missing in the output, users can check `errors` element of the output to find the corresponding error messages.
#'
#' Parallelization is supported via \pkg{future}. Progress bars are customized using \pkg{progressr}.
#'
#' @returns A list with the following elements:
#'  \itemize{
#'    \item{\code{summary}}: a data frame of descriptive statistics for ligand/receptor gene expressions in sender/receiver cell types.
#'    \item{\code{test}}: a data frame containing the results of differential cell-cell communication analysis, including effect sizes, p-values, adjusted p-values, etc.
#'    \item(\code{errors}): a list of communication pairs for which model fitting failed, along with corresponding error messages.
#'    \item{\code{warnings}}: a list of communication pairs for which warnings are issued during model fitting, along with corresponding warning messages.
#'    \item{\code{messages}}: a list of communication pairs for which diagnostic messages are generated during model fitting, along with corresponding messages.
#'  }
#'
#' @import data.table
#' @importFrom future.apply future_lapply
#' @importFrom lme4 lmer fixef
#' @importFrom GLMMadaptive mixed_model marginal_coefs
#' @import progressr
#' @importFrom utils data
#' @importFrom stats as.formula binomial dlogis glm lm p.adjust pchisq plogis pnorm qnorm rgamma rnbinom rnorm runif vcov coef
#'
#' @examples
#' \dontrun{
#' ## Run with example data
#' data(data.sim, package = "ccc")
#' data(lr.sim, package = "ccc")
#' expression_matrix <- log_normalize(data.sim$counts)
#' metadata <- data.sim$metadata
#' 
#' # Run sequentially
#' a <- ccc_diff(
#'   expression_matrix = expression_matrix, metadata = metadata,
#'   id_col = "sample", lr = lr.sim, sender = "CT1", receiver = c("CT2", "CT3"),
#'   contrast = c(grp2 = 1, grp1 = -1), lmm_re = TRUE, logmm_re = TRUE
#' )
#' head(a$summary)
#' head(a$test)
#'
#' # Run in parallel
#' library(future)
#' oplan <- plan(multisession, workers = 4L)
#' a <- ccc_diff(
#'   expression_matrix = expression_matrix, metadata = metadata,
#'   id_col = "sample", lr = lr.sim, sender = "CT1", receiver = c("CT2", "CT3"),
#'   contrast = c(grp2 = 1, grp1 = -1), lmm_re = TRUE, logmm_re = TRUE
#' )
#' plan(oplan)
#' head(a$summary)
#' head(a$test)
#' }
#'
#' @export

ccc_diff <- function(expression_matrix, metadata, contrast,
                     cell_id_col = "cell_id", cell_type_col = "cell_type",
                     group_col = "group", covar_col = NULL, cdr = TRUE,
                     id_col = "id", lmm_re = TRUE, logmm_re = TRUE,
                     sender = NULL, receiver = NULL,
                     lr = "omnipathr", multi_sub = "minimum",
                     sandwich = FALSE, verbose = TRUE, min_cell = 10,
                     min_pct = 0.01, large_n = 2, min_total_pct = 0,
                     threshold = 0, sep_detection = TRUE, sep_prop = 0, sep_n = 0,
                     padj_method = "BH", cell_type_padj = TRUE,
                     control_logm = list(),
                     control_lmm = lme4::lmerControl(), control_logmm = list(),
                     chunk_size = 10) {
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

  assertthat::assert_that(assertthat::is.flag(cdr))
  assertthat::assert_that(assertthat::is.flag(lmm_re))
  assertthat::assert_that(assertthat::is.flag(logmm_re))
  assertthat::assert_that(assertthat::is.flag(sandwich))
  assertthat::assert_that(assertthat::is.flag(sep_detection))
  assertthat::assert_that(assertthat::is.flag(cell_type_padj))

  required_args <- c("expression_matrix", "metadata", "contrast")
  passed_args <- names(match.call())[-1]
  missing_args <- setdiff(required_args, passed_args)
  if (length(missing_args) > 0) {
    stop("Missing required arguments: ", paste(missing_args, collapse = ", "))
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

  if (is.null(rownames(expression_matrix)) || is.null(colnames(expression_matrix))) {
    stop("'expression_matrix' must have rownames (genes) and colnames (cell ids).")
  }
  # contrast <- contrast[, colSums(abs(contrast) > 0)]
  if (any(is.na(expression_matrix))) {
    stop("Missing values not allowed in expression matrix, please remove NA values.")
  }
  if (any(is.na(metadata))) {
    stop("Missing values not allowed in metadata, please remove NA values.")
  }

  multi_sub <- match.arg(multi_sub, choices = c("minimum", "arithmetic_mean", "geometric_mean", "min_avg_gene", "min_rate_gene"))

  err_msg <- " in 'covar_col'. Please use a different argument or rename it."
  if ("id" %in% covar_col) {
    stop(paste0("\"id\"", err_msg))
  }
  if ("group" %in% covar_col) {
    stop(paste0("\"group\"", err_msg))
  }
  if ("cdr" %in% covar_col && isTRUE(cdr)) {
    stop(paste0("\"cdr\"", err_msg))
  }

  if ("y" %in% colnames(metadata)) {
    stop("\"y\" in column names of 'metadata'. Please rename it.")
  }
  if ("z" %in% colnames(metadata)) {
    stop("\"z\" in column names of 'metadata'. Please rename it.")
  }

  covar_exists <- covar_col %in% colnames(metadata)
  if (!is.null(covar_col) && !all(covar_exists)) {
    stop(paste0("The following columns specified in 'covar_col' do not exist in 'metadata': ", paste(covar_col[!covar_exists], collapse = ", ")))
  }
  if (isFALSE(group_col %in% colnames(metadata))) {
    stop(paste0("'group_col' \"", group_col, "\" does not exist in 'metadata'."))
  }
  if (isFALSE(cell_id_col %in% colnames(metadata))) {
    stop(paste0("'cell_id_col' \"", cell_id_col, "\" does not exist in 'metadata'."))
  }
  if (isFALSE(cell_type_col %in% colnames(metadata))) {
    stop(paste0("'cell_type_col' \"", cell_type_col, "\" does not exist in 'metadata'."))
  }

  if (!all(colnames(contrast) %in% unique(metadata[[group_col]]))) {
    stop(paste0("'contrast' contains group levels that are not present in the \"", group_col, "\" column of 'metadata'."))
  }
  metadata <- metadata[metadata[[group_col]] %in% colnames(contrast), ]

  if (isTRUE(lmm_re) || isTRUE(logmm_re)) {
    if (is.null(id_col)) {
      stop("'id_col' is not specified, while 'lmm_re' or 'logmm_re' is TRUE")
    }
    if (isFALSE(id_col %in% colnames(metadata))) {
      stop(paste0("'id_col' \"", id_col, "\" does not exist in 'metadata'."))
    }
  }
  if (isFALSE(lmm_re) && isFALSE(logmm_re) && !is.null(id_col)) {
    message("'id_col' is not NULL. This input will be ignored, because 'lmm_re' and 'logmm_re' are FALSE")
    message("To suppress this message, set 'id_col = NULL'.")
    id_col <- NULL
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
    id_col = id_col, group_col = group_col
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
  filtered_obj <- filter_cell_type(
    metadata = metadata, sender = sender, receiver = receiver,
    min_cell = min_cell, contrast = contrast
  )
  metadata_subset <- filtered_obj$metadata_subset
  sender <- filtered_obj$sender
  receiver <- filtered_obj$receiver

  rm(metadata, filtered_obj)
  gc()


  ### compute cdr if specified by the user
  if (isTRUE(cdr)) {
    metadata_subset <- compute_cdr(
      expression_matrix = expression_matrix,
      metadata_subset = metadata_subset, threshold = threshold
    )
  }

  ### filter lr; fit models; conduct tests
  # Create all possible combinations of sender and receiver
  sender_receiver_combinations <- expand.grid(sender = sender, receiver = receiver)

  pairs4analysis <- base::merge(sender_receiver_combinations, lr_table, by = NULL)

  setDT(pairs4analysis)
  npairs <- nrow(pairs4analysis)
  unique_ids <- unique(metadata_subset[, id])
  i_s <- seq(1L, nrow(pairs4analysis), by = chunk_size)
  if (verbose) {
    p <- progressr::progressor(along = i_s)
    message("Starting statistical analysis...")
  }

  run_analysis <- function(i) {
    chunk <- pairs4analysis[i:min(i + chunk_size - 1L, npairs), ]


    results.summary <- results.test <- results.error <- results.warning <- results.message <- list()
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
      data_sender_ligand[, y := compute_expression_value(ligand_expr_values, multi_sub, threshold)]
      data_receiver_receptor[, y := compute_expression_value(receptor_expr_values, multi_sub, threshold)]

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

      #######################
      # another filtering step
      total_pct_ligand <- data_sender_ligand[, mean(y > threshold)]
      total_pct_receptor <- data_receiver_receptor[, mean(y > threshold)]
      if (min(total_pct_ligand, total_pct_receptor) < min_total_pct) {
        next
      }

      ##
      # Add indicator column
      data_sender_ligand[, z := ifelse(y > threshold, 1, 0)]
      data_receiver_receptor[, z := ifelse(y > threshold, 1, 0)]

      # Subset data
      data_sender_ligand_1 <- data_sender_ligand[z == 1]
      data_receiver_receptor_1 <- data_receiver_receptor[z == 1]

      # descriptive statistics summary
      dt.summary.ligand <- compute_group_stats(dt = data_sender_ligand, prefix = "ligand.")
      dt.summary.receptor <- compute_group_stats(dt = data_receiver_receptor, prefix = "receptor.")
      dt.summary <- merge(dt.summary.ligand, dt.summary.receptor, by = "group")
      dt.summary[, c("sender", "receiver", "ligand", "receptor") := list(sender, receiver, ligand, receptor)]
      setcolorder(dt.summary, c("sender", "receiver", "ligand", "receptor", setdiff(names(dt.summary), c("sender", "receiver", "ligand", "receptor"))))
      results.summary[[length(results.summary) + 1L]] <- dt.summary

      # Define covariates
      if (is.null(covar_col) && isFALSE(cdr)) {
        covariates <- "group"
      } else {
        covar <- c(covar_col, if (isTRUE(cdr)) "cdr")

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

      fit_model <- function(part, data, formula, name) {
        stopifnot(part == "linear" || part == "logistic")
        warnings. <- list()
        messages. <- list()
        fit <- tryCatch(
          withCallingHandlers(
            {
              if (part == "linear") {
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
              } else if (part == "logistic") {
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
            error_messages[[name]] <<- e$message
            return(NULL)
          }
        )
        if (length(warnings.) > 0) {
          warning_messages[[name]] <<- warnings.
        }
        if (length(messages.) > 0) {
          the_messages[[name]] <<- messages.
        }
        fit
      }

      ##
      warning_messages <- the_messages <- error_messages <- list()
      
      fit.l.linear <- fit_model(part = "linear", data = data_sender_ligand_1, formula = formula_linear, name = "ligand.linear")
      fit.r.linear <- fit_model(part = "linear", data = data_receiver_receptor_1, formula = formula_linear, name = "receptor.linear")
      fit.l.logistic <- fit_model(part = "logistic", data = data_sender_ligand, formula = formula_logistic, name = "ligand.logistic")
      fit.r.logistic <- fit_model(part = "logistic", data = data_receiver_receptor, formula = formula_logistic, name = "receptor.logistic")

      if (length(warning_messages) > 0) {
        results.warning[[paste(sender, receiver, ligand, receptor, sep = "-")]] <- warning_messages
      }
      if (length(the_messages) > 0) {
        results.message[[paste(sender, receiver, ligand, receptor, sep = "-")]] <- the_messages
      }
      if (length(error_messages) > 0) {
        results.error[[paste(sender, receiver, ligand, receptor, sep = "-")]] <- error_messages
      }

      results.test[[length(results.test) + 1L]] <- ccc_test(
        fit.l.linear = fit.l.linear, fit.l.logistic = fit.l.logistic,
        fit.r.linear = fit.r.linear, fit.r.logistic = fit.r.logistic,
        contrast = contrast, lmm_re = lmm_re, logmm_re = logmm_re,
        sandwich = sandwich,
        sender = sender, receiver = receiver,
        ligand = ligand, receptor = receptor
      )
    }
    if (verbose) {
      p()
    }
    list(
      descriptive_stats = rbindlist(results.summary),
      test_results = rbindlist(results.test, fill = TRUE),
      errors = results.error,
      warnings = results.warning,
      messages = results.message
    )
  }

  setDTthreads(threads = 1)
  results_obj <- future_lapply(i_s, FUN = run_analysis, future.seed = TRUE)
  setDTthreads(threads = old_nthreads)

  list.descriptive_stats <- lapply(results_obj, \(x) x$descriptive_stats)
  list.test_results <- lapply(results_obj, \(x) x$test_results)
  list.errors <- lapply(results_obj, \(x) x$errors)
  list.warnings <- lapply(results_obj, \(x) x$warnings)
  list.messages <- lapply(results_obj, \(x) x$messages)

  dt.test.all <- rbindlist(list.test_results, fill = TRUE)
  # cell_type_padj

  # pval_cols <- c("pvalue_linear", "pvalue_logistic", "pvalue_hurdle", "pvalue_2part")
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

  list(
    summary = as.data.frame(rbindlist(list.descriptive_stats)),
    test = as.data.frame(dt.test.all),
    errors = do.call(c, list.errors),
    warnings = do.call(c, list.warnings),
    messages = do.call(c, list.messages)
  )
}
