#' Log Normalize Data
#'
#' Log normalize the count matrix. Gene counts for each cell are divided by the total counts for that cell and multiplied by `scale.factor`. This is incremented by 1 and then log transformed.
#'
#' @param count_matrix a numeric matrix of raw counts, with rows corresponding
#'   to genes and columns corresponding to cells.
#' @param num_scale_factor numeric. The scale factor for cell-level
#'   normalization. Defaults to NULL. If this parameter is specified, it
#'   overrides `scale_factor`.
#' @param scale_factor character. The method for cell-level normalization.
#'   This argument only works when `num_scale_factor` is NULL.
#'   \itemize{
#'     \item \dQuote{\code{geometric}}: (the default) uses the geometric mean of library sizes.
#'     \item \dQuote{\code{arithmetic}}: uses the arithmetic mean of library sizes.
#'     \item \dQuote{\code{median}}: uses the median of library sizes.
#'   }
#' @param base numeric. The base of the log transformation. Defaults to `2`.
#' @param trim_lower,trim_upper numeric. The fraction of observations to be trimmed
#'   from the lower and upper ends of the library size distribution before its
#'   center is computed. These parameters only work when `num_scale_factor` is NULL.
#'   Each value should be in \eqn{[0, 0.45]}. Both default to 0.1.
#' @returns A numeric matrix of normalized counts.
#' @examples
#' set.seed(123)
#' raw <- matrix(rnbinom(100, mu = 1.5, size = 2), nrow = 20)
#' normalized <- log_normalize(raw)
#' normalized
#'
#' @export
log_normalize <- function(count_matrix,
                          num_scale_factor = NULL,
                          scale_factor = "geometric",
                          base = 2,
                          trim_lower = 0.1,
                          trim_upper = 0.1) {
  lib_size <- Matrix::colSums(count_matrix)
  if (is.null(num_scale_factor)) {
    if (trim_lower < 0 | trim_lower > 0.45 | trim_upper < 0 | trim_upper > 0.45) {
      stop("'trim_lower' and 'trim_upper' must be between 0 and 0.45.")
    }
    lib_size4nsf <- lib_size
    if (trim_lower > 0 | trim_upper > 0) {
      lib_size4nsf <- sort(lib_size4nsf)
      n <- length(lib_size4nsf)
      lo <- floor(n * trim_lower) + 1
      hi <- n - floor(n * trim_upper)
      lib_size4nsf <- lib_size4nsf[lo:hi]
    }
    scale_factor <- match.arg(scale_factor, c("geometric", "arithmetic", "median"))
    num_scale_factor <- switch(scale_factor,
      "geometric" = exp(mean(log(lib_size4nsf))),
      "arithmetic" = mean(lib_size4nsf),
      "median" = median(lib_size4nsf)
    )
  }
  expression_matrix <- logb(sweep(count_matrix, MARGIN = 2L, STATS = lib_size, FUN = "/") * num_scale_factor + 1, base = base)
  return(expression_matrix)
}
