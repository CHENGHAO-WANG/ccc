#' Log Normalize Data
#' 
#' Log normalize the count matrix. Gene counts for each cell are divided by the total counts for that cell and multiplied by `scale.factor`. This is incremented by 1 and then log transformed.
#' 
#' @param count_matrix a numeric matrix of raw counts, with rows corresponding to genes and columns corresponding to cells.
#' @param scale_factor numeric. The scale factor fro cell-level normalization. Defaults to 10000.
#' @param base numeric. The base of the log transformation. Defaults to `exp(1)`.
#' @returns a numeric matrix of normalized counts.
#' @examples
#' set.seed(123)
#' raw <- matrix(rnbinom(100, mu = 1.5, size = 2), nrow = 20)
#' normalized <- log_normalize(raw)
#' normalized
#' 
#' @export

# log_normalize <- function(count_matrix, MARGIN = 2L, scale_factor = 10000, base = exp(1)) {
#   
#   lib_size <- switch(as.character(MARGIN),
#                      "1" = rowSums(count_matrix),
#                      "2" = colSums(count_matrix),
#                      stop("'MARGIN' must be 1 or 2"))
#   expression_matrix <- logb(sweep(count_matrix, MARGIN = MARGIN, STATS = lib_size, FUN = "/") * scale_factor + 1, base = base)
#   return(expression_matrix)
# }
log_normalize <- function(count_matrix, scale_factor = 10000, base = exp(1)) {
  
  lib_size <- Matrix::colSums(count_matrix)
  expression_matrix <- logb(sweep(count_matrix, MARGIN = 2L, STATS = lib_size, FUN = "/") * scale_factor + 1, base = base)
  return(expression_matrix)
}
