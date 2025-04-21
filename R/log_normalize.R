

#' @export
log_normalize <- function(count_matrix, MARGIN = 2L, scale_factor = 10000, base = exp(1)) {
  
  lib_size <- switch(as.character(MARGIN),
                     "1" = rowSums(count_matrix),
                     "2" = colSums(count_matrix),
                     stop("'MARGIN' must be 1 or 2"))
  expression_matrix <- logb(sweep(count_matrix, MARGIN = MARGIN, STATS = lib_size, FUN = "/") * scale_factor + 1, base = base)
  return(expression_matrix)
}

