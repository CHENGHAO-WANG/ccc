#' 
#' 
#' 
#' 
#' @param multi_sub it 
#'  \itemize{
#'    \item
#'    \item
#'  }

ccc_analysis <- function(expression_matrix, metadata, covar_to_test, covar = NULL, cdr = TRUE,
                         id = NULL, lmm_re = TRUE, logmm_re = TRUE,
                         lr = c("omnipathr","ramilowski"),
                         multi_sub = c("minimum","arithmetic_mean","geometric_mean","min_avg_gene","min_rate_gene"),
                         contrast,
                         verbose = TRUE,
                         filter
) {
  
  
}


#' 
#' 
lr_prep <- function(lr) {
  if (is.character(lr)) {
    lr_name <- lr[1L]
  } else if (is.data.frame(lr)) {
    stopifnot(colnames(lr) == c("ligand","receptor"))
    lr_name <- "user"
    lr_user <- lr
  } 
  
  lr.omnipathr <- function() {data("omnipathr", envir = environment())
    omnipathr}
  lr.ramilowski <- function() {data("ramilowski", envir = environment())
    ramilowski}

  lr_table <- switch(EXPR = lr_name,
                     "omnipathr" = lr.omnipathr(),
                     "ramilowski" = lr.ramilowski(),
                     "user" = lr_user,
                     stop("'lr' should be \"omnipathr\" or \"ramilowski\" or a data.frame of ligand-receptor pairs")
  )
  
  return(lr_table)
}



# intercell_network <- OmnipathR::intercell_network(ligand_receptor = TRUE, high_confidence = TRUE, simplify = TRUE)
# omnipathr <- data.frame(ligand=intercell_network$source_genesymbol,receptor=intercell_network$target_genesymbol)
# save(omnipathr, file = "omnipathr.rda")



oplan <- plan(multisession, workers = 2)
on.exit(plan(oplan), add = TRUE)
