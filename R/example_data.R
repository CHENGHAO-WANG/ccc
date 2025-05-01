#' Example Dataset
#' 
#' Example scRNA-seq data from simulation.
#' 
#' @format ## `data.sim`
#' A list with the following elements:
#' \describe{
#'   \item{\code{counts}}{a matrix of counts with 1000 rows (genes) and 1200 columns (cells)}
#'   \item{\code{metadata}}{a data frame of cell-level metadata with 1200 rows and 4 columns (`cell_id`, `sample`, `group`, and `cell_type`)}
#'   \item{\code{gene_info}}{a data frame with 1000 rows and 9 columns:
#'   - `gene_id`: gene symbol
#'   - `baseline_log_mu`: baseline log mean count
#'   - `dispersion`: dispersion
#'   - `within_sample_correlation`: does this gene have a sample-level random intercept?
#'   - `CT1_logfc`, `CT2_logfc`, `CT3_logfc`: log fold change of cell-type-specific genes in "CT1", "CT2", and "CT3" cell types
#'   - `grp2_logfc`, `grp3_logfc`: log fold change comparing "grp2" and "grp3" to "grp1" group
#'   }
#' }
#' 
#' @source Created using [sim_count()].
"data.sim"

#' Example Ligand-Receptor Pairs
#' 
#' Example ligand-receptor pairs database from simulation. This is designed to be used in conjunction with [data.sim].
#' 
#' @format ## `lr.sim`
#' A data frame with 20 rows and 2 columns: `ligand` and `receptor`.
#' 
#' @source Created using [sim_lr()].
"lr.sim"