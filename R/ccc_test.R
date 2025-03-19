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
#'  

ccc_analysis <- function(expression_matrix, metadata,
                         cell_id = "cell_id", cell_type = "cell_type", group = "group", covar = NULL, cdr = TRUE,
                         id = NULL, lmm_re = TRUE, logmm_re = TRUE,
                         sender = NULL, receiver = NULL,
                         lr = c("omnipathr","ramilowski"),
                         multi_sub = c("minimum","arithmetic_mean","geometric_mean","min_avg_gene","min_rate_gene"),
                         contrast,
                         verbose = TRUE,
                         min_pct = 0.01, large_n = 2, min_avg_pct = 0,
                         min_cell = 10,
                         cutoff = 0,
                         padj_method = "BH", sr_adj = TRUE,
                         control_lm, control_logmm,
                         ...
) {
  
  
}

#'
#' @import data.table



rename_metadata <- function(metadata, cell_id, id, group, cell_type) {
  metadata <- as.data.table(metadata)
  if (is.null(id)) {
    setnames(metadata, old = c(cell_id, group, cell_type), new = c("cell_id", "group", "cell_type"))
    metadata[, id := group]
  } else {
    setnames(metadata, old = c(cell_id, id, group, cell_type), new = c("cell_id", "id", "group", "cell_type"))
  }
  return(metadata)
}



#' 
#' 
prep_lr <- function(lr) {
  if (is.character(lr)) {
    lr_name <- lr[1L]
  } else if (is.data.frame(lr)) {
    stopifnot(colnames(lr) == c("ligand","receptor"))
    lr_name <- "user"
    lr_user <- lr
  } 
  
  lr.omnipathr <- function() {
    data("omnipathr", envir = environment())
    omnipathr
    }
  lr.ramilowski <- function() {
    data("ramilowski", envir = environment())
    ramilowski
    }

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


# 
# oplan <- plan(multisession, workers = 2)
# on.exit(plan(oplan), add = TRUE)

#' 
#' @import data.table
#' 
filter_cell_type <- function(metadata, sender, receiver, min_cell, contrast) {
  
  if (is.null(sender)) {
    sender <- unique(metadata$cell_type)
  }
  if (is.null(receiver)) {
    receiver <- unique(metadata$cell_type)
  }
  
  # Apply contrast filtering
  if (is.vector(contrast)) {
    zero_groups <- names(contrast)[contrast == 0]
  } else if (is.matrix(contrast)) {
    zero_groups <- colnames(contrast)[colSums(contrast) == 0]
  } else {
    stop("'contrast' must be either a named vector or a matrix with column names.")
  }
  
  if (length(zero_groups) > 0) {
    metadata <- metadata[!group %in% zero_groups]
  }
  
  # # Check if any rows remain after contrast filtering
  # if (nrow(metadata) == 0) {
  #   stop("No rows remain after applying contrast filtering.")
  # }
  
  # Find rows where cell_type appears in either sender OR receiver
  valid_cell_types <- union(sender, receiver)
  metadata_subset <- metadata[cell_type %in% valid_cell_types]
  
  # # Check if any rows remain
  # if (nrow(metadata_subset) == 0) {
  #   stop("No rows remain after filtering by sender or receiver cell types.")
  # }
  
  # Count occurrences of each id-cell_type combination
  counts <- metadata_subset[, .N, by = .(id, cell_type)]
  
  # Keep only cell_type where count >= min_cell for each id
  metadata_subset <- metadata_subset[counts[N >= min_cell], on = .(id, cell_type)]
  
  # # Check again if any rows remain after applying min_cell filter
  # if (nrow(metadata_subset) == 0) {
  #   stop("No rows remain after applying min_cell filter.")
  # }
  
  # Ensure each id has the same unique set of cell_type values
  valid_cell_types_by_id <- metadata_subset[, .(unique_cell_types = list(unique(cell_type))), by = id]
  
  # Find the intersection of cell_types across all ids
  common_cell_types <- Reduce(intersect, valid_cell_types_by_id$unique_cell_types)
  
  # Filter metadata_subset so that only these common cell_types remain for each id
  metadata_subset <- metadata_subset[cell_type %in% common_cell_types]
  
  # Check if any rows remain after ensuring consistent cell_types across ids
  if (nrow(metadata_subset) == 0) {
    stop("No cell types remain after 'min_cell' filtering.")
  }
  
  # Check if sender or receiver contain elements not in the final subset
  remaining_cell_types <- unique(metadata_subset$cell_type)
  
  missing_sender <- setdiff(sender, remaining_cell_types)
  missing_receiver <- setdiff(receiver, remaining_cell_types)
  
  if (length(missing_sender) > 0 || length(missing_receiver) > 0) {
    warning_msg <- "Some cell types in 'sender' or 'receiver' do not appear in the final subset."
    if (length(missing_sender) > 0) {
      warning_msg <- paste0(warning_msg, "\nMissing in sender: ", paste(missing_sender, collapse = ", "))
    }
    if (length(missing_receiver) > 0) {
      warning_msg <- paste0(warning_msg, "\nMissing in receiver: ", paste(missing_receiver, collapse = ", "))
    }
    warning(warning_msg)
  }
  
  sender <- intersect(sender, remaining_cell_types)
  receiver <- intersect(receiver, remaining_cell_types)
  
  if (length(sender) == 0) {
    stop("No sender cell types remain after 'min_cell' filtering.")
  }
  if (length(receiver) == 0) {
    stop("No receiver cell types remain after 'min_cell' filtering.")
  }
  
  return(list(metadata_subset = metadata_subset, sender = sender, receiver = receiver))
}


#' 
#' @import data.table
#' 
filter_lr <- function(expression_matrix, metadata_subset,
                      cdr,
                      sender, receiver,
                      lr_table,
                      multi_sub = c("minimum","arithmetic_mean","geometric_mean","min_avg_gene","min_rate_gene"),
                      contrast,
                      verbose = TRUE,
                      min_pct = 0.01, large_n = 2, min_avg_pct = 0,
                      min_cell = 10,
                      cutoff = 0,
                      ...) {
  if (is.null(id)) id <- group
  
}


run_analysis <- function(){
  # Create all possible combinations of sender and receiver
  sender_receiver_combinations <- expand.grid(sender = sender, receiver = receiver)
  
  # Merge with lr_table to get all possible combinations
  result <- merge(sender_receiver_combinations, lr_table, by = NULL)
  
  # Convert to data.table
  result_dt <- data.table(result)
}




