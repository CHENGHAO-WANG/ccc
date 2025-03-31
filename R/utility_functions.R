

center_covar <- function(dt, covar) {
  dt_copy <- copy(dt) # Avoid modifying original data.table
  
  categorical_covar <- covar[sapply(dt_copy[, covar, with = FALSE], function(col) is.factor(col) || is.character(col))]
  
  # Convert categorical covar columns to numeric (one-hot encoding)
  for (col_name in categorical_covar) {
    col <- dt_copy[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt_copy[, (new_col_name) := as.integer(col == level)]
    }
    dt_copy[, (col_name) := NULL] # Remove original categorical column
    #Add new one hot encoded columns to the covar list, and remove the original column.
    covar <- c(covar[covar != col_name], names(dt_copy)[startsWith(names(dt_copy), paste0(col_name, "_"))])
  }
  
  # Center specified columns
  for (col_name in covar) {
    dt_copy[, (col_name) := scale(dt_copy[[col_name]], center = TRUE, scale = FALSE)]
  }
  
  return(dt_copy)
}

preprocess_dt_center_covar <- function(dt, covar, make_copy = TRUE) {
  if (make_copy) {
    dt_copy <- copy(dt)
  } else {
    dt_copy <- dt
  }
  
  covar_exists <- covar %in% names(dt_copy)
  if (!all(covar_exists)) {
    stop(paste0("The following columns specified in covar do not exist in the data.table: ", paste(covar[!covar_exists], collapse = ", ")))
  }
  
  categorical_covar <- covar[sapply(dt_copy[, covar, with = FALSE], function(col) is.factor(col) || is.character(col))]
  
  for (col_name in categorical_covar) {
    col <- dt_copy[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt_copy[, (new_col_name) := as.integer(col == level)]
    }
    dt_copy[, (col_name) := NULL]
    covar <- c(covar[covar != col_name], names(dt_copy)[startsWith(names(dt_copy), paste0(col_name, "_"))])
  }
  
  for (col_name in covar) {
    dt_copy[, (col_name) := scale(dt_copy[[col_name]], center = TRUE, scale = FALSE)]
  }
  
  return(dt_copy)
}

# Example Usage:
dt <- data.table(
  A = factor(c("a", "b", "a")),
  B = c(1, 2, 3),
  C = c("x", "y", "z"),
  D = 4:7,
  E = 8:11
)

covar_cols <- c("A", "B", "E")

# Default: make_copy = TRUE
preprocessed_dt_copy <- preprocess_dt_center_covar(dt, covar_cols)
print(preprocessed_dt_copy)
print(dt) # Original 'dt' is unchanged

# make_copy = FALSE
preprocessed_dt_no_copy <- preprocess_dt_center_covar(dt, covar_cols, make_copy = FALSE)
print(preprocessed_dt_no_copy)
print(dt) # Original 'dt' is now modified

library(data.table)

#' Preprocess a data.table: center specified columns (numerical and converted categorical) in place.
#'
#' @param dt A data.table (will be modified directly).
#' @param covar A vector of column names to be centered.
#' @return A character vector containing the names of the columns that were centered.
#'
#' @examples
#' dt <- data.table(
#'   A = factor(c("a", "b", "a")),
#'   B = c(1, 2, 3),
#'   C = c("x", "y", "z"),
#'   D = 4:7,
#'   E = 8:11
#' )
#' covar_cols <- c("A", "B", "E")
#' centered_cols <- preprocess_dt_center_covar_inplace(dt, covar_cols)
#' print(centered_cols)
#' print(dt) # The original 'dt' is now modified
center_covar<- function(dt, covar) {
  # Ensure dt is a data.table
  setDT(dt)
  
  covar_exists <- covar %in% names(dt)
  # if (!all(covar_exists)) {
  #   stop(paste0("The following columns specified in covar do not exist in the data.table: ", paste(covar[!covar_exists], collapse = ", ")))
  # }
  
  categorical_covar <- covar[sapply(dt[, covar, with = FALSE], function(col) is.factor(col) || is.character(col))]
  centered_column_names <- character(0) # Initialize an empty vector to store centered column names
  
  for (col_name in categorical_covar) {
    col <- dt[[col_name]]
    levels_col <- unique(col)
    for (level in levels_col) {
      new_col_name <- paste0(col_name, "_", level)
      dt[, (new_col_name) := as.integer(col == level)]
      centered_column_names <- c(centered_column_names, new_col_name) # Add new column name
    }
    dt[, (col_name) := NULL]
    covar <- c(covar[covar != col_name], names(dt)[startsWith(names(dt), paste0(col_name, "_"))])
  }
  
  for (col_name in covar) {
    dt[, (col_name) := scale(dt[[col_name]], center = TRUE, scale = FALSE)]
    centered_column_names <- c(centered_column_names, col_name) # Add centered column name
  }
  
  return(unique(centered_column_names)) # Return the unique names of centered columns
}

# Example Usage:
dt <- data.table(
  A = factor(c("a", "b", "a")),
  B = c(1, 2, 3),
  C = c("x", "y", "z"),
  D = 4:7,
  E = 8:11
)
original_dt_address <- address(dt) # Get the memory address of the original dt
covar_cols <- c("A", "B", "E")
centered_cols <- preprocess_dt_center_covar_inplace(dt, covar_cols)
print("Centered Columns:", centered_cols)
print("Modified Data Table:")
print(dt)
print(paste("Original data.table was modified:", address(dt) == original_dt_address))

dt2 <- data.table(
  F = c(1,2,3,4),
  G = factor(c("red", "blue", "red", "green")),
  H = 10:13
)
original_dt2_address <- address(dt2)
covar_cols2 <- c("G", "H")
centered_cols2 <- preprocess_dt_center_covar_inplace(dt2, covar_cols2)
print("Centered Columns 2:", centered_cols2)
print("Modified Data Table 2:")
print(dt2)
print(paste("Original data.table 2 was modified:", address(dt2) == original_dt2_address))



detect_re_separation <- function(dt, z_col, id_col, num_ids = NULL, sep_prop = 0, sep_n = 0) {
  # Ensure dt is a data.table
  setDT(dt)
  
  # Count the number of distinct id's
  if (is.null(num_ids)) {
    num_ids <- dt[, uniqueN(get(id_col))]
  }
  
  # Count the number of ids with all z values as either 0 or 1
  count_sep_ids <- dt[, all(get(z_col) == 0) || all(get(z_col) == 1), by = get(id_col)][, sum(V1)]
  
  # Check if the proportion exceeds the threshold
  detection <- count_sep_ids > prop * num_ids && count_sep_ids > sep_n
  
  return(detection)
}

detect_all_zeros <- function(dt, z_col, id_col) {
  setDT(dt)
  dt[, all(get(z_col) == 0), by = get(id_col)][, any(V1)]
}

dt1 <- data.table(
    id = c(1, 1, 2, 2, 3, 3),
    z = c(0, 0, 1, 1, 0, 0)
  )
dt2 <- data.table(
    id = c(1, 1, 2, 2, 3, 3),
    z = c(1, 0, 1, 1, 0, 1)
  )

detect_all_zeros <- function(dt, id_col, id) {
  setDT(dt)
  
  unique_id_col <- unique(dt[[id_col]]) # Get unique values from id_col
  !all(id %in% unique_id_col)
}


vcov_lm <- function(fit, group_names) {
  
}


vcov_logm <- function(fit, group_names)