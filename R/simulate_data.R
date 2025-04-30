#' Simulate Data
#'
#' Simulate multi-sample multi-group single-cell RNA-seq data using negative binomial distribution.
#' 
#' @param seed random seed to use.
#' @param n_genes number of genes to simulate.
#' @param n_cell_types number of cell types to simulate.
#' @param cell_type_prob a vector of length `n_cell_types` containing weights for simulating number of cells in each cell type. Need not sum to 1.
#' @param sample_group a vector specifying the number of samples in each group. The length of this vector specifies the total number of groups. The first group is considered as the reference group.
#' @param cells_per_sample a vector specifying the number of cells in each sample.
#' @param re_prop numeric specifies the proportion of genes with sample-level random intercepts.
#' @param re_sigma_shape,re_sigma_rate Simulate gene-wise standard deviation for sample-level random intercepts from a inverse Gamma distribution with provided shape and rate. 
#' @param grp_specific_prop numeric specifies the proportions of differentially expressed (DE) genes between each non-reference group and the reference group.
#' @param grp_similarity numeric in [0, 1] specifies the overlap across the DE genes in each non-reference group.
#' @param avg_logfc a vector of length `length(sample_group) - 1`. Specifies the average log fold change for each non-reference group compared to the reference group.
#' @param up_grp_prob a vector of length `length(sample_group) - 1`. Specifies the probabilities of up-regulation in each non-reference group.
#' @param ct_specific_prop numeric specifies the proportions of cell-type-specific genes in each cell type. 
#' @param logfc_ct a vector of length `n_cell_types`. Specifies the average log fold change for each cell type.
#' @param up_ct_prob a vector of length `n_cell_types`. Specifies the probabilities of up-regulation in each cell type.
#' @param cell_cont_prop,cell_disc_prop,sample_cont_prop,sample_disc_prop proportions of genes associated with cell-level continuous covariate, cell-level discrete covariate, sample-level continuous covariate, and sample-level discrete covariate respectively. Default to 0.
#' @param is_cc_confounder,is_cd_confounder,is_sc_confounder,is_sd_confounder logical specifying whether each covariate and group effects are confounded. Default to FALSE.
#' @param baseline_log_mu_mean,baseline_log_mu_sd  Simulate baseline log mean for each gene from a normal distribution with provided mean and sd.
#' @param phi_shape,phi_rate Simulate the reciprocal of gene-wise dispersion parameter from a Gamma distribution with provided shape and rate.
#' @param lib_size_factor_shape,lib_size_factor_rate Simulate library size factors for each cell from a Gamma distribution with provided shape and rate. Recommend these be equal. 
#' 
#' @returns descriptions
#' 
#' @examples
#' # example code
#' 
#' @export


sim.count <- function(seed = NULL, n_genes = 1000, n_cell_types = 3, cell_type_prob = rep(1, n_cell_types),
                      sample_group = c(4, 4, 4), cells_per_sample = rep(100, sum(sample_group)),
                      re_prop = 0.8, re_sigma_shape = 3, re_sigma_rate = 2,
                      grp_specific_prop = 0.2, grp_similarity = 0.5, avg_logfc = rep(1, length(sample_group) - 1), up_grp_prob = rep(0.5, length(sample_group) - 1),
                      ct_specific_prop = 0.05, logfc_ct = rep(1, n_cell_types), up_ct_prob = rep(1, n_cell_types),
                      cell_cont_prop = 0, cell_disc_prop = 0, sample_cont_prop = 0, sample_disc_prop = 0,
                      is_cc_confounder = FALSE, is_cd_confounder = FALSE, is_sc_confounder = FALSE, is_sd_confounder = FALSE,
                      baseline_log_mu_mean = 1, baseline_log_mu_sd = 0.3,
                      phi_shape = 2, phi_rate = 1,
                      lib_size_factor_shape = 3, lib_size_factor_rate = 3) {
  
  # -----------------------------
  # Setup
  # -----------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  # n_genes <- 1000
  all_genes <- paste0("G", seq_len(n_genes))
  
  gene_info <- data.table(
    gene_id = all_genes
  )
  
  # n_cell_types <- 3
  # cell_type_prob <- NULL
  if (length(cell_type_prob) != n_cell_types) {
    cell_type_prob <- rep(1, n_cell_types)
    message("'cell_type_prob' of incorrect length. Use the default instead.")
  }
  
  # sample_group <- c(4, 4, 4)
  n_groups <- length(sample_group)
  #n_samples_per_group <- 4
  # cells_per_sample <- rep(30, sum(sample_group))
  
  # avg_logfc <- rep(1, n_groups - 1)  # user-defined average log fold change for group effects
  # up_grp_prob <- rep(0.5, n_groups - 1)
  
  # Optional user-provided vector for cell type-specific logFC and positive proportions
  # logfc_ct <- NULL
  # up_ct_prob <- NULL
  
  # Covariate effect proportions
  # cell_cont_prop <- 0
  # cell_disc_prop <- 0
  # sample_cont_prop <- 0
  # sample_disc_prop <- 0
  
  # Covariate confounder flags
  # is_cc_confounder <- FALSE
  # is_cd_confounder <- FALSE
  # is_sc_confounder <- FALSE
  # is_sd_confounder <- FALSE
  
  groups <- paste0("grp", seq_len(n_groups))
  group_labels <- rep(groups, times = sample_group)
  sample_ids <- paste0("S", seq_len(length(group_labels)))
  total_samples <- length(sample_ids)
  
  # Check cells_per_sample
  stopifnot(length(cells_per_sample) == total_samples)
  
  # Cell type labels
  cell_types <- paste0("CT", seq_len(n_cell_types))
  
  # Metadata: assign cells to samples and cell types
  metadata <- data.table(
    cell_id = paste0("cell", seq_len(sum(cells_per_sample))),
    sample = rep(sample_ids, times = cells_per_sample),
    group = rep(group_labels, times = cells_per_sample),
    cell_type = sample(cell_types, size = sum(cells_per_sample), replace = TRUE, prob = cell_type_prob)
  )
  
  n_cells <- nrow(metadata)
  
  # -----------------------------
  # Library size factor
  # -----------------------------
  # lib_size_factor_shape <- 3
  # lib_size_factor_rate <- 3
  if (lib_size_factor_shape != lib_size_factor_rate) {
    message("It is recommended that 'lib_size_factor_shape' and 'lib_size_factor_rate' be equal.")
  }
  lib_size <- rgamma(n_cells, shape = lib_size_factor_shape, rate = lib_size_factor_rate)  # mean = 1, variance = 1/3
  
  # -----------------------------
  # Simulate Covariates
  # -----------------------------
  include_cell_cont <- cell_cont_prop > 0
  include_cell_disc <- cell_disc_prop > 0
  include_sample_cont <- sample_cont_prop > 0
  include_sample_disc <- sample_disc_prop > 0
  
  # Cell-level continuous
  if (include_cell_cont) {
    if (is_cc_confounder) {
      mu_per_group <- rnorm(n_groups, mean = 0, sd = 1)
      var_per_group <- rinvgamma(n_groups, shape = 3, scale = 2)
      metadata[, cell_continuous := rnorm(.N, mean = mu_per_group[match(group, groups)], sd = sqrt(var_per_group[match(group, groups)]))]
    } else {
      metadata[, cell_continuous := rnorm(.N, mean = 0, sd = 1)]
    }
  }
  # if (include_cell_cont) {
  #   metadata$cell_cont <- rnorm(n_cells, mean = 0, sd = 1)
  # }
  # Cell-level discrete
  if (include_cell_disc) {
    levels <- c("A", "B", "C")
    if (is_cd_confounder) {
      probs_list <- lapply(seq_len(n_groups), function(i) {
        p <- runif(3)
        p / sum(p)
      })
      metadata[, cell_discrete := mapply(function(grp) sample(levels, 1, prob = probs_list[[match(grp, groups)]]), group)]
    } else {
      metadata[, cell_discrete := sample(levels, .N, replace = TRUE, prob = c(0.45, 0.35, 0.2))]
    }
  }
  # if (include_cell_disc) {
  #   metadata$cell_disc <- sample(c("A", "B", "C"), size = n_cells, replace = TRUE, prob = c(0.45, 0.35, 0.2))
  # }
  
  # Sample-level continuous
  if (include_sample_cont) {
    if (is_sc_confounder) {
      mu_per_group <- rnorm(n_groups, mean = 0, sd = 1)
      var_per_group <- rinvgamma(n_groups, shape = 3, scale = 2)
      sample_covariate <- data.table(sample = sample_ids, group = group_labels)
      sample_covariate[, sample_continuous := rnorm(.N, mean = mu_per_group[match(group, groups)], sd = sqrt(var_per_group[match(group, groups)]))]
      metadata <- merge(metadata, sample_covariate[, .(sample, sample_continuous)], by = "sample")
    } else {
      sample_covariate <- data.table(sample = sample_ids)
      sample_covariate[, sample_continuous := rnorm(.N, mean = 0, sd = 1)]
      metadata <- merge(metadata, sample_covariate, by = "sample")
    }
  }
  # Sample-level discrete
  if (include_sample_disc) {
    levels <- c("X", "Y")
    if (is_sd_confounder) {
      probs_list <- lapply(seq_len(n_groups), function(i) {
        p <- runif(2)
        p / sum(p)
      })
      sample_covariate <- data.table(sample = sample_ids, group = group_labels)
      sample_covariate[, sample_discrete := mapply(function(grp) sample(levels, 1, prob = probs_list[[match(grp, groups)]]), group)]
      metadata <- merge(metadata, sample_covariate[, .(sample, sample_discrete)], by = "sample")
    } else {
      sample_covariate <- data.table(sample = sample_ids)
      sample_covariate[, sample_discrete := sample(levels, .N, replace = TRUE, prob = c(0.5, 0.5))]
      metadata <- merge(metadata, sample_covariate, by = "sample")
    }
  }
  # if (include_sample_cont || include_sample_disc) {
  #   sample_covs <- data.table(sample = sample_ids)
  #   if (include_sample_cont) {
  #     sample_covs$sample_cont <- rnorm(total_samples, mean = 0, sd = 1)
  #   }
  #   if (include_sample_disc) {
  #     sample_covs$sample_disc <- sample(c("X", "Y"), size = total_samples, replace = TRUE, prob = c(0.5, 0.5))
  #   }
  #   metadata <- merge(metadata, sample_covs, by = "sample", sort = FALSE)
  # }
  
  # -----------------------------
  # Gene categories
  # -----------------------------
  # ct_specific_prop <- 0.2
  # similarity <- 0.5
  # grp_specific_prop <- 0.2
  # grp_similarity <- 0.5
  # re_prop <- 0.8  # genes with random sample effect
  
  # n_ct <- round(n_genes * ct_specific_prop)
  # n_grp <- round(n_genes * grp_specific_prop)
  n_sample_eff <- round(n_genes * re_prop)
  
  # stopifnot(n_ct + n_grp + n_sample_eff <= n_genes)
  
  # available_genes <- seq_len(n_genes)
  
  # ct_specific_genes <- sample(available_genes, n_ct)
  # available_genes <- setdiff(available_genes, ct_specific_genes)
  
  # group_genes <- sample(all_genes, n_grp)
  # available_genes <- setdiff(available_genes, group_genes)
  
  sample_effect_genes <- sample(all_genes, n_sample_eff)
  
  # -----------------------------
  # Baseline expression and phi
  # -----------------------------
  # baseline_log_mu_mean <- 1
  # baseline_log_mu_sd <- 0.3
  # phi_shape <- 2
  # phi_rate <- 1
  baseline_log_mu <- rnorm(n_genes, mean = baseline_log_mu_mean, sd = baseline_log_mu_sd)
  phi <- rgamma(n_genes, shape = phi_shape, rate = phi_rate)  # gene-specific phi
  
  gene_info[, dispersion := 1/phi, baseline_log_mu := baseline_log_mu]
  # -----------------------------
  # Sample random intercepts for sample-effect genes
  # -----------------------------
  random_intercepts <- matrix(0, n_genes, total_samples)
  rownames(random_intercepts) <- all_genes
  colnames(random_intercepts) <- sample_ids
  
  # re_sigma_shape <- 3
  # re_sigma_rate <- 2
  random_intercepts_sigma <- sqrt(rinvgamma(n = length(sample_effect_genes), shape = re_sigma_shape, rate = re_sigma_rate))
  
  # random_intercepts[sample_effect_genes, ] <- matrix(
  #   rnorm(length(sample_effect_genes) * total_samples, mean = 0, sd = 0.3),
  #   nrow = length(sample_effect_genes),
  #   ncol = total_samples
  # )
  random_intercepts[sample_effect_genes, ] <- matrix(rnorm(n_sample_eff * total_samples), nrow = n_sample_eff, ncol = total_samples) * random_intercepts_sigma 
  
  gene_info[, within_sample_correlation := FALSE]  # Initialize to FALSE
  gene_info[gene_id %in% sample_effect_genes, within_sample_correlation := TRUE]
  
  # -----------------------------
  # Cell type effects
  # -----------------------------
  n_ct_specific <- round(n_genes * ct_specific_prop)
  
  # Check if the total number of cell type-specific genes exceeds the total number of available genes
  total_cell_type_specific <- n_ct_specific * n_celltypes #- similarity * (n_celltypes - 1)
  
  if (total_cell_type_specific > n_genes) {
    stop("The total number of cell type-specific genes exceeds the total number of available genes. Please provide a different value for 'ct_specific_prop' or 'n_celltypes'.")
  }
  
  logfc_ct = rep(1, n_cell_types)
  # if (is.null(logfc_ct) || length(logfc_ct) != n_cell_types) {
  #   logfc_ct <- rgamma(n_cell_types, shape = 4, rate = 4)
  #   message("logfc_ct not provided or incorrect length. Simulated instead.")
  # }
  # if (is.null(up_ct_prob) || length(up_ct_prob) != n_cell_types) {
  #   up_ct_prob <- rbeta(n_cell_types, shape1 = 4, shape2 = 4)
  #   message("up_ct_prob not provided or incorrect length. Simulated instead.")
  # }
  if (length(logfc_ct) != n_cell_types) {
    logfc_ct <- rep(1, n_cell_types)
    message("'logfc_ct' of incorrect length. Use the default instead.")
  }
  if (length(up_ct_prob) != n_cell_types) {
    up_ct_prob <- rep(1, n_cell_types)
    message("'up_ct_prob' of incorrect length. Use the default instead.")
  }
  
  cell_type_effects <- matrix(0, nrow = n_genes, ncol = n_cell_types)
  colnames(cell_type_effects) <- cell_types
  rownames(cell_type_effects) <- all_genes
  
  available_genes <- all_genes
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    # genes_for_ct <- sample(ct_specific_genes, round(length(ct_specific_genes) / 2))
    genes_for_ct <- sample(available_genes, n_ct_specific)
    available_genes <- setdiff(available_genes, genes_for_ct)
    signs <- sample(c(-1, 1), size = length(genes_for_ct), replace = TRUE,
                    prob = c(1 - up_ct_prob[i], up_ct_prob[i]))
    effect_sizes <- rgamma(length(genes_for_ct), shape = 4, rate = 4 / logfc_ct[i])
    cell_type_effects[genes_for_ct, ct] <- signs * effect_sizes
    
    colname <- paste0(ct, "_logfc")
    gene_info[, (colname) := 0]  # Initialize to 0
    gene_info[gene_id %in% genes_for_ct, (colname) := (signs * effect_sizes)[match(gene_id, genes_for_ct)]]
  }
  
  # if (n_ct_specific > 0) {
  #   is_ct_specific <- matrix(FALSE, nrow = n_genes, ncol = n_celltypes)
  #   colnames(is_ct_specific) <- cell_types
  #   
  #   n_common <- round(n_ct_specific * similarity)
  #   n_real_per_ct <- n_ct_specific - n_common
  #   
  #   common_genes <- sample(all_genes, n_common)
  #   
  #   available_genes <- setdiff(all_genes, common_genes)
  #   cell_type_genes_list <- vector("list", n_celltypes)
  #   for (i in seq_len(n_celltypes)) {
  #     selected <- sample(available_genes, n_real_per_ct)
  #     cell_type_genes_list[[i]] <- c(selected, common_genes)
  #     available_genes <- setdiff(available_genes, selected)
  #   }
  #   ##
  #   for (i in seq_len(n_celltypes)) {
  #     is_ct_specific[all_genes %in% cell_type_genes_list[[i]], i] <- TRUE
  #   }
  #   
  #   cell_type_effects <- matrix(0, nrow = n_genes, ncol = n_cell_types)
  #   colnames(cell_type_effects) <- cell_types
  #   rownames(cell_type_effects) <- all_genes
  #   
  #   for (i in seq_along(cell_types)) {
  #     ct <- cell_types[i]
  #     # genes_for_ct <- sample(ct_specific_genes, round(length(ct_specific_genes) / 2))
  #     genes_for_ct <- cell_type_genes_list[[i]]
  #     signs <- sample(c(-1, 1), size = length(genes_for_ct), replace = TRUE,
  #                     prob = c(1 - up_ct_prob[i], up_ct_prob[i]))
  #     effect_sizes <- rgamma(length(genes_for_ct), shape = 4, rate = 4 / logfc_ct[i])
  #     cell_type_effects[genes_for_ct, ct] <- signs * effect_sizes
  #   }
  # }
  
  # -----------------------------
  # group effects (gamma * random sign)
  # -----------------------------
  n_grp_specific <- round(n_genes * grp_specific_prop)
  # group_genes <- sample(all_genes, n_grp)
  # Number of common and true specific genes
  n_common_group <- round(grp_similarity * n_grp_specific)
  n_true_specific_group <- n_group_specific - n_common_group
  
  total_needed <- n_common_group + n_true_specific_group * (n_groups - 1)
  
  if (total_needed > n_genes) {
    stop("The total number of group-specific genes exceeds the total number of available genes. Please provide a different value for 'group_effect_specific', 'similarity', or 'sample_group'.")
  }
  
  common_group_genes <- sample(all_genes, n_common_group)
  available_genes <- setdiff(all_genes, common_group_genes)
  
  group_genes_list <- vector("list", n_groups - 1)
  
  for (i in seq_len(n_groups - 1)) {
    selected_genes <- sample(available_genes, n_true_specific_group)
    group_genes_list[[i]] <- c(selected_genes, common_group_genes)
    available_genes <- setdiff(available_genes, selected_genes)
  }
  
  reference_group <- groups[1]
  group_effects <- matrix(0, nrow = n_genes, ncol = n_groups)
  colnames(group_effects) <- groups
  rownames(group_effects) <- all_genes
  
  # if (is.null(avg_logfc) || length(avg_logfc) != (n_groups - 1)) {
  #   avg_logfc <- rgamma(n_groups - 1, shape = 4, rate = 4)
  #   message("avg_logfc not provided or incorrect length. Simulated instead.")
  # }
  # if (is.null(up_grp_prob) || length(up_grp_prob) != (n_groups - 1)) {
  #   up_grp_prob <- rbeta(n_groups - 1, shape1 = 4, shape2 = 4)
  #   message("up_grp_prob incorrect length. Simulated instead.")
  # }
  if (length(avg_logfc) != (n_groups - 1)) {
    avg_logfc <- rep(1, n_groups - 1)
    message("'avg_logfc' of incorrect length. Use the default instead.")
  }
  if (is.null(up_grp_prob) || length(up_grp_prob) != (n_groups - 1)) {
    up_grp_prob <- rep(0.5, n_groups - 1)
    message("'up_grp_prob' of incorrect length. Use the default instead.")
  }
  
  
  for (i in seq_along(groups)) {
    grp <- groups[i]
    group_genes <- group_genes_list[[i-1]]
    if (i == 1) {
      # group_effects[group_genes, grp] <- 0
      next
    }
    signs <- sample(c(-1, 1), size = n_grp_specific, replace = TRUE,
                    prob = c(1 - up_grp_prob[i-1], up_grp_prob[i-1]))
    effect_sizes <- rgamma(n_grp_specific, shape = 4, rate = 4 / avg_logfc[i-1])
    group_effects[group_genes, grp] <- signs * effect_sizes
    
    colname <- paste0(grp, "_logfc")
    gene_info[, (colname) := 0]  # Initialize to 0
    gene_info[gene_id %in% group_genes, (colname) := (signs * effect_sizes)[match(gene_id, group_genes)]]
  }
  
  # -----------------------------
  # Covariate effects on genes
  # -----------------------------
  
  covariate_effects <- list()
  
  if (include_cell_cont) {
    idx <- sample(seq_len(n_genes), round(n_genes * cell_cont_prop))
    covariate_effects$cell_cont <- numeric(n_genes)
    covariate_effects$cell_cont[idx] <- rnorm(length(idx), 0, 0.3)
    colname <- "cell_level_continuous_covariate_association"
    gene_info[, (colname) := FALSE]
    gene_info[idx, (colname) := TRUE]
  }
  if (include_cell_disc) {
    idx <- sample(seq_len(n_genes), round(n_genes * cell_disc_prop))
    covariate_effects$cell_disc <- numeric(n_genes)
    covariate_effects$cell_disc[idx] <- rnorm(length(idx), 0, 0.3)
    colname <- "cell_level_discrete_covariate_association"
    gene_info[, (colname) := FALSE]
    gene_info[idx, (colname) := TRUE]
  }
  if (include_sample_cont) {
    idx <- sample(seq_len(n_genes), round(n_genes * sample_cont_prop))
    covariate_effects$sample_cont <- numeric(n_genes)
    covariate_effects$sample_cont[idx] <- rnorm(length(idx), 0, 0.3)
    colname <- "sample_level_continuous_covariate_association"
    gene_info[, (colname) := FALSE]
    gene_info[idx, (colname) := TRUE]
  }
  if (include_sample_disc) {
    idx <- sample(seq_len(n_genes), round(n_genes * sample_disc_prop))
    covariate_effects$sample_disc <- numeric(n_genes)
    covariate_effects$sample_disc[idx] <- rnorm(length(idx), 0, 0.3)
    colname <- "sample_level_discrete_covariate_association"
    gene_info[, (colname) := FALSE]
    gene_info[idx, (colname) := TRUE]
  }
  
  
  
  # -----------------------------
  # Simulate expression
  # -----------------------------
  expr_mat <- matrix(0, nrow = n_genes, ncol = n_cells)
  rownames(expr_mat) <- all_genes
  colnames(expr_mat) <- metadata$cell_id
  
  for (j in seq_len(n_cells)) {
    cell <- metadata[j]
    sample_idx <- which(sample_ids == cell$sample)
    ct_idx <- which(cell_types == cell$cell_type)
    
    log_mu <- baseline_log_mu +
      random_intercepts[, sample_idx] +
      cell_type_effects[, ct_idx]
    
    if (cell$group != reference_group) {
      grp_idx <- which(colnames(group_effects) == cell$group)
      log_mu <- log_mu + group_effects[, grp_idx]
    }
    
    if (include_cell_cont) {
      log_mu <- log_mu + covariate_effects$cell_cont * cell$cell_cont
    }
    if (include_cell_disc) {
      log_mu <- log_mu + covariate_effects$cell_disc * as.numeric(cell$cell_disc == "B")
    }
    if (include_sample_cont) {
      log_mu <- log_mu + covariate_effects$sample_cont * cell$sample_cont
    }
    if (include_sample_disc) {
      log_mu <- log_mu + covariate_effects$sample_disc * as.numeric(cell$sample_disc == "Y")
    }
    
    mu <- exp(log_mu) * lib_size[j]
    expr_mat[, j] <- rnbinom(n_genes, mu = mu, size = 1 / phi)
  }
  
  
  
  if () {
    gene_info <- data.table(
      is_cell_type_specific = seq_len(n_genes) %in% ct_specific_genes,
      #is_group_specific = seq_len(n_genes) %in% group_genes,
      is_sample_effect = seq_len(n_genes) %in% sample_effect_genes
    )
    
  }
  
  if (n_groups >= 2) {
    
  }
  rownames(expr_mat) <- paste0(prefix.g, all_genes, suffix.g)
  gene_info[, gene_id := paste0(prefix.g, gene_id, suffix.g)]
  gene_info <- as.data.frame(gene_info)
  # -----------------------------
  # Output
  # -----------------------------
  #expr_mat <- Matrix(expr_mat, sparse = TRUE)
  list(
    counts = expr_mat,
    metadata = metadata,
    gene_info = gene_info
  )
  
}


