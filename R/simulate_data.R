#'
#'
#'

sim.count <- function() {
  
}


# Simulate count matrix (genes x cells) with metadata
# Each gene can be non-DE, cell type-specific DE, or group-specific DE
# Count values come from a negative binomial distribution with gene-specific dispersion

library(data.table)
# library(Matrix)


# -----------------------------
# User-defined parameters
# -----------------------------
seed <- 123
if (!is.null(seed)) {
  set.seed(seed)
}
n_genes <- 1000
all_genes <- paste0("G", seq_len(n_genes))

gene_info <- data.table(
  gene_id = all_genes
)

n_cell_types <- 3
cell_type_prop <- NULL
if (is.null(cell_type_prop) || length(cell_type_prop) != n_cell_types) {
  cell_type_prop <- rep(1, times = n_cell_types)
  message("'cell_type_prop' not provided or incorrect length. Simulate cell type labels with equal probabilities.")
}

sample_group <- c(4, 4, 4)
n_groups <- length(sample_group)  # user-defined number of groups
#n_samples_per_group <- 4
cells_per_sample <- rep(30, sum(sample_group))

average_logfc <- rep(1, n_groups - 1)  # user-defined average log fold change for group effects
prop_positive_group_effects <- rep(0.5, n_groups - 1)

# Optional user-provided vector for cell type-specific logFC and positive proportions
logfc_ct <- NULL
prop_pos_ct <- NULL

# Covariate effect proportions
prop_cell_continuous <- 0
prop_cell_discrete <- 0
prop_sample_continuous <- 0
prop_sample_discrete <- 0

# Covariate confounder flags
is_cell_continuous_confounder <- FALSE
is_cell_discrete_confounder <- FALSE
is_sample_continuous_confounder <- FALSE
is_sample_discrete_confounder <- FALSE

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
  cell_type = sample(cell_types, size = sum(cells_per_sample), replace = TRUE, prob = cell_type_prop)
)

n_cells <- nrow(metadata)

# -----------------------------
# Library size factor
# -----------------------------
lib_size_factor_shape <- 3
lib_size_factor_rate <- 3
if (lib_size_factor_shape != lib_size_factor_rate) {
  message("It is recommended that 'lib_size_factor_shape' and 'lib_size_factor_rate' be equal.")
}
lib_size <- rgamma(n_cells, shape = lib_size_factor_shape, rate = lib_size_factor_rate)  # mean = 1, variance = 1/3

# -----------------------------
# Simulate Covariates
# -----------------------------
include_cell_cont <- prop_cell_continuous > 0
include_cell_disc <- prop_cell_discrete > 0
include_sample_cont <- prop_sample_continuous > 0
include_sample_disc <- prop_sample_discrete > 0

# Cell-level continuous
if (include_cell_cont) {
  if (is_cell_continuous_confounder) {
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
  if (is_cell_discrete_confounder) {
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
  if (is_sample_continuous_confounder) {
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
  if (is_sample_discrete_confounder) {
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
prop_cell_type_specific <- 0.2
similarity <- 0.5
prop_group_specific <- 0.2
group_effect_similarity <- 0.5
prop_sample_effect_genes <- 0.8  # genes with random sample effect

# n_ct <- round(n_genes * prop_cell_type_specific)
# n_grp <- round(n_genes * prop_group_specific)
n_sample_eff <- round(n_genes * prop_sample_effect_genes)

# stopifnot(n_ct + n_grp + n_sample_eff <= n_genes)

# available_genes <- seq_len(n_genes)

# ct_specific_genes <- sample(available_genes, n_ct)
# available_genes <- setdiff(available_genes, ct_specific_genes)

# group_genes <- sample(all_genes, n_grp)
# available_genes <- setdiff(available_genes, group_genes)

sample_effect_genes <- sample(all_genes, n_sample_eff)

# -----------------------------
# Baseline expression and dispersion
# -----------------------------
baseline_log_mu_mean <- 1
baseline_log_mu_sd <- 0.3
dispersion_shape <- 2
dispersion_rate <- 1
baseline_log_mu <- rnorm(n_genes, mean = baseline_log_mu_mean, sd = baseline_log_mu_sd)
dispersion <- rgamma(n_genes, shape = dispersion_shape, rate = dispersion_rate)  # gene-specific dispersion

gene_info[, disperion := dispersion, baseline_log_mu := baseline_log_mu]
# -----------------------------
# Sample random intercepts for sample-effect genes
# -----------------------------
random_intercepts <- matrix(0, n_genes, total_samples)
rownames(random_intercepts) <- all_genes
colnames(random_intercepts) <- sample_ids

random_intercepts_sigma_shape <- 3
random_intercepts_sigma_rate <- 2
random_intercepts_sigma <- sqrt(rinvgamma(n = length(sample_effect_genes), shape = random_intercepts_sigma_shape, rate = random_intercepts_sigma_rate))

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
n_ct_specific <- round(n_genes * prop_cell_type_specific)

# Check if the total number of cell type-specific genes exceeds the total number of available genes
total_cell_type_specific <- n_ct_specific * n_celltypes #- similarity * (n_celltypes - 1)

if (total_cell_type_specific > n_genes) {
  stop("The total number of cell type-specific genes exceeds the total number of available genes. Please provide a different value for 'prop_cell_type_specific' or 'n_celltypes'.")
}

if (is.null(logfc_ct) || length(logfc_ct) != n_cell_types) {
  logfc_ct <- rgamma(n_cell_types, shape = 4, rate = 4)
  message("logfc_ct not provided or incorrect length. Simulated instead.")
}
if (is.null(prop_pos_ct) || length(prop_pos_ct) != n_cell_types) {
  prop_pos_ct <- rbeta(n_cell_types, shape1 = 4, shape2 = 4)
  message("prop_pos_ct not provided or incorrect length. Simulated instead.")
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
                  prob = c(1 - prop_pos_ct[i], prop_pos_ct[i]))
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
#                     prob = c(1 - prop_pos_ct[i], prop_pos_ct[i]))
#     effect_sizes <- rgamma(length(genes_for_ct), shape = 4, rate = 4 / logfc_ct[i])
#     cell_type_effects[genes_for_ct, ct] <- signs * effect_sizes
#   }
# }

# -----------------------------
# group effects (gamma * random sign)
# -----------------------------
n_grp_specific <- round(n_genes * prop_group_specific)
# group_genes <- sample(all_genes, n_grp)
# Number of common and true specific genes
n_common_group <- round(group_effect_similarity * n_grp_specific)
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

if (is.null(average_logfc) || length(average_logfc) != (n_groups - 1)) {
  average_logfc <- rgamma(n_groups - 1, shape = 4, rate = 4)
  message("average_logfc not provided or incorrect length. Simulated instead.")
}
if (is.null(prop_positive_group_effects) || length(prop_positive_group_effects) != (n_groups - 1)) {
  prop_positive_group_effects <- rbeta(n_groups - 1, shape1 = 4, shape2 = 4)
  message("prop_positive_group_effects incorrect length. Simulated instead.")
}

for (i in seq_along(groups)) {
  grp <- groups[i]
  group_genes <- group_genes_list[[i-1]]
  if (i == 1) {
    # group_effects[group_genes, grp] <- 0
    next
  }
  signs <- sample(c(-1, 1), size = n_grp_specific, replace = TRUE,
                  prob = c(1 - prop_positive_group_effects[i-1], prop_positive_group_effects[i-1]))
  effect_sizes <- rgamma(n_grp_specific, shape = 4, rate = 4 / average_logfc[i-1])
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
  idx <- sample(seq_len(n_genes), round(n_genes * prop_cell_continuous))
  covariate_effects$cell_cont <- numeric(n_genes)
  covariate_effects$cell_cont[idx] <- rnorm(length(idx), 0, 0.3)
  colname <- "cell_level_continuous_covariate_association"
  gene_info[, (colname) := FALSE]
  gene_info[idx, (colname) := TRUE]
}
if (include_cell_disc) {
  idx <- sample(seq_len(n_genes), round(n_genes * prop_cell_discrete))
  covariate_effects$cell_disc <- numeric(n_genes)
  covariate_effects$cell_disc[idx] <- rnorm(length(idx), 0, 0.3)
  colname <- "cell_level_discrete_covariate_association"
  gene_info[, (colname) := FALSE]
  gene_info[idx, (colname) := TRUE]
}
if (include_sample_cont) {
  idx <- sample(seq_len(n_genes), round(n_genes * prop_sample_continuous))
  covariate_effects$sample_cont <- numeric(n_genes)
  covariate_effects$sample_cont[idx] <- rnorm(length(idx), 0, 0.3)
  colname <- "sample_level_continuous_covariate_association"
  gene_info[, (colname) := FALSE]
  gene_info[idx, (colname) := TRUE]
}
if (include_sample_disc) {
  idx <- sample(seq_len(n_genes), round(n_genes * prop_sample_discrete))
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
  expr_mat[, j] <- rnbinom(n_genes, mu = mu, size = 1 / dispersion)
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
