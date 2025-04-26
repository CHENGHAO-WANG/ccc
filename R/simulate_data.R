


# Simulate count matrix (genes x cells) with metadata
# Each gene can be non-DE, cell type-specific DE, or condition-specific DE
# Count values come from a negative binomial distribution with gene-specific dispersion

library(data.table)
# library(Matrix)


# -----------------------------
# User-defined parameters
# -----------------------------
set.seed(123)
n_genes <- 1000
n_cell_types <- 3
n_conditions <- 3  # user-defined number of conditions
n_samples_per_condition <- 4
cells_per_sample <- rep(30, n_conditions * n_samples_per_condition)

average_logfc <- rep(1, n_conditions - 1)  # user-defined average log fold change for condition effects
prop_positive_condition_effects <- rep(0.5, n_conditions - 1)

# Optional user-provided vector for cell type-specific logFC and positive proportions
logfc_ct <- NULL
prop_pos_ct <- NULL

# Covariate effect proportions
prop_cell_continuous <- 0
prop_cell_discrete <- 0
prop_sample_continuous <- 0
prop_sample_discrete <- 0

conditions <- paste0("cond", seq_len(n_conditions))
condition_labels <- rep(conditions, each = n_samples_per_condition)
sample_ids <- paste0("S", seq_len(length(condition_labels)))
total_samples <- length(sample_ids)

# Check cells_per_sample
stopifnot(length(cells_per_sample) == total_samples)

# Cell type labels
cell_types <- paste0("CT", seq_len(n_cell_types))

# Metadata: assign cells to samples and cell types
metadata <- data.table(
  cell_id = paste0("cell", seq_len(sum(cells_per_sample))),
  sample = rep(sample_ids, times = cells_per_sample),
  condition = rep(condition_labels, times = cells_per_sample),
  cell_type = sample(cell_types, size = sum(cells_per_sample), replace = TRUE)
)

n_cells <- nrow(metadata)

# -----------------------------
# Library size factor
# -----------------------------
lib_size <- rgamma(n_cells, shape = 10, rate = 10)  # mean = 1, variance = 0.1

# -----------------------------
# Simulate Covariates
# -----------------------------
include_cell_cont <- prop_cell_continuous > 0
include_cell_disc <- prop_cell_discrete > 0
include_sample_cont <- prop_sample_continuous > 0
include_sample_disc <- prop_sample_discrete > 0

if (include_cell_cont) {
  metadata$cell_cont <- rnorm(n_cells, mean = 0, sd = 1)
}
if (include_cell_disc) {
  metadata$cell_disc <- sample(c("A", "B", "C"), size = n_cells, replace = TRUE, prob = c(0.45, 0.35, 0.2))
}
if (include_sample_cont || include_sample_disc) {
  sample_covs <- data.table(sample = sample_ids)
  if (include_sample_cont) {
    sample_covs$sample_cont <- rnorm(total_samples, mean = 0, sd = 1)
  }
  if (include_sample_disc) {
    sample_covs$sample_disc <- sample(c("X", "Y"), size = total_samples, replace = TRUE, prob = c(0.5, 0.5))
  }
  metadata <- merge(metadata, sample_covs, by = "sample", sort = FALSE)
}



# -----------------------------
# Gene categories (mutually exclusive)
# -----------------------------
prop_cell_type_specific <- 0.2
prop_condition_specific <- 0.2
prop_sample_effect_genes <- 0.2  # genes with random sample effect

n_ct <- round(n_genes * prop_cell_type_specific)
n_cond <- round(n_genes * prop_condition_specific)
n_sample_eff <- round(n_genes * prop_sample_effect_genes)

stopifnot(n_ct + n_cond + n_sample_eff <= n_genes)

available_genes <- seq_len(n_genes)

ct_specific_genes <- sample(available_genes, n_ct)
available_genes <- setdiff(available_genes, ct_specific_genes)

condition_genes <- sample(available_genes, n_cond)
available_genes <- setdiff(available_genes, condition_genes)

sample_effect_genes <- sample(available_genes, n_sample_eff)

# -----------------------------
# Baseline expression and dispersion
# -----------------------------
baseline_log_mu <- rnorm(n_genes, mean = 1, sd = 0.3)
dispersion <- rgamma(n_genes, shape = 2, rate = 1)  # gene-specific dispersion

# -----------------------------
# Sample random intercepts for sample-effect genes
# -----------------------------
random_intercepts <- matrix(0, n_genes, total_samples)
rownames(random_intercepts) <- paste0("G", seq_len(n_genes))
colnames(random_intercepts) <- sample_ids

random_intercepts[sample_effect_genes, ] <- matrix(
  rnorm(length(sample_effect_genes) * total_samples, mean = 0, sd = 0.3),
  nrow = length(sample_effect_genes),
  ncol = total_samples
)

# -----------------------------
# Cell type effects
# -----------------------------
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
rownames(cell_type_effects) <- paste0("G", seq_len(n_genes))

for (i in seq_along(cell_types)) {
  ct <- cell_types[i]
  genes_for_ct <- sample(ct_specific_genes, round(length(ct_specific_genes) / 2))
  signs <- sample(c(-1, 1), size = length(genes_for_ct), replace = TRUE,
                  prob = c(1 - prop_pos_ct[i], prop_pos_ct[i]))
  effect_sizes <- rgamma(length(genes_for_ct), shape = 4, rate = 4 / logfc_ct[i])
  cell_type_effects[genes_for_ct, ct] <- signs * effect_sizes
}

# -----------------------------
# Condition effects (gamma * random sign)
# -----------------------------
reference_condition <- conditions[1]
condition_effects <- matrix(0, nrow = n_genes, ncol = n_conditions - 1)
colnames(condition_effects) <- conditions[-1]
rownames(condition_effects) <- paste0("G", seq_len(n_genes))

if (is.null(average_logfc) || length(average_logfc) != (n_conditions - 1)) {
  average_logfc <- rgamma(n_conditions - 1, shape = 4, rate = 4)
  message("average_logfc not provided or incorrect length. Simulated instead.")
}
if (length(prop_positive_condition_effects) != (n_conditions - 1)) {
  prop_positive_condition_effects <- rbeta(n_conditions - 1, shape1 = 4, shape2 = 4)
  message("prop_positive_condition_effects incorrect length. Simulated instead.")
}

for (i in seq_along(conditions[-1])) {
  cond <- conditions[-1][i]
  signs <- sample(c(-1, 1), size = length(condition_genes), replace = TRUE,
                  prob = c(1 - prop_positive_condition_effects[i], prop_positive_condition_effects[i]))
  effect_sizes <- rgamma(length(condition_genes), shape = 4, rate = 4 / average_logfc[i])
  condition_effects[condition_genes, cond] <- signs * effect_sizes
}

# -----------------------------
# Covariate effects on genes
# -----------------------------

covariate_effects <- list()

if (include_cell_cont) {
  idx <- sample(seq_len(n_genes), round(n_genes * prop_cell_continuous))
  covariate_effects$cell_cont <- numeric(n_genes)
  covariate_effects$cell_cont[idx] <- rnorm(length(idx), 0, 0.3)
}
if (include_cell_disc) {
  idx <- sample(seq_len(n_genes), round(n_genes * prop_cell_discrete))
  covariate_effects$cell_disc <- numeric(n_genes)
  covariate_effects$cell_disc[idx] <- rnorm(length(idx), 0, 0.3)
}
if (include_sample_cont) {
  idx <- sample(seq_len(n_genes), round(n_genes * prop_sample_continuous))
  covariate_effects$sample_cont <- numeric(n_genes)
  covariate_effects$sample_cont[idx] <- rnorm(length(idx), 0, 0.3)
}
if (include_sample_disc) {
  idx <- sample(seq_len(n_genes), round(n_genes * prop_sample_discrete))
  covariate_effects$sample_disc <- numeric(n_genes)
  covariate_effects$sample_disc[idx] <- rnorm(length(idx), 0, 0.3)
}



# -----------------------------
# Simulate expression
# -----------------------------
expr_mat <- matrix(0, nrow = n_genes, ncol = n_cells)
rownames(expr_mat) <- paste0("G", seq_len(n_genes))
colnames(expr_mat) <- metadata$cell_id

for (j in seq_len(n_cells)) {
  cell <- metadata[j]
  sample_idx <- which(sample_ids == cell$sample)
  ct_idx <- which(cell_types == cell$cell_type)

  log_mu <- baseline_log_mu +
    random_intercepts[, sample_idx] +
    cell_type_effects[, ct_idx]

  if (cell$condition != reference_condition) {
    cond_idx <- which(colnames(condition_effects) == cell$condition)
    log_mu <- log_mu + condition_effects[, cond_idx]
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

# -----------------------------
# Output
# -----------------------------
#expr_mat <- Matrix(expr_mat, sparse = TRUE)
list(
  counts = expr_mat,
  metadata = metadata,
  gene_info = data.table(
    gene_id = paste0("G", seq_len(n_genes)),
    is_cell_type_specific = seq_len(n_genes) %in% ct_specific_genes,
    is_condition_specific = seq_len(n_genes) %in% condition_genes,
    is_sample_effect = seq_len(n_genes) %in% sample_effect_genes
  )
)
