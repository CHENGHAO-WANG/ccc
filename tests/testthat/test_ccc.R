#context("ccc.diff")

# Simulate count matrix (genes x cells) with metadata
# Each gene can be non-DE, cell type-specific DE, or condition-specific DE
# Count values come from a negative binomial distribution with gene-specific dispersion

library(data.table)
library(Matrix)

# -----------------------------
# User-defined parameters
# -----------------------------
set.seed(123)
n_genes <- 1000
n_cell_types <- 3
n_conditions <- 2
n_samples_per_condition <- 4
cells_per_sample <- 100

total_samples <- n_conditions * n_samples_per_condition
sample_ids <- paste0("S", seq_len(total_samples))

# Cell type labels
cell_types <- paste0("CT", seq_len(n_cell_types))

# Metadata: assign cells to samples and cell types
metadata <- data.table(
  cell_id = paste0("cell", seq_len(cells_per_sample * total_samples)),
  sample = rep(sample_ids, each = cells_per_sample),
  condition = rep(rep(c("A", "B"), each = n_samples_per_condition), each = cells_per_sample),
  cell_type = sample(cell_types, size = cells_per_sample * total_samples, replace = TRUE)
)

n_cells <- nrow(metadata)

# -----------------------------
# Gene categories
# -----------------------------
prop_cell_type_specific <- 0.2
prop_condition_specific <- 0.2
prop_sample_effect_genes <- 0.2  # genes with random sample effect

n_ct <- round(n_genes * prop_cell_type_specific)
n_cond <- round(n_genes * prop_condition_specific)
n_sample_eff <- round(n_genes * prop_sample_effect_genes)

ct_specific_genes <- sample(seq_len(n_genes), n_ct)
remaining_genes <- setdiff(seq_len(n_genes), ct_specific_genes)
condition_genes <- sample(remaining_genes, n_cond)
remaining_genes <- setdiff(remaining_genes, condition_genes)
sample_effect_genes <- sample(seq_len(n_genes), n_sample_eff)  # from all genes

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
cell_type_effects <- matrix(0, nrow = n_genes, ncol = n_cell_types)
colnames(cell_type_effects) <- cell_types
rownames(cell_type_effects) <- paste0("G", seq_len(n_genes))

for (ct in cell_types) {
  genes_for_ct <- sample(ct_specific_genes, round(length(ct_specific_genes) / 2))
  cell_type_effects[genes_for_ct, ct] <- rnorm(length(genes_for_ct), mean = 1, sd = 0.3)
}

# -----------------------------
# Condition effects (gamma * random sign)
# -----------------------------
condition_effect <- rep(0, n_genes)
prop_positive_condition_effects <- 0.5
n_condition <- length(condition_genes)
signs <- sample(c(-1, 1), size = n_condition, replace = TRUE,
                prob = c(1 - prop_positive_condition_effects, prop_positive_condition_effects))
effect_sizes <- rgamma(n_condition, shape = 2, scale = 0.5)
condition_effect[condition_genes] <- signs * effect_sizes

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
    cell_type_effects[, ct_idx] +
    ifelse(cell$condition == "B", condition_effect, 0)
  
  mu <- exp(log_mu)
  
  expr_mat[, j] <- rnbinom(n_genes, mu = mu, size = 1 / dispersion)
}

# -----------------------------
# Output
# -----------------------------
expr_mat <- Matrix(expr_mat, sparse = TRUE)
sim.ge <- list(
  counts = expr_mat,
  metadata = metadata,
  gene_info = data.table(
    gene_id = paste0("G", seq_len(n_genes)),
    is_cell_type_specific = seq_len(n_genes) %in% ct_specific_genes,
    is_condition_specific = seq_len(n_genes) %in% condition_genes,
    is_sample_effect = seq_len(n_genes) %in% sample_effect_genes
  )
)


expression_matrix <- log.normalize(count_matrix = as.matrix(expr_mat))
#expression_matrix <- Matrix(expression_matrix, sparse = TRUE)

set.seed(42)

# Step 1: Create gene IDs
all_genes <- paste0("G", 1:1000)

# Step 2: Split into two disjoint subsets
subset_size <- 50
gene1_pool <- sample(all_genes, subset_size)
gene2_pool <- setdiff(all_genes, gene1_pool)
gene2_pool <- sample(gene2_pool, subset_size)

# Step 3: Function to generate complexes or single genes
make_complex <- function(genes) {
  used_genes <- character()  # Track genes that are part of complexes or are single
  
  # Function to create a complex or single gene
  create_complex_or_single <- function() {
    remaining_genes <- setdiff(genes, used_genes)
    
    prob <- runif(1)
    
    if (length(remaining_genes) == 2) {
      prob <- runif(1, min = 0.25, max = 1)
    } else if (length(remaining_genes) == 1) {
      return(remaining_genes)
    }
    
    # 25% chance for a 2-gene complex, 10% for a 3-gene complex
    if (prob < 0.25) {
      complex_genes <- sample(remaining_genes, 2)
      used_genes <<- union(used_genes, complex_genes)
      return(paste(sort(complex_genes), collapse = "_"))
    } else if (prob < 0.35) {
      complex_genes <- sample(remaining_genes, 3)
      used_genes <<- union(used_genes, complex_genes)
      return(paste(sort(complex_genes), collapse = "_"))
    } else {
      single_gene <- sample(remaining_genes, 1)
      used_genes <<- union(used_genes, single_gene)
      return(single_gene)
    }
  }
  
  # Generate complexes or single genes
  complexes_or_singles <- sapply(genes, function(g) create_complex_or_single())
  
  return(complexes_or_singles)
}

# Step 4: Create gene sets with complexes
gene1_set <- unique(make_complex(gene1_pool))
gene2_set <- unique(make_complex(gene2_pool))

# Step 5: Form all possible cross-pairs
all_combos <- expand.grid(gene1 = gene1_set, gene2 = gene2_set, stringsAsFactors = FALSE)

# Step 6: Randomly sample desired number of interactions
interaction_df <- all_combos[sample(nrow(all_combos), 100), ]

# View a few rows
colnames(interaction_df) <- c("ligand", "receptor")

library(future)
oplan <- plan(multisession, workers = 4L)
start.time <- Sys.time()
expression_matrix <- log.normalize(count_matrix = as.matrix(expr_mat))
a <- ccc.diff(expression_matrix = expression_matrix, metadata = metadata,
                  group_col = "condition", id_col = "sample", lr = interaction_df[3:6,],
                  contrast = c(A = 1, B = -1), lmm_re = F, logmm_re = F, verbose = T)
Sys.time() - start.time
plan(oplan)
# expression_matrix <- log.normalize(count_matrix = expr_mat)
# a <- ccc.diff(expression_matrix = expression_matrix, metadata = metadata,
#                   group_col = "condition", id_col = "sample", lr = interaction_df[3:6,],
#                   contrast = c(A = 1, B = -1), lmm_re = F, logmm_re = F, verbose = T)
# 
# 
# start.time <- Sys.time()
# b <- ccc.diff(expression_matrix = expression_matrix, metadata = metadata,
#                   group_col = "condition", id_col = "sample", lr = interaction_df[3:10,],
#                   contrast = c(A = 1, B = -1), logmm_re = F)
# Sys.time() - start.time
