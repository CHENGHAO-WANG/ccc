
dat.sim <- sim_count(seed = 123)

lrdb.sim <- sim_lr(seed = 456, n_lr = 10)

expression_matrix <- log_normalize(dat.sim$counts)
metadata <- dat.sim$metadata

a <- ccc_diff(expression_matrix = expression_matrix, metadata = metadata, contrast = c(grp1 = -1, grp2 = 1),
              id_col = "sample", lr = lrdb.sim, sender = "CT1", receiver = "CT3", logmm_re = FALSE)

a1 <- ccc_diff(expression_matrix = expression_matrix, metadata = metadata, contrast = c(grp1 = -1, grp2 = 1/2, grp3 = 1/2),
              id_col = "sample", lr = lrdb.sim, sender = "CT1", receiver = "CT3", logmm_re = FALSE)
