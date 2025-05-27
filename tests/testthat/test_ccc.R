
dat.sim <- sim_count(seed = 123)

lrdb.sim <- sim_lr(seed = 456, n_lr = 10)

expression_matrix <- log_normalize(dat.sim$counts)
metadata <- dat.sim$metadata

test_that("errors and messages", {
  a <- ccc_diff(expression_matrix = expression_matrix, metadata = metadata,
                id_col = "sample", lr = lrdb.sim, sender = "CT1", receiver = "CT3", 
                logmm_re = FALSE, verbose = FALSE)
  a.result <- ccc_test(ccc_obj = a, contrast = c(grp1 = -1, grp2 = 1/2, grp3 = 1/2),
                       verbose = FALSE)
  b <- ccc_enrich(expression_matrix = expression_matrix, metadata = metadata,
                  id_col = "sample", lr = lrdb.sim, sender = "CT1", receiver = "CT3",
                  logmm_re = FALSE, verbose = FALSE)
  b.result <- ccc_test(ccc_obj = b, verbose = FALSE)
  expect_error(ccc_test(ccc_obj = a, verbose = FALSE))
  expect_error(ccc_test(ccc_obj = a, contrast = c(-1, 1), verbose = FALSE))
  expect_message(ccc_test(ccc_obj = a, contrast = c(grp1 = -1, grp2 = 1), test_type = "chisq", ha = "greater", verbose = FALSE))
  expect_message(ccc_test(ccc_obj = b, contrast = c(grp1 = -1, grp2 = 1), verbose = FALSE))
  expect_error(ccc_test(ccc_obj = b, test_type = "z", ha == "greater.abs", c_linear = -1, verbose = FALSE))
  expect_message(ccc_test(ccc_obj = b, test_type = "chisq", c_linear = 1, verbose = FALSE))
})
