
data.sim <- sim_count(seed = 100)

usethis::use_data(data.sim, overwrite = TRUE)

lr.sim <- sim_lr(seed = 200, n_lr = 20)

usethis::use_data(lr.sim)