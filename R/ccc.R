

difflr <- function(geData.list, ctData.list, celltypes.L, celltypes.R, group.vec, ref,
                   lrdb = NULL,
                   re.cont = TRUE, re.disc = TRUE,
                   min.prop = 0.10, filter.all = FALSE,
                   verbose = TRUE) {

  env <- environment()

  data("lrdb", envir=env)

  precheck()
  Rename()
  preprocess()

}


