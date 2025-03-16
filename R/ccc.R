
# not export

difflr0 <- function(geData.list, ctData.list, celltypes.L, celltypes.R, group.vec, ref,
                   lrdb = NULL,
                   re.cont = TRUE, re.disc = TRUE,
                   min.prop = 0.10, filter.all = FALSE,
                   verbose = TRUE) {

  env <- environment()

  # if (is.null(lrdb)) {
  #   data("lrdb", envir=env)
  # }
  #
  # precheck(geData.list = geData.list, ctData.list = ctData.list, group.vec = group.vec, ref = ref)
  #
  # REname(geData.list = geData.list, ctData.list = ctData.list, lrdb = lrdb, env = env)

  preprocess(geData.list = geData.list, ctData.list = ctData.list, lrdb = lrdb, celltypes.L = celltypes.L,
             celltypes.R = celltypes.R, env = env)

  existing.lr <- prefilter(lrdb = lrdb, min.prop = min.prop, filter.all = filter.all)



}


# main function
# geData.list is gene expression matrices; row for genes; cols for cells; rownames are gene names
# ctData.list is cell type information
# group.vec is a vector indicating sample group assignment
# celltypes.L is the cell type of interest expressing ligand
# celltypes.R is the cell type of interest expressing receptor
# lrdb is the database for ligand-receptor pairs

difflr <- function(geData.list, ctData.list, celltypes.L = NULL, celltypes.R = NULL,
           group.vec, ref, test.all = TRUE,
           lrdb = NULL,
           re.cont = TRUE, re.disc = TRUE,
           min.prop = 0.10, filter.all = FALSE,
           verbose = TRUE) {

  env <- environment()

  if (is.null(lrdb)) {
    data("lrdb", envir=env)
  }
  #browser()
  precheck(geData.list = geData.list, ctData.list = ctData.list, group.vec = group.vec, ref = ref, env = env)

  REname(geData.list = geData.list, ctData.list = ctData.list, lrdb = lrdb, env = env)

  if (is.null(celltypes.L) & is.null(celltypes.R)) {

    celltype.pairs <- get_celltype_pairs(ctData.list)

    # ......

  } else if (!is.null(celltypes.L) & !is.null(celltypes.R)) {

    preprocess(geData.list = geData.list, ctData.list = ctData.list, lrdb = lrdb, celltypes.L = celltypes.L, celltypes.R = celltypes.R, env = env)

    existing.lr <- prefilter(lrdb = lrdb, min.prop = min.prop, filter.all = filter.all, geData.L = geData.L, geData.R = geData.R)

    if (verbose) {
      stat_analysis
    } else {
      suppressMessages( stat_analysis() )
    }





  } else {
    stop("Please specify celltypes.L and celltypes.R, or set both as NULL to test all cell type pairs")
  }






}
