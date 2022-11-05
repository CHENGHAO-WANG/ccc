
# check the input

precheck <- function(geData.list, ctData.list, group.vec, ref) {
  
  message("Checking data dimensions ... \n")
  
  if (!is.vector(group.vec) & !is.factor(group.vec)) stop("group.vec must be a vector or a factor")
  
  n.grp <- length(unique(group.vec))
  
  if (n.grp > 6) warning(paste0("Too many group levels. <= 6 is recommended. Currently ",n.grp))
  if (n.grp < 2) stop("Less than 2 group levels. Check your group.vec input or try diffLR2")
  
  if (!is.list(geData.list) | !is.list(ctData.list)) {
    if ()
  }
  
  if (!all(sapply(list(length(geData.list),length(ctData.list)), function(x) x == nrow(design_matrix) ))) {
    stop("geData.list, ctData.list, group.vec indicate different number of samples")
  }
  
  
}

prefilter() <- function() {

}
