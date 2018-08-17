#' Interpolate FD values between slices
#' @param origMaxNoSlices: No. of slices listed in XLS
#' @param interpNoSlices: No. of slices to be interpolated to
interpfunc <- function(x, origMaxNoSlices, interpNoSlices){
    tmp <- spline(1:origMaxNoSlices, x, interpNoSlices)
    t(tmp$y)
}

interpSlices <- function(x, origMaxNoSlices=12, interpNoSlices=10){
    tmp <- apply(x, 1, function(x) {
        interpfunc(x, origMaxNoSlices, interpNoSlices)
    })
    t(tmp)
}
