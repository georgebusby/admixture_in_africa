#' this is an internal function that splits the results of a MALDER analysis
#'
#' Converts a Z-score to a p-value.
#' @param x a Z-score
#' @keywords Busby_bespoke
#' @export
#' @examples
#' z2p(3.4)

z2p <- function(x){
    p <- 2*pnorm(-abs(x))
}
