#' this is an internal function that generates time in generations from a date
#'
#' This function converts a date to time in generations.
#' @param x a time in generations
#' @param year0 the year to start counting back from. Defaults to 1950
#' @param gen_length generation length in years. Defaults to 29 years
#' @keywords Busby_bespoke
#' @export
#' @examples
#' makeGenFromDate(20, year0=2000, gen_length=30)

makeGenFromDate <- function(x, year0=1950, gen_length=29){
    y <- round((((year0-x)/gen_length)-1),0)
    return(y)
}
