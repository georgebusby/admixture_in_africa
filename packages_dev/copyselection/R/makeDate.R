#' this is an internal function that turns an FS type matrix into rows averaged across clusters/groups
#'
#' This function converts a date in generations to a date in time.
#' @param x a time in generations
#' @param year0 the year to start counting back from. Defaults to 1950
#' @param gen_length generation length in years. Defaults to 29 years
#' @param add_BCE if TRUE converts negative years (i.e. before 0CE) to year in BCE. Defaults to TRUE
#' @keywords Busby_bespoke
#' @export
#' @examples
#' makeDate(20, gen_length=30,)

makeDate <- function(x, year0=1950, gen_length=29, add_BCE=TRUE){
    y <- round(year0-((x+1)*gen_length),0)
    if(is.na(y)){
        y <- NA
    } else   
        if(add_BCE == TRUE & y < 0) {
            y <- paste0(-y,"B")
        }
    return(y)
}