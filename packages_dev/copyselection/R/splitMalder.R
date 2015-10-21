#' this is an internal function that splits the results of a MALDER analysis
#'
#' This function splits the results of a MALDER analysis into useful components.
#' @param x a line from a MALDER output file
#' @keywords Busby_bespoke
#' @export
#' @examples
#' splitmalder(...)

splitMalder <- function(x){
    
    x <- strsplit(x,split="\\+\\/\\-")[[1]]
    a <- as.numeric(x[1])
    aci <- as.numeric(strsplit(x[2],split="\\(")[[1]][1])
    z <- as.numeric(gsub("\\)","",strsplit(x[2],split="\\(Z\\=")[[1]][2]))
    return(c(a,aci,z))
}
