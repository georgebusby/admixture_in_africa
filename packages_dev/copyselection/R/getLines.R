#' internal function for pulling specific lines out of a GLOBETROTTER output file
#'
#' This function pulls specified lines out of a GLOBETROTTER output file and is used by getGlobetrotter.R
#' @param x line number containing headers
#' @param y line number containing values
#' @param tmp output file that have been scanned in as a variable
#' @keywords Busby_bespoke
#' @return a two row table with the header and values
#' @export
#' @examples
#' getLines(12,13,tmp)

getLines <- function(x,y,tmp){
    tt <- rbind(unlist(strsplit(tmp[x],split="\\ ")),unlist(strsplit(tmp[y],split="\\ ")))
    return(tt)
}

