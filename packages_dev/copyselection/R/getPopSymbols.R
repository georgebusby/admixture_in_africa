#' find colour, pch, and rim colour for an ethnic group
#'
#' This function finds various plotting parameters for a given ethnic group, as outlined in the leginfo object
#' @param pop_vec a vector (can be length 1) of ethnic groups
#' @param leginfo an object containing data to match to the ethnic group, must have the following columns: "EthnicGroup" "Country" "Region" "Colour" "pch" "rim", where rim defines whether the rim is plotted as the same colour as the symbol
#' @keywords Busby_bespoke
#' @return a dataframe containing the ethnic group, colour, pch, and rim colour for points.
#' @export
#' @examples
#' getPopSymbols("FULA",leginfo)



getPopSymbols <- function(pop_vec,leginfo)
{
    col2plot <- pch2plot <- rim2plot <- c()
    for(i in pop_vec)
    {
        col2plot <- c(col2plot,as.character(leginfo$Colour[leginfo$EthnicGroup==i]))
        pch2plot <- c(pch2plot,as.numeric(leginfo$poppch[leginfo$EthnicGroup==i]))
        rim2plot <- c(rim2plot,as.numeric(leginfo$rim[leginfo$EthnicGroup==i]))
    }
    rim2plot[rim2plot==0] <- "#000000"
    rim2plot[rim2plot==1] <- as.character(col2plot)[rim2plot==1]
    out_mat <- data.frame(cbind(pop_vec,col2plot,as.numeric(pch2plot),rim2plot),stringsAsFactors=F)
    colnames(out_mat)[3] <- "pch2plot"
    return(out_mat)
}
