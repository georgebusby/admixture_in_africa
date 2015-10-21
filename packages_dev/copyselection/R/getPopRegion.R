#' this is an internal function that converts a list of population names to regions
#'
#' This function finds the region that an ethnic group.
#' @param pops a list of names of ethnic groups
#' @param popkey a table with at least two columns: Ethnic_Group which is the population that you want to look up and RegionM which is the region that the Ethnic_Group
#' @param subname should the subnames of populations like FULAI be used to query the table, or should FULA be used? Defaults to TRUE
#' @keywords Busby_bespoke
#' @return returns a vector of regions
#' @export
#' @examples
#' getPopRegion("ESOMALI", popkey)

getPopRegion <- function(pops, popkey, subname=TRUE)
{
    regions <- c()
    for(i in pops)
    {
        ## need to change the colour of Karretjie -- which here represents all of South Africa -- to be same as Malawi
        if(i == "KARRETJIE") i  <- "MALAWI"
        regt <- popkey$RegionM[tidyNames(popkey$Ethnic_Group, fula=subname)==i]
        if(length(regt)==0) regt <- "Eurasia"
        regions <- c(regions,unique(regt))
    }
    return(regions)
}
