#' this is an internal function that converts a list of population names to regions
#'
#' This function finds the region that an ethnic group.
#' @param pops a list of names of ethnic groups
#' @param popkey a table with at least two columns: Ethnic_Group which is the population that you want to look up and RegionM which is the region that the Ethnic_Group
#' @keywords Busby_bespoke
#' @export returns a vector of regions
#' @examples
#' findPopsRegion("ESOMALI", popkey)

findPopRegion <- function(pops, popkey)
{
    regions <- c()
    for(i in pops)
    {
        ## need to change the colour of Karretjie -- which here represents all of South Africa -- to be same as Malawi
        if(i == "KARRETJIE") i  <- "MALAWI"
        regt <- popkey$RegionM[tidyNames(popkey$Ethnic_Group)==i]
        if(length(regt)==0) regt <- "Eurasia"
        regions <- c(regions,unique(regt))
    }
    return(regions)
}
