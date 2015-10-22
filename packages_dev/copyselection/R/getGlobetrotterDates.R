#' extract the relevant information from globetrotter's date bootstraps output files.
#'
#' This function extracts dates from GLOBETROTTER's date output files. It assumes that the dating output files have the form POP.0.boot.txt, or in the case of a null bootstrap, it can either be POP.0.null.boot.txt, in which case use null=T, or POP.null.0.txt, in which case use null2=T. It also assumed that 20 bootstraps have been performed labelled 0 to 19. 
#' @param inroot the root of the file before the bootstrap number, e.g. DIRECTORY/POP.globetrotter.main.bt
#' @param num_boots how many different bootstrap files have been generated? Defaults to 20
#' @param suffix everything after the BTSTRAPNUMBER in the file name
#' @keywords Busby_bespoke
#' @return a matrix of dates
#' @export
#' @examples
#' getGlobetrotterDates(...)


getGlobetrotterDates <- function(inroot,num_boots=20,suffix=".txt")
{
    datemat <- matrix(0,nrow=0,ncol=4)
    for(i in 0:num_boots)
    {
        infile <- paste0(inroot,i,suffix)
        if(file.exists(infile)) 
        {
            ld <- read.table(infile, header=T)
            datemat <- rbind(datemat,ld)
        }
    }
    return(datemat)
}


