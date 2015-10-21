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

rowsAsMapClusts <- function(x,y1,stat=mean){# x is clusts as MapState, y1 is matrix
    yMat <- c()
    for(i in 1:length(x)){
        clust.name <- names(x[i])
        x1 <- as.vector(unlist(x[i]))
        x1.length <- length(x1)
        y.clust <- y1[x1,]
        if(x1.length>1){
            y2 <- apply(y.clust,2,stat)
        } else {
            y2 <- y.clust
        }
        yMat <- rbind(yMat,y2)
        rownames(yMat)[length(yMat[,1])] <- clust.name
    }
    return(yMat)
}