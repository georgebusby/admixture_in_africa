#' this is an internal function that turns an FS type matrix into rows averaged across clusters/groups
#'
#' This function converts a an fineSTRUCTURE type matrix into rows averaged across clusters/groups.
#' @param x a list of clusters, as output as a MapState in D. Lawson's FinestructureLibrary.R
#' @param y1 a matrix of copying vectors, where each row is an individual labelled as in the clusters in x
#' @param stat should rows be averaged or summed together? Defaults to mean average
#' @keywords Busby_bespoke
#' @export
#' @examples
#' rowsAsMapCluste(mapstatelist, fs_matrix)

rowsAsMapClusts <- function(x,y1,stat=mean){
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