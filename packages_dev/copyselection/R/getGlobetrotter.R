#' extract the relevant information from globetrotter's output files.
#'
#' This function extracts information on admixture events from a GLOBETROTTER output file
#' @param infile the file name of the GLOBETROTTER output file
#' @param donors a matrix of donor population copying vectors which are used to generate copying vectors for the admixture sources from the source composition betas. Default to NULL, in which case no source copying vectors are generated
#' @keywords Busby_bespoke
#' @return a two element list: the first element is the full results, and the second are the admixture sources as a matrix (should 'donors' be defined above)
#' @export
#' @examples
#' getGlobetrotter("FULAInolocal.props.out",donors=NULL)


getGlobetrotter <- function(infile,donors=NULL)
{
    
    tmp <- scan(infile,what="char",flush=T,fill=T,sep=";",quiet=T)
    res <- gsub("\\)","",strsplit(tmp[1],split="\\ ")[[1]][8])
    date1 <- getLines(4,5,tmp)
    date2 <- getLines(8,9,tmp)
    ## SOURCES
    date1pc1a <- getLines(12,13,tmp)
    date1pc1b <- getLines(14,15,tmp)
    date1pc2a <- getLines(18,19,tmp)
    date1pc2b <- getLines(20,21,tmp)
    date2pc1a <- getLines(24,25,tmp)
    date2pc1b <- getLines(26,27,tmp)
    date2pc2a <- getLines(30,31,tmp)
    date2pc2b <- getLines(32,33,tmp)
    
    ## STORE RESULTS
    fullres <- cbind(matrix(c("result",res),2,1),date1,date2)
    colnames(fullres) <- fullres[1,]
    fullres <- fullres[2,]
    if(res == "Unknown")
    {
        fullres["proportion.source1"] <- date1pc1a[2,1]
        fullres["bestmatch.event1.source1"] <- date1pc1a[1,1]
    }
    
    ## GENERATE A MATRIX WITH ADMIXTURE BETAS AND PROPORTION
    if(!is.null(donors) & res != "Unknown")
    {
        sourcemat <- matrix(0,nrow=8,ncol=length(donors)+1)
        colnames(sourcemat) <- c("proportion",donors)
        sourcemat[1,date1pc1a[1,]] <- date1pc1a[2,]
        sourcemat[2,date1pc1b[1,]] <- date1pc1b[2,]
        sourcemat[3,date1pc2a[1,]] <- date1pc2a[2,]
        sourcemat[4,date1pc2b[1,]] <- date1pc2b[2,]
        sourcemat[5,date2pc1a[1,]] <- date2pc1a[2,]
        sourcemat[6,date2pc1b[1,]] <- date2pc1b[2,]
        sourcemat[7,date2pc2a[1,]] <- date2pc2a[2,]
        sourcemat[8,date2pc2b[1,]] <- date2pc2b[2,]
        
        rownames(sourcemat) <- c("1date.pc1a","1date.pc1b","1date.pc2a","1date.pc2b",
                                 "2date.pc1a","2date.pc1b","2date.pc2a","2date.pc2b")
        sourcemat <- apply(sourcemat,2,as.double)
        return(list(fullres,sourcemat))
    }  else 
    {
        return(list(fullres,NULL))
    }
}
