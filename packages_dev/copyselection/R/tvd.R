#' this is an internal function that generates a table of TVD estimates bases on a matrix of copying vectors
#'
#' This function generates a table of TVD estimates bases on a matrix of copying vectors.
#'
#' @param x a matrix of copying vectors
#' @keywords Busby_bespoke
#' @return returns a table of pairwise TVD for each pair of rows of matrix 'x'
#' @export
#' @examples
#' tvd(x)

tvd <- function(x){
    tvd <- c()
    pop1 <- c()
    pop2 <- c()
    for(i in 1:nrow(x))
    {
        for(j in 2:nrow(x))
        {
            if(j > i)
            {
                tvd <- c(tvd,0.5*sum(abs(x[i,]-x[j,])))   
                pop1 <- c(pop1,rownames(x)[i])
                pop2 <- c(pop2,rownames(x)[j])
            }  
        }
    }
    d <- data.frame(cbind(pop1,pop2,tvd))
    colnames(d) <- c("pop1","pop2","tvd")
    d
}
