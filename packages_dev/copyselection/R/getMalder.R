#' this is an internal function that generates a date from a time in generations
#'
#' This function takes MALDER output and tabulates
#' @param file_name name of file to pull MALDER results from.
#' @param pop1 name of population.
#' @keywords Busby_bespoke
#' @export
#' @examples
#' getMalder(...)

getMalder <- function(file_name, pop1)
{
    all_res2 <- c()
    
    ######################################
    ## FIND NUMBER OF ADMIXTURE EVENTS
    res <- system(paste0("grep ^RESULT ",file_name," | awk '{print $1}' | sort | uniq "),intern=T)
    if(length(res)>0)
    {
        if(length(res)>1)
        {
            res <- res[length(res)-1]
            num_events <- as.numeric(strsplit(res,split="\\_")[[1]][2])
        } else
        {
            res <- ""
            num_events <- 0
        }
        ## GET DATA
        #if(num_events == 1 ) one_event <- c(one_event,pop1)
        #if(num_events > 2) multi_events <- c(multi_events,pop1)
        print(c(pop1,num_events))
        if(num_events > 2) num_events <- 2
        if(res == "RESULT_3") res <- "RESULT_2"
        if(num_events <= 1)
        {
            res <- "RESULT_1"
            tmp_tab <- system(paste0("grep ",res," ",file_name),intern=T)
            res_tab <- c()
            for(i in 1:length(tmp_tab))
            {
                i <- strsplit(tmp_tab[i],split="\t")[[1]]
                res_tab <- rbind(res_tab,i)
            }
            colnames(res_tab) <- res_tab[1,]
            res_tab <- res_tab[-1,]
            res_tab2 <- c()
            for(i in 1:nrow(res_tab))
            {
                a1 <- splitMalder(res_tab[i,"amp0"])
                t1 <- splitMalder(res_tab[i,"time0"])
                res_tab2 <- rbind(res_tab2,c(a1,t1,rep(0,6)))
            }
            colnames(res_tab2) <- c("amp0","amp0.ci","amp0.Z",
                                    "time0","time0.ci","time0.Z",
                                    "amp1","amp1.ci","amp1.Z",
                                    "time1","time1.ci","time1.Z")
            
            res_tab <- cbind(pop1,num_events,res_tab[,"refpops"],res_tab2)
            all_res2 <- rbind(all_res2,res_tab)
        }
        if(num_events == 2)
        {
            tmp_tab <- system(paste0("grep ",res," ",file_name),intern=T)
            res_tab <- c()
            for(i in 1:length(tmp_tab))
            {
                i <- strsplit(tmp_tab[i],split="\t")[[1]]
                res_tab <- rbind(res_tab,i)
            }
            colnames(res_tab) <- res_tab[1,]
            res_tab <- res_tab[-1,]
            res_tab2 <- c()
            for(i in 1:nrow(res_tab))
            {
                a1 <- splitMalder(res_tab[i,"amp0"])
                t1 <- splitMalder(res_tab[i,"time0"])
                a2 <- splitMalder(res_tab[i,"amp1"])
                t2 <- splitMalder(res_tab[i,"time1"])
                res_tab2 <- rbind(res_tab2,c(a1,t1,a2,t2))
            }
            colnames(res_tab2) <- c("amp0","amp0.ci","amp0.Z",
                                    "time0","time0.ci","time0.Z",
                                    "amp1","amp1.ci","amp1.Z",
                                    "time1","time1.ci","time1.Z")
            
            res_tab <- cbind(pop1,num_events,res_tab[,"refpops"],res_tab2)
            all_res2 <- rbind(all_res2,res_tab)
        }
        if(num_events == 3)
        {
            tmp_tab <- system(paste0("grep ",res," ",file_name),intern=T)
            res_tab <- c()
            for(i in 1:length(tmp_tab))
            {
                i <- strsplit(tmp_tab[i],split="\t")[[1]]
                res_tab <- rbind(res_tab,i)
            }
            colnames(res_tab) <- res_tab[1,]
            res_tab <- res_tab[-1,]
            res_tab2 <- c()
            for(i in 1:nrow(res_tab))
            {
                a1 <- splitMalder(res_tab[i,"amp0"])
                t1 <- splitMalder(res_tab[i,"time0"])
                a2 <- splitMalder(res_tab[i,"amp1"])
                t2 <- splitMalder(res_tab[i,"time1"])
                a3 <- splitMalder(res_tab[i,"amp2"])
                t3 <- splitMalder(res_tab[i,"time1"])
                res_tab2 <- rbind(res_tab2,c(a1,t1,a2,t2,a3,t3))
            }
            colnames(res_tab2) <- c("amp0","amp0.ci","amp0.Z",
                                    "time0","time0.ci","time0.Z",
                                    "amp1","amp1.ci","amp1.Z",
                                    "time1","time1.ci","time1.Z",
                                    "amp2","amp2.ci","amp2.Z",
                                    "time2","time2.ci","time2.Z")
            
            res_tab <- cbind(pop1,num_events,res_tab[,"refpops"],res_tab2)
            all_res2 <- rbind(all_res2,res_tab)
        }
    }
    colnames(all_res2)[2:3] <- c("result","test.pops")
    return(all_res2)
}
