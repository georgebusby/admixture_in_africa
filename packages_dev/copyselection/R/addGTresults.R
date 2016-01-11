#' Pull out data from a main GLOBETROTTER results table depending on what event-type is called
#'
#' This function extracts information on admixture events from a GLOBETROTTER results table and applies them some data frames that can be used for plotting.
#' @param pops a vector of populations in "pltable" (the default name for the GT results table)
#' rev_pops a vector of populations where we want to define the main source as the one that contribute least to the admixture event. This is useful when we want to line up populations to have the the source coming from the same place together on a plot, even if the source from one is 'major' and the other'minor'
#' @keywords Busby_bespoke
#' @return this is a rather unwieldy function that adds to the various results tables defined in the script prior to running this.
#' @export
#' @examples
#' addGTresults(pop)


addGTresults <- function(pltable, rev_pops){
    all_dates <- all_plot_mat <- all_src_mat <- dboots <- allsources <- c()
    for(pop in as.character(pltable$Cluster))
    {    
        print(pop)
        res <- pltable$Result[pltable$Cluster==pop]
        if(res == "1D")
        {
            a <- pltable$Analysis[pltable$Cluster==pop]
            prop <- as.numeric(pltable$alpha[pltable$Cluster==pop])
            gen <- pltable$date.1D[pltable$Cluster==pop]
            revg <- grep("B",gen)
            gen <- gsub("B","",gen)
            if(length(revg>0))
            {
                gen[revg] <- -as.numeric(gen[revg])
            } else
            {
                gen <- as.numeric(gen)
            }
            gens <- pltable$Date.CI[pltable$Cluster==pop]
            gens <- strsplit(gsub("\\)","",gsub("\\(","",gens)),split="\\-")[[1]]
            revg <- grep("B",gens)
            gens <- gsub("B","",gens)
            gens <- as.numeric(gens)
            if(length(revg>0)) gens[revg] <- -gens[revg]
            gen1 <- as.numeric(gens[1])
            gen2 <- as.numeric(gens[2])
            if(a == "main.null") a <- "null"
            src1 <- admixturesources2[paste(pop,a,1,sep="."),]
            src2 <- admixturesources2[paste(pop,a,2,sep="."),]
            bsrc1 <- pltable$best.source1[pltable$Cluster==pop]
            bsrc2 <- pltable$best.source2[pltable$Cluster==pop]
            
            if(pop %in% rev_pops)
            {
                prop <- 1-prop
                src1 <- admixturesources2[paste(pop,a,2,sep="."),]
                src2 <- admixturesources2[paste(pop,a,1,sep="."),]
                bsrc1 <- pltable$best.source2[pltable$Cluster==pop]
                bsrc2 <- pltable$best.source1[pltable$Cluster==pop]
                
            }
            pkey <- pop #paste(pop,a,sep=".")
            tds <- dateboots[dateboots$pop==pkey&dateboots[,2]==a,4]
            test <- length(tds)<100
            if(test) tds <- c(tds,rep(NA,(100-length(tds))))
            dboots <- rbind(dboots,c(pop,a,res,bsrc1,bsrc2,tds))
            allsources <- rbind(allsources,c(pop,res,1,prop,bsrc1,src1),c(pop,res,2,(1-prop),bsrc2,src2))
        }
        
        if(res == "1MW")
        {
            #res <- pltable$Result[pltable$Cluster==pop]
            a <- pltable$Analysis[pltable$Cluster==pop]
            prop <- as.numeric(pltable$alpha[pltable$Cluster==pop])
            prop2 <- as.numeric(pltable$alpha2[pltable$Cluster==pop])
            gen <- as.numeric(pltable$date.1D[pltable$Cluster==pop])
            revg <- grep("B",gen)
            gen <- gsub("B","",gen)
            if(length(revg>0))
            {
                gen[revg] <- -as.numeric(gen[revg])
            } else
            {
                gen <- as.numeric(gen)
            }
            gens <- pltable$Date.CI[pltable$Cluster==pop]
            gens <- strsplit(gsub("\\)","",gsub("\\(","",gens)),split="\\-")[[1]]
            revg <- grep("B",gens)
            gens <- gsub("B","",gens)
            gens <- as.numeric(gens)
            if(length(revg>0)) gens[revg] <- -gens[revg]
            gen1 <- as.numeric(gens[1])
            gen2 <- as.numeric(gens[2])
            if(a == "main.null") a <- "null"
            src1 <- admixturesources2[paste(pop,a,1,sep="."),]
            src2 <- admixturesources2[paste(pop,a,2,sep="."),]
            src3 <- admixturesources2[paste(pop,a,3,sep="."),]
            src4 <- admixturesources2[paste(pop,a,4,sep="."),]
            bsrc1 <- pltable$best.source1[pltable$Cluster==pop]
            bsrc2 <- pltable$best.source2[pltable$Cluster==pop]
            bsrc3 <- pltable$best.source1.ev2[pltable$Cluster==pop]
            bsrc4 <- pltable$best.source2.ev2[pltable$Cluster==pop]
            
            if(pop %in% rev_pops)
            {
                prop <- 1-prop
                src1 <- admixturesources2[paste(pop,a,2,sep="."),]
                src2 <- admixturesources2[paste(pop,a,1,sep="."),]
                src3 <- admixturesources2[paste(pop,a,4,sep="."),]
                src4 <- admixturesources2[paste(pop,a,3,sep="."),]
                bsrc1 <- pltable$best.source2[pltable$Cluster==pop]
                bsrc2 <- pltable$best.source1[pltable$Cluster==pop]
                bsrc3 <- pltable$best.source2.ev2[pltable$Cluster==pop]
                bsrc4 <- pltable$best.source1.ev2[pltable$Cluster==pop]
            }
            
            pkey <- pop # paste(pop,a,sep=".")
            tds <- dateboots[dateboots$pop==pkey&dateboots[,2]==a,4]
            test <- length(tds)<100
            if( test ) tds <- c(tds,rep(NA,(100-length(tds))))
            dboots <- rbind(dboots,c(pop,a,res,bsrc1,bsrc2,tds))
            dboots <- rbind(dboots,c(pop,a,res,bsrc3,bsrc4,tds))
            
            allsources <- rbind(allsources,
                                c(pop,res,1,prop,bsrc1,src1),
                                c(pop,res,2,(1-prop),bsrc2,src2),
                                c(pop,res,3,prop2,bsrc3,src3),
                                c(pop,res,4,(1-prop2),bsrc4,src4))
            
        }
        
        if(res == "2D")
        {
            a <- pltable$Analysis[pltable$Cluster==pop]
            prop <- as.numeric(pltable$alpha2.date1[pltable$Cluster==pop])
            prop2 <- as.numeric(pltable$alpha2.date2[pltable$Cluster==pop])
            gen <- as.numeric(pltable$date.2D.1[pltable$Cluster==pop])
            revg <- grep("B",gen)
            gen <- gsub("B","",gen)
            if(length(revg>0))
            {
                gen[revg] <- -as.numeric(gen[revg])
            } else
            {
                gen <- as.numeric(gen)
            }
            gens <- pltable$Date2a.CI[pltable$Cluster==pop]
            if(!is.null(gens))
            {
                gens <- strsplit(gsub("\\)","",gsub("\\(","",gens)),split="\\-")[[1]]
                revg <- grep("B",gens)
                gens <- gsub("B","",gens)
                gens <- as.numeric(gens)
                if(length(revg>0)) gens[revg] <- -gens[revg]
                gen1 <- as.numeric(gens[1])
                gen2 <- as.numeric(gens[2])
            } else
            {
                gen1 <- NA
                gen2 <- NA
            }
            if(a == "main.null") a <- "null"
            src1 <- admixturesources2[paste(pop,a,5,sep="."),]
            src2 <- admixturesources2[paste(pop,a,6,sep="."),]
            src3 <- admixturesources2[paste(pop,a,7,sep="."),]
            src4 <- admixturesources2[paste(pop,a,8,sep="."),]
            
            bsrc1 <- pltable$best.source1.date1[pltable$Cluster==pop]
            bsrc2 <- pltable$best.source2.date1[pltable$Cluster==pop]
            bsrc3 <- pltable$best.source1.date2[pltable$Cluster==pop]
            bsrc4 <- pltable$best.source2.date2[pltable$Cluster==pop]
            
            #       if(pop %in% c("SEBANTU","AMAXHOSA","SWBANTU"))
            #       {
            #         prop <- 1-prop
            #         src1 <- admixturesources2[paste(pop,a,6,sep="."),]
            #         src2 <- admixturesources2[paste(pop,a,5,sep="."),]
            #         src3 <- admixturesources2[paste(pop,a,8,sep="."),]
            #         src4 <- admixturesources2[paste(pop,a,7,sep="."),]
            #         bsrc1 <- pltable$best.source2.date1[pltable$Cluster==pop]
            #         bsrc2 <- pltable$best.source1.date1[pltable$Cluster==pop]
            #         bsrc3 <- pltable$best.source2.date2[pltable$Cluster==pop]
            #         bsrc4 <- pltable$best.source1.date2[pltable$Cluster==pop]
            #         
            #       }
            
            pkey <- pop #paste(pop,a,sep=".")
            tds <- date2boots[date2boots$pop==pkey,4]
            test <- length(tds)<100
            if( test ) tds <- c(tds,rep(NA,(100-length(tds))))
            dboots <- rbind(dboots,c(pop,a,res,bsrc1,bsrc2,tds))
            tds <- date2boots[date2boots$pop==pkey,5]
            test <- length(tds)<100
            if( test ) tds <- c(tds,rep(NA,(100-length(tds))))
            dboots <- rbind(dboots,c(pop,a,res,bsrc3,bsrc4,tds))
            
            allsources <- rbind(allsources,
                                c(pop,res,1,prop,bsrc1,src1),
                                c(pop,res,2,(1-prop),bsrc2,src2),
                                c(pop,res,3,prop2,bsrc3,src3),
                                c(pop,res,4,(1-prop2),bsrc4,src4))      
        }
        
        ## THESE ARE ALL THE FIRST EVENTS
        if(noprop == F)
        {
            plot_mat <- matrix(c(src1[levels(popplot)]*prop,0.1,
                                 src2[levels(popplot)]*(1-prop)))
        }
        if(noprop == T)
        {
            plot_mat <- matrix(c(src1[levels(popplot)]*0.5,0.1,
                                 src2[levels(popplot)]*0.5))
        }
        all_plot_mat <- cbind(all_plot_mat,plot_mat)
        all_dates <- rbind(all_dates,c(gen,gen1,gen2,prop,res))
        rownames(all_dates)[nrow(all_dates)] <- pop
        colnames(all_plot_mat)[ncol(all_plot_mat)] <- pop
        ## NOW GET THE SECOND EVENTS
        if(res %in% c("1MW","2D"))
        {
            if(res == "1MW")
            {
                prop <- as.numeric(pltable$alpha2[pltable$Cluster==pop])
                gen <- as.numeric(pltable$date.1D[pltable$Cluster==pop])
                revg <- grep("B",gen)
                gen <- gsub("B","",gen)
                if(length(revg>0))
                {
                    gen[revg] <- -as.numeric(gen[revg])
                } else
                {
                    gen <- as.numeric(gen)
                }
                gens <- pltable$Date.CI[pltable$Cluster==pop]
                if(!is.null(gens))
                {
                    gens <- strsplit(gsub("\\)","",gsub("\\(","",gens)),split="\\-")[[1]]
                    revg <- grep("B",gens)
                    gens <- gsub("B","",gens)
                    gens <- as.numeric(gens)
                    if(length(revg>0)) gens[revg] <- -gens[revg]
                    gen1 <- as.numeric(gens[1])
                    gen2 <- as.numeric(gens[2])
                } else
                {
                    gen1 <- NA
                    gen2 <- NA
                }
            }
            
            if(res == "2D")
            {
                prop <- as.numeric(pltable$alpha2.date2[pltable$Cluster==pop])
                gen <- pltable$date.2D.2[pltable$Cluster==pop]
                revg <- grep("B",gen)
                gen <- gsub("B","",gen)
                if(length(revg>0))
                {
                    gen[revg] <- -as.numeric(gen[revg])
                } else
                {
                    gen <- as.numeric(gen)
                }
                gens <- pltable$Date2b.CI[pltable$Cluster==pop]
                if(!is.null(gens))
                {
                    gens <- strsplit(gsub("\\)","",gsub("\\(","",gens)),split="\\-")[[1]]
                    revg <- grep("B",gens)
                    gens <- gsub("B","",gens)
                    gens <- as.numeric(gens)
                    if(length(revg>0)) gens[revg] <- -gens[revg]
                    gen1 <- as.numeric(gens[1])
                    gen2 <- as.numeric(gens[2])
                } else
                {
                    gen1 <- NA
                    gen2 <- NA
                }
                
                if(pop %in% rev_pops) prop <- 1-prop
            }
            
            
            if(noprop == F)
            {
                plot_mat <- matrix(c(src3[levels(popplot)]*prop,0.1,
                                     src4[levels(popplot)]*(1-prop)))
            }
            if(noprop == T)
            {
                plot_mat <- matrix(c(src3[levels(popplot)]*0.5,0.1,
                                     src4[levels(popplot)]*0.5))
            }
            all_plot_mat <- cbind(all_plot_mat,plot_mat)
            all_dates <- rbind(all_dates,c(gen,gen1,gen2,prop,res))
            rownames(all_dates)[nrow(all_dates)] <- paste(pop,"_a",sep="")
            colnames(all_plot_mat)[ncol(all_plot_mat)] <- paste(pop,"_a",sep="")
        }
    }
    return(list(all_dates,all_plot_mat,all_src_mat,dboots,allsources))
}

