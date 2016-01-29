### SCRIPT TO PLOT MALDER RESULTS FOR PAPER -- CAN USE RESULTS FROM DIFFERENT GENETIC MAPS ###
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
library(xtable)
setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
## LOAD POPORDER FILE ##
popplot <- scan("/mnt/kwiat/home/popgen/scripts/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplotorder <- popplot

## CHOOSE ONLY THE AFRICAN POPS FROM POPPLOT ORDER
popplotorder <- popplotorder[c(1:16,18:49)]
poplist <- factor(popplotorder,levels=popplotorder)
poplist<- as.matrix(poplist)
n_pops <- nrow(poplist)
in_dir <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/alder/output/"
analys <- c("HAPMAP","CeuMap","AfrMap")
malder_file <- "~/R/Copy/Rprojects/AfricaPOPGEN/manuscript/f3tables/AllPopsMalderFinalALLmindis.txt"
malder_plot <- read.table(malder_file,header=T, as.is=T)
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
leginfo <- read.table(leginfo_file,header=T,comment.char="", as.is = T)
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)

### DEFINE COLOURS ##
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )


###################################################################
## DO THE DIFFERENT MAPS INFER THE SAME NUMBER OF EVENTS? ##
for(i in unique(malder_plot$EthnicGroup))
{
    tmp <- malder_plot[malder_plot$EthnicGroup==i,]
    if(!nrow(tmp) %in% c(3,6)) print(tmp)
}
###################################################################

pop_vec <- unique(malder_plot$EthnicGroup)
n_pops <- length(pop_vec)

### FOR A GIVEN ANALYSIS, SUBSET THE MAIN DATABASE AND PLOT

for(analy in c("HAPMAP","AfrMap","CeuMap"))
{
    ## ORDER pop_vec BY THE DATES AND ANCESTRIES IN THE HAPMAP ##
    order_analy <- analy
    malder <- malder_plot[malder_plot$Rec.Map==order_analy,]
    dup_pops <- malder$EthnicGroup[duplicated(malder$EthnicGroup)]
    ## ORDER ON THE FIRST EVENT SO MAKE SURE THIS IS THE FIRST ONE IN THE TABLE
    for(i in dup_pops)
    {
        ii <- malder[malder$EthnicGroup==i,]
        ii <- ii[order(ii$Date.Gens),]
        malder[malder$EthnicGroup==i,] <- ii
    }
    malder <- malder[!duplicated(malder$EthnicGroup),]
    malder$pop1 <- sapply(malder$Test.Pops,function(x)strsplit(x,split="\\;")[[1]][1])
    malder$pop2 <- sapply(malder$Test.Pops,function(x)strsplit(x,split="\\;")[[1]][2])
    malder$reg1 <- factor(getPopRegion(tidyNames(malder$pop1,fula=T),popkey),levels=rev(ancreg_list))
    malder$reg2 <- factor(getPopRegion(tidyNames(malder$pop2,fula=T),popkey),levels=rev(ancreg_list))
    malder <- malder[order(malder$reg1,malder$reg2,malder$Date.Gens,decreasing=T),]
    pop_vec <- popplot[c(1:16,18:49)]#malder$EthnicGroup
    
    malder <- malder_plot[malder_plot$Rec.Map==analy,]
    malder[malder$N.evs==0,2:ncol(malder)] <- 0
    malder[malder$Main.Anc.p>0.05,2:ncol(malder)] <- 0
    
    ##################################################################
    ## PLOT DATES AND MAP COMPARISONS
    pdf(paste0("figures/AllPopsMalderTimesMapCompsMinDis",analy,".pdf"),height=9,width=9)
        layout(matrix(c(1,1,1,2,2,2,3,3,3,4,5,6),3,4),widths=c(4,1,1,3))
        par(mar=c(4,12,4,0.5))
        
        comp_ax <- c(0,400)
        x_labs3 <- c(2000,0,-2000,-5000,-10000)    
        x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE")
        #comp_ax <- c(0,240)
        #x_labs3 <- c(2000,0,-2000,-5000)
        #x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE")
    
        x_max <- min(x_labs3)
        plot(0,0,xlim=rev(range(x_labs3)),
             ylim=c(n_pops,1),type="n",axes=F,xlab="",ylab="")
        text(poplabpos+2000,-1,labels=LETTERS[1],adj=0,las=1,cex=2,lwd=3,xpd=T)
        axis(1,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5)
        for(j in x_labs3) abline(v=j,lty=2)
        mtext("Date of Admixture",1,line=3)
        
        ## COLOURS FOR EACH GROUP
        y_ax_cols <- c()
        for(i in pop_vec) 
        {
            if(i=="SEMI-BANTU") i <- "SEMI-BANTU"
            ii <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
            y_ax_cols <- c(y_ax_cols,ii)
        }
        
        ## SET THE SCENE
        poplabpos <- 4800
        ev1pos1 <- 4550
        ev1pos2 <- 3900
        ev2pos1 <- 3250
        ev2pos2 <- 2600
#         poplabpos <- 4600
#         ev1pos1 <- 4100
#         ev1pos2 <- 3600
#         ev2pos1 <- 3100
#         ev2pos2 <- 2600
    
        
        for(i in 1:n_pops) axis(2,pos=poplabpos,at=i,
                                labels=tidyNames(pop_vec[i],fula=T,khoesan=T,tig=T),
                                col.axis=y_ax_cols[i],las=2,tck=0,lwd=0,line=-0.5)
        
        for(i in seq(0.5,(n_pops+1),1)) abline(h=i,lty=3,lwd=0.5)
        lines(x=rep(median(c(ev2pos1,ev1pos2)),2),y=c(0,n_pops+1),lwd=1,xpd=T)
        
        text(y=-1,x=median(c(ev1pos1,ev1pos2)),labels="1ST EVENT",xpd=T,srt=35,adj=0)
        text(y=-1,x=median(c(ev2pos1,ev2pos2)),labels="2ND EVENT",xpd=T,srt=35,adj=0)
        
        ## NOW PLOT DATES AND ANCESTRIES ##
        
    
        onedate_mat <- twodate_mat <- matrix(NA,nr=length(pop_vec),nc=(length(ancreg_list)*2))
        tpopcols <- c(grep("tpop1",colnames(malder),value=T)[2:9],
                      grep("tpop2",colnames(malder),value=T)[2:9])
        colnames(onedate_mat) <- colnames(twodate_mat) <- tpopcols
        for(i in (1:length(pop_vec)))
        {
            pop1 <- as.character(pop_vec)[i]
            tt <- malder[malder$EthnicGroup==pop1,]
            tt <- tt[order(tt$Date.Gens),]
            res_tmp <- nrow(tt)
            if(nrow(tt)>0)
            {
                if(tt$Date.Gens>0)
                {
                    for(j in 1:nrow(tt))
                    {
                        date <- tt$Date.Gens[j]
                        dateci <- tt$Date.Gens.CI[j]
                        donpop1 <- strsplit(tt$Test.Pops[j],split="\\;")[[1]][1]
                        donanc1 <- tt$tpop1[j]
                        donpop2 <- strsplit(tt$Test.Pops[j],split="\\;")[[1]][2]
                        donanc2 <- tt$tpop2[j]
                        mainanc <- tt$Main.Anc[j]
                        mainancp <- tt$Main.Anc.p[j]
                    
                        ## PLOT DATE
                        dh <- makeDate(date-dateci,add_BCE=F)
                        dl <- makeDate(date+dateci,add_BCE=F)
                        if(!is.na(dl) & dl > x_labs3[length(x_labs3)])
                        {
                            arrows(dh,i,dl,i,angle=90,length=0.025,code=3)
                        } else
                        {
                            arrows(dh,i,dl,i,angle=90,length=0.025,code=3)
                            dl <- x_labs3[length(x_labs3)]
                            xrange <- range(x_labs3)[2]-range(x_labs3)[1]
                            dl <- dl - 0.04*xrange
                            arrows(dh,i,dl,i,angle=45,length=0.025,code=2)
                        }
                        
                        
                        pcol1 <- pcolshex[ancreg_list==donanc1]
                        pcol2 <- pcolshex[ancreg_list==donanc2]
                        pnt_bg <- "black"
                        #if(mainancp>0.001) pnt_bg <- "grey"
                        #if(mainancp>0.05) pnt_bg <- "white"
                        points(makeDate(date,add_BCE=F),i,pch=21,bg=pnt_bg, cex=2)
                    
                        ## PLOT EVENT ANCESTRIES
                        if(j == 1)
                        {
                            points(ev1pos1,i,pch=15,col=pcol1,xpd=T,cex=2)
                            points(ev1pos2,i,pch=15,col=pcol2,xpd=T,cex=2)
                            if(donanc1 == mainanc) points(ev1pos1,i,pch=22,xpd=T,cex=2,lwd=1.5)
                            if(donanc2 == mainanc) points(ev1pos2,i,pch=22,xpd=T,cex=2,lwd=1.5)
                            ## ADD VALUES TO ONE-DATE MATRIX
                            onedate_mat[i,] <- as.numeric(tt[j,tpopcols])
                        }
                        
                        if(j == 2)
                        {    
                            points(ev2pos1,i,pch=15,col=pcol1,xpd=T,cex=2)
                            points(ev2pos2,i,pch=15,col=pcol2,xpd=T,cex=2)
                            if(donanc1 == mainanc) points(ev2pos1,i,pch=22,xpd=T,cex=2,lwd=1.5)
                            if(donanc2 == mainanc) points(ev2pos2,i,pch=22,xpd=T,cex=2,lwd=1.5)
                            ## ADD VALUES TO ONE-DATE MATRIX
                            twodate_mat[i,] <- as.numeric(tt[j,tpopcols]) 
                        }     
                    }
                }
            }
        }

        ## PLOT SOME
        for(mat in c(1,2))
        {
            if(mat == 1) tmpmat <- onedate_mat
            if(mat == 2) tmpmat <- twodate_mat
            tmpmat[tmpmat>2] <- NA
            tmpmat[is.na(tmpmat)]  <- "white"
        
            for(j in 1:ncol(tmpmat))
            {
                tmpmat[tmpmat[,j]!="white",j] <- rep(pcolshex,2)[j]
            }
        
            ## EMPTY PLOT
            par(mar=c(4,0.1,4,0))
            plot(0,0,xlim=c(0,ncol(tmpmat)),
                ylim=c(n_pops,1),type="n",axes=F,xlab="",ylab="")
            xleft <- rep(1:ncol(tmpmat),each=nrow(tmpmat))    
            ybottom <- rep(1:nrow(tmpmat), times=ncol(tmpmat))
            rect(xleft-1, ybottom-0.5, xleft, ybottom+0.5, col=tmpmat,border=NA)
            abline(v=8,lwd=1,lty=2)
            if(mat == 2) abline(v=-0.5,xpd=F)
            text(y=0,x=0*(ncol(tmpmat)/4),labels="SOURCE 1",xpd=T,srt=35,adj=0)
            text(y=0,x=2*(ncol(tmpmat)/4),labels="SOURCE 2",xpd=T,srt=35,adj=0)
            if(mat == 1) text(y=-3,x=2*(ncol(tmpmat)/4),labels="1ST EVENT",xpd=T,srt=0,adj=0.5)
            if(mat == 2) text(y=-3,x=2*(ncol(tmpmat)/4),labels="2ND EVENT",xpd=T,srt=0,adj=0.5)
        }
    
        ## PLOT LEGEND 
        par(mar=c(1,1.5,2,0))
        plot(0,0,axes=F,xlab="",ylab="",type="n")
        legend_text <- c("West African Niger-Congo",
                         "Central West African Niger-Congo",
                         "East African Niger-Congo",
                         "South African Niger-Congo",
                         "East African Nilo-Saharan",
                         "East African Afroasiatic",
                         "KhoeSan",
                         "Eurasia",
                         "main event ancestry")
        l <- legend("top",legend=legend_text,pch=c(rep(15,8),22),
                    col=c(pcolshex[c(1:3,6,4,5,7,8)],"black"),bty="n",
                    pt.bg=c(pcolshex[c(1:3,6,4,5,7,8)],"white"),
                    ncol=1,xpd=T,pt.cex=3,x.intersp=1,y.intersp=1.25,
                    pt.lwd=2,cex=1.25, title="Ancestry Region")
        par(mar=c(5,5,2,1))
    
        ## PLOT DATE COMPARISONS
        ## USE ONLY THOSE POPULATIONS WHERE THE SAME EVENTS ARE INFERRED IN ALL THREE MAP ANALYSES ##
        test <- malder_plot$Date.Gens>0 & malder_plot$N.evs>0
        dates2plot <- malder_plot[test,]
        matching_events <- c()
        for(i in unique(dates2plot$EthnicGroup))
        {
            num_runs <- nrow(dates2plot[dates2plot$EthnicGroup==i,])
            num_ev_run <- unique(dates2plot$N.evs[dates2plot$EthnicGroup==i])
            if(num_runs %in% c(3,6) & length(num_ev_run) == 1)
            {
                matching_events <- c(matching_events,i)
            }
        }    
    
        dates2plot <- dates2plot[dates2plot$EthnicGroup%in%matching_events,]    
        tmp <- dates2plot[dates2plot$N.evs<2,]
        dates2plot<- dates2plot[dates2plot$N.evs==2,]
        for(i in unique(dates2plot$EthnicGroup))
        {
            for(j in unique(dates2plot$Rec.Map))
            {
                test <- dates2plot$EthnicGroup==i & dates2plot$Rec.Map == j
                dates2plot[test,] <- dates2plot[test,][order(dates2plot$Date.Gens[test]),]    
            }
        }
        dates2plot <- rbind(tmp,dates2plot)        
        dates2plot <- dates2plot[order(dates2plot$EthnicGroup),]
    
        xplot <- "AfrMap"
        xplotlab <- "YRI"
        yplot <- "HAPMAP"
    
        for(j in 1:2)
        {
            if(j == 1)
            {
                xplot <- "CeuMap"
                xplotlab <- "CEU"
            }
            if(j == 2)
            {
                xplot <- "AfrMap"
                xplotlab <- "YRI"
            }
            plot(0,0,
                 xlab=paste("Date with",xplotlab, "genetic map (Hinch et al 2011)"),
                 ylab="Date with HAPMAP worldwide genetic map",
                 xlim=sapply(comp_ax,function(x)makeDate(x,add_BCE=F)),
                 ylim=sapply(rev(comp_ax),function(x)makeDate(x,add_BCE=F)),
                 yaxt="n",type="n",xaxt="n")
            xaxis <- c(2000,0,-2000,-5000,-10000)
            xaxislabs <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE")
            axis(1,at=xaxis,labels=xaxislabs,tick=F,padj=0.5,line=0)
            axis(2,at=xaxis,labels=xaxislabs,las=2,tick=F,hadj=0.5,line=0.25)
            for(k in xaxis) abline(h=k,lty=2)
            for(k in xaxis) abline(v=k,lty=2)
            abline(a=0,b=1,lwd=2,col="red")
            mtext(3,line=0,text=LETTERS[j+1],adj=0,las=1,cex=1.5,lwd=3,xpd=T)
        
            plotx <- dates2plot$Rec.Map == xplot
            ploty <- dates2plot$Rec.Map == yplot
            plot_popsx <- dates2plot$EthnicGroup[plotx] 
            plot_popsy <- dates2plot$EthnicGroup[ploty]
            plot_pops <- intersect(plot_popsx,plot_popsy)
            plotx <- dates2plot$Rec.Map == xplot & dates2plot$EthnicGroup %in% plot_pops
            ploty <- dates2plot$Rec.Map == yplot & dates2plot$EthnicGroup %in% plot_pops
            
            ## PLOT CIs
            for(i in (1:sum(plotx)))
            {
                xs <- c(dates2plot$Date.Gens[plotx][i]+dates2plot$Date.Gens.CI[plotx][i],
                        dates2plot$Date.Gens[plotx][i]-dates2plot$Date.Gens.CI[plotx][i])
                ys <- c(dates2plot$Date.Gens[ploty][i],dates2plot$Date.Gens[ploty][i])
                lines(sapply(xs,function(x)makeDate(x,add_BCE=F)),
                      sapply(ys,function(x)makeDate(x,add_BCE=F)))
                xs <- c(dates2plot$Date.Gens[plotx][i],dates2plot$Date.Gens[plotx][i])
                ys <- c(dates2plot$Date.Gens[ploty][i]+dates2plot$Date.Gens.CI[ploty][i],
                        dates2plot$Date.Gens[ploty][i]-dates2plot$Date.Gens.CI[ploty][i])
                lines(sapply(xs,function(x)makeDate(x,add_BCE=F)),
                      sapply(ys,function(x)makeDate(x,add_BCE=F)))
            }
            ## GET SYMBOLS, COLS, ETC FOR EACH POP AND PLOT
            pop_pnts <- getPopSymbols(tidyNames(dates2plot$EthnicGroup[plotx],fula=F),leginfo)
            points(sapply(dates2plot$Date.Gens[plotx],function(x){makeDate(x,add_BCE=F)}),
                   sapply(dates2plot$Date.Gens[ploty],function(x){makeDate(x,add_BCE=F)}),
                   pch=as.numeric(pop_pnts$pch2plot),
                   col=pop_pnts$rim2plot,
                   bg=pop_pnts$col2plot)
                   
        }
    dev.off()                
}

