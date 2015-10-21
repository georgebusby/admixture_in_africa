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
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"

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
malder_file <- "~/R/Copy/Rprojects/AfricaPOPGEN/manuscript/f3tables/AllPopsMalderFinalALL.txt"
malder_plot <- read.table(malder_file,header=T, as.is=T)
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
leginfo <- read.table(leginfo_file,header=T,comment.char="", as.is = T)
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)

### DEFINE COLOURS ##
pcolshex <- c("#0000CD", "#03B4CC", "#A65628", "#FF7F00", "#984EA3", "#4DAF4A", "#CCCC00")[c(1,2,4,5,3,6,7)]
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan",
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
    pop_vec <- malder$EthnicGroup
    
    malder <- malder_plot[malder_plot$Rec.Map==analy,]
    ##################################################################
    ## PLOT DATES AND MAP COMPARISONS
    pdf(paste0("figures/AllPopsMalderTimesMapComps",analy,".pdf"),height=9,width=9)
        layout(matrix(c(1,1,1,2,3,4),3,2),widths=c(6,3))
        par(mar=c(4,14,4,1))
        x_labs3 <- c(2000,0,-2000,-5000,-10000)
        x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE")
        x_max <- min(x_labs3)
        plot(0,0,xlim=rev(range(x_labs3)),
             ylim=c(n_pops,1),type="n",axes=F,xlab="",ylab="")
        text(5000,-1,labels=LETTERS[1],adj=0,las=1,cex=2,lwd=3,xpd=T)
        x_labs3 <- c(2000,0,-2000,-5000,-10000)
        x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE")
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
        for(i in 1:n_pops) axis(2,pos=4000,at=i,
                                labels=tidyNames(pop_vec[i],fula=T),
                                col.axis=y_ax_cols[i],las=2,tck=0,lwd=0,line=-0.5)
        
        for(i in seq(0.5,(n_pops+1),1)) abline(h=i,lty=3,lwd=0.5)
        lines(x=c(3200,3200),y=c(0,n_pops+1),lwd=1,xpd=T)
        
        text(y=-1,x=3600,labels="1ST EVENT",xpd=T,srt=35,adj=0)
        text(y=-1,x=2800,labels="2ND EVENT",xpd=T,srt=35,adj=0)
        
        ## NOW PLOT DATES AND ANCESTRIES ##
        
        ev1posa <- 3800
        ev1posb <- 3400
        ev2posa <- 3000
        ev2posb <- 2600
        
        for(i in (1:length(pop_vec)))
        {
            pop1 <- as.character(pop_vec)[i]
            tt <- malder[malder$EthnicGroup==pop1,]
            tt <- tt[order(tt$Date.Gens),]
            res_tmp <- nrow(tt)
            if(nrow(tt)>0)
            {
                for(j in 1:nrow(tt))
                {
                    date <- tt$Date.Gens[j]
                    dateci <- tt$Date.Gens.CI[j]
                    donpop1 <- strsplit(tt$Test.Pops[j],split="\\;")[[1]][1]
                    donanc1 <- getPopRegion(tidyNames(donpop1,fula=T),popkey)
                    donpop2 <- strsplit(tt$Test.Pops[j],split="\\;")[[1]][2]
                    donanc2 <- getPopRegion(tidyNames(donpop2,fula=T),popkey)
                    mainanc <- tt$Main.Anc[j]
                    mainancp <- tt$Main.Anc.p[j]
                
                    ## PLOT DATE
                    arrows(makeDate(date-dateci,add_BCE=F),i,
                           makeDate(date+dateci,add_BCE=F),i,
                           angle=90,length=0.025,code=3)
                    pcol1 <- pcolshex[ancreg_list==donanc1]
                    pcol2 <- pcolshex[ancreg_list==donanc2]
                    pnt_bg <- "black"
                    if(mainancp>0.001) pnt_bg <- "grey"
                    if(mainancp>0.05) pnt_bg <- "white"
                    points(makeDate(date,add_BCE=F),i,pch=21,bg=pnt_bg, cex=2)
                
                    ## PLOT EVENT ANCESTRIES
                    if(j == 1)
                    {
                        points(ev1posa,i,pch=15,col=pcol1,xpd=T,cex=2)
                        points(ev1posb,i,pch=15,col=pcol2,xpd=T,cex=2)
                        if(donanc1 == mainanc) points(ev1posa,i,pch=22,xpd=T,cex=2,lwd=1.5)
                        if(donanc2 == mainanc) points(ev1posb,i,pch=22,xpd=T,cex=2,lwd=1.5)
                    }
                    
                    if(j == 2)
                    {    points(ev2posa,i,pch=15,col=pcol1,xpd=T,cex=2)
                         points(ev2posb,i,pch=15,col=pcol2,xpd=T,cex=2)
                         if(donanc1 == mainanc) points(ev2posa,i,pch=22,xpd=T,cex=2,lwd=1.5)
                         if(donanc2 == mainanc) points(ev2posb,i,pch=22,xpd=T,cex=2,lwd=1.5)
                    }     
                }
            }
        }
    
        ## PLOT LEGEND 
        par(mar=c(1,1.5,3,0))
        plot(0,0,axes=F,xlab="",ylab="",type="n")
        legend_text <- c("West African Niger-Congo",
                         "Central West African Niger-Congo",
                         "East African Niger-Congo",
                         "South African Niger-Congo",
                         "Nilo-Saharan / Afro-Asiatic",
                         "KhoeSan",
                         "Eurasia",
                         "main event ancestry",
                         "high confidence date (P<0.001)",
                         "low confidence date (P<0.05)",
                         "non-significant date (P>0.05)")
        l <- legend("top",legend=legend_text,pch=c(rep(15,7),22,21,21,21),
                    col=c(pcolshex[c(1:3,5,4,6,7)],"black","black","black","black"),bty="n",
                    pt.bg=c(pcolshex[c(1:3,5,4,6,7)],"white","black","grey","white"),
                    ncol=1,xpd=T,pt.cex=3,x.intersp=1,y.intersp=1.25,
                    pt.lwd=2,cex=1.25, title="Ancestry Region")
        par(mar=c(5,5,2,1))
    
        ## PLOT DATE COMPARISONS
        ## USE ONLY THOSE POPULATIONS WHERE THE SAME EVENTS ARE INFERRED IN ALL THREE MAP ANALYSES ##
        matching_events <- c()
        for(i in unique(malder_plot$EthnicGroup))
        {
            num_runs <- nrow(malder_plot[malder_plot$EthnicGroup==i,])
            if(num_runs %in% c(3,6))
            {
                matching_events <- c(matching_events,i)
            }
        }
        
        dates2plot <- malder_plot[malder_plot$EthnicGroup%in%matching_events,]
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
                 xlim=sapply(c(0,400),function(x)makeDate(x,add_BCE=F)),
                 ylim=sapply(c(400,0),function(x)makeDate(x,add_BCE=F)),
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
            plot_pops <- dates2plot$EthnicGroup[plotx] 
            
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
            pop_pnts <- getPopSymbols(tidyNames(plot_pops,fula=F),leginfo)
            points(sapply(dates2plot$Date.Gens[plotx],function(x)makeDate(x,add_BCE=F)),
                   sapply(dates2plot$Date.Gens[ploty],function(x)makeDate(x,add_BCE=F)),
                   pch=as.numeric(pop_pnts$pch2plot),
                   col=pop_pnts$rim2plot,
                   bg=pop_pnts$col2plot)
                   
        }
    dev.off()                
}
