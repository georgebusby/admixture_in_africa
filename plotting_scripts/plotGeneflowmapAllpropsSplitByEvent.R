### SCRIPT TO PLOT GLOBETROTTER RESULTS ON MAP ###
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
library(xtable)
library("rworldmap")
library("maptools")
#source("functions/convertfactorstostrings.R")
#source("functions/makeTransparent.R")

setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
## GET DATES AS WELL
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"

leginfo <-read.table(leginfo_file,header=T,comment.char="")
newi <- plotPopPoints(latlong_file,leginfo_file,poppos_file,pt_cex=2)
pointsonly <- unique(cbind(newi$Lat,newi$Long))

## SET PLOTTING PARAMETERS
plotSep <- FALSE ## plot individually
plotAdmixture <- FALSE ## plot an admixture plot?
plotEquator <- FALSE ## should equator and tropics be plotted
map.xlim.2plot<-c(-20,50) ## X-AXIS LIMITS OF MAP
map.ylim.2plot<-c(-32,28) ## Y-AXIS LIMITS OF MAP
## POSITION OF EACH OF THE LEGENDS
legposx <- c(map.xlim.2plot[1]-5,map.xlim.2plot[2]-5,map.xlim.2plot[1]+8)
legposy <- c(5,15,-5)
pt_cex <- 1.5



eventlist <- c("Western Bantu","Eastern Bantu","West Africa","Eurasia")

## make some bespoke xs and ys
europe <- c("CEU","GBR","FIN","IBS","TSI")
asia <- c("CDX","KHV","JPT","GIH","CHS","CHB")
americas <-  c("PELII","PEL")

euroxy <- cbind(c(europe,asia,americas),
                c(seq(map.ylim.2plot[1]+10,map.ylim.2plot[2]-10,length=length(europe)-2),
                  c(map.ylim.2plot[2],map.ylim.2plot[2]),
                  seq(map.ylim.2plot[1]+10,map.ylim.2plot[2]-10,length=length(asia)),c(0,0)),
                c(rep(map.xlim.2plot[1]-0,length(europe)-2),
                  c(median(map.xlim.2plot)-5,median(map.xlim.2plot)+5),
                  rep(map.xlim.2plot[2]+0,length(asia)),
                  c(0,0)),
                c(rep("#eeee13",3),rep("#CCCC00",2),
                  "#8B8878","#787833","#8B8878","#787833",rep("#8B8878",4)))
                
euroxy <- data.frame(euroxy,stringsAsFactors = FALSE)
colnames(euroxy) <- c("pop","y","x","cols")

propslessthan <- 0.025

######################################################
## PLOT FIGURE 6 FROM MAIN PAPER

pdf("figures/AfricadmixtureEventsMapAllPropsOver025SplitByEvent.pdf",height=10,width=10)
layout(matrix(1:4,2,2,byrow=T))
for(k in 1:length(eventlist))
{
    date2plot <- c(-10000,2000)
    ## plot map
    par(mar=c(1,1,1,1))
    mapCountryData(newmap2,nameColumnToPlot="colour",catMethod="categorical",
                   oceanCol="paleturquoise1",missingCountryCol="white",addLegend=F,
                   mapTitle="", colourPalette=rep("palegoldenrod",length(levels(mapinfo$colour))),
                   xlim=map.xlim.2plot,
                   ylim=map.ylim.2plot)
    box()
    tmpdates <- c()
    ## plot populatoin points
    for(i in 1:nrow(pointsonly)) points(pointsonly[i,2],pointsonly[i,1],pch=20,col="black",cex=1)
    linewds <- c(1,2.5,4,6,6,8)
    admixturesources2 <- read.table("data/MalariaGenGlobetrotterAdmixtureSources3.txt",header=T,row.names=1,as.is=T)
    admixturesources2 <- admixturesources2[,levels(popplot)]
    final.res2plot <- read.table("data/MalariaGenGlobetrotter2plot.txt",header=T,row.names=1,as.is=T)
    dateboots <- read.table("data/MalariaGenGlobetrotterOneDateBootstraps.txt",header=T,row.names=1,as.is=T)
    date2boots <- read.table("data/MalariaGenGlobetrotterTwoDateBootstraps.txt",header=T,row.names=1,as.is=T)
    best_ald <- read.table("data/MalariaGenBestAlder.txt",header=T,as.is=T,comment.char="",fill=T)

    for(i in 1:nrow(final.res2plot))
    {
        recip <- final.res2plot$Cluster[i]
        donor1 <- c(final.res2plot$best.source1[i],
                    final.res2plot$best.source1.ev2[i],
                    final.res2plot$best.source1.date1[i],
                    final.res2plot$best.source1.date2[i])
        donor2 <- c(final.res2plot$best.source2[i],
                    final.res2plot$best.source2.ev2[i],
                    final.res2plot$best.source2.date1[i],
                    final.res2plot$best.source2.date2[i])
        dates <- c(final.res2plot$date.1D[i],
                   final.res2plot$date.2D.1[i],
                   final.res2plot$date.2D.2[i])
        props <- c(final.res2plot$alpha[i],
                   final.res2plot$alpha2[i],
                   final.res2plot$alpha2.date1[i],
                   final.res2plot$alpha2.date2[i])
        event <- final.res2plot$Result[i]
        if(event %in% c("1D","1D(2D)"))
        {
            don1 <- donor1[1]
            don2 <- donor2[1]
            donprops1 <- paste(recip,"main",1,sep=".")
            donprops1 <- admixturesources2[donprops1,]
            donprops2 <- paste(recip,"main",2,sep=".")
            donprops2 <- admixturesources2[donprops2,]
            date <- dates[1]
            if(gsub("[0-9]","",date) != "")
            {
                date <- -as.numeric(gsub("B","",date))
            } else
            {
                date <- as.numeric(date)
            }
            prop1 <- round((as.numeric(props[1])*100)/20)
            if(prop1 > 5) prop1 <- 5
            prop1 <- linewds[prop1]
            prop2 <- round(((1-as.numeric(props[1]))*100)/20)
            if(prop2 > 5) prop2 <- 5
            prop2 <- linewds[prop2]
            x1 <- newi$Long[newi$EthnicGroup==tidyNames(recip)]
            y1 <- newi$Lat[newi$EthnicGroup==tidyNames(recip)]
            
            if(date >= date2plot[1] & date < date2plot[2])
            {
                bothdons <- c(don1,don2)
                bothdonprops <- c(donprops1*as.numeric(props[1]),donprops2*(1-as.numeric(props[1])))
            } else 
            {
                bothdons <- c()
                bothdonprops <- c()
            }
            } else if(event == "1MW")
            {
                don1 <- donor1[1]
                don2 <- donor2[1]
                don3 <- donor1[2]
                don4 <- donor2[2]
                donprops1 <- paste(recip,"main",1,sep=".")
                donprops1 <- admixturesources2[donprops1,]
                donprops2 <- paste(recip,"main",2,sep=".")
                donprops2 <- admixturesources2[donprops2,]
                donprops3 <- paste(recip,"main",3,sep=".")
                donprops3 <- admixturesources2[donprops3,]
                donprops4 <- paste(recip,"main",4,sep=".")
                donprops4 <- admixturesources2[donprops4,]
                date <- dates[1]
                if(gsub("[0-9]","",date) != "")
                {
                    date <- -as.numeric(gsub("B","",date))
                } else
                {
                    date <- as.numeric(date)
                }
                prop1 <- round((as.numeric(props[1])*100)/20)
                if(prop1 > 5) prop1 <- 5
                prop1 <- linewds[prop1]
                prop2 <- round(((1-as.numeric(props[1]))*100)/20)
                if(prop2 > 5) prop2 <- 5
                prop2 <- linewds[prop2]
                prop3 <- round((as.numeric(props[2])*100)/20)
                if(prop3 > 5) prop3 <- 5
                prop3 <- linewds[prop3]
                prop4 <- round(((1-as.numeric(props[2]))*100)/20)
                if(prop4 > 5) prop4 <- 5
                prop4 <- linewds[prop4]
                
                x1 <- newi$Long[newi$EthnicGroup==tidyNames(recip)]
                y1 <- newi$Lat[newi$EthnicGroup==tidyNames(recip)]
                
                if(date >= date2plot[1] & date < date2plot[2])
                {
                    bothdons <- c(don1,don2,don3,don4)
                    bothdonprops <- c(donprops1*as.numeric(props[1]),donprops2*(1-as.numeric(props[1])),
                                      donprops3*as.numeric(props[2]),donprops4*(1-as.numeric(props[2])))
                } else 
                {
                    bothdons <- c()
                    bothdonprops <- c()
                }
                } else if(event == "2D")
                {
                    don1 <- donor1[3]
                    don2 <- donor2[3]
                    don3 <- donor1[4]
                    don4 <- donor2[4]
                    donprops1 <- paste(recip,"main",5,sep=".")
                    donprops1 <- admixturesources2[donprops1,]
                    donprops2 <- paste(recip,"main",6,sep=".")
                    donprops2 <- admixturesources2[donprops2,]
                    donprops3 <- paste(recip,"main",7,sep=".")
                    donprops3 <- admixturesources2[donprops3,]
                    donprops4 <- paste(recip,"main",8,sep=".")
                    donprops4 <- as.vector(admixturesources2[donprops4,])
                    date <- dates[c(2,3)]
                    
                    for(j in 1:2)
                    {
                        if(gsub("[0-9]","",date[j]) != "")
                        {
                            date[j] <- -as.numeric(gsub("B","",date[j]))
                        } else
                        {
                            date[j] <- as.numeric(date[j])
                        }
                    }
                    
                    prop1 <- round((as.numeric(props[3])*100)/20)
                    if(prop1 > 5) prop1 <- 5
                    prop1 <- linewds[prop1]
                    prop2 <- round(((1-as.numeric(props[3]))*100)/20)
                    if(prop2 > 5) prop2 <- 5
                    prop2 <- linewds[prop2]
                    prop3 <- round((as.numeric(props[4])*100)/20)
                    if(prop3 > 5) prop3 <- 5
                    prop3 <- linewds[prop3]
                    prop4 <- round(((1-as.numeric(props[4]))*100)/20)
                    if(prop4 > 5) prop4 <- 5
                    prop4 <- linewds[prop4]
                    
                    x1 <- newi$Long[newi$EthnicGroup==tidyNames(recip)]
                    y1 <- newi$Lat[newi$EthnicGroup==tidyNames(recip)]
                    if(date[1] >= date2plot[1] & date[1] < date2plot[2] & date[2] >= date2plot[1] & date[2] < date2plot[2])
                    {
                        bothdons <- c(don1,don2,don3,don4)
                        bothdonprops <- c(donprops1*as.numeric(props[3]),donprops2*(1-as.numeric(props[3])),
                                          donprops3*as.numeric(props[4]),donprops4*(1-as.numeric(props[4])))
                    } else if(date[1] >= date2plot[1] & date[1] < date2plot[2])
                    {
                        bothdons <- c(don1,don2)
                        bothdonprops <- c(donprops1*as.numeric(props[3]),donprops2*(1-as.numeric(props[3])))
                    } else if(date[2] >= date2plot[1] & date[2] < date2plot[2])
                    {
                        bothdons <- c(don3,don4)
                        bothdonprops <- c(donprops3*as.numeric(props[4]),donprops4*(1-as.numeric(props[4])))
                    } else
                    {
                        bothdons <- c()
                        bothdonprops <- c()
                    }
                }
            ## sum across pops from same region  
            bothdonpropscnt <- sapply(names(bothdonprops),
                                      function(x){as.character(popkey$Country[popkey$Ethnic_Group==x])})
            bothdonprops <- as.vector(unlist(bothdonprops))
            bothdonprops[bothdonprops<propslessthan] <- 0
            bothdonpropscnt <- bothdonpropscnt[bothdonprops>=propslessthan]
            bothdonprops <- bothdonprops[bothdonprops>=propslessthan]
            
            for(don1 in unique(bothdonpropscnt))
            {
                don <- as.character(popkey$Ethnic_Group[popkey$Country==don1][1])
                if(don %in% europe)
                {
                    x2 <- as.numeric(euroxy$x[euroxy$pop==don])
                    y2 <- as.numeric(euroxy$y[euroxy$pop==don])
                    plotcol <- as.character(euroxy$cols[euroxy$pop==don])
                } else if(don %in% asia)
                {
                    x2 <- as.numeric(euroxy$x[euroxy$pop==don])
                    y2 <- as.numeric(euroxy$y[euroxy$pop==don])
                    plotcol <- as.character(euroxy$cols[euroxy$pop==don])
                } else if(don %in% americas)
                {
                    x2 <- as.numeric(euroxy$x[euroxy$pop==don])
                    y2 <- as.numeric(euroxy$y[euroxy$pop==don])
                    plotcol <- as.character(euroxy$cols[euroxy$pop==don])
                } else
                {
                    x2 <- newi$Long[newi$EthnicGroup==tidyNames(don)]
                    y2 <- newi$Lat[newi$EthnicGroup==tidyNames(don)]
                    
                    donnew <- don
                    if(donnew == "SEMI-BANTU") donnew <- "SEMI.BANTU"
                    if(donnew == "YORUBA") donnew <- "YRI"
                    if(donnew == "/GUI//GHANA_KGAL") donnew <- "GUIGHANAKGAL"
                    if(donnew == "JU/HOANSI") donnew <- "JUHOAN"
                    if(donnew == "MANDINKA") donnew <- "MANDINKAI"
                    #plotcol <- pcolshex[regions==as.character(popkey$RegionM[popkey$Ethnic_Group==donnew])]
                    plotcol <- as.character(leginfo$Colour[leginfo$EthnicGroup==tidyNames(don)])
                    
                    print(paste(plotcol,tidyNames(don)))
                }
                plotlwd <- sum(bothdonprops[bothdonpropscnt==don1])
                plotlwd <- round((as.numeric(plotlwd)*10))
                if(plotlwd > 5) plotlwd <- 5
                plotlwd <- plotlwd + 1
                plotlwd <- linewds[plotlwd]
                
            ## nw only, plot certin donors
                if(eventlist[k]== "Eurasia")
                {
                    if(don1 %in% as.character(popkey$Country[popkey$RegionM=="Eurasia"]))
                    {
                        arrows(x2,y2,x1,y1,lwd=plotlwd,length=0.1,lty=1,
                               col=makeTransparent(plotcol,200),lend=2,angle=45)
                        tmpdates <- c(tmpdates,date)
                    }
                }
                if(eventlist[k]== "Western Bantu")
                {
                    if(don1 %in% c("Cameroon"))
                    {
                        arrows(x2,y2,x1,y1,lwd=plotlwd,length=0.1,lty=1,
                               col=makeTransparent(plotcol,200),lend=2,angle=45)
                        tmpdates <- c(tmpdates,date)
                    }
                }
                if(eventlist[k]== "Eastern Bantu")
                {
                    if(don1 %in% c("Malawi","Tanzania","Kenya","SouthAfrica"))
                    {
                        arrows(x2,y2,x1,y1,lwd=plotlwd,length=0.1,lty=1,
                               col=makeTransparent(plotcol,200),lend=2,angle=45)
                        tmpdates <- c(tmpdates,date)
                    }
                }
                if(eventlist[k]== "West Africa")
                {
                    if(don1 %in% c("Gambia","BurkinaFaso","Ghana","Mali",
                                   "Nigeria","Ethiopia","Sudan","Somalia"))
                    {
                        arrows(x2,y2,x1,y1,lwd=plotlwd,length=0.1,lty=1,
                               col=makeTransparent(plotcol,200),lend=2,angle=45)
                        tmpdates <- c(tmpdates,date)
                    }
                }
            }
    }
    ## ADD PANEL LETTERS
    mtext(3,line=-2,text=paste("",LETTERS[k]),cex=2,adj=0)

    ## ADD GENE-FLOW MODEL POINTS
    if(eventlist[k]== "Eurasia")
    {
        text(as.numeric(euroxy$x)[1:9],
             as.numeric(euroxy$y)[1:9],
             labels=euroxy$pop[1:9])
        ## colonial era Khoesan
        x3 <- newi$Long[newi$EthnicGroup=="=KHOMANI"]
        y3 <- newi$Lat[newi$EthnicGroup=="=KHOMANI"]-5
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="1")
        ## medieval swahili
        x3 <- newi$Long[newi$EthnicGroup=="KAUMA"]
        y3 <- newi$Lat[newi$EthnicGroup=="KAUMA"]-5
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="3")
        ## eurasian sahara
        x3 <- newi$Long[newi$EthnicGroup=="BAMBARA"]
        y3 <- newi$Lat[newi$EthnicGroup=="BAMBARA"]+5
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="4")
        ## eurasian ethiopia
        x3 <- newi$Long[newi$EthnicGroup=="AMHARA"]
        y3 <- newi$Lat[newi$EthnicGroup=="AMHARA"]+5
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="5")
        legend("top",legend="Eurasian gene-flow into Africa",bty="n")
    }
    if(eventlist[k]== "Western Bantu")
    {
        ## western bantu
        x3 <- newi$Long[newi$EthnicGroup=="XUN"]-5
        y3 <- newi$Lat[newi$EthnicGroup=="XUN"]
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="2")
        legend("bottomleft",legend=c("0-5%","5-10%","10-20%","20-50%",">50%"),lwd=linewds,
               title="Admixture Proportions",bg="white",border="white")
        legend("top",legend="Recent Western Bantu gene-flow",bty="n")
    }
    if(eventlist[k]== "Eastern Bantu")
    {
        ## eastern bantu
        x3 <- newi$Long[newi$EthnicGroup=="MALAWI"]+5
        y3 <- newi$Lat[newi$EthnicGroup=="MALAWI"]
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="6")
        legend("top",legend="Eastern Bantu gene-flow",bty="n")
    }

    if(eventlist[k]== "West Africa")
    {
        ## east to south
        x3 <- newi$Long[newi$EthnicGroup=="JU/HOANSI"]
        y3 <- newi$Lat[newi$EthnicGroup=="JU/HOANSI"]-5
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="7")
        ## east to west
        x3 <- newi$Long[newi$EthnicGroup=="BANTU"]+10
        y3 <- newi$Lat[newi$EthnicGroup=="BANTU"]+10
        points(x3,y3,pch=21,bg="white",cex=3)
        text(x3,y3,labels="8")
        legend("top",legend="East / West gene-flow",bty="n")
    }
    print(tmpdates)
}

dev.off()

###########################################################################
### PLOT BY TIME PERIOD
pdf(paste("figures/AfricadmixtureEventsMapAllPropsOver",
          gsub("0\\.","",propslessthan),"AncestryRegions.pdf",sep=""),width=10,height=10)
layout(matrix(c(1:4),2,2,byrow=T))
for(k in 1:4)
{
    if(k == 1)
    {
        date2plot <- c(-10000,-50)
        legtext <- "before 0CE"
    }
    if(k == 2)
    {
        date2plot <- c(-50,1000)
        legtext <- "0-1000CE"
    }
    if(k == 3)
    {
        date2plot <- c(1000,1500)
        legtext <- "1000-1500CE"
    }
    if(k == 4)
    {
        date2plot <- c(1500,2000)
        legtext <- "1500-2000CE"
    } 
        ## plot map
        par(mar=c(1,1,1,1))
        mapCountryData(newmap2,nameColumnToPlot="colour",catMethod="categorical",
                       oceanCol="paleturquoise1",missingCountryCol="white",addLegend=F,
                       mapTitle="", colourPalette=rep("palegoldenrod",length(levels(mapinfo$colour))),
                       xlim=map.xlim.2plot,
                       ylim=map.ylim.2plot)
        box()
        tmpdates <- c()
        ## plot populatoin points
        for(i in 1:nrow(pointsonly)) points(pointsonly[i,2],pointsonly[i,1],pch=20,col="black",cex=1)
        linewds <- c(1,2.5,4,6,6,8)
        admixturesources2 <- read.table("data/MalariaGenGlobetrotterAdmixtureSources3.txt",header=T,row.names=1,as.is=T)
        admixturesources2 <- admixturesources2[,levels(popplot)]
        final.res2plot <- read.table("data/MalariaGenGlobetrotter2plot.txt",header=T,row.names=1,as.is=T)
        dateboots <- read.table("data/MalariaGenGlobetrotterOneDateBootstraps.txt",header=T,row.names=1,as.is=T)
        date2boots <- read.table("data/MalariaGenGlobetrotterTwoDateBootstraps.txt",header=T,row.names=1,as.is=T)
        best_ald <- read.table("data/MalariaGenBestAlder.txt",header=T,as.is=T,comment.char="",fill=T)
        
        for(i in 1:nrow(final.res2plot))
        {
            recip <- final.res2plot$Cluster[i]
            donor1 <- c(final.res2plot$best.source1[i],
                        final.res2plot$best.source1.ev2[i],
                        final.res2plot$best.source1.date1[i],
                        final.res2plot$best.source1.date2[i])
            donor2 <- c(final.res2plot$best.source2[i],
                        final.res2plot$best.source2.ev2[i],
                        final.res2plot$best.source2.date1[i],
                        final.res2plot$best.source2.date2[i])
            dates <- c(final.res2plot$date.1D[i],
                       final.res2plot$date.2D.1[i],
                       final.res2plot$date.2D.2[i])
            props <- c(final.res2plot$alpha[i],
                       final.res2plot$alpha2[i],
                       final.res2plot$alpha2.date1[i],
                       final.res2plot$alpha2.date2[i])
            event <- final.res2plot$Result[i]
            if(event %in% c("1D","1D(2D)"))
            {
                don1 <- donor1[1]
                don2 <- donor2[1]
                donprops1 <- paste(recip,"main",1,sep=".")
                donprops1 <- admixturesources2[donprops1,]
                donprops2 <- paste(recip,"main",2,sep=".")
                donprops2 <- admixturesources2[donprops2,]
                date <- dates[1]
                if(gsub("[0-9]","",date) != "")
                {
                    date <- -as.numeric(gsub("B","",date))
                } else
                {
                    date <- as.numeric(date)
                }
                prop1 <- round((as.numeric(props[1])*100)/20)
                if(prop1 > 5) prop1 <- 5
                prop1 <- linewds[prop1]
                prop2 <- round(((1-as.numeric(props[1]))*100)/20)
                if(prop2 > 5) prop2 <- 5
                prop2 <- linewds[prop2]
                x1 <- newi$Long[newi$EthnicGroup==tidyNames(recip)]
                y1 <- newi$Lat[newi$EthnicGroup==tidyNames(recip)]
                
                if(date >= date2plot[1] & date < date2plot[2])
                {
                    bothdons <- c(don1,don2)
                    bothdonprops <- c(donprops1*as.numeric(props[1]),donprops2*(1-as.numeric(props[1])))
                } else 
                {
                    bothdons <- c()
                    bothdonprops <- c()
                }
                } else if(event == "1MW")
                {
                    don1 <- donor1[1]
                    don2 <- donor2[1]
                    don3 <- donor1[2]
                    don4 <- donor2[2]
                    donprops1 <- paste(recip,"main",1,sep=".")
                    donprops1 <- admixturesources2[donprops1,]
                    donprops2 <- paste(recip,"main",2,sep=".")
                    donprops2 <- admixturesources2[donprops2,]
                    donprops3 <- paste(recip,"main",3,sep=".")
                    donprops3 <- admixturesources2[donprops3,]
                    donprops4 <- paste(recip,"main",4,sep=".")
                    donprops4 <- admixturesources2[donprops4,]
                    date <- dates[1]
                    if(gsub("[0-9]","",date) != "")
                    {
                        date <- -as.numeric(gsub("B","",date))
                    } else
                    {
                        date <- as.numeric(date)
                    }
                    prop1 <- round((as.numeric(props[1])*100)/20)
                    if(prop1 > 5) prop1 <- 5
                    prop1 <- linewds[prop1]
                    prop2 <- round(((1-as.numeric(props[1]))*100)/20)
                    if(prop2 > 5) prop2 <- 5
                    prop2 <- linewds[prop2]
                    prop3 <- round((as.numeric(props[2])*100)/20)
                    if(prop3 > 5) prop3 <- 5
                    prop3 <- linewds[prop3]
                    prop4 <- round(((1-as.numeric(props[2]))*100)/20)
                    if(prop4 > 5) prop4 <- 5
                    prop4 <- linewds[prop4]
                    
                    x1 <- newi$Long[newi$EthnicGroup==tidyNames(recip)]
                    y1 <- newi$Lat[newi$EthnicGroup==tidyNames(recip)]
                    
                    if(date >= date2plot[1] & date < date2plot[2])
                    {
                        bothdons <- c(don1,don2,don3,don4)
                        bothdonprops <- c(donprops1*as.numeric(props[1]),donprops2*(1-as.numeric(props[1])),
                                          donprops3*as.numeric(props[2]),donprops4*(1-as.numeric(props[2])))
                    } else 
                    {
                        bothdons <- c()
                        bothdonprops <- c()
                    }
                } else if(event == "2D")
                {
                    don1 <- donor1[3]
                    don2 <- donor2[3]
                    don3 <- donor1[4]
                    don4 <- donor2[4]
                    donprops1 <- paste(recip,"main",5,sep=".")
                    donprops1 <- admixturesources2[donprops1,]
                    donprops2 <- paste(recip,"main",6,sep=".")
                    donprops2 <- admixturesources2[donprops2,]
                    donprops3 <- paste(recip,"main",7,sep=".")
                    donprops3 <- admixturesources2[donprops3,]
                    donprops4 <- paste(recip,"main",8,sep=".")
                    donprops4 <- as.vector(admixturesources2[donprops4,])
                    date <- dates[c(2,3)]
                    
                    tmp <- date
                    tmp[grep("B",tmp)] <- paste("-",tmp[grep("B",tmp)],sep="")
                    tmp <- as.vector(sapply(gsub("B","",tmp),as.numeric))
                    date <- tmp
                    
                    prop1 <- round((as.numeric(props[3])*100)/20)
                    if(prop1 > 5) prop1 <- 5
                    prop1 <- linewds[prop1]
                    prop2 <- round(((1-as.numeric(props[3]))*100)/20)
                    if(prop2 > 5) prop2 <- 5
                    prop2 <- linewds[prop2]
                    prop3 <- round((as.numeric(props[4])*100)/20)
                    if(prop3 > 5) prop3 <- 5
                    prop3 <- linewds[prop3]
                    prop4 <- round(((1-as.numeric(props[4]))*100)/20)
                    if(prop4 > 5) prop4 <- 5
                    prop4 <- linewds[prop4]
                    
                    x1 <- newi$Long[newi$EthnicGroup==tidyNames(recip)]
                    y1 <- newi$Lat[newi$EthnicGroup==tidyNames(recip)]
                    if(date[1] >= date2plot[1] & date[1] < date2plot[2] & date[2] >= date2plot[1] & date[2] < date2plot[2])
                    {
                        bothdons <- c(don1,don2,don3,don4)
                        bothdonprops <- c(donprops1*as.numeric(props[3]),donprops2*(1-as.numeric(props[3])),
                                          donprops3*as.numeric(props[4]),donprops4*(1-as.numeric(props[4])))
                    } else if(date[1] >= date2plot[1] & date[1] < date2plot[2])
                    {
                        bothdons <- c(don1,don2)
                        bothdonprops <- c(donprops1*as.numeric(props[3]),donprops2*(1-as.numeric(props[3])))
                    } else if(date[2] >= date2plot[1] & date[2] < date2plot[2])
                    {
                        bothdons <- c(don3,don4)
                        bothdonprops <- c(donprops3*as.numeric(props[4]),donprops4*(1-as.numeric(props[4])))
                    } else
                    {
                        bothdons <- c()
                        bothdonprops <- c()
                    }
                }
                ## sum across pops from same region  
                bothdonpropscnt <- sapply(names(bothdonprops),
                                          function(x){as.character(popkey$Country[popkey$Ethnic_Group==x])})
                bothdonprops <- as.vector(unlist(bothdonprops))
                bothdonprops[bothdonprops<propslessthan] <- 0
                bothdonpropscnt <- bothdonpropscnt[bothdonprops>=propslessthan]
                bothdonprops <- bothdonprops[bothdonprops>=propslessthan]
                
                for(don1 in unique(bothdonpropscnt))
                {
                    don <- as.character(popkey$Ethnic_Group[popkey$Country==don1][1])
                    if(don %in% europe)
                    {
                        x2 <- as.numeric(euroxy$x[euroxy$pop==don])
                        y2 <- as.numeric(euroxy$y[euroxy$pop==don])
                        plotcol <- as.character(euroxy$cols[euroxy$pop==don])
                    } else if(don %in% asia)
                    {
                        x2 <- as.numeric(euroxy$x[euroxy$pop==don])
                        y2 <- as.numeric(euroxy$y[euroxy$pop==don])
                        plotcol <- as.character(euroxy$cols[euroxy$pop==don])
                    } else if(don %in% americas)
                    {
                        x2 <- as.numeric(euroxy$x[euroxy$pop==don])
                        y2 <- as.numeric(euroxy$y[euroxy$pop==don])
                        plotcol <- as.character(euroxy$cols[euroxy$pop==don])
                    } else
                    {
                        x2 <- newi$Long[newi$EthnicGroup==tidyNames(don)]
                        y2 <- newi$Lat[newi$EthnicGroup==tidyNames(don)]
                        
                        donnew <- don
                        if(donnew == "SEMI-BANTU") donnew <- "SEMI.BANTU"
                        if(donnew == "YORUBA") donnew <- "YRI"
                        if(donnew == "/GUI//GHANA_KGAL") donnew <- "GUIGHANAKGAL"
                        if(donnew == "JU/HOANSI") donnew <- "JUHOAN"
                        if(donnew == "MANDINKA") donnew <- "MANDINKAI"
                        #plotcol <- pcolshex[regions==as.character(popkey$RegionM[popkey$Ethnic_Group==donnew])]
                        plotcol <- as.character(leginfo$Colour[leginfo$EthnicGroup==tidyNames(don)])
                        
                        print(paste(plotcol,tidyNames(don)))
                    }
                    plotlwd <- sum(bothdonprops[bothdonpropscnt==don1])
                    plotlwd <- round((as.numeric(plotlwd)*10))
                    if(plotlwd > 5) plotlwd <- 5
                    plotlwd <- plotlwd + 1
                    plotlwd <- linewds[plotlwd]
                    
                    ## nw only, plot certin donors
                    arrows(x2,y2,x1,y1,lwd=plotlwd,length=0.1,lty=1,
                           col=makeTransparent(plotcol,200),lend=2,angle=45)
                    tmpdates <- c(tmpdates,date)
                }
            }
    ## ADD EURASIA LABELS
    text(as.numeric(euroxy$x)[c(1,2,4,5,7,9)],
         as.numeric(euroxy$y)[c(1,2,4,5,7,9)],
         labels=euroxy$pop[c(1,2,4,5,7,9)])
    ## ADD PANEL LETTERS
    mtext(3,line=-2,text=paste("",LETTERS[k]),cex=2,adj=0)
    legend("top",legend=legtext,bty="n")
    
    ## ADD LINE LEGEND
    if(k == 1)
    {
        legend("bottomleft",legend=c("0-5%","5-10%","10-20%","20-50%",">50%"),lwd=linewds,
               title="Admixture Proportions",bg="white",border="white")
    }
    
    
}
        
dev.off()
    
