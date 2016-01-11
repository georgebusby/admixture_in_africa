############################################################
## SCRIPT TO PLOT NICE MAP OF AFRICA WITH MALARIAGEN POPS ##
## HIGHLIGHTED                                            ##
############################################################
## LOAD MAP LIBRARIES ##
library("rworldmap")
library("maptools")
###########################################################
## LOAD DATA
## file outlining population provenance, colours, and symbols
mapinfo <- read.table(mapinfo_file,header=T,comment.char="", as.is = T)

###########################################################
## SET UP THE MAPINFO DATAFRAME
## make the country colour transparent?
if(transparent_cols == T) mapinfo$colour <- makeTransparent(mapinfo$colour,alpha=175)
mapinfo <- data.frame(mapinfo)
## set colour to factor
mapinfo$colour <- factor(mapinfo$colour,levels=unique(mapinfo$colour))
## generate a vector of regions for each ethnic group
tmpgrp <- gsub("MANDINKA","MANDINKAII",mapinfo$EthnicGroup)
newcols <- getPopRegion(tmpgrp,popkey)
newcols <- factor(newcols, levels=regions)
mapinfo <- cbind(mapinfo,newcols)
colnames(mapinfo)[ncol(mapinfo)] <- "region"
## convert to colours
newcols <- pcolshex[newcols]
## add colour column to mapinfo
mapinfo <- cbind(mapinfo,newcols)
## BRING IN A MAP FROM THE rworldmaps PACKAGE
newmap <- getMap(resolution="low")

###########################################################
###########################################################
###########################################################
## PLOT MAP WITH COUNTRY COLOURS
if(plot_ancestry_regions == FALSE)
{
    newmap2 <- joinCountryData2Map(mapinfo,joinCode="ISO_A3",nameJoinColumn="pop",verbose=FALSE)
    mapCountryData(newmap2,nameColumnToPlot="colour",catMethod="categorical",
                   oceanCol=map_ocean_col,missingCountryCol=map_miss_col,addLegend=F,
                   mapTitle="", colourPalette=makeTransparent(levels(mapinfo$colour)),
                   xlim=map_xlim,ylim=map_ylim)
}

###########################################################
### OR PLOT MAP 7 REGION ANCESTRIES
if(plot_ancestry_regions == TRUE)
{
    newmap2 <- joinCountryData2Map(mapinfo,joinCode="ISO_A3",nameJoinColumn="pop",verbose=FALSE)
    mapCountryData(newmap2,nameColumnToPlot="newcols",catMethod="categorical",
                   oceanCol=map_ocean_col,missingCountryCol=map_miss_col,addLegend=F,
                   mapTitle="",colourPalette=levels(mapinfo$newcols),
                   xlim=map.xlim.2plot,
                   ylim=map.ylim.2plot)
}

## PLOT EQUATOR AND TROPICS?
if(plot_equator==T)
{
    abline(h=0,lwd=1)
    abline(h=23, lty=2,lwd=1)
    abline(h=-23, lty=2,lwd=1)
}

box()

if(plot_letter == TRUE) legend("bottomright",legend="A",cex=2,bty="n")

## PLOT POINTS, ALLOWING FOR MULTIPLE POPS FROM THE SAME COUNTRY
if(plot_points == TRUE) newi <- plotPopPoints(latlong_file,leginfo_file,poppos_file,pt_cex=2)

###########################################################
## PLOT THREE SEPARATE LEGENDS
if(plot_legends == TRUE)
{
    afr_cnts <- list(c("Gambia","Mali","BurkinaFaso","Ghana","Nigeria","Cameroon"),
                     c("Sudan","Ethiopia","Somalia","Kenya"),
                     c("Angola","Namibia","Botswana","SouthAfrica","Tanzania","Malawi"))
    names(afr_cnts) <- c("WEST/CENTRAL","EASTERN","SOUTHERN")
    malgen_cnts <- c("Gambia","Mali","BurkinaFaso","Ghana","Cameroon","Kenya","Tanzania","Malawi")
    for(i in 1:length(afr_cnts))
    {
        leginfo <-read.table(leginfo_file,header=T,comment.char="")
        legnumbers <- read.table(popnums_file,header=F)
        legnumbers[,1] <- tidyNames(legnumbers[,1])
        pop_nums <- c()
        for(j in as.character(leginfo$EthnicGroup))
        {
            pn <- legnumbers[legnumbers$V1==j,2]
            if(length(pn)>0)
            {
                pop_nums <- c(pop_nums,pn)    
            }  else {
                pop_nums <- c(pop_nums,"X")
            }
        }
        leginfo <- cbind(leginfo,pop_nums)
        leginfo <- leginfo[leginfo$Country%in%unlist(afr_cnts[i]),]
        leginfo$Country <- factor(leginfo$Country,levels=unlist(afr_cnts[i]))
        leginfo <- leginfo[order(leginfo$Country,leginfo$Colour,leginfo$EthnicGroup),]
        if(i==3)
        {
            leginfo$Country <- factor(leginfo$Country,levels=c("Tanzania","Malawi","SouthAfrica","Namibia","Botswana","Angola"))
            leginfo <- leginfo[order(leginfo$Country,decreasing=F),]
        }
        legrim <- rep("#000000",nrow(leginfo))
        legrim[leginfo$rim==1]  <- as.character(leginfo$Colour)[leginfo$rim==1]
        
        par(mar=c(0,0,0,0))
        plot(0,0,type="n",axes=F,xlab="",ylab="")
        legend("topleft",ncol=1,
               legend=c(paste0(tidyNames(as.character(leginfo$EthnicGroup),khoesan=T),
                                               " [",as.character(leginfo$pop_nums),"]")),
               pch=as.numeric(leginfo$poppch),
               pt.bg=as.character(leginfo$Colour),
               pt.cex=c(rep(pt_cex,length(leginfo$poppch)),0),
               pt.lwd=0.5,y.intersp = 1.5,
               col=legrim,cex=0.75,
               text.col=as.character(leginfo$Colour),xpd=T,
               bty="n",
               title.col="black",title.adj=0)
    }
}    