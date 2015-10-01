#' this is an internal function that finds new lat/longs for groups where there are multiple populations from the same country and plots onto a map
#'
#' This function re-positions points on a map using three different input files, for use with plotting samples on a map, for example of Africa
#' @param latlong_file a CSV file containing all ethnic groups in the analysis, must have the following columns: "EthnicGroup" "Country" "Lat" "Long"
#' @param leginfo_file a file matching the ethnic groups to plotting parameters, must have the following columns: "EthnicGroup" "Country" "Region" "Colour" "pch" "rim", where rim defines whether the rim is plotted as the same colour as the symbol
#' @param poppos_file a file outlining for particular countries, whether latitude and/or longitude should be shifted up or down from the original country latitude or longitude
#' @param pt_cex the size of the symbols on the map
#' @param pt_lwd the width of the outline to the symbols
#' @keywords Busby_bespoke
#' @return plots new points and returns a dataframe with plotting parameters (helpful for legends)
#' @export
#' @examples
#' plotPopPoints("data/MalariaGenAdmixturePopulationKeyLatLongs.txt",
#'               "data/MalariaGenAdmixturePopulationKey.txt",
#'               "data/MalariaGenAdmixturePopulationKeyMapPositions.txt",
#'               pt_cex=2,pt_lwd=0.75)




plotPopPoints <- function(latlong_file,leginfo_file,poppos_file,pt_cex=1,pt_lwd=0.5)
{
    ll <- read.csv(latlong_file,header=T)
    ll$EthnicGroup <- toupper(ll$EthnicGroup)
    leginfo <-read.table(leginfo_file,header=T,comment.char="")
    leginfo$EthnicGroup <- toupper(leginfo$EthnicGroup)
    ## READ IN CONTROLS FOR POSITIONING OF MULTIPLE POINTS FROM SAME LAT/LONG
    poppos <- read.table(poppos_file,header=T)
    ## GENERATE A NEW SET OF LAT/LONGS FOR POINTS, ALLOWING FOR MULTIPLE POINTS FROM THE SAME POSITION
    newi <- c()
    for(i in levels(ll$Country))
    {
        ii <- ll[ll$Country==i,]
        ii <- cbind(ii,ii$Lat,ii$Long)
        colnames(ii)[c((ncol(ii)-1),(ncol(ii)))] <- c("newlat","newlong")
        if(nrow(ii)>1)
        {
            dups <- ii[duplicated(ii[,c("Lat","Long")]),"EthnicGroup"]
            nd <- length(dups)
            if(nd>0)
            {
                splitlat <- as.character(poppos[poppos$Country==i,"latlong"])
                if(splitlat=="newlat") splitlong <- "newlong"
                if(splitlat=="newlong") splitlong <- "newlat"
                updown <- poppos[poppos$Country==i,"updown"] * 2
                for(j in 1:nd)
                {
                    if(j%%2==1)
                    {
                        ii[duplicated(ii[,c("Lat","Long")]),splitlat][j] <- ii[duplicated(ii[,c("Lat","Long")]),splitlat][j]-(((j/2)+0.5)*3)
                    }
                    if(j%%2==0)
                    {
                        ii[duplicated(ii[,c("Lat","Long")]),splitlat][j] <- ii[duplicated(ii[,c("Lat","Long")]),splitlat][j]+((j/2)*3)
                    }
                }
                dups <- unique(ii[duplicated(ii[,c("Lat","Long")]),splitlong])
                ii[ii[,splitlong]==dups,splitlong] <- dups+updown
                ## sort points so that they are alphabetical
                ii[ii[,splitlong]==(dups+updown),splitlat] <- sort(ii[ii[,splitlong]==(dups+updown),splitlat],decreasing=T)
                ## add points and lines
                points(ii$Long,ii$Lat,pch=20,col="black",cex=1)
                for(k in 1:length(ii$Long))
                {
                    lines(x=c(ii$Long[k],ii$newlong[k]),y=c(ii$Lat[k],ii$newlat[k]),col="black",lwd=0.5)
                }
            }
        }
        newi <- rbind(newi,ii)
    }
    ## NOW GET POINT INFO FROM LEGINFO FILE AND PLOT WITH POSITIONS IN newi
    pops <- newi$EthnicGroup
    pntinfo <- getPopSymbols(pops,leginfo)
    pntinfo <- cbind(ll,pntinfo)
    points(newi$newlong,newi$newlat,
           pch=as.numeric(as.character(pntinfo$pch2plot)),
           bg=as.character(pntinfo$col2plot),
           col=as.character(pntinfo$rim2plot),
           cex=pt_cex,lwd=pt_lwd)
    return(pntinfo)
    
}
