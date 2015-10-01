#' this is an internal function that finds new lat/longs for groups where there are multiple populations from the same country
#'
#' This function re-positions points on a map using three different input files
#' @param latlong_file a file containing all ethnic groups in the analysis, must have the following columns: "EthnicGroup" "Country" "Lat" "Long"
#' @param leginfo_file a file matching the ethnic groups to plotting parameters, must have the following columns: "EthnicGroup" "Country" "Region" "Colour" "pch" "rim", where rim defines whether the rim is plotted as the same colour as the symbol
#' @param poppos_file a file outlining for particular countries, whether latitude and/or longitude should be shifted up or down from the original country latitude or longitude
#' @keywords Busby_bespoke
#' @export returns a dataframe with plotting parameters
#' @examples
#' plotPopPoints(latlong_file,leginfo_file,poppos_file)




plotPopPoints <- function(latlong_file,leginfo_file,poppos_file)
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
    p <- r <- co <- c()
    for(i in newi$EthnicGroup)
    {
        ii <- leginfo[leginfo$EthnicGroup==i,]
        p <- c(p,ii$poppch)
        r <- c(r,ii$rim)
        co <- c(co,as.character(ii$Colour))
    }
    
    pntinfo <- cbind(ll,co,p,r)
    pntcol <- rep("#000000",nrow(ll))
    pntcol[pntinfo$r==1]  <- as.character(pntinfo$co)[pntinfo$r==1]
    pntinfo <- cbind(pntinfo,pntcol)
    points(newi$newlong,newi$newlat,pch=pntinfo$p,cex=pt_cex,lwd=0.5,bg=as.character(pntinfo$co),col=pntcol)
    return(pntinfo)
    
}
