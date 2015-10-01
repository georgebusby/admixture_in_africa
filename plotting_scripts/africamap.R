############################################################
## SCRIPT TO PLOT NICE MAP OF AFRICA WITH MALARIAGEN POPS ##
## HIGHLIGHTED                                            ##
############################################################
## LOAD MAP LIBRARIES ##
library("rworldmap")
library("maptools")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"
poppos_file <- "data/MalariaGenAdmixturePopulationKeyMapPositions.txt"
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"
###########################################################
## LOAD DATA
## file outlining population provenance, colours, and symbols
mapinfo <- read.table(mapinfo_file,header=T,comment.char="", as.is = T)
## population latitudes and longitudes

## controls how pops with same lat/long are plotted

## define the region that a population belongs to 
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
###########################################################
## SET PLOTTING PARAMETERS
transparent_cols <- FALSE ## use transparent colours for malariagen countries
plot_equator <- FALSE ## should equator and tropics be plotted
plot_ancestry_regions <- FALSE ## should countries be coloured by ancestry region or individually?
map_xlim<-c(-20,50) ## X-AXIS LIMITS OF MAP
map_ylim<-c(-32,28) ## Y-AXIS LIMITS OF MAP
regions <- unique(popkey$RegionM) ## regions
pcolshex <- c("#0000CD", "#03B4CC", "#A65628", "#FF7F00", "#984EA3", "#4DAF4A", "#CCCC00") ## colours
map_ocean_col <- "paleturquoise1" ## colour of ocean in map
map_miss_col <- "white" ## colour of non-focal countries in map

###########################################################
## SET UP THE MAPINFO DATAFRAME
## make the country colour transparent?
if(transparent_cols == T) mapinfo$colour <- makeTransparent(mapinfo$colour,alpha=175)
mapinfo <- data.frame(mapinfo)
## set colour to factor
mapinfo$colour <- factor(mapinfo$colour,levels=unique(mapinfo$colour))
## generate a vector of regions for each ethnic group
newcols <- findPopRegion(mapinfo$EthnicGroup,popkey)
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

## PLOT POINTS, ALLOWING FOR MULTIPLE POPS FROM THE SAME COUNTRY
newi <- plotPopPoints(latlong_file,leginfo_file,poppos_file,pt_cex=2)


### newi can now be used to make legends ..
 
