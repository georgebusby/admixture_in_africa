### SCRIPT TO CALL THE VARIOUS PLOTTING SCRIPTS TO PRODUCE A SINGLE PLOT ###
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
fsanalyname <- 'MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalForce.A'
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"
poppos_file <- "data/MalariaGenAdmixturePopulationKeyMapPositions.txt"
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"
popnums_file <- "data/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"
pca_file <- "data/Africa300KPCS.txt"
tree_file <- paste0("data/",fsanalyname,".mcmc.tree.xml")
mat_file <- paste0("data/",fsanalyname,"CoAncestry.txt")

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)

# pdf("figures/Figure1newLS.pdf",width=10,height=5)
#     par(mar=c(0.5,0.5,0.5,0.5))
#     layout(matrix(c(1,2,3,4,11,10,
#                     1,2,3,4,7,8,
#                     5,6,6,6,7,8,
#                     5,6,6,6,12,9),4,6,byrow=T),
#            width=c(2.5,rep((2.5/3),3),1,4),
#            height=c(1,1.5,2,0.5))
pdf("figures/Figure1new.pdf",width=6,height=12)
    par(mar=c(0.5,0.5,0.5,0.5))
    layout(matrix(c(1,1,2,3,4,
                    5,5,6,6,6,
                    11,10,10,10,10,
                    7,8,8,8,8,
                    12,9,9,9,9),5,5,byrow=T),
           width=c(1,2,1,1,1),
           height=c(3,3,0.5,5,0.5))


    ###########################################################
    ## 00 MAP AND LEGENDS
    transparent_cols <- FALSE ## use transparent colours for malariagen countries
    plot_equator <- FALSE ## should equator and tropics be plotted
    plot_ancestry_regions <- FALSE ## should countries be coloured by ancestry region or individually?
    map_xlim<-c(-20,50) ## X-AXIS LIMITS OF MAP
    map_ylim<-c(-32,28) ## Y-AXIS LIMITS OF MAP
    regions <- unique(popkey$RegionM) ## regions
    pcolshex <- c("#0000CD", "#03B4CC", "#A65628", "#FF7F00", "#984EA3", "#4DAF4A", "#CCCC00") ## colours
    map_ocean_col <- "paleturquoise1" ## colour of ocean in map
    map_miss_col <- "white" ## colour of non-focal countries in map
    pt_cex <- 2
    pt_lwd <- 0.75
    plot_legends <- TRUE
    plot_letter <- TRUE
    plot_points <- TRUE
    source("plotting_scripts/africamap.R")
    ###########################################################
    ## 01 PC 1 v 2
    xplot <- 2 ## x-axis pc to plot
    yplot <- 1 ## y-axis pc to plot
    revX <- T ## should x-axis be reversed
    revY <- T ## should y-axis be reversed
    panel_letter <- "B"
    par(mar=c(0.5,0.5,0.5,0.5))
    source("plotting_scripts/africapca.R")
    ###########################################################
    ## 02 PC 1 v 3
    xplot <- 3 ## x-axis pc to plot
    yplot <- 1 ## y-axis pc to plot
    revX <- F ## should x-axis be reversed
    revY <- T ## should y-axis be reversed
    panel_letter <- "C"
    source("plotting_scripts/africapca.R")
    ###########################################################
    ## 03 FS TREE
    plot_vert <- TRUE
    plot_tree_labels <- FALSE
#     rect_splits <- c(0,12,237,910,1176,1895,2198)
#     rect_cols <- c("#CCCC00","#32CD32","#FF7F00","#8B0A50","#0000CD","#03B4CC")
#     rect_labels <- c("","SOUTH_AFRICA_[KS / NC]","EAST_AFRICA_[NC]","EAST_AFRICA_[AA / NS]",
#                      "WEST_AFRICA_[NC]","CENTRAL_WEST_AFRICA_[NC]")
    ## THIS IS FOR ALL EIGHT ANCESTRY REGIONS    
    rect_splits <- c(0,12,195,339,910,1080,1176,1895,2198)
    rect_cols <- c("#CCCC00","#32CD32","#A65628","#FF7F00","#FF69B4","#984EA3","#0000CD","#03B4CC")
    rect_labels <- c("","SOUTH_AFRICA_[KS]","SOUTH_AFRICA_[NC]",
                     "EAST_AFRICA_[NC]","EAST AFRICA_[AA]","EAST AFRICA_[NS]",
                     "WEST_AFRICA_[NC]","CENTRAL_WEST_AFRICA_[NC]")

    rect_split_x <- 1250 ## where to make the rectangles more transparent
    par(mar=c(0.5,0.5,0,0))
    source("plotting_scripts/africatree.R")
    ##########################################################
    ## 04 FS COPYING MATRIX
    tmatmax <- 30 # cap the heatmap
    tmatmin <- 0
    plot_scale <- TRUE
    par(mar=c(0.5,0,0,0.5))
    source("plotting_scripts/africacopyingmatrix.R")
    ###########################################################
    ## 05 FS TREE_2
    plot_vert <- FALSE
    par(mar=c(0,0,1,0.5))
    source("plotting_scripts/africatree.R")
    ###########################################################
    ## 06 PLOT PANEL LETTER TOP LEFT ABOVE TREE
    par(mar=c(0,0,1,0))    
    plot(0,0,type="n",axes=F,xlab="",ylab="")
    legend("topleft",legend="D",cex=2,bty="n")


dev.off()
