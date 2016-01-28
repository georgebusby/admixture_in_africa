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
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
popnums_file <- "data/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"
pca_file <- "data/Africa300KPCS.txt"
tree_file <- paste0("data/",fsanalyname,".mcmc.tree.xml")
mat_file <- paste0("data/",fsanalyname,"CoAncestry.txt")

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)

pdf("figures/LocalAncestryOverview.pdf",width=6,height=6)
    par(mar=c(0.5,0.5,0.5,0.5))
    layout(matrix(c(1,5,
                    2,2,
                    3,3,
                    4,4),4,2,byrow=T),
           heights=c(3,0.75,0.75,1.5),widths=c(3,3))
           
    ###########################################################
    ## 00 MAP AND LEGENDS
    transparent_cols <- FALSE ## use transparent colours for malariagen countries
    plot_equator <- FALSE ## should equator and tropics be plotted
    plot_ancestry_regions <- TRUE ## should countries be coloured by ancestry region or individually?
    map_xlim<-c(-20,50) ## X-AXIS LIMITS OF MAP
    map_ylim<-c(-32,28) ## Y-AXIS LIMITS OF MAP
    #regions <- unique(popkey$RegionM) ## regions
    pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
    regions <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
    map_ocean_col <- "paleturquoise1" ## colour of ocean in map
    map_miss_col <- "white" ## colour of non-focal countries in map
    pt_cex <- 2
    pt_lwd <- 0.75
    plot_legends <- FALSE
    plot_letter <- FALSE
    plot_points <- FALSE
    source("plotting_scripts/africamap.R")
    ## ADD POINTS FOR COUNTRY SAMPLING POSITIONS
    latlongs <- read.csv(latlong_file,header=T)
    points(latlongs$Long,latlongs$Lat,pch=20,cex=1)    
    ## ADD A LEGEND
    legend("bottomleft",legend=gsub("Africa ","",gsub("\\_"," ",regions)),
           fill=pcolshex,border=pcolshex,bty="n",title="Ancestry Regions")

    #########################################################
    ## 01 PLOT A PAINTED CHROMOSOME
    library("rhdf5")
    ## FUNCTIONS
    hap2sampleindex <- function(hap,nsamps=10){
        ## finds the first sample index for a haplotype
        sample <- (hap*nsamps)-(nsamps-1)
        return(sample)
    }

    datafile <- '/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
    ## CHROMOSOME 2
    chrom <- "02"
    ## GET MAP AND POSITION INFO
    map <- data.frame(t(h5read(datafile,paste0("/paintings/chrom",chrom,"/map"))))
    H5close()
    colnames(map) <- c("position","recrate")
    ## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
    psamples <- t(h5read(datafile,paste0("/paintings/samples/individuals")))
    H5close()
    colnames(psamples) <- c("ind","region","X")
    psamplesind <- (1:nrow(psamples))[psamples[,"ind"] == "FULA_S3"]
    ## 2 haps per sample!!
    psampleshap <- hap2sampleindex(psamplesind,2)
    psamplesindsamp <- hap2sampleindex(psampleshap)
    
    for(analysis in c("local","nonlocal"))
    {
        tmp <- psamplesindsamp:(psamplesindsamp+9)
        paintedchrom <- h5read(datafile,paste0("/paintings/chrom",chrom,"/",analysis),
                               index=list(tmp,1:nrow(map)))
        H5close()
        ## SWITCH DONORS TO REGIONS
        happops <- c()
        for(i in 1:nrow(psamples)) happops <- c(happops,rep(as.character(psamples[i,"region"]),2))
        happops <- gsub("SEMI.BANTU","SEMI-BANTU",happops)
        hapregs <- c()
        for(i in happops) hapregs <- c(hapregs,as.character(popkey$RegionM[popkey$Ethnic_Group==i]))
        hapregs <- factor(hapregs,levels=ancreg_list)
        paintedchromreg <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
        for(i in 1:nrow(paintedchrom)) paintedchromreg[i,] <- hapregs[paintedchrom[i,]]
        ## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
        paintedchromregprop <- matrix(0,nr=length(ancreg_list),nc=ncol(paintedchromreg))
        for(i in 1:ncol(paintedchromreg))
        {
            print(i)
            tmp <- table(paintedchromreg[,i])
            paintedchromregprop[as.numeric(names(tmp)),i] <- tmp
        }
    
        ## BAR WIDTHS
        chromlength <- as.numeric(as.character(map$recrate))
        chrompos <- as.numeric(as.character(map$position))
        chromposI <- c(diff(chrompos),0)
        chrom2plot <- paintedchromregprop[8:1,]
       #chrom2plot <- t(apply(chrom2plot,2,which.max))
    
        if(analysis == "local") par(mar=c(0.5,1,1.5,1))
        if(analysis == "nonlocal") par(mar=c(0.5,1,1.5,1))
        barplot(chrom2plot,
                width=chromposI,
                col=rev(pcolshex),xaxs="i",yaxs="i",
                space=0,axes=F,xaxt="n",border=NA,
                xlim=c(0,sum(chromposI)),
                xlab="")
        mtext(3,text=paste(analysis,"painting of a Fulani individual: one haplotype x 10 samples"),
              line=0.5,cex=0.75)
    }

    ## NOW SHOW PAINTINGS ACROSS ALL FULAI
        psamplesind <- (1:nrow(psamples))[psamples[,"region"] == "FULAI"]
        ## 2 haps per sample!!
        psampleshap <- hap2sampleindex(psamplesind,2)
        psampleshap <- sort(c(psampleshap,psampleshap+1))
        psamplesindsamp <- hap2sampleindex(psampleshap)
        tmp <- c()
        for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))
        paintedchrom <- h5read(datafile,paste0("/paintings/chrom",chrom,"/nonlocal"),
                               index=list(tmp,1:nrow(map)))
        H5close()
        ## SWITCH DONORS TO REGIONS
        happops <- c()
        for(i in 1:nrow(psamples)) happops <- c(happops,rep(as.character(psamples[i,"region"]),2))
        happops <- gsub("SEMI.BANTU","SEMI-BANTU",happops)
        hapregs <- c()
        for(i in happops) hapregs <- c(hapregs,as.character(popkey$RegionM[popkey$Ethnic_Group==i]))
        hapregs <- factor(hapregs,levels=ancreg_list)
        paintedchromreg <- matrix(nc=ncol(paintedchrom),nr=nrow(paintedchrom))
        for(i in 1:nrow(paintedchrom)) paintedchromreg[i,] <- hapregs[paintedchrom[i,]]
        ## NOW GET PROPOTION OF ANCESTRY AT EACH SNP
        paintedchromregprop <- matrix(0,nr=length(ancreg_list),nc=ncol(paintedchromreg))
        for(i in 1:ncol(paintedchromreg))
        {
            print(i)
            tmp <- table(paintedchromreg[,i])
            paintedchromregprop[as.numeric(names(tmp)),i] <- tmp
        }
        
        ## BAR WIDTHS
        chromlength <- as.numeric(as.character(map$recrate))
        chrompos <- as.numeric(as.character(map$position))
        chromposI <- c(diff(chrompos),0)
        chrom2plot <- paintedchromregprop[8:1,]
        #chrom2plot <- t(apply(chrom2plot,2,which.max))
        
        par(mar=c(3,1,2,1))
        barplot(chrom2plot,
                width=chromposI,
                col=rev(pcolshex),xaxs="i",yaxs="i",
                space=0,axes=F,xaxt="n",border=NA,
                xlim=c(0,sum(chromposI)),
                xlab="")
        mtext(3,text=paste("all nonlocal paintings of 72 Fulani individuals"),line=0.5,cex=0.75)
        xat <- pretty(chrompos)
        xat[length(xat)] <- sum(chromposI)
        axis(1,at=xat,labels=round(xat/1e6))
        mtext(1,text=paste("position on chromosome",as.numeric(chrom)),line=2,cex=0.75)

        ## FILLER
        par(mar=c(0,0.5,0,0))
        plot(0,0,xlab="",ylab="",type="n",axes=F)
        mtext(2,text="ARG?",cex=5,las=1,adj=0)

dev.off()