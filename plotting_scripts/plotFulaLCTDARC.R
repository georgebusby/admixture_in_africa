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


duffy <- 159e6
lct <- 136e6

duffy_reg <- c()
lct_reg <- c()

pdf("figures/FULAIancestryLCTDARC.pdf",width=10,height=6)
    layout(matrix(c(1,5,
                    3,7,
                    2,6,
                    4,8),4,2,byrow=T),
           heights=c(3,1,0.75,1.25))
    leftmar <- 3
    ## 01 PLOT A PAINTED CHROMOSOME
    library("rhdf5")
    ## FUNCTIONS
    hap2sampleindex <- function(hap,nsamps=10){
        ## finds the first sample index for a haplotype
        sample <- (hap*nsamps)-(nsamps-1)
        return(sample)
    }
    datafile <- '/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
    for(region in c(duffy,lct))
    {
        if(region == lct) chrom <- "02"
        if(region == duffy) chrom <- "01"
        ## GET MAP AND POSITION INFO
        map <- data.frame(t(h5read(datafile,paste0("/paintings/chrom",chrom,"/map"))))
        H5close()
        colnames(map) <- c("position","recrate")
        ## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
        psamples <- t(h5read(datafile,paste0("/paintings/samples/individuals")))
        H5close()
        colnames(psamples) <- c("ind","region","X")
        psamplesind <- (1:nrow(psamples))[psamples[,"region"] == "FULAI"]
        ## 2 haps per sample!!
        psampleshap <- hap2sampleindex(psamplesind,2)
        psampleshap <- sort(c(psampleshap,psampleshap+1))
        psamplesindsamp <- hap2sampleindex(psampleshap)
        tmp <- c()
        for(i in psamplesindsamp) tmp <- c(tmp,i:(i+9))
        mappos <- map$position
        mappos <- as.numeric(as.character(mappos))
        tmp1 <- (1:nrow(map))[mappos>=(region-5e6)&mappos<=(region+5e6)]
        
        paintedchrom <- h5read(datafile,paste0("/paintings/chrom",chrom,"/nonlocal"),
                               index=list(tmp,tmp1))
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
        
        ###########################################################
        ## PLOT ACTUAL HAPLOTYPES
        chromlength <- as.numeric(as.character(map$recrate)[tmp1])
        chrompos <- as.numeric(as.character(map$position)[tmp1])
        chromposI <- c(diff(chrompos),0)
        chromplot <- paintedchromreg
        chromcols <- pcolshex[sort(unique(unlist(apply(paintedchromreg,2,unique))))]
        chromplot[chromplot!=8] <- 2
        chromcols[1] <- "black"
        par(mar=c(0,leftmar,3,1))
        image(chrompos,
              1:nrow(chromplot),
              t(chromplot),
              col=chromcols,
              axes=F,xlab="",ylab="")
        xat <- pretty(chrompos)
        xat[length(xat)] <- max(chrompos)
        axis(3,at=xat,labels=round(xat/1e6))
        mtext(3,text=paste("position on chromosome",as.numeric(chrom),"(cM)"),line=2,cex=0.75)
        
        
        ## FIND MIN SNPS
        propeur <- apply(chromplot,2,function(x){sum(x==8)}/length(x))
        if(region == duffy) snps <- chrompos[propeur==min(propeur)]
        if(region == lct) snps <- chrompos[propeur==max(propeur)]
        innerregion <- range(snps) + c(-2.5e5,2.5e5)
        
        abline(v=innerregion[1],col="red",lwd=2)
        abline(v=innerregion[2],col="red",lwd=2)
        
        ###########################################################
        ## PLOT PROPORTIONS
        chrom2plot <- paintedchromregprop[8:1,]
        par(mar=c(0,leftmar,0.5,1))
        barplot(chrom2plot,
                width=chromposI,
                col=rev(pcolshex),xaxs="i",yaxs="i",
                space=0,axes=F,xaxt="n",border=NA,
                xlim=c(0,sum(chromposI)),
                xlab="")
        abline(v=innerregion[1]-chrompos[1],col="red",lwd=2,xpd=T)
        abline(v=innerregion[2]-chrompos[1],col="red",lwd=2,xpd=T)        
        ## POSITIOS FOR THE LINES BELOW
        x.tmp1<- grconvertX(innerregion[1]-chrompos[1],from='user',to='ndc')
        x.tmp2<- grconvertX(innerregion[2]-chrompos[1],from='user',to='ndc')
        y.tmp <- grconvertY(0,from='user',to='ndc')
        
        ###########################################################
        ## 02 PLOT RECRATES
        chr_rec<-read.table(paste0("/mnt/kwiat/data/galton/malariagen/human/reference/mathgen.stats.ox.ac.uk/impute//2013-02-14_ALL_1000G_phase1integrated_v3/genetic_map_chr",as.numeric(chrom),"_combined_b37.txt"),header=T)
        par(mar=c(1,leftmar,1,1))
        outer.region1 <- c(chrompos[1],chrompos[length(chrompos)])
        rec_plot_x <- chr_rec$position[chr_rec$position>=outer.region1[1]&chr_rec$position<=outer.region1[2]]
        rec_plot_y<- chr_rec$COMBINED_rate.cM.Mb.[chr_rec$position>=outer.region1[1]&chr_rec$position<=outer.region1[2]]
        rec_plot_y <- aggregate(rec_plot_y,by=list(round(rec_plot_x,-4)),mean)
        
        plot(rec_plot_y[,1],rec_plot_y[,2],
             xlim=outer.region1,ylim=c(0,30),type="S",
             ylab="",
             xlab="",xaxt="n",
             cex.lab=2,yaxt="n",bty="7",xaxs="i")
        x_at=seq(round(outer.region1[1]-1e6,-6),
                 round(outer.region1[2]+1e6,-6),1e6) ## change for more/less ticks
        #axis(1,at=x_at,labels=x_at/1e6,cex.axis=1.5)
        y_at=c(0,15,30)
        axis(2,at=y_at,labels=y_at,cex.axis=0.75,las=1)
        mtext(2,text="cM/bp",adj=0.5,line=1.5,cex=0.75)
        box()
        abline(v=innerregion[1],col="red",lwd=2,xpd=T)
        abline(v=innerregion[2],col="red",lwd=2,xpd=T)
        ###########################################################
        
        ###########################################################
        ## 03 PLOT GENES
        source("~/R/Copy/Rprojects/AfricaPOPGEN/functions/hitplots.R")
        genes = load.genes()
        outer.region <- innerregion
        chrom2 <- chrom
        genes = genes[ which( genes$chromosome == chrom2 & genes$txEnd >= outer.region[1] & genes$txStart <= outer.region[2] ),, drop = FALSE ]
        par(mar=c(3,leftmar,1,1))
        if(nrow(genes)>0)
        {
            plot.genes(outer.region, genes, 1 )
        } else {
            plot(0,0,type="n",axes=F,xlab="",ylab="")
        }
        ## ADD LINES
        segments(outer.region[1],3,
                 grconvertX(x.tmp1,  from='ndc'),grconvertY(y.tmp, from='ndc'),
                 lty=1,col="red",lwd=2,xpd=T)
        segments(outer.region[1],3,
                 outer.region[1],0, lwd=2 ,lty=1,col="red")
        segments(outer.region[2],3,
                 grconvertX(x.tmp2,  from='ndc'),grconvertY(y.tmp, from='ndc'),
                 lty=1,col="red",lwd=2,xpd=T)
        segments(outer.region[2],3,
                 outer.region[2],0, lwd=2, lty=1,col="red")
        xat <- pretty(innerregion)
        axis(1,at=xat,labels=xat/1e6)
        mtext(1,text=paste("position on chromosome",as.numeric(chrom),"(cM)"),line=2,cex=0.75)
    }



dev.off()