#########################################################################
### PLOT OF TVD v F_st FOR PAPER 
## NOTE THAT I FIRST USED SMARTPCA TO GENERATE F_st OUTPUT
setwd("~/repos/popgen/")
############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
source(paste0(main_dir,"popgen/packages_ext/FinestructureLibrary.R"))
setwd("~/repos/popgen/")
###########################################################
## DEFINE DATAFILES
fst_file <- "data/AllPopsPairwiseFstResults.txt"
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
cv_file <- "data/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalALL2.chunklengths.out"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)

###########################################################
## GET FULL LIST OF POPULATIONS IN THE CORRECT PLOTTING ORDER ##
popplot <- scan("data/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplot <- popplot[popplot!="SEMI.BANTU"]
## DEFINE COLOURS FOR REGIONS
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

###########################################################
## LOAD FST DATAFILE
fstdata <- read.table(fst_file,header=F,as.is=T)
colnames(fstdata) <- c("pop1","pop2","fst","se")

## GENERATE A MATRIX OF PAIRWISE FST VALUES
fstmat <- matrix(NA,ncol=length(popplot),nr=length(popplot))
rownames(fstmat) <- popplot
colnames(fstmat) <- rev(popplot)
for(i in 1:nrow(fstdata))
{
    p1 <- fstdata$pop1[i]
    p2 <- fstdata$pop2[i]
    fst <- fstdata$fst[i]
    ser <- fstdata$se[i]
    rownum <- (1:nrow(fstmat))[colnames(fstmat)==p1]
    colnum <- (1:ncol(fstmat))[colnames(fstmat)==p2]
    fstmat[p1,p2] <- fst
    fstmat[p2,p1] <- fst
}

## MAKE SOME HEAT MAP COLOURS
cols<-MakeColorYRP(final=c(0.2,0.2,0.2))
ignorebelow <- 0
ignoreabove <- 0.25
colscale<-c(ignorebelow,ignoreabove)
fstmat[fstmat>ignoreabove] <- ignoreabove

###########################################################
## PLOT FST MATRIX
pdf("figures/AllPopsPairwiseFst.pdf",height=10,width=11)
    layout(matrix(c(1,2),1,2),widths=c(10,1))
    par(mar=c(0.5,6,6,0.5))
    image(1:dim(fstmat)[1],
          1:dim(fstmat)[2],fstmat,
          xaxt="n",yaxt="n",xlab="",ylab="",
          col=cols,zlim=colscale,xaxs="i",yaxs="i",
          useRaster=T)
    for(i in popplot) 
    {
        ii <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
        axpos1 <- (1:length(popplot))[i==popplot]
        axpos2 <- (length(popplot):1)[i==popplot]
        mtext(2,text=tidyNames(i,fula=T,khoesan=T),at=axpos2,col=ii,las=2,line=0.5)
        mtext(3,text=tidyNames(i,fula=T,khoesan=T),at=axpos1,col=ii,las=2,line=0.5)
    }
    ## SCALE
    scalesplit <- 10
    colindex<-t(matrix(seq(min(colscale),max(colscale),length.out=scalesplit),
                       nrow=1,ncol=scalesplit)) # colour scale
    par(mar=c(5,0.5,15,4))
    image(1,1:scalesplit,t(colindex),
          xaxt="n",yaxt="n",xlab="",ylab="",
          col=cols,zlim=range(colindex),main="",
          useRaster=T)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=11)
    scalelocs <- round(scalelocs,3)
    #scalelocs[1] <- paste0("<",scalelocs[1])
    scalelocs[length(scalelocs)] <- paste0(scalelocs[length(scalelocs)],"+")
    scalephysicalpos <- seq(1,scalesplit,length.out=11)
    axis(4,at=scalephysicalpos,labels=scalelocs,cex.axis=1,las=1)
    box()

    mtext(3,text=expression(F[ST]))

dev.off()

###########################################################
### NOW MAKE NICE TABLE FOR SUPPLEWMENT
for(i in 1:nrow(fstdata))
{
    p1 <- fstdata$pop1[i]
    p2 <- fstdata$pop2[i]
    fst <- fstdata$fst[i]
    ser <- fstdata$se[i]
    
    rownum <- (1:nrow(fstmat))[colnames(fstmat)==p1]
    colnum <- (1:ncol(fstmat))[colnames(fstmat)==p2]
    
    if(rownum < colnum) fstmat[p1,p2] <- ser
    if(rownum > colnum) fstmat[p2,p1] <- ser
}

fstmat <- fstmat[,rownames(fstmat)]
rownames(fstmat) <- tidyNames(rownames(fstmat),fula=T)
colnames(fstmat) <- tidyNames(colnames(fstmat),fula=T)

library(xtable)
newtab <- xtable(fstmat[,1:48]*1000)

### print all events
newlines <- c()
newlineord <- as.character(sapply(popplot,function(x){
    ancreg_list[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==x])]}))
for(i in (2:length(newlineord))) if(newlineord[i]!=newlineord[(i-1)]) newlines <- c(newlines,i)

addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(-1)
addtorow$pos[[2]] <- c(newlines-1)

addtorow$command[[1]]  <- "fst results table"
addtorow$command[[2]] <- "\\hline \n"
print(newtab, floating=FALSE,
      tabular.environment="longtable", 
      comment=FALSE,
      include.colnames=TRUE,
      include.rownames=TRUE,
      caption.placement="top",file="figures/AllPopsFstAfrica.tex",
      booktabs=TRUE,add.to.row=addtorow)#,

newtab <- xtable(fstmat[,49:60]*1000)

### print all events
newlines <- c()
newlineord <- as.character(sapply(popplot,function(x){
    ancreg_list[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==x])]}))
for(i in (2:length(newlineord))) if(newlineord[i]!=newlineord[(i-1)]) newlines <- c(newlines,i)

addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(-1)
addtorow$pos[[2]] <- c(newlines-1)

addtorow$command[[1]]  <- "fst results table"
addtorow$command[[2]] <- "\\hline \n"
print(newtab, floating=FALSE,
      tabular.environment="longtable", 
      comment=FALSE,
      include.colnames=TRUE,
      include.rownames=TRUE,
      caption.placement="top",file="figures/AllPopsFstEurasia.tex",
      booktabs=TRUE,add.to.row=addtorow)#,

###########################################################
#### ESIMATES OF FST VALUES FOR PARTICULAR POPS
# popa <- "JUHOAN"
# popa <- "CHONYI"
# tmp <- fstdata[fstdata$pop1==popa | fstdata$pop2 == popa,]
# nopops <- c(asia,americas) #europe
# tmp <- tmp[!tmp$pop1 %in% nopops,]
# tmp <- tmp[!tmp$pop2 %in% nopops,]
# mean(tmp$fst)
# 
# popa <- nopops
# tmp <- fstdata
# tmp <- tmp[tmp$pop1 %in% nopops | tmp$pop2 %in% nopops,]
# notboth <- tmp$pop1 %in% nopops & tmp$pop2 %in% nopops
# tmp <- tmp[!notboth,]
# 
# sd(tmp$fst)
###########################################################
## LOAD COPYING VECTORS
lengths <- read.table(cv_file, header=T,row.names=1)
x <- lengths
colnames(x) <- rownames(x)

## LOAD INDIVIDUAL CLUSTER ASSIGNMENT:: NEEDS SERVER ACCESS
final_clusts <- vector("list",length(popplot))
for(i in 1:length(popplot))
{
    ii <- as.character(popplot[i])
    names(final_clusts)[i] <- ii
    if(ii =="SEMI.BANTU") ii <- "SEMI-BANTU"
    iinds <- scan(paste0("/mnt/kwiat/home/popgen/scripts/finalpoplists/",ii,"finalinds.txt"),what="char")
    final_clusts[[i]] <- iinds
}
paintedinds <- read.table("/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/samplelists/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCP.inds",sep=" ")

## CHANGE ACTUAL SAMPLE IDS TO CP LABELS
final_clusts2<- lapply(final_clusts,function(x){as.character(paintedinds[match(unlist(x),paintedinds[,2]),1])})

## NOW GENERATE POPULATION AVERAGE COPYING VECTORS ##
x <- rowsAsMapClusts(final_clusts2,x,mean)
x <- t(rowsAsMapClusts(final_clusts2,t(x),sum))
x <- x/rowSums(x)

## ESTIMATE PAIRWISE TVD LIKE LESLIE ET AL 2015; VAN DORP ET AL 2015
tvdall <- tvd(x)

## MAKE A SIMILAR MATRIX TO FST
tvdmat <- matrix(NA,ncol=length(popplot),nr=length(popplot))
rownames(tvdmat) <- popplot
colnames(tvdmat) <- rev(popplot)
for(i in 1:nrow(tvdall))
{
    p1 <- as.character(tvdall$pop1[i])
    p2 <- as.character(tvdall$pop2[i])
    tvdi <- as.numeric(as.character(tvdall$tvd[i]))
    rownum <- (1:nrow(tvdmat))[colnames(tvdmat)==p1]
    colnum <- (1:ncol(tvdmat))[colnames(tvdmat)==p2]
    tvdmat[p1,p2] <- tvdi
    tvdmat[p2,p1] <- tvdi
}


ignorebelow <- 0
ignoreabove <- 1
colscale<-c(ignorebelow,ignoreabove)
tvdmat[tvdmat>ignoreabove] <- ignoreabove

pdf("figures/AllPopsPairwiseTVD.pdf",height=10,width=11)
    layout(matrix(c(1,2),1,2),widths=c(10,1))
    par(mar=c(0.5,6,6,0.5))
    image(1:dim(tvdmat)[1],
          1:dim(tvdmat)[2],tvdmat,
          xaxt="n",yaxt="n",xlab="",ylab="",
          col=cols,zlim=colscale,xaxs="i",yaxs="i",
          useRaster=T)
    for(i in popplot) 
    {
        ii <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
        axpos1 <- (1:length(popplot))[i==popplot]
        axpos2 <- (length(popplot):1)[i==popplot]
        mtext(2,text=tidyNames(i,fula=T,khoesan=T),at=axpos2,col=ii,las=2,line=0.5)
        mtext(3,text=tidyNames(i,fula=T,khoesan=T),at=axpos1,col=ii,las=2,line=0.5)
    }
    ## SCALE
    scalesplit <- 10
    colindex<-t(matrix(seq(min(colscale),max(colscale),length.out=scalesplit),
                       nrow=1,ncol=scalesplit)) # colour scale
    par(mar=c(5,0.5,15,4))
    image(1,1:scalesplit,t(colindex),
          xaxt="n",yaxt="n",xlab="",ylab="",
          col=cols,zlim=range(colindex),main="",
          useRaster=T)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=11)
    scalelocs <- round(scalelocs,3)
    scalelocs[length(scalelocs)] <- paste0(scalelocs[length(scalelocs)],"+")
    scalephysicalpos <- seq(1,scalesplit,length.out=11)
    axis(4,at=scalephysicalpos,labels=scalelocs,cex.axis=1,las=1)
    box()
    
    mtext(3,text=expression(TVD))

dev.off()


######################################################################
## COMBINE BOTH OF THE DATAFRAMES TO MAKE A SINGLE MATRIX AND PLOT
## put them both in the same matrix
bothmat <- tvdmat
ignorebelow <- 0
ignoreabove <- 0.75
bothmat[bothmat>ignoreabove] <- ignoreabove

## now re-scale so that the units are same as fstmat (max(fstmat) = 0.25)
bothmat <- bothmat/3
bothmat <- bothmat[rev(rownames(bothmat)),]
repmat <- t(fstmat[rev(rownames(fstmat)),]) #rev(colnames(fstmat))
bothmat[upper.tri(bothmat)] <- repmat[upper.tri(repmat)]
bothmat <- bothmat[rev(rownames(bothmat)),]

ignorebelow <- 0
ignoreabove <- 0.25
colscale<-c(ignorebelow,ignoreabove)

## combine into same dataframe
alldata <- tvdall
alldata$fst <- 0
for(i in 1:nrow(alldata))
{
    p1 <- as.character(alldata$pop1[i])
    p2 <- as.character(alldata$pop2[i])
    fstrow <- (as.character(fstdata$pop1) == p1 & as.character(fstdata$pop2) == p2) | (as.character(fstdata$pop1) == p2 & as.character(fstdata$pop2) == p1)  
    alldata$fst[i] <- fstdata$fst[fstrow]
}

atvd <- as.numeric(as.character(alldata$tvd))
afst <- as.numeric(as.character(alldata$fst))

pdf("figures/AllPopsPairwiseFSTTVD.pdf",height=10,width=12)
    layout(matrix(c(1,2,1,3),byrow=T,2,2),widths=c(10,2),heights=c(8,2))
    par(mar=c(0.5,8,8,0.5))
    image(1:dim(bothmat)[1],
          1:dim(bothmat)[2],bothmat,
          xaxt="n",yaxt="n",xlab="",ylab="",
          col=cols,zlim=colscale,xaxs="i",yaxs="i",
          useRaster=T)
    for(i in popplot) 
    {
        ii <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
        axpos1 <- (1:length(popplot))[i==popplot]
        axpos2 <- (length(popplot):1)[i==popplot]
        mtext(2,text=tidyNames(i,fula=T,khoesan=T),at=axpos2,col=ii,las=2,line=0.5)
        mtext(3,text=tidyNames(i,fula=T,khoesan=T),at=axpos1,col=ii,las=2,line=0.5)
    }
    box()
    mtext(3,text="A",adj=0,line=4,cex=2,xpd=T,at=-5)
    
    ## SCALE
    scalesplit <- 11
    colindex<-t(matrix(seq(min(colscale),max(colscale),length.out=scalesplit),
                       nrow=1,ncol=scalesplit)) # colour scale
    par(mar=c(5,5,15,5))
    image(1,1:scalesplit,t(colindex),
          xaxt="n",yaxt="n",xlab="",ylab="",
          col=cols,zlim=range(colindex),main="",
          useRaster=T)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=11)
    scalelocs <- round(scalelocs,3)
    scalelocs2 <- scalelocs*3
    scalelocs[1] <- paste0(">",scalelocs[1])
    scalelocs2[1] <- paste0(">",scalelocs2[1])
    scalelocs[length(scalelocs)] <- paste0(scalelocs[length(scalelocs)],"+")
    scalelocs2[length(scalelocs2)] <- paste0(scalelocs2[length(scalelocs2)],"+")
    scalephysicalpos <- seq(1,scalesplit,length.out=11)
    axis(4,at=scalephysicalpos,labels=scalelocs,cex.axis=1,las=1,adj=0.5,lwd=0)
    axis(2,at=scalephysicalpos,labels=scalelocs2,cex.axis=1,las=1,adj=0.5,lwd=0)
    box()
    
    title1 <- as.list(expression(italic(TVD), "(lower)"))
    title2 <- as.list(expression(italic(F[ST]), "(upper)"))
    
    mtext(do.call(expression, title1 ),side=3,at=0.2,line = c(1,0))
    mtext(do.call(expression, title2 ),side=3,at=1.8,line = c(1,0))                  
    
    par(mar=c(4,4,0.5,0.5))
    plot(atvd,afst,pch=20,xlab="",ylab="",yaxt="n",xaxt="n")
    yax <- c(0,0.1,0.2,0.3)
    xax <- c(0,0.3,0.6,0.9)
    axis(2,las=2,at=yax,labels=yax)
    axis(1,at=xax,labels=xax)
    ctest <- cor.test(atvd,afst)
    r2 <-round(ctest$estimate,2)
    
    mtext(3,text="B",adj=0,line=0,cex=2,xpd=T,at=-0.5)
    mtext(3,text=expression(" " ~ italic(R^2) == 0.79),line=-2,adj=0,cex=0.75)
    mtext(3,text=expression(" " ~ italic(P) < 0.0001),line=-3,adj=0,cex=0.75)
    
    mtext(expression(italic(TVD)),side=1,line=c(2))
    mtext(expression(italic(F[ST])),side=2,line=c(2))

dev.off()

###########################################################################
## MAKE A NICE TABLE FOR THE SUPPLEMENT
library(xtable)
temptab <- round(tvdmat[,60:13]*1000)
newtab <- xtable(temptab, digits=0)

### print all events
newlines <- c()
newlineord <- as.character(sapply(popplot,function(x){
    ancreg_list[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==x])]}))
for(i in (2:length(newlineord))) if(newlineord[i]!=newlineord[(i-1)]) newlines <- c(newlines,i)

addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(-1)
addtorow$pos[[2]] <- c(newlines-1)

addtorow$command[[1]]  <- "tvd results table"
addtorow$command[[2]] <- "\\hline \n"
print(newtab, floating=FALSE,
      tabular.environment="longtable", 
      comment=FALSE,digits=0,
      include.colnames=TRUE,
      include.rownames=TRUE,
      caption.placement="top",file="figures/AllPopsTVDAfrica.tex",
      booktabs=TRUE,add.to.row=addtorow)#,

temptab <- round(tvdmat[,12:1]*1000)
newtab <- xtable(temptab, digits=0)

### print all events
newlines <- c()
newlineord <- as.character(sapply(popplot,function(x){
    ancreg_list[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==x])]}))
for(i in (2:length(newlineord))) if(newlineord[i]!=newlineord[(i-1)]) newlines <- c(newlines,i)

addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(-1)
addtorow$pos[[2]] <- c(newlines-1)

addtorow$command[[1]]  <- "TVD results table"
addtorow$command[[2]] <- "\\hline \n"
print(newtab, floating=FALSE,
      tabular.environment="longtable", 
      comment=FALSE,
      include.colnames=TRUE,
      include.rownames=TRUE,
      caption.placement="top",file="figures/AllPopsTVDEurasia.tex",
      booktabs=TRUE,add.to.row=addtorow)#,
