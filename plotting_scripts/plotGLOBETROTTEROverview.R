### PLOT OVERVIEW OF GLOBETROTTER

#library(extrafont)
#loadfonts()
## to get some nice fonts
#plot_font <-  "Ubuntu"
xAxCex <- 0.75
lab_text_cex <- 0.75

#######
pdf("figures/OverviewOfPAINTINGandGLOBETROTTER.pdf",height=5,width=10)
####### PLOT PAINTED CHROMOSOMES
par(mar=c(0,0.5,2,0.5))
plot(0,0,type="n",axes=F,xlab="",ylab="",
     xlim=c(0,10),ylim=c(0.5,4))

pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
ancreg_list2 <- c("WNC","CNC","ENC","ENS","EAA","SNC","SKS","EUR")
pcoltable <- cbind(ancreg_list,ancreg_list2,pcolshex)

## CHROMOSOME START
left1 <- 0.6
chromlength <- 4.125
chromheight <- 0.15
bottom1 <- 3.5
linelwd <- 1.5


plotChromo <- function(y,col1,pname)
{
    rect(left1,bottom1,
         left1+chromlength,bottom1+chromheight,
         col=col1,border=col1)
    text(left1,bottom1+(0.5*chromheight),labels=pname,pos=2)
}

snpposs <- seq(from=left1+0.025,to=left1+chromlength-0.025,length=20)
chunkposs<- c()
for(i in 1:nrow(pcoltable))
{
    plotChromo(bottom1,pcoltable[i,3],pcoltable[i,2])
    chunkposs <- c(chunkposs,bottom1+(0.5*chromheight))
    for(j in snpposs)
    {
        text(x=j,y=chunkposs[length(chunkposs)],label="0")
    }
    bottom1 <- bottom1 - 0.25
}

bottom1 <- bottom1 - 0.25
plotChromo(bottom1,pcoltable[1,3],"YRI")
## ADD SOME ADMIXTURE CHUNKS
chunks <- snpposs[c(3,4,7,9,10,13,15,16,18,19,20)]

chunks2 <- c(1,3,1,3,6,1,3,4,7,1,8)
chunkcols <- pcoltable[chunks2,3]
chunkpops <- pcoltable[chunks2,1]

x1 <- left1
distbetweensnps <- 0.5*unique(diff(snpposs))[1]
for(i in 1:(length(chunks)-1))
{
    chunkpos <- chunks[i]
    chunkcol <- chunkcols[i]
    chunkpop <- chunkpops[i]
    yline <- chunkposs[ancreg_list==chunkpop]
     ## ADD LINES TO TOP PAINTINGS
    if(i < length(chunks))
    {
        yline1 <- yline
        yline2 <- chunkposs[ancreg_list==chunkpops[(i+1)]]
        segments(chunkpos-distbetweensnps,yline1,
                 chunkpos-distbetweensnps,yline2,lwd=linelwd)    
    }
    ## ADD PAINTED CHUNKS
    x2 <- chunkpos-distbetweensnps
    segments(x1,yline,x2,yline,lwd=linelwd)    
    rect(x1,bottom1,x2,bottom1+chromheight,col=chunkcol,border=chunkcol)
    x1 <- x2
    if(i == length(chunks)-1)
    {
        chunkcol <- chunkcols[length(chunks)]
        chunkpop <- chunkpops[length(chunks)]
        chunkpos <- chunks[length(chunks)-1]
        yline1 <- yline
        yline2 <- chunkposs[ancreg_list==chunkpop]
        segments(chunkpos-distbetweensnps,yline2,
                 left1+(1*chromlength),yline2,lwd=linelwd)
        x2 <- left1+chromlength
        rect(x1,bottom1,x2,bottom1+chromheight,col=chunkcol,border=chunkcol)
    }
}
## put in the final line
abline(h=bottom1+0.5-(0.5*chromheight),lty=3)
## some additional text
text(0,median(chunkposs),labels="PAINTING DONORS",srt=90)
## add some SNPs to recipient
for(j in snpposs)   text(x=j,y=bottom1+(0.5*chromheight),label="0")

###### 2nd panel == copying vectors
## add connector

cvmat <- matrix(0.125,nc=length(ancreg_list),nr=length(ancreg_list))
colnames(cvmat) <- rownames(cvmat) <- ancreg_list2
cvmat[,2] <- 0

cvbeg <- 5
for(i in 1:length(chunkposs))
{
    x1 <- x2 <- cvbeg
    for(j in 1:nrow(cvmat))
    {
        x2 <- x2+cvmat[i,j]
        rect(x1,chunkposs[i]-chromheight*0.5,
             x2,chunkposs[i]+chromheight*0.5,
             col=pcoltable[j,3],border=pcoltable[j,3])
        x1 <- x2
    }
}
text(x=mean(c(cvbeg,x2)),y=chunkposs[1]+1.5*chromheight,label="COPYING\nVECTORS")
#abline(v=x2+0.25)


####### 3rd panel == beta coefficients
cvbeg <- x2+0.5
for(i in 2)
{
    x1 <- x2 <- cvbeg
    for(j in 1:nrow(cvmat))
    {
        x2 <- x2+cvmat[i,j]
        rect(x1,chunkposs[i]-chromheight*0.5,
             x2,chunkposs[i]+chromheight*0.5,
             col=pcoltable[j,3],border=pcoltable[j,3])
        x1 <- x2
    }
}
text(x=mean(c(cvbeg,x2)),y=chunkposs[1]+1.5*chromheight,label="MIXTURE\nMODEL")

## now make yoruba a mixture of different groups

## plot YORUBA
x <- 2
x1 <- cvbeg-1
x2 <- x1+0.2
y1 <- bottom1+3*chromheight
for(j in 1:ncol(cvmat))
{
    y2 <- y1-cvmat[i,j]
    rect(x1,y1,x2,y2,
         col=pcoltable[j,3],border=pcoltable[j,3])
    y1 <- y2
}

##
text(x=x2+0.25,y=bottom1,label="=",cex=2)
x1 <- cvbeg
popcv <- c("WNC","ENC","SKS")
popmix <- c(0.6,0.3,0.1)
cnt<-1
for(x in c(1,3,5))
{
    x2 <- x1+0.2
    y1 <- bottom1+3*chromheight
    for(j in 1:ncol(cvmat))
    {
        y2 <- y1-cvmat[i,j]
        rect(x1,y1,x2,y2,
             col=pcoltable[j,3],border=pcoltable[j,3])
        y1 <- y2
    }
    ppop <- popcv[cnt]
    pmix <- popmix[cnt]
    text(x=x2-0.25,y=bottom1,label=paste(pmix,"x",sep=""),cex=1,adj=1)
    text(mean(c(x1,x2)),y=y2,label=ppop,adj=c(0.5,1))
    x1 <- x2+0.75
    cnt <- cnt+1
}



dev.off()



#########################################################
## 01 PLOT A PAINTED YORUBA CHROMOSOME
library("rhdf5")
## FUNCTIONS
hap2sampleindex <- function(hap,nsamps=10){
    ## finds the first sample index for a haplotype
    sample <- (hap*nsamps)-(nsamps-1)
    return(sample)
}

datafile <- '/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
## CHROMOSOME 2
chrom <- "22"
## GET MAP AND POSITION INFO
map <- data.frame(t(h5read(datafile,paste0("/paintings/chrom",chrom,"/map"))))
H5close()
colnames(map) <- c("position","recrate")
## GET 10X SAMPLES OF A SINGLE PAINTED CHROMOSOME
psamples <- t(h5read(datafile,paste0("/paintings/samples/individuals")))
H5close()
colnames(psamples) <- c("ind","region","X")
psamplesind <- (1:nrow(psamples))[psamples[,"ind"] == "YRI1"]
## 2 haps per sample!!
psampleshap <- hap2sampleindex(psamplesind,2)
psamplesindsamp <- hap2sampleindex(psampleshap)

for(analysis in c("nonlocal"))
{
    tmp <- psamplesindsamp :(psamplesindsamp+9)
    paintedchrom <- h5read(datafile,paste0("/paintings/chrom",chrom,"/",analysis),
                           index=list(tmp,1:nrow(map)))
    H5close()
    ## SWITCH DONORS TO REGIONS
    happops <- c()
    for(i in 1:nrow(psamples)) happops <- c(happops,rep(as.character(psamples[i,"region"]),2))
    #happops <- gsub("SEMI.BANTU","SEMI-BANTU",happops)
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
    chrom2plot <- paintedchromregprop
    
    if(analysis == "local") par(mar=c(0.5,1,1.5,1))
    if(analysis == "nonlocal") par(mar=c(0.5,1,1.5,1))
    barplot(chrom2plot,col=pcolshex,xaxs="i",yaxs="i",
            space=0,axes=F,xaxt="n",border=NA,
            width=chromposI,
            xlab="",horiz=F)
    
    
    image(t(t(paintedchromreg[3,])),col=pcolshex,axes=F)

    image(t(paintedchromreg),col=pcolshex)
    
}
#########################################################

yricv <- table(paintedchromreg)
yricv <- yricv/sum(yricv)

barplot(paintedchromreg,pols)
