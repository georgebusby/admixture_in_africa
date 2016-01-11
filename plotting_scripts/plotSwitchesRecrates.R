## ALL DATA ARE STORED IN AN hdf5 FILE

#datafile <- '/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
datafile <- '/well/malariagen/malariagen/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'

library("rhdf5")
library("bigmemory")


## FIRST PULL IN SWITCH DATA FROM HDF5 FILE
switchdata <- c()


chroms <- c(1:22)
for(chrom in chroms)
{
    if (chrom < 10 ) chrom <- paste0("0",chrom)
    chrom <- as.character(chrom)
    print(paste0("getting switches for chromosome: ", chrom))
    ## GET SWITCH INFORMATION
    switches <- data.frame(t(h5read(datafile,paste0("/lengths/chrom",chrom,"/switches"))),
                           stringsAsFactors = FALSE)
    H5close()
    colnames(switches) <- c("counts","counts1st5")
    ## GET MAP AND POSITION INFO
    map <- data.frame(t(h5read(datafile,paste0("/paintings/chrom",chrom,"/map"))))
    H5close()
    colnames(map) <- c("position","recrate")
    tmpout <- cbind(chrom,map,switches)
    switchdata <- rbind(switchdata,tmpout)    
}    

switchdata <- apply(switchdata,2,as.numeric)
## COMPUTE THE PHYSICAL DISTANCE BETWEEN SNPS
posdiff <- c()
for(chrom in unique(switchdata[,"chrom"]))
{
    tpos <- switchdata[switchdata[,"chrom"]==chrom,"position"]
    tposdiff <- c(diff(tpos),0)
    posdiff <- c(posdiff,tposdiff)
}    

## NOW COPMUTE THE PERBP SWITCH RATE AND COMPARE TO RECRATE
xs <- switchdata[,"recrate"]
ys <- switchdata[,"counts"]/posdiff

############################################################
## PLOT
png("figures/SwitchVersusRecrateComparison.png",height=800,width=800)
par(mar=c(4,4,1,1),cex = 1.2)
## A: recrate comparisons
xLim <- c(0,10e-7)
yLim <- c(0,200)
plot(xs,ys,pch = 21,col = "black", bg= "red",
     xlab="recombination rate (Morgans per bp * 10e-7)",
     ylab="painting switch rate (per bp)",
     axes = FALSE, 
     xlim = xLim, ylim = yLim)

## x-axis
xlimat <- seq(xLim[1],xLim[2],length=5)
axis(1,at = xlimat, labels = xlimat/1e-7)
## y-axis
ylimat <- seq(yLim[1],yLim[2],length=5)
axis(2,at = ylimat, labels = ylimat,las= 2)
box()
dev.off() 

############################################################
## PROPORTION VERSUS LENGTHS
regions <- data.frame(t(h5read(datafile,paste0("/paintings/samples/regions"))))
colnames(regions) <- c("pop","country","region")
samps <- data.frame(t(h5read(datafile,paste0("/paintings/samples/individuals"))))
colnames(samps) <- c("ind","pop","x")


## LET'S LOOK AT EURASIAN IN FULAI ANCESTRY
## NB THIS IS RUN IN A SEPARATE PYTHON PROGRAM

pop  <- 'FULAI'
popfile <- paste0('/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings',pop,'.hdf5')
region <- 'Eurasia'
proplengths <- c()
for(chrom in chroms[c(1,2,4,5,7:20,22)]) ## missing 3, 6, 21
{
    if(chrom < 10) chrom <- paste0('0',chrom)
    print(paste0("getting props data for chromosome: ", chrom))
    temp <- t(h5read(popfile,paste0('/lengths/chrom',chrom,'/',pop,'/',region,'/switches')))
    proplengths <- rbind(proplengths,cbind(chrom,temp))
}

proplengths <- data.frame(proplengths,stringsAsFactors=F)
colnames(proplengths)<- c("chrom","freq","switchlength","reclength")
proplengths <- apply(proplengths,2,as.numeric)


png('figures/FULAIswitchlengthVproportion.png',height=1600,width=1600)
layout(matrix(c(1:4),2,2, byrow=T))
par(mar=c(4,4,1,1), cex=1.5)
plot(proplengths[,"switchlength"],proplengths[,"freq"],
     xlab = expression(italic(x_f)),ylab=expression(italic(f)),
     pch = 21,col = "black", bg= "red",
     main='switch lengths',
     axes=F)
axis(1)
axis(2,las=2)
box()
## FIND AND PLOT LACTASE AND DUFFY REGION
#snps <- read.table("/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/AllPops330Kphased.legend")
snps <- read.table("/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/AllPops330Kphased.legend")
colnames(snps)<- c("chrom","rsid","pos","a0","a1")

lct <- (1:nrow(snps))[snps$chrom==2 & snps$pos>136.5e6 & snps$pos < 136.75e6]
lct <- proplengths[lct,]
for(i in 1:nrow(lct))
{
    points(lct[i,"switchlength"],lct[i,"freq"],
           pch = 21,col = "black", bg= "blue")
}
## DUFFY
dar <- (1:nrow(snps))[snps$chrom==1 & snps$pos>159e6 & snps$pos < 159.5e6]
dar <- proplengths[dar,]
for(i in 1:nrow(dar))
{
    points(dar[i,"switchlength"],dar[i,"freq"],
           pch = 21,col = "black", bg= "green")
}

## RECOMBINATION LENGTHS
plot(proplengths[,"reclength"],proplengths[,"freq"],
     xlab = expression(italic(x_f)),ylab=expression(italic(f)),
     pch = 21,col = "black", bg= "red",
     main='recombination rate lengths',
     axes=F)
axis(1)
axis(2,las=2)
box()
## LACTASE
for(i in 1:nrow(lct))
{
    points(lct[i,"reclength"],lct[i,"freq"],
           pch = 21,col = "black", bg= "blue")
}
## DUFFY
dar <- (1:nrow(snps))[snps$chrom==1 & snps$pos>159e6 & snps$pos < 159.5e6]
dar <- proplengths[dar,]
for(i in 1:nrow(dar))
{
    points(dar[i,"reclength"],dar[i,"freq"],
           pch = 21,col = "black", bg= "green")
}



## RECOMBINATION LENGTHS v SWITCH LENGTHS
plot(proplengths[,"reclength"],proplengths[,"switchlength"],
     xlab = "length based on recombination map",ylab="length based on switches",
     pch = 21,col = "black", bg= "red",
     main='',
     axes=F)
axis(1)
axis(2,las=2)
box()


dev.off()


