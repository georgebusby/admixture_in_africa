##########################################################
## PLOT FINESTRUCTURE MATRIX
##########################################################

library(ape)
source(paste0(main_dir,"popgen/packages_ext/FinestructureLibrary.R"))
source(paste0(main_dir,"popgen/packages_ext/FinestructureLibrary_GB.R"))

##########################################################
## LOAD DATA
tmpmat <- read.table(mat_file,header=T,row.names=1)
tmpmat <- as.matrix(tmpmat)
tmpmat[tmpmat>tmatmax] <- tmatmax  
tmpmat[tmpmat<tmatmin] <- tmatmin  
## MAKE SOME HEATMAPPY COLOURS
cols <- MakeColorYRP(final=c(0.2,0.2,0.2))
ignorebelow <- tmatmin
colscale<-c(max(ignorebelow,min(tmpmat[tmpmat>ignorebelow])),max(tmpmat))

##########################################################
## PLOT HEATMAP
image(1:dim(tmpmat)[1],
      1:dim(tmpmat)[2],tmpmat,
      xaxt="n",yaxt="n",xlab="",ylab="",
      col=cols,zlim=colscale,xaxs="i",yaxs="i",
      useRaster=T)

##########################################################
## SCALE
if(plot_scale == TRUE)
{
    scalesplit <- 25
    colindex<-t(matrix(seq(min(colscale),max(colscale),length.out=scalesplit),
                       nrow=1,ncol=scalesplit)) # colour scale
    par(mar=c(2,1,1,2))
    image(1:scalesplit,1,(colindex),
          xaxt="n",yaxt="n",xlab="",ylab="",
          col=cols,zlim=range(colindex),main="NUMBER OF CHUNKS COPIED",
          useRaster=T)
    scalelocs <- min(colindex)+(max(colindex)-min(colindex))*seq(0,1,length.out=7)
    scalelocs <- round(scalelocs,0)
    scalelocs[length(scalelocs)] <- paste0(scalelocs[length(scalelocs)],"+")
    scalephysicalpos <- seq(1,scalesplit,length.out=7)
    axis(1,at=scalephysicalpos,labels=scalelocs,cex.axis=1,las=1)
}