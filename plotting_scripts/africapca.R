#############################################################
## PCA - need to run the actual pcas elsewhere, perhaps    ##
##       using snpRelate package
#############################################################
## LOAD FILES
pca_table <- read.table(pca_file)
leginfo <-read.table(leginfo_file,header=T,comment.char="")

###########################################################
## GET EIGEN VECTORS AND EIGENVALS
eigvecs <- pca_table[,5:(ncol(pca_table)-1)]
eigvals <- pca_table[,ncol(pca_table)]

## DEFINE PCH AND COLS
pop_vec <- as.character(pca_table[,2])
pop_pnts <- getPopSymbols(pop_vec,leginfo)

## GET POINTS AND REVERSE AS NECESSARY
xpnts <- eigvecs[,xplot]
ypnts <- eigvecs[,yplot]
## REVERSE AXES IF NECESSARY
if(revX==T) xpnts <- -(xpnts)
if(revY==T) ypnts <- -(ypnts)

## PLOT 
plot(xpnts,ypnts,main="",xlab="",ylab="",axes=T,type="n",xaxt="n",yaxt="n")
abline(v=0,h=0)
points(xpnts,ypnts,bg=pop_pnts$col2plot,
       pch=as.numeric(pop_pnts$pch2plot),
       col=pop_pnts$rim2plot,
       cex=pt_cex,lwd=pt_lwd)
text(x=min(xpnts),y=-(max(ypnts)-min(ypnts))/30,adj=0,cex=1,
     label=paste0("PC",xplot, " (", signif(eigvals[xplot]*100,3),"%)"))
text(x=(max(xpnts)-min(xpnts))/30,y=min(ypnts),adj=0,srt=90,cex=1,
     label=paste0("PC",yplot, " (", signif(eigvals[yplot]*100,3),"%)"))