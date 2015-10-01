#############
## GB finestructure library
## various functins that help with analysising 
## and plotting finestructure output

copydir <- "~/R/Copy"
source(paste0(copydir,"/Rprojects/GLOBETROTTER/libraries/FinestructureDendrogram.R"))
source(paste0(copydir,"/Rprojects/GLOBETROTTER/libraries/FinestructureLibrary.R"))
## MATHEMATICAL FUNCTIONS
## standard error
se <- function(x){sd(x)/sqrt(length(x)-1)}
## calculate mean of n-1
mean2 <- function(x){(sum(x)/(length(x)-1))}



### COLUMN MANIPULATION FUNCTIONS
### average columns into clusters
## x =  list of clusters in finsestructure format NB - function to get these first?
## y = chromopainter matrix

colsAsClusters <- function(x,y,stat=mean){
  y1 <- c()
  for(i in 1:length(x)){
    clust.name <- unlist(strsplit(x,split="\\(")[[i]][1])
    x1 <- unlist(strsplit(sub("\\)","",strsplit(x,split="\\(")[[i]][2]),split=","))
    x1.pops <- sort(unique(gsub("[0-9]","",x1)))
    x1.length <- length(x1)
    y.clust <- y[,x1]
    if(x1.length>1){
      yNew <- apply(y.clust,1,stat)
    } else {
      yNew <- y.clust
    }
    y1 <- cbind(y1,yNew)
    colnames(y1)[length(y1[1,])] <- clust.name
  }
  return(y1)
}

### rows as clusters
## x = finestructure clusters
## y1 = chromopainter matrix with optional cols as clusters

rowsAsClusters <- function(x,y1,stat=mean){
  yMat <- c()
  for(i in 1:length(x)){
    clust.name <- unlist(strsplit(x,split="\\(")[[i]][1])
    x1 <- unlist(strsplit(sub("\\)","",strsplit(x,split="\\(")[[i]][2]),split=","))
    x1.pops <- sort(unique(gsub("[0-9]","",x1)))
    x1.length <- length(x1)
    y.clust <- y1[x1,]
    if(x1.length>1){
      y2 <- apply(y.clust,2,stat)
    } else {
      y2 <- y.clust
    }
    yMat <- rbind(yMat,y2)
    rownames(yMat)[length(yMat[,1])] <- clust.name
  }
  return(yMat)
}


## as above but uses the fineSTRUCTURE map state format for cluster labels
rowsAsMapClusts <- function(x,y1,stat=mean){# x is clusts as MapState, y1 is matrix
  yMat <- c()
  for(i in 1:length(x)){
    clust.name <- names(x[i])
    x1 <- as.vector(unlist(x[i]))
    x1.length <- length(x1)
    y.clust <- y1[x1,]
    if(x1.length>1){
      y2 <- apply(y.clust,2,stat)
    } else {
      y2 <- y.clust
    }
    yMat <- rbind(yMat,y2)
    rownames(yMat)[length(yMat[,1])] <- clust.name
  }
  return(yMat)
}


## makes a matrix with all inds as recips and superinds only as donors
makeSuperindmatrix <- function(superinds,dataraw){ # superinds is a list of superinds in FS format, dataraw is a full chunkcounts matrix
  superindmatrix <- c()
  for(pop in superinds){
    popname <- strsplit(pop,split="\\(")[[1]][1]
    popinds <- strsplit(pop,split="\\(")[[1]][2]
    popinds <- unlist(strsplit(sub("\\)","",popinds),split=","))
    popcopyvec <- apply(dataraw[,popinds],1,mean)
    superindmatrix <- cbind(superindmatrix,popcopyvec)
    colnames(superindmatrix)[ncol(superindmatrix)] <- popname
  }
  return(superindmatrix)
}



### collapse an ind v ind matrix into pop rows and cols
## x is the matrix
## note: think about whether the rows/cols should be summed or MEANED

popsAsCols <- function(x,stat=mean){
  names <- colnames(x)
  indPops <- gsub("[0-9]","",names)
  pops <- unique(indPops)
  y <- c()
  for(i in pops){
    popName <- indPops %in% grep(i,indPops,value=TRUE)
    print(i)
    x1 <- x[,popName]
    x1[x1==0] <- NA
    if(!is.null(dim(x1))){
    x2 <- apply(x1,1,stat, na.rm=TRUE)
    } else {
    x2 <- x1
    }
    y <- cbind(y,x2)
    colnames(y)[length(y[1,])] <- i
  }
  return(y)
}


popsAsRows <- function(x,stat=mean){
  names <- rownames(x)
  indPops <- sapply(names,popIn)
  pops <- unique(indPops)
  y <- c()
  for(i in pops){
    popName <- indPops %in% grep(i,indPops,value=TRUE)
    x1 <- x[popName, ]
    x1[x1==0] <- NA
    if(!is.null(dim(x1))){
      x2 <- apply(x1,2,stat, na.rm=TRUE)
    } else {
      x2 <- x1
    }
    y <- rbind(y,x2)
    rownames(y)[length(y[,1])] <- i
  }
  return(y)
}



matAsPops <- function(x){
  y <- popsAsCols(x)
  y <- popsAsRows(y)
  return(y)
}


## get orderof individuals from a finestructure analysis
## note that there is no structure within clusters

getTreeOrder <- function(x,tree){#x is a list of FS clusters
  clustTable <- matrix(nrow=length(x),ncol=2)
  cnt <- 1
  for(i in x){
    clustTable[cnt,1] <- unlist(strsplit(i,split="\\("))[1]
    clustTable[cnt,2] <- gsub("\\)","",unlist(strsplit(i,split="\\("))[2])
    cnt<- cnt+1
  }
  indOrder <- c()
  for(i in tree$tip.label){
    indOrder <- paste(indOrder,clustTable[clustTable[,1]==i,2],sep=",")
  }
  indOrder <- unlist(strsplit(sub("\\,","",indOrder),split=","))
  return(indOrder)
}

## function to average chunkvectors across superinds
## clustsVec= postions of clusters in cluster file here superinds are found
## y=chunk matrix
collapseSuperInds <- function(clustsVec,y){
  y1 <- c()
  inds2drop <- c()
  for(i in 1:length(clustsVec)){
    clust.name <- unlist(strsplit(clustsVec,split="\\(")[[i]][1])
    clust.name <- gsub("[0-9]","",clust.name)
    clust.name <- gsub("\\[","",clust.name)
    clust.name <- gsub("\\]","",clust.name)
    x1 <- unlist(strsplit(sub("\\)","",strsplit(clustsVec,split="\\(")[[i]][2]),split=","))
    x1.pops <- sort(unique(gsub("[0-9]","",x1)))
    x1.length <- length(x1)
    y.clust <- y[,x1]
    if(x1.length>1){
      y2 <- apply(y.clust,1,mean)
    } else {
      y2 <- y.clust   
    }
    y1 <- cbind(y1,y2)
    colnames(y1)[ncol(y1)] <- clust.name
    inds2drop<- c(inds2drop,x1)
  }
  y1Dropped <- y[,!colnames(y)%in%inds2drop]
  y1 <- cbind(y1Dropped,y1)
  ## now rows
  yMat <- c()
  inds2drop <- c()
  for(i in 1:length(clustsVec)){
    clust.name <- unlist(strsplit(clustsVec,split="\\(")[[i]][1])
    clust.name <- gsub("[0-9]","",clust.name)
    clust.name <- gsub("\\[","",clust.name)
    clust.name <- gsub("\\]","",clust.name)
    x1 <- unlist(strsplit(sub("\\)","",strsplit(clustsVec,split="\\(")[[i]][2]),split=","))
    x1.pops <- sort(unique(gsub("[0-9]","",x1)))
    x1.length <- length(x1)
    y.clust <- y1[x1,]
    if(x1.length>1){
      y2 <- apply(y.clust,2,mean)
    } else {
      y2 <- y.clust   
    }
    yMat <- rbind(yMat,y2)
    rownames(yMat)[length(yMat[,1])] <- clust.name
    inds2drop <- c(inds2drop,x1)
  }
  yMatDropped <- y1[!rownames(y1)%in%inds2drop,]
  yMat <- rbind(yMatDropped,yMat)
  return(yMat)
}


### TEXT MANIPULATION FUNCTIONS ###
splitVec <- function(x){
  y <- unlist(strsplit(x,split=" "))
  return(y)
}

uniqueClusterNames <- function(allClusts){
  clustVec <- strsplit(allClusts, split=",")
  clustKeep <- c()
  clustNames <- c()
  newTips <- list()
  for(i in 1:length(clustVec)){
    cl <- gsub("[0-9]","",clustVec[i][[1]])
    numPops <- length(unname(table(cl)))
    pops <- names(table(cl))
    nums <- unname(table(cl))
    nameVec <- c()
    for(j in 1:numPops){
      nameVec <- paste(nameVec,paste(pops[j],"[",nums[j],"]",sep=""),sep="_")
    }
    newTips <- c(newTips,nameVec)
    newTips <- sub("_","",newTips)
  }
  for(i in newTips){
    test <- sum(newTips%in%i)
    if(test>1){
      pos2change <- (1:length(newTips))[newTips%in%i]
      for(j in 2:test){
        newTips[pos2change[j]] <- paste(newTips[pos2change[j]],letters[j-1],sep="")
      }
    }
  }
  return(newTips)
}


NameExpandGB<-function(x){
  ret<-list()
  for(i in x){
    tlist<-NULL
    j<-strsplit(i,"_")[[1]]
    for(k in j) {
      tlist<-c(tlist,getIndivsFromSummary(k))
    }
    ret[[i]]<-tlist
  }
  ret
}


## function generates additional files required for plotting
## using finestructure mcmc file
## assumes that files are labelled 
## <maindir><indir><analysis><root>.chunkcounts.out
## <maindir><outdir><analysis><root><replicate>.xml
## <maindir><outdir><analysis><root><replicate><treetype>.tree.xml
generateFSfiles <- function(maindir,analy,replicate,treetype){
  
  ## define files
  chunkfile <- paste0(maindir,"input/",analy,".chunkcounts.out")
  mcmcfile <- paste0(maindir,"output/",analy,replicate,".xml")
  treefile <- paste0(maindir,"output/",analy,replicate,treetype,".tree.xml")
  
  ## generate mapstate file
  mappopchunkfile <- paste0(maindir,"output/",analy,replicate,".mapstate.csv") 
  if(!file.exists(mappopchunkfile)){
    system( paste("finestructure -X -Y -e X2",chunkfile,treefile,mappopchunkfile) )
  }
  ## generate mean coincidence file
  meancoincidencefile <- paste0(dir,"output/",analysis,replicate,".meancoincidence.csv")
  if(!file.exists(meancoincidencefile)){
    system( paste("finestructure -X -Y -e meancoincidence",chunkfile,mcmcfile,meancoincidencefile) )
  }
}

getMapstate <- function(maindir, analy, replicate, treetype){
  ## define files
  mcmcfile <- paste0(maindir,"output/",analy,replicate,".xml")
  treefile <- paste0(maindir,"output/",analy,replicate,treetype,".tree.xml")
  ## read in xml
  mcmcxml <- xmlTreeParse(mcmcfile)
  ## convert this into a data frame
  mcmcdata <- as.data.frame.myres(mcmcxml) 
  ## read the tree as xml format
  treexml <- xmlTreeParse(treefile)
  ## extract the tree into ape's phylo format
  ttree <- extractTree(treexml) 
  ## convert to dendrogram format
  tdend <- myapetodend(ttree,factor=1, tol=1) 
  ## Now we work on the MAP state
  ## map state as a finestructure clustering
  mapstate <- extractValue(treexml,"Pop")
  ## .. and as a list of individuals in populations
  mapstatelist <- popAsList(mapstate)
  ## population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
  #popnames <- lapply(mapstatelist,clusterNames) 
  #names(popnames)<-popnames
  mapstatelist <- makeUniqueClusterNames(mapstatelist)
  
  #names(mapstatelist) <- popnames
  mapstatelist
}

getDend <- function(maindir, analy, replicate, treetype){
  ## define files
  mcmcfile <- paste0(maindir,"output/",analy,replicate,".xml")
  treefile <- paste0(maindir,"output/",analy,replicate,treetype,".tree.xml")
  ## read in xml
  mcmcxml <- xmlTreeParse(mcmcfile)
  ## convert this into a data frame
  mcmcdata <- as.data.frame.myres(mcmcxml) 
  ## read the tree as xml format
  treexml <- xmlTreeParse(treefile)
  ## extract the tree into ape's phylo format
  ttree <- extractTree(treexml) 
  ## convert to dendrogram format
  tdend <- myapetodend(ttree,factor=1, tol=1) 
  tdend
}



#########################################################################
### ADMIXTURE FUNCTIONS ###
## get info from Myers results.fit object output by ADMIXprops

getInfo <- function(r,al,va){
  ## function gets info from resuls database
  # r is usually results.fit
  # al is usually avrelmse
  # va is the variable we want
  tmp <- sort(unlist(r[names(r)==va][which(unlist(r[names(r)==al])==min(unlist(r[names(r)==al]),na.rm=T))]))
  return(tmp)
}

getInfoGiveAlpha <- function(r,alph,va){
  ## funtion like getInfo but allows results to be got at specific alpha
  # r is usually results.fit
  # alph is level of alpha
  # va is the variable we want
  tmp <- sort(unlist(r[names(r)==va][which(unlist(r[names(r)=="alpha"])==(alph))]))
  return(tmp)
}

## make source vectors from GH Admixture output
makeSrcVec <- function(x,src){
  # x=legnth matrix with rows as clusts
  # src = list of sources
  # returns the source copying vector
  names(src) <- gsub("pop[0-9]fit.","",names(src))
  sourceAll <- rep(0,length(colnames(x)))  
  for(i in 1:length(names(src))){
    tmpSrc <- names(src[i])
    tmpProp <- as.numeric(src[i])
    sourceAll <- sourceAll+(x[tmpSrc,]*tmpProp)
  }
  return(sourceAll)
}


## plot source heatmaps
plotSources <- function(x,scalemax,order.pop,col,colindex){    
  x[x>scalemax] <- scalemax
  x <- x[,order(x[order.pop,])]
  image(1:nrow(x),
        1:ncol(x),
        x,
        xaxt ="n",yaxt="n",
        xlab="", ylab="",zlim=range(colindex),
        col=some.colorsEnd, main="")
  axis(3,at=1:dim(x)[1],
       labels=rownames(x),
       lwd=0,las=2,cex.axis=0.75)
}

## plot sources without labels
plotSources2 <- function(x,scalemax,col,colindex,cors,main.title=""){    
  x[x>scalemax] <- scalemax
  image(1:nrow(x),
        1:ncol(x),
        x,
        xaxt ="n",yaxt="n",
        xlab="", ylab="",zlim=range(colindex),
        col=some.colorsEnd,main=main.title,frame.plot=FALSE)
  
}

### PLOTTING FUNCTIONS ###
## make a color scale coutesy of dan lawson
MakeColorYRP<-function(colby=0.05,final=NULL){
  ## makes yellow/red/purple colour scheme, adding rgb(final) if not null
  tmp<-c(rgb(1,seq(1,colby,-colby),colby),rgb(1,colby,seq(colby,1,colby)),rgb(seq(1,colby,-colby),colby,1.0))
  if(is.null(final)) return(tmp)
  c(tmp,rgb(final[1],final[2],final[3]))
}


## slightly update color.legend function to control border of the legend box
color.legend2=function (xl, yb, xr, yt, legend, rect.col, cex = 1, align = "lt", 
                        gradient = "x", border, ...) {
  oldcex <- par("cex")
  par(xpd = TRUE, cex = cex)
  gradient.rect(xl, yb, xr, yt, col = rect.col, nslices = length(rect.col), 
                gradient = gradient, border = border)
  if (gradient == "x") {
    xsqueeze <- (xr - xl)/(2 * length(rect.col))
    textx <- seq(xl + xsqueeze, xr - xsqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      texty <- yb - 0.2 * strheight("O")
      textadj <- c(0.5, 1)
    }
    else {
      texty <- yt + 0.2 * strheight("O")
      textadj <- c(0.5, 0)
    }
  }
  else {
    ysqueeze <- (yt - yb)/(2 * length(rect.col))
    texty <- seq(yb + ysqueeze, yt - ysqueeze, length.out = length(legend))
    if (match(align, "rb", 0)) {
      textx <- xr + 0.2 * strwidth("O")
      textadj <- c(0, 0.5)
    }
    else {
      textx <- xl - 0.2 * strwidth("O")
      textadj <- c(1, 0.5)
    }
  }
  text(textx, texty, labels = legend, adj = textadj, ...)
  par(xpd = FALSE, cex = oldcex)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################





################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

### FUNCTIONS TO PLOT COANCESTRY MATRICES FROM CHROMOPAINTER ###
## do eigen analysis first
## this way the plot of the eigen object can be manipulated
## risk of signs changing if eigen called within function


finePCA=function(mat2pca, ## covariance matrix to be PCA'd
                 eigY, ## eigen object of cov matrix(mat2plot) 
                 eig1=1, ## eigen vector 1 to plot
                 eig2=2, ## eigen vector 2 to plot
                 x2plot=2, ## eigen on x-axis
                 y2plot=1, ## eigen on y axis
                 revXax=F, ## should x-axis be reversed
                 revYax=T, ## should y-axis be reversed
                 pch.labcex=1.5, ## pch size for plot labels
                 labcex=1.5, ## cex of axis labels
                 plotClusts=F, ## if TRUE, gub rownames for pop-level legend
                 axisDig=3, ## number of signif dig on axis labels
                 mainPlotWdth=0.75, ## width of main plot
                 multiPlot=F, ## is this plot part of many layout plots? if T, then no layout calls
                 plotColTable=NULL, ## table with pop / col / pch columns
                 plotBespokeCols=F, ## if T then use plotColTable for colours and pch
                 plotLegend=T, ## should a legend be plotted?
                 plotColours=T, ## should different coloured points be plotted?
                 plotAxes=T, ## should axes labels be plotted?
                 legCols=1, ## number of cols in legend
                 legCex=1,## cex of legend text
                 pchVec=c(20,15,4,18,17,8,7,6,12,13,16,1,2,3,19,14)){ #vector of psh symbols
  
  
  ################################################################################
  ## prepare colours and labels for plot
  if(plotClusts==T){
    pop.labels = gsub("[0-9]","",rownames(mat2pca))
    pops=sort(unique(pop.labels))
    popColLabs=c()
    count=1
    for(i in 1:length(pops))
    {
      popColLabs=c(popColLabs,i)
      count=count+1
    }
    ## prepare vector of labels for points
    
    popPchLabs=c()
    for(pchType in pchVec)
    {
      popPchLabs=c(popPchLabs,rep(pchType,8))
    }
    pop.cols.labels=cbind(pops,popColLabs,popPchLabs[1:length(pops)])
    
    col.labels=c()
    for(i in 1:length(pop.labels))
    {
      col2add=pop.cols.labels[pop.cols.labels[,1]==pop.labels[i],2]
      col.labels=c(col.labels,as.integer(col2add))
    }
    
    pch.labels=c()
    for(i in 1:length(pop.labels))
    {
      pch2add=pop.cols.labels[pop.cols.labels[,1]==pop.labels[i],3]
      pch.labels=c(pch.labels,as.integer(pch2add))
    } 
  } else {
    pop.labels=rownames(mat2pca)
    pops=unique(pop.labels)
    
    if(plotBespokeCols==T){
      popColLabs=c()
      count=1
      for(i in 1:length(pops))
      {
        popCol=plotColTable[plotColTable[,1]==pops[i],2]
        popCols=rep(popCol,length((1:length(pop.labels))[pop.labels%in%pops[i]]))
        popColLabs=c(popColLabs,popCols)
        count=count+1
      }
      popPchLabs=rep(20,length(popColLabs))
    }
    else{
      popColLabs=c()
      count=1
      for(i in 1:length(pops))
      {
        popCol=rep(i,length((1:length(pop.labels))[pop.labels%in%pops[i]]))
        popColLabs=c(popColLabs,popCol)
        count=count+1
      }
      ## prepare vector of labels for points
      popPchLabs=c()
      count=1
      for(i in 1:length(pops))
      {
        popPch=rep(i,length((1:length(pop.labels))[pop.labels%in%pops[i]]))
        popPchLabs=c(popPchLabs,popPch)
        count=count+1
      }
      #for(pchType in pchVec)
      #{
      #popPchLabs=c(popPchLabs,rep(pchType,8))
      #}
    }
    pop.cols.labels=cbind(pops,unique(popColLabs),unique(popPchLabs))#[1:length(pops)])
    col.labels=popColLabs
    pch.labels=popPchLabs
  }
  ################################################################################
  
  ## get eig values from eigY object
  eigVal1=round((Mod(eigY$values[eig1])/sum(Mod(eigY$values)))*100,digits=axisDig)
  eigVal2=round((Mod(eigY$values[eig2])/sum(Mod(eigY$values)))*100,digits=axisDig)
  
  eigen1=Re(eigY$vectors[,eig1])
  eigen2=Re(eigY$vectors[,eig2])
  
  xlim.plot=c(min(get(paste("eigen",x2plot,sep=""))),max(get(paste("eigen",x2plot,sep=""))))
  ylim.plot=c(min(get(paste("eigen",y2plot,sep=""))),max(get(paste("eigen",y2plot,sep=""))))
  
  
  ## prepare axes to plot ##
  if(revXax==TRUE){xlim.plot=rev(xlim.plot)}
  if(revYax==TRUE){ylim.plot=rev(ylim.plot)}
  
  if(x2plot==1){
    xlabel=paste("PC",eig1," (",eigVal1,"%)",sep="")
    ylabel=paste("PC",eig2," (",eigVal2,"%)",sep="")
    xText=eigen1
    yText=eigen2
  }else{
    xlabel=paste("PC",eig2," (",eigVal2,"%)",sep="")
    ylabel=paste("PC",eig1," (",eigVal1,"%)",sep="")
    xText=eigen2
    yText=eigen1
  }
  
  ## sort out layout for two plots ##
  if(plotLegend==T&multiPlot==F){layout(matrix(c(1,2),1,2,byrow=T),widths=c(mainPlotWdth,(1-mainPlotWdth)))}
  
  ## make empty plot 
  par(mar=c(0,0,0,0))
  if(plotAxes==F){
    xlabel=""
    ylabel=""
  }
  plot(0,0,type='n',main="",xlab=xlabel,ylab=ylabel,cex.lab=labcex,
       xlim=xlim.plot,ylim=ylim.plot,
       frame.plot=FALSE,axes=FALSE)
  ## add points
  if(plotColours==F)
  {
    col.labels="grey"
    pch.labels=20
  }
  points(xText,yText,col=col.labels,cex=pch.labcex, pch=pch.labels)
  
  if(plotAxes==T){
    abline(h=0)
    abline(v=0)
  }
  if(plotLegend==T){
    ## have separate plot for legend ##
    par(mar=c(0,0,0,0))
    plot(0,0,type="n", frame.plot=F,axes=F)
    legend("left", legend=pop.cols.labels[,1],
           cex=legCex,ncol=legCols,
           col=pop.cols.labels[,2],
           pch=as.integer(pop.cols.labels[,3]), bty="n")
  }
}

#############

finePCAsources=function(mat2pca, ## covariance matrix to be PCA'd
                        eigY, ## eigen object of cov matrix(mat2plot) 
                        eig1=1, ## eigen vector 1 to plot
                        eig2=2, ## eigen vector 2 to plot
                        x2plot=2, ## eigen on x-axis
                        y2plot=1, ## eigen on y axis
                        revXax=F, ## should x-axis be reversed
                        revYax=T, ## should y-axis be reversed
                        pch.labcex=1.5, ## pch size for plot labels
                        plotClusts=F, ## if TRUE, gub rownames for pop-level legend
                        axisDig=3, ## number of signif dig on axis labels
                        mainPlotWdth=0.75, ## width of main plot
                        multiPlot=F, ## is this plot part of many layout plots? if T, then no layout calls
                        plotLegend=F, ## should a legend be plotted?
                        plotKerns=F, ## should ellipses be plot
                        plotColours=F, ## should different coloured points be plotted?
                        plotCrosses=F, ## should cross lines at 0 be plotted?
                        margins=c(4,4,2,2),
                        plotAxes=T, ## should axes labels be plotted?
                        xAxPos=1,
                        yAxPos=2,
                        plotFrame=F,
                        giveAxLims=NULL,
                        xlim.plot=NULL,
                        ylim.plot=NULL,
                        plotSourcePoints=F, ## should the source points be plotted using the vectors below
                        pointColours=NULL, ## vector of colours for points
                        pointREG=NULL, ## if plotKerns==T; this is the factors in the plot
                        pointPCH=NULL, ## vector of pch for points
                        pointCEX=NULL, ## vector of point sizes for points
                        pointRIM=NULL,
                        pointLWD=1,
                        pointTEXT=NULL, ## vector of text to put in points
                        pointTEXTcol="black", ## colour of text in points
                        pointTEXTcex=1, ## size of text in points
                        ellipLabels=F, ## should labels be plotted of the ellips?
                        ellipColsTab=NULL, ## table of colours and regions
                        cPoint=0, ## if 0 then no points drawn, if greater, then pch of points
                        legendColTable=NULL, ## optional table for legend
                        plotLegSamePlot=F, ## optional plot legend on same plot?
                        legCols=1, ## number of cols in legend
                        legCex=1,## cex of legend text
                        pchVec=c(20,15,4,18,17,8,7,6,12,13,16,1,2,3,19,14), #vector of psh symbols
                        clusterlist=run1mapstate, ## list of clusters for kern plotting
                        vertHR=80, ## home range plot size
                        si=si, ## list of superinds
                        euroReg=clustercols[,c(5,1)], ## list of colours for each world rgion
                        indnames=rownames(clengths3), ## list of ind names to reference eigenY
                        verEXTENT=2)## the Kern plottinh
                        {
  
  
  ## prepare colours and labels for plot
  if(plotClusts==T){
    pop.labels = gsub("[0-9]","",rownames(mat2pca))
    pops=sort(unique(pop.labels))
    popColLabs=c()
    count=1
    for(i in 1:length(pops)){
      popColLabs=c(popColLabs,i)
      count=count+1
    }
    ## prepare vector of labels for points
    popPchLabs=c()
    for(pchType in pchVec){
      popPchLabs=c(popPchLabs,rep(pchType,8))
    }
    pop.cols.labels=cbind(pops,popColLabs,popPchLabs[1:length(pops)])
    
    col.labels=c()
    for(i in 1:length(pop.labels)){
      col2add=pop.cols.labels[pop.cols.labels[,1]==pop.labels[i],2]
      col.labels=c(col.labels,as.integer(col2add))
    }
    
    pch.labels=c()
    for(i in 1:length(pop.labels)){
      pch2add=pop.cols.labels[pop.cols.labels[,1]==pop.labels[i],3]
      pch.labels=c(pch.labels,as.integer(pch2add))
    } 
  } else {
    pop.labels=rownames(mat2pca)
    pops=unique(pop.labels)
    popColLabs=c()
    count=1
    for(i in 1:length(pops)){
      popCol=rep(i,length((1:length(pop.labels))[pop.labels%in%pops[i]]))
      popColLabs=c(popColLabs,popCol)
      count=count+1
    }
    ## prepare vector of labels for points
    popPchLabs=c()
    count=1
    for(i in 1:length(pops)){
      popPch=rep(i,length((1:length(pop.labels))[pop.labels%in%pops[i]]))
      popPchLabs=c(popPchLabs,popPch)
      count=count+1
    }
    #for(pchType in pchVec)
    #{
    #popPchLabs=c(popPchLabs,rep(pchType,8))
    #}
    pop.cols.labels=cbind(pops,unique(popColLabs),unique(popPchLabs))#[1:length(pops)])
    col.labels=popColLabs
    pch.labels=popPchLabs
  }
  
  ## get eig values from eigY object
  eigVal1=round((Mod(eigY$values[eig1])/sum(Mod(eigY$values)))*100,digits=axisDig)
  eigVal2=round((Mod(eigY$values[eig2])/sum(Mod(eigY$values)))*100,digits=axisDig)
  
  eigen1=Re(eigY$vectors[,eig1])
  eigen2=Re(eigY$vectors[,eig2])
  
  if(is.null(giveAxLims)){
    xlim.plot=c(min(get(paste("eigen",x2plot,sep=""))),max(get(paste("eigen",x2plot,sep=""))))
    ylim.plot=c(min(get(paste("eigen",y2plot,sep=""))),max(get(paste("eigen",y2plot,sep=""))))}
  ## prepare axes to plot ##
  if(revXax==TRUE){xlim.plot=rev(xlim.plot)}
  if(revYax==TRUE){ylim.plot=rev(ylim.plot)}
  
  if(x2plot==1){
    xlabel=paste("PC",eig1," (",eigVal1,"%)",sep="")
    ylabel=paste("PC",eig2," (",eigVal2,"%)",sep="")
    xText=eigen1
    yText=eigen2
  }else{
    xlabel=paste("PC",eig2," (",eigVal2,"%)",sep="")
    ylabel=paste("PC",eig1," (",eigVal1,"%)",sep="")
    xText=eigen2
    yText=eigen1
  }
  
  ## sort out layout for two plots ##
  if(plotLegend==T&multiPlot==F){layout(matrix(c(1,2),1,2,byrow=T),widths=c(mainPlotWdth,(1-mainPlotWdth)))}
  ## make empty plot
  par(mar=margins)
  
  if(plotAxes==F){
    xlabel=""
    ylabel=""
  }
  plot(0,0,type='n',main="",xlab="",ylab="",cex.lab=pch.labcex,
       xlim=xlim.plot,ylim=ylim.plot,
       frame.plot=plotFrame,axes=F)
  if(plotAxes==T){
    axis(xAxPos,cex.lab=pch.labcex)
    mtext(xlabel,side=xAxPos,cex=pch.labcex,line=2.5)
    axis(yAxPos,cex.lab=pch.labcex)
    mtext(ylabel,side=yAxPos,cex=pch.labcex,line=2.5)}
  
  ## add elliptical plots
  if(plotKerns==T){
## get region colours
    for(kernReg in names(clusterlist)){ ## changed from clustNameVec
      print(kernReg)
      ii <- run1mapstate[[kernReg]]
      if(length(ii)>6){
        iireg <- names(si)[lapply(si,function(x){x[x%in%ii[1]]})==ii[1]]
        kernCol <- as.character(euroReg[euroReg[,2]==iireg,1])
       
        xy=t(rbind(xText[1:nrow(mat2pca)],yText[1:nrow(mat2pca)]))
        
        clusts2plotKern <- indnames%in%ii
        rownames(xy) = NULL
        id <- clusts2plotKern
        xy <- data.frame(xy)
        xyid <- xy[id,]  
        
        coordinates(xyid) <- colnames(xyid)[1:2]
        
        ## Run the kernel density analysis, this uses the least-square cross-validation (LSCV) and 80% threshhold
        hr=kernelUD(xyid, h="LSCV",extent=verEXTENT)
        ver=getverticeshr(hr, vertHR,standardise=T)
        
        plot(ver,add=T,xlab="",ylab="",xlim=xlim.plot,ylim=ylim.plot,
             border=kernCol,pbg="transparent",lwd=1,col=kernCol) #"transparent"
      }
    }  
  }
    
    
    
#     colVec=c() 
#     for(i in levels(as.factor(pointREG))){colVec=c(colVec,as.character(ellipColsTab[,2][ellipColsTab[,1]==i]))}
#     
#     eigenY2=cbind(xText,yText)[!is.na(pointColours),]
#     s.class(eigenY2[,c(1,2)],
#             as.factor(pointREG),
#             col=colVec,
#             grid=F,addaxes=F,add.plot=T,
#             xlim=xlim.plot,ylim=ylim.plot,
#             label=ellipLabels,cpoint=cPoint,
#             axesell=F,cstar=0,cellipse=1.5)
#   }
  
  ## add points
  if(plotColours==F){
    col.labels="grey"
    pch.labels=20
  }
  if(plotSourcePoints==T){
    col.labels=pointColours
    pch.labels=pointPCH
    pch.labcex=pointCEX
    rim.labels=pointRIM
  }
  points(xText,yText,col=col.labels,cex=pch.labcex, pch=pch.labels,lwd=pointLWD,bg=rim.labels)
  text(xText,yText,col=pointTEXTcol,labels=pointTEXT,cex=pointTEXTcex)
  
  if(plotCrosses==T) {abline(h=0)}
  if(plotCrosses==T) {abline(v=0)}
  
  if(plotLegend==T){
    ## have separate plot for legend ##
    par(mar=c(0,0,0,0))
    if(is.null(legendColTable)&plotLegSamePlot==F){
      plot(0,0,type="n", frame.plot=F,axes=F)
      legend("left", legend=pop.cols.labels[,1],
             cex=legCex,ncol=legCols,
             col=pop.cols.labels[,2],
             pch=as.integer(pop.cols.labels[,3]), bty="n")}
    else if(plotLegSamePlot==F){
      plot(0,0,type="n", frame.plot=F,axes=F)
      legend("left", legend=legendColTable[,1],
             cex=legCex,ncol=legCols,
             col=legendColTable[,2],
             pch=pch.labels, bty="n")}
    else {
      legend("bottomright",legend=legendColTable[,1],
             cex=legCex,ncol=legCols,
             col=legendColTable[,2],
             pch=pch.labels, bty="n")}
    
  }
} 

###########


elegantTips=function(clustVec){    
  newTips=list()
  for(i in 1:length(clustVec)){
    cl=gsub("[0-9]","",clustVec[i])
    cl=strsplit(cl,split=",")
    numPops=length(unname(table(cl)))
    pops=names(table(cl))
    nums=unname(table(cl))
    nameVec=c()
    for(j in 1:numPops){
      nameVec=paste(nameVec,paste(pops[j],"[",nums[j],"]",sep=""),sep="_")
    }
    newTips=c(newTips,nameVec)
    newTips=sub("_","",newTips)
  }     
  return(newTips)
}


## remove clusters with > 5 inds
#clust2keep=clusterInfo[clusterInfo$n>=5,"clustName"]
## for some reason the cluster names are different so need to be rearranged

reOrderClusterNamesAplh=function(x){
  x2=c()
  for(i in x){
    j=paste(sort(unlist(strsplit(i,split="\\;"))),sep="",collapse="")
    x2=c(x2,j)
  }
  return(x2)
}


## function to group clusters based on the non-European region(s) that they copy most from 
getClusterSourceRegion=function(nonEuroSourceRegionList, ## list of clusters and world regions that they copy most from
                                regCol=1, ## column of nonEuroSourceRegion List with most copied world region
                                regVec=c("Yoruba","Mandenka","NorthAfrica"), ## vector contains regions of interest
                                regName="African", ## name of the major region
                                newColNames=c("cluster","source")){ ## colnames of output
  x=rownames(nonEuroSourceRegionList)[nonEuroSourceRegionList[,regCol]%in%regVec]
  x=gsub("\\.[0-9]","",x)
  x=data.frame(x,regName)
  colnames(x)=newColNames
  return(x)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

### other necessary functionts

colLab <- function(n) {
  if(is.leaf(n)) {
    a <- attributes(n)
    # clusMember - a vector designating leaf grouping
    # labelColors - a vector of colors for the above grouping
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
  }
  n
}

### GB CODE


## makes list out of clusters
makeMapstate=function(clusterlist)
{
  worldreglist=list()
  worldregnames=c()
  for(i in clusterlist)
  {
    temp=popAsList(i)
    worldreglist=c(worldreglist,temp[2])
    worldregnames=c(worldregnames,temp[1])
    #worldregnames=c(worldregnames,as.character(unlist(temp[1])))
  }
  names(worldreglist)=worldregnames
  return(worldreglist)
}

clusterNames <- function(cl) { #renames a list of clusters by the inds within
  #n <- c()
  #for(i in cl){
    tb <- sort(table(gsub("[0-9]","",unlist(cl))),decreasing=T)
    nn <- gsub(" ","",paste(names(tb),tb,collapse="_"))
  #  n <- c(n,nn)
  #}
  nn
}

### function used in findFStrees ans elsewhere to make unique cluster names

makeUniqueClusterNames<-function(x){
  
  clusternames <- unlist(lapply(x,clusterNames))
  names(x) <- clusternames
  numpops <- length(unname(table(names(x))))
  pops <- names(table(names(x)))
  oldnames <- sort(names(x))
  nums <- unname(table(names(x)))
  namevec <- c()
  for(j in 1:numpops){
    if(nums[j]==1){
      namevec <- c(namevec,pops[j])
    } else {
      n<-pops[j]
      for(k in 2:nums[j]){
         if(k<27){
            n<-c(n,paste0(pops[j],letters[k-1]))
         } else {
           n<-c(n,paste0(pops[j],letters[k-26],1))
         }
      }
      namevec <- c(namevec,n)
    }
  }
      
  x <- x[order(names(x))]
  names(x) <- namevec
  
  return(x)
}


findFSsubtrees <- function(tdend, nGroups){
  
  ## now split the tree into nGroups ...
  h <- attr(tdend,"height")
  l <- 1
  x <- nGroups
  while(l!=x){
    dendtmp <- cut(tdend,h=h)
    l <- length(dendtmp$lower)
    h <- max(unlist(lapply(dendtmp$lower,function(x){attr(x,"height")})))-1
   }
  
  clustergroups <- list()
  for(i in 1:x){
    clustergroups[[i]] <- labels(dendtmp$lower[[i]])
  }
  

  clustergroups<- makeUniqueClusterNames(clustergroups)
  
  dendtmp <- makemydend(tdend, clustergroups)
  dendtmp <- fixMidpointsComplete(dendtmp)
  #names(clustergroups) <- lapply(clustergroups,clusterNames)
  
  return(list(dendtmp,clustergroups))
}



## fnction to collapse inds from same node into a
## single lesaf of the tree

collapseTreesClusters=function(ttree,mapstatelist)
{
  
  newTips=c()
  for(i in 1:length(ttree$tip.label)){
    oldTip=ttree$tip.label[i]
    for(j in 1:length(mapstatelist)){
      clust=as.character(unlist(mapstatelist[[j]]))
      if(oldTip %in% clust==TRUE){newTips=c(newTips,j)} 
    }
  }
  
  ttree2=ttree
  ttree2$tip.label=newTips
  tree3=ttree2
  
  ###
  ind=1
  while(ind < length(tree3$tip.label)){
    if(tree3$tip.label[ind]==tree3$tip.label[(ind+1)])
    {tree3=drop.tip(tree3,(ind+1))
     ind=ind}
    else {ind=ind+1}
  }
  
  inds2keep=ttree$tip.label
  
  ## MAKE A VECTOR OF THE FIRST INDIVIDUALS IN EACH CLUSTER ##
  xFirstInd=strsplit(names(mapstatelist),split=",", fixed=TRUE)
  f1=function(x) unlist(strsplit(unlist(x[1]),split=" ")[1])
  xFirstInd=unlist(lapply(xFirstInd,f1))
  
  
  ## FIND THE TIP OF THE TREE THAT CONTAINS THE FIRST IND AND REPLACE WITH CLUSTER INDS ##
  
  test=ttree$tip.label
  
  for(i in 1:length(xFirstInd)){
    indSwap=inds2keep %in% xFirstInd[i]
    if(length(mapstatelist[[i]])>1){
      test[indSwap]=names(mapstatelist[i])
    } else {
      test[indSwap]=paste(names(mapstatelist[i]),
                          names(mapstatelist[i]),sep=",")
    }
  }
  
  
  ## REMOVE TIPS WITH ONE INDIVIDUAL ##
  test2=ttree$tip.label
  vec=c()
  clustVec=c()
  for(i in 1:length(test2)){
    if(length(strsplit(test[[i]],split=",")[[1]])==1){
      vec=c(vec,i)
    } else {
      clustVec=c(clustVec,test[[i]])
    }
  }
  tree4=drop.tip(ttree,ttree$tip.label[vec])
  tree4$tip.label=clustVec
  
  
  ## sort tip labels aphabetically ##
  for(i in 1:length(tree4$tip.label))
  {
    x=strsplit(tree4$tip.label[i], split=",")
    x=mixedsort(unique(x[[1]]))
    z=paste(x, collapse=",")
    tree4$tip.label[i]=z
  }
  
  ## remove NAs brought in during the sorting
  tree4$tip.label=gsub(",NA","",tree4$tip.label)
  tree4$tip.label=gsub("NA,","",tree4$tip.label)
  
  
  ## make elegant tree labels ##
  allClusts=tree4$tip.label
  clustVec=strsplit(allClusts, split=",")
  clustKeep=c()
  clustNames=c()
  
  newTips=list()
  for(i in 1:length(clustVec)){
    cl=gsub("[0-9]","",clustVec[i][[1]])
    numPops=length(unname(table(cl)))
    pops=names(table(cl))
    nums=unname(table(cl))
    nameVec=c()
    for(j in 1:numPops){
      nameVec=paste(nameVec,paste(pops[j],"[",nums[j],"]",sep=""),sep="_")
    }
    newTips=c(newTips,nameVec)
    newTips=sub("_","",newTips)
  }
  
  
  for(i in newTips)
  {
    test=sum(newTips%in%i)
    if(test>1){
      pos2change=(1:length(newTips))[newTips%in%i]
      for(j in 2:test)
      {
        newTips[pos2change[j]]=paste(newTips[pos2change[j]],letters[j-1],sep="")
      }
    }
  }
  
  clustIndsList=paste(newTips,"(",tree4$tip.label,")",sep="")
  
  tree4$tip.label=newTips
  
  return(tree4)
  
}




#################################################################
#################################################################
###
### plot.phylo with extra abilities
## cex.node: changes the node cex
## node.label.offset: node label offset
## font.node: font of the node
## node.label.vert.offset: movement of node label above(+), or below (-) edge

plot.phylo.GB <- function (x, type = "phylogram", use.edge.length = TRUE, node.pos = NULL, 
                           show.tip.label = TRUE, show.node.label = FALSE, edge.color = "black", 
                           edge.width = 1, edge.lty = 1, font = 3, cex = par("cex"), 
                           adj = NULL, srt = 0, no.margin = FALSE, root.edge = FALSE, 
                           label.offset = 0, underscore = FALSE, x.lim = NULL, y.lim = NULL, 
                           direction = "rightwards", lab4ut = "horizontal", tip.color = "black", 
                           plot = TRUE,
                           cex.node=cex, node.label.offset=label.offset,
                           node.label.vert.offset=0, font.node=font, ...) 
{
  Ntip <- length(x$tip.label)
  if (Ntip == 1) {
    warning("found only one tip in the tree")
    return(NULL)
  }
  if (any(tabulate(x$edge[, 1]) == 1)) 
    stop("there are single (non-splitting) nodes in your tree; you may need to use collapse.singles()")
  .nodeHeight <- function(Ntip, Nnode, edge, Nedge, yy) .C("node_height", 
                                                           as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
                                                                                                                1]), as.integer(edge[, 2]), as.integer(Nedge), as.double(yy), 
                                                           DUP = FALSE, PACKAGE = "ape")[[6]]
  .nodeDepth <- function(Ntip, Nnode, edge, Nedge) .C("node_depth", 
                                                      as.integer(Ntip), as.integer(Nnode), as.integer(edge[, 
                                                                                                           1]), as.integer(edge[, 2]), as.integer(Nedge), double(Ntip + 
                                                                                                                                                                   Nnode), DUP = FALSE, PACKAGE = "ape")[[6]]
  .nodeDepthEdgelength <- function(Ntip, Nnode, edge, Nedge, 
                                   edge.length) .C("node_depth_edgelength", as.integer(Ntip), 
                                                   as.integer(Nnode), as.integer(edge[, 1]), as.integer(edge[, 
                                                                                                             2]), as.integer(Nedge), as.double(edge.length), double(Ntip + 
                                                                                                                                                                      Nnode), DUP = FALSE, PACKAGE = "ape")[[7]]
  Nedge <- dim(x$edge)[1]
  Nnode <- x$Nnode
  ROOT <- Ntip + 1
  type <- match.arg(type, c("phylogram", "cladogram", "fan", 
                            "unrooted", "radial"))
  direction <- match.arg(direction, c("rightwards", "leftwards", 
                                      "upwards", "downwards"))
  if (is.null(x$edge.length)) 
    use.edge.length <- FALSE
  if (type %in% c("unrooted", "radial") || !use.edge.length || 
        is.null(x$root.edge) || !x$root.edge) 
    root.edge <- FALSE
  if (type == "fan" && root.edge) {
    warning("drawing root edge with type = 'fan' is not yet supported")
    root.edge <- FALSE
  }
  phyloORclado <- type %in% c("phylogram", "cladogram")
  horizontal <- direction %in% c("rightwards", "leftwards")
  xe <- x$edge
  if (phyloORclado) {
    phyOrder <- attr(x, "order")
    if (is.null(phyOrder) || phyOrder != "cladewise") {
      x <- reorder(x)
      if (!identical(x$edge, xe)) {
        ereorder <- match(x$edge[, 2], xe[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
      }
    }
    yy <- numeric(Ntip + Nnode)
    TIPS <- x$edge[x$edge[, 2] <= Ntip, 2]
    yy[TIPS] <- 1:Ntip
  }
  z <- reorder(x, order = "pruningwise")
  if (phyloORclado) {
    if (is.null(node.pos)) {
      node.pos <- 1
      if (type == "cladogram" && !use.edge.length) 
        node.pos <- 2
    }
    if (node.pos == 1) 
      yy <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
    else {
      ans <- .C("node_height_clado", as.integer(Ntip), 
                as.integer(Nnode), as.integer(z$edge[, 1]), as.integer(z$edge[, 
                                                                              2]), as.integer(Nedge), double(Ntip + Nnode), 
                as.double(yy), DUP = FALSE, PACKAGE = "ape")
      xx <- ans[[6]] - 1
      yy <- ans[[7]]
    }
    if (!use.edge.length) {
      if (node.pos != 2) 
        xx <- .nodeDepth(Ntip, Nnode, z$edge, Nedge) - 
        1
      xx <- max(xx) - xx
    }
    else {
      xx <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                                 z$edge.length)
    }
  }
  else switch(type, fan = {
    TIPS <- x$edge[which(x$edge[, 2] <= Ntip), 2]
    xx <- seq(0, 2 * pi * (1 - 1/Ntip), 2 * pi/Ntip)
    theta <- double(Ntip)
    theta[TIPS] <- xx
    theta <- c(theta, numeric(Nnode))
    theta <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, theta)
    if (use.edge.length) {
      r <- .nodeDepthEdgelength(Ntip, Nnode, z$edge, Nedge, 
                                z$edge.length)
    }
    else {
      r <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
      r <- 1/r
    }
    xx <- r * cos(theta)
    yy <- r * sin(theta)
  }, unrooted = {
    nb.sp <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
    XY <- if (use.edge.length) 
      unrooted.xy(Ntip, Nnode, z$edge, z$edge.length, nb.sp)
    else unrooted.xy(Ntip, Nnode, z$edge, rep(1, Nedge), 
                     nb.sp)
    xx <- XY$M[, 1] - min(XY$M[, 1])
    yy <- XY$M[, 2] - min(XY$M[, 2])
  }, radial = {
    X <- .nodeDepth(Ntip, Nnode, z$edge, Nedge)
    X[X == 1] <- 0
    X <- 1 - X/Ntip
    yy <- c((1:Ntip) * 2 * pi/Ntip, rep(0, Nnode))
    Y <- .nodeHeight(Ntip, Nnode, z$edge, Nedge, yy)
    xx <- X * cos(Y)
    yy <- X * sin(Y)
  })
  if (phyloORclado) {
    if (!horizontal) {
      tmp <- yy
      yy <- xx
      xx <- tmp - min(tmp) + 1
    }
    if (root.edge) {
      if (direction == "rightwards") 
        xx <- xx + x$root.edge
      if (direction == "upwards") 
        yy <- yy + x$root.edge
    }
  }
  if (no.margin) 
    par(mai = rep(0, 4))
  if (is.null(x.lim)) {
    if (phyloORclado) {
      if (horizontal) {
        x.lim <- c(0, NA)
        pin1 <- par("pin")[1]
        strWi <- strwidth(x$tip.label, "inches")
        xx.tips <- xx[1:Ntip] * 1.04
        alp <- try(uniroot(function(a) max(a * xx.tips + 
                                             strWi) - pin1, c(0, 1e+06))$root, silent = TRUE)
        if (is.character(alp)) 
          tmp <- max(xx.tips) * 1.5
        else {
          tmp <- if (show.tip.label) 
            max(xx.tips + strWi/alp)
          else max(xx.tips)
        }
        x.lim[2] <- tmp
      }
      else x.lim <- c(1, Ntip)
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                        cex)
        x.lim <- c(min(xx) - offset, max(xx) + offset)
      }
      else x.lim <- c(min(xx), max(xx))
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                        cex)
        x.lim <- c(0 - offset, max(xx) + offset)
      }
      else x.lim <- c(0, max(xx))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar(x$tip.label) * 0.03 * cex)
        x.lim <- c(-1 - offset, 1 + offset)
      }
      else x.lim <- c(-1, 1)
    })
  }
  else if (length(x.lim) == 1) {
    x.lim <- c(0, x.lim)
    if (phyloORclado && !horizontal) 
      x.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      x.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                         cex)
    if (type == "radial") 
      x.lim[1] <- if (show.tip.label) 
        -1 - max(nchar(x$tip.label) * 0.03 * cex)
    else -1
  }
  if (phyloORclado && direction == "leftwards") 
    xx <- x.lim[2] - xx
  if (is.null(y.lim)) {
    if (phyloORclado) {
      if (horizontal) 
        y.lim <- c(1, Ntip)
      else {
        y.lim <- c(0, NA)
        pin2 <- par("pin")[2]
        strWi <- strwidth(x$tip.label, "inches")
        yy.tips <- yy[1:Ntip] * 1.04
        alp <- try(uniroot(function(a) max(a * yy.tips + 
                                             strWi) - pin2, c(0, 1e+06))$root, silent = TRUE)
        if (is.character(alp)) 
          tmp <- max(yy.tips) * 1.5
        else {
          tmp <- if (show.tip.label) 
            max(yy.tips + strWi/alp)
          else max(yy.tips)
        }
        y.lim[2] <- tmp
      }
    }
    else switch(type, fan = {
      if (show.tip.label) {
        offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                        cex)
        y.lim <- c(min(yy) - offset, max(yy) + offset)
      }
      else y.lim <- c(min(yy), max(yy))
    }, unrooted = {
      if (show.tip.label) {
        offset <- max(nchar(x$tip.label) * 0.018 * max(yy) * 
                        cex)
        y.lim <- c(0 - offset, max(yy) + offset)
      }
      else y.lim <- c(0, max(yy))
    }, radial = {
      if (show.tip.label) {
        offset <- max(nchar(x$tip.label) * 0.03 * cex)
        y.lim <- c(-1 - offset, 1 + offset)
      }
      else y.lim <- c(-1, 1)
    })
  }
  else if (length(y.lim) == 1) {
    y.lim <- c(0, y.lim)
    if (phyloORclado && horizontal) 
      y.lim[1] <- 1
    if (type %in% c("fan", "unrooted") && show.tip.label) 
      y.lim[1] <- -max(nchar(x$tip.label) * 0.018 * max(yy) * 
                         cex)
    if (type == "radial") 
      y.lim[1] <- if (show.tip.label) 
        -1 - max(nchar(x$tip.label) * 0.018 * max(yy) * 
                   cex)
    else -1
  }
  if (phyloORclado && direction == "downwards") 
    yy <- y.lim[2] - yy
  if (phyloORclado && root.edge) {
    if (direction == "leftwards") 
      x.lim[2] <- x.lim[2] + x$root.edge
    if (direction == "downwards") 
      y.lim[2] <- y.lim[2] + x$root.edge
  }
  asp <- if (type %in% c("fan", "radial", "unrooted")) 
    1
  else NA
  plot(0, type = "n", xlim = x.lim, ylim = y.lim, ann = FALSE, 
       axes = FALSE, asp = asp, ...)
  if (plot) {
    if (is.null(adj)) 
      adj <- if (phyloORclado && direction == "leftwards") 
        1
    else 0
    if (phyloORclado && show.tip.label) {
      MAXSTRING <- max(strwidth(x$tip.label, cex = cex))
      loy <- 0
      if (direction == "rightwards") {
        lox <- label.offset + MAXSTRING * 1.05 * adj
      }
      if (direction == "leftwards") {
        lox <- -label.offset - MAXSTRING * 1.05 * (1 - 
                                                     adj)
      }
      if (!horizontal) {
        psr <- par("usr")
        MAXSTRING <- MAXSTRING * 1.09 * (psr[4] - psr[3])/(psr[2] - 
                                                             psr[1])
        loy <- label.offset + MAXSTRING * 1.05 * adj
        lox <- 0
        srt <- 90 + srt
        if (direction == "downwards") {
          loy <- -loy
          srt <- 180 + srt
        }
      }
    }
    if (type == "phylogram") {
      phylogram.plot(x$edge, Ntip, Nnode, xx, yy, horizontal, 
                     edge.color, edge.width, edge.lty)
    }
    else {
      if (type == "fan") {
        ereorder <- match(z$edge[, 2], x$edge[, 2])
        if (length(edge.color) > 1) {
          edge.color <- rep(edge.color, length.out = Nedge)
          edge.color <- edge.color[ereorder]
        }
        if (length(edge.width) > 1) {
          edge.width <- rep(edge.width, length.out = Nedge)
          edge.width <- edge.width[ereorder]
        }
        if (length(edge.lty) > 1) {
          edge.lty <- rep(edge.lty, length.out = Nedge)
          edge.lty <- edge.lty[ereorder]
        }
        circular.plot(z$edge, Ntip, Nnode, xx, yy, theta, 
                      r, edge.color, edge.width, edge.lty)
      }
      else cladogram.plot(x$edge, xx, yy, edge.color, edge.width, 
                          edge.lty)
    }
    if (root.edge) 
      switch(direction, rightwards = segments(0, yy[ROOT], 
                                              x$root.edge, yy[ROOT]), leftwards = segments(xx[ROOT], 
                                                                                           yy[ROOT], xx[ROOT] + x$root.edge, yy[ROOT]), 
             upwards = segments(xx[ROOT], 0, xx[ROOT], x$root.edge), 
             downwards = segments(xx[ROOT], yy[ROOT], xx[ROOT], 
                                  yy[ROOT] + x$root.edge))
    if (show.tip.label) {
      if (is.expression(x$tip.label)) 
        underscore <- TRUE
      if (!underscore) 
        x$tip.label <- gsub("_", " ", x$tip.label)
      if (phyloORclado) 
        text(xx[1:Ntip] + lox, yy[1:Ntip] + loy, x$tip.label, 
             adj = adj, font = font, srt = srt, cex = cex, 
             col = tip.color)
      if (type == "unrooted") {
        if (lab4ut == "horizontal") {
          y.adj <- x.adj <- numeric(Ntip)
          sel <- abs(XY$axe) > 0.75 * pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            1.05
          sel <- abs(XY$axe) > pi/4 & abs(XY$axe) < 0.75 * 
            pi
          x.adj[sel] <- -strwidth(x$tip.label)[sel] * 
            (2 * abs(XY$axe)[sel]/pi - 0.5)
          sel <- XY$axe > pi/4 & XY$axe < 0.75 * pi
          y.adj[sel] <- strheight(x$tip.label)[sel]/2
          sel <- XY$axe < -pi/4 & XY$axe > -0.75 * pi
          y.adj[sel] <- -strheight(x$tip.label)[sel] * 
            0.75
          text(xx[1:Ntip] + x.adj * cex, yy[1:Ntip] + 
                 y.adj * cex, x$tip.label, adj = c(adj, 0), 
               font = font, srt = srt, cex = cex, col = tip.color)
        }
        else {
          adj <- abs(XY$axe) > pi/2
          srt <- 180 * XY$axe/pi
          srt[adj] <- srt[adj] - 180
          adj <- as.numeric(adj)
          xx.tips <- xx[1:Ntip]
          yy.tips <- yy[1:Ntip]
          if (label.offset) {
            xx.tips <- xx.tips + label.offset * cos(XY$axe)
            yy.tips <- yy.tips + label.offset * sin(XY$axe)
          }
          font <- rep(font, length.out = Ntip)
          tip.color <- rep(tip.color, length.out = Ntip)
          cex <- rep(cex, length.out = Ntip)
          for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                                 cex = cex[i], x$tip.label[i], adj = adj[i], 
                                 font = font[i], srt = srt[i], col = tip.color[i])
        }
      }
      if (type %in% c("fan", "radial")) {
        xx.tips <- xx[1:Ntip]
        yy.tips <- yy[1:Ntip]
        angle <- atan2(yy.tips, xx.tips)
        if (label.offset) {
          xx.tips <- xx.tips + label.offset * cos(angle)
          yy.tips <- yy.tips + label.offset * sin(angle)
        }
        s <- xx.tips < 0
        angle <- angle * 180/pi
        angle[s] <- angle[s] + 180
        adj <- as.numeric(s)
        font <- rep(font, length.out = Ntip)
        tip.color <- rep(tip.color, length.out = Ntip)
        cex <- rep(cex, length.out = Ntip)
        for (i in 1:Ntip) text(xx.tips[i], yy.tips[i], 
                               x$tip.label[i], font = font[i], cex = cex[i], 
                               srt = angle[i], adj = adj[i], col = tip.color[i])
      }
    }
    if (show.node.label) 
      text(xx[ROOT:length(xx)] + node.label.offset, yy[ROOT:length(yy)] + node.label.vert.offset, 
           x$node.label, adj = adj, font = font.node, srt = srt, 
           cex = cex.node)
  }
  L <- list(type = type, use.edge.length = use.edge.length, 
            node.pos = node.pos, show.tip.label = show.tip.label, 
            show.node.label = show.node.label, font = font, cex = cex, 
            adj = adj, srt = srt, no.margin = no.margin, label.offset = label.offset, 
            x.lim = x.lim, y.lim = y.lim, direction = direction, 
            tip.color = tip.color, Ntip = Ntip, Nnode = Nnode)
  assign("last_plot.phylo", c(L, list(edge = xe, xx = xx, yy = yy)), 
         envir = .PlotPhyloEnv)
  invisible(L)
}
