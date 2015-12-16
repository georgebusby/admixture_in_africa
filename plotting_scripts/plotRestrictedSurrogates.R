### SCRIPT NO_LOCAL SURROGATE ANALYSIS ###
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
library(xtable)
library("rworldmap")
library("maptools")
#source("functions/convertfactorstostrings.R")
#source("functions/makeTransparent.R")

setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
## GET DATES AS WELL
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"

leginfo <-read.table(leginfo_file,header=T,comment.char="")



library("rworldmap")
library("maptools")
source("functions/convertfactorstostrings.R")
source("functions/makeTransparent.R")
source("functions/plottingFunctions.R")
info3 <- read.table("admixture/MalariaGenAdmixturePopulationKey2.txt",header=T,comment.char="")
sifile <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/forcefiles/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinal.force"
superinds <- scan(sifile,what="char",comment.char="#")

###########################################################
## SET PLOTTING PARAMETERS
plotSep <- FALSE ## plot individually
plotAdmixture <- FALSE ## plot an admixture plot?
plotEquator <- FALSE ## should equator and tropics be plotted
map.xlim.2plot<-c(-20,50) ## X-AXIS LIMITS OF MAP
map.ylim.2plot<-c(-32,28) ## Y-AXIS LIMITS OF MAP
## POSITION OF EACH OF THE LEGENDS
legposx <- c(map.xlim.2plot[1]-5,map.xlim.2plot[2]-5,map.xlim.2plot[1]+8)
legposy <- c(5,15,-5)
pt_cex <- 1.5

## DEFINE GROUPS OF COUNTRIES AS LIST 
afr_cnts <- list(c("Gambia","Mali","BurkinaFaso","Ghana","Nigeria","Cameroon"),
                 c("Sudan","Ethiopia","Somalia","Kenya"),
                 c("Angola","Namibia","Botswana","SouthAfrica","Tanzania","Malawi"))
names(afr_cnts) <- c("WEST/CENTRAL","EASTERN","SOUTHERN")
malgen_cnts <- c("Gambia","Mali","BurkinaFaso","Ghana","Cameroon","Kenya","Tanzania","Malawi")
## PATHS TO FILES WITH POPULATION INFO AND LAT LONGS
leginfo_file <- "admixture/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "admixture/MalariaGenAdmixturePopulationKeyLatLongs.txt"
## controls how pops with same lat/long are plotted
poppos_file <- "admixture/MalariaGenAdmixturePopulationKeyMapPositions.txt" 
pop_nums_file <- "admixture/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"

##########
## I WANT TO PLOT HOW THE MAJOR SOURCE CONTIRBUTION CHANGES WHEN
## GROUPS ARE REMOVED FROM THE ANALYSIS
## NOW GET POINT INFO FROM LEGINFO FILE AND PLOT WITH POSITIONS IN newi

#pop_list <- scan("/mnt/kwiat/home/popgen/scripts/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
pop_list <- scan("~/Backup/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
pop_list <- pop_list[pop_list!="SEMI-BANTU"]
pop_list[pop_list=="SEMI.BANTU"] <- "SEMI-BANTU"
pops2 <- pop_list
pops3 <- pops2
pops3[pops2=="SEMI-BANTU"] <- "SEMI.BANTU"
p <- c()
r <- c()
co <- c()
for(i in pops2){
  ii <- info3[info3$EthnicGroup==tidyNames(i),]
  p <- c(p,ii$poppch)
  r <- c(r,ii$rim)
  co <- c(co,as.character(ii$popcols.info.Country.))
}
pntinfo <- cbind(pops2,co,p,r)
pntinfo <- data.frame(pntinfo)
pntcol <- rep("#000000",length(pops2))
pntcol[pntinfo$r==1]  <- as.character(pntinfo$co)[pntinfo$r==1]

## empty plot
##########
date2plot <- c(-10000,2000)
propslessthan <- 0.05
linewds <- c(0.5,1,2,3,4,5)
minmajor <- "major"  
date2plot <- c(-10000,2000)

##########
## NEED TO PULL IN ADMIXTURESOURECES 3 AND 4 

admixturesources5 <- rbind(admixturesources3,admixturesources4)
srckey <- gsub("\\.[0-9]","",rownames(admixturesources5))
plotruns <- unique(srckey)
altruns <- c("main","nolocal", "nolocalmalawi",
             "nolocalwest","nolocaleast", "nolocalsouth")


plot_mat1 <- plot_mat2 <- matrix(0,nrow=length(plotruns),ncol=60)
rownames(plot_mat1) <- rownames(plot_mat2) <- plotruns
colnames(plot_mat1) <- colnames(plot_mat2) <- colnames(admixturesources5)[21:80]
date_mat <- c()

for(i in plotruns)
{
  recip <- as.character(sapply(i,function(x){strsplit(x,split="\\.")[[1]][1]}))
  arun <- as.character(sapply(i,function(x){strsplit(x,split="\\.")[[1]][2]}))
  if(recip == "SEMI") recip <- "SEMI.BANTU"
  if(arun == "BANTU") arun <- as.character(sapply(i,function(x){strsplit(x,split="\\.")[[1]][3]}))
  doncols <- c("bestmatch.event1.source1","bestmatch.event1.source2",
                "bestmatch.event2.source1","bestmatch.event2.source2",
                "bestmatch.date1.source1","bestmatch.date1.source2",
                "bestmatch.date2.source1","bestmatch.date2.source2")
  donor1 <- admixturesources5[srckey==i,doncols[c(1,3,5,7)]][1,]
  donor2 <- admixturesources5[srckey==i,doncols[c(2,4,6,8)]][1,]
  datecols <- c("gen.1date","gen.2dates.date1","gen.2dates.date2")
  dates <- admixturesources5[srckey==i,datecols][1,]
  propcols <- c("proportion.source1","proportion.event2.source1",
                 "proportion.date1.source1","proportion.date2.source1")
  props <-  admixturesources5[srckey==i,propcols][1,]
  #event <- admixturesources5[srckey==i,"result"][1]
  event <- admixturesources5[srckey==paste(recip,"main",sep="."),"result"][1]
  if(event == "one-date") event <- "1D"
  if(event == "one-date-multiway") event <- "1MW"
  if(event == "multiple-dates") event <- "2D"
  if(event %in% c("1D","1D(2D)"))
  {
      don1 <- donor1[1]
      don2 <- donor2[1]
      donprops1 <- paste(recip,arun,1,sep=".")
      donprops1 <- admixturesources5[donprops1,21:80]
      donprops1 <- donprops1[as.numeric(donprops1)>0]
      #donprops1 <- donprops1/sum(donprops1)
      donprops2 <- paste(recip,arun,2,sep=".")
      donprops2 <- admixturesources5[donprops2,21:80]
      donprops2 <- donprops2[as.numeric(donprops2)>0]
      #donprops2 <- donprops2/sum(donprops2)
      
      date <- makeDate(as.numeric(dates[1]))
      if(gsub("[0-9]","",date) != "")
      {
        date <- -as.numeric(gsub("B","",date))
      } else
      {
        date <- as.numeric(date)
      }
      
      prop1 <- round((as.numeric(props[1])*100)/20)
      if(prop1 > 5) prop1 <- 5
      prop1 <- linewds[prop1]
      
      prop2 <- round(((1-as.numeric(props[1]))*100)/20)
      if(prop2 > 5) prop2 <- 5
      prop2 <- linewds[prop2]
            
      x1 <- recipcolumn
      y1 <- (1:length(pops3))[pops3==recip]
      
      bothdons1 <- c(don2)
      bothdonprops1 <- c(as.numeric(donprops2)*(1-as.numeric(props[1])))
      names(bothdonprops1) <- c(names(donprops2))
      bothdons2 <- c(don1)
      bothdonprops2 <- c(as.numeric(donprops1)*(as.numeric(props[1])))
      names(bothdonprops2) <- c(names(donprops1))
      
      plot_mat1[i,names(bothdonprops1)] <- bothdonprops1
      plot_mat2[i,names(bothdonprops2)] <- bothdonprops2
      
      date_mat <- rbind(date_mat,c(i,event,date,0))
      
  } 
  if(event == "1MW")
  {
      don1 <- donor1[1]
      don2 <- donor2[1]
      don3 <- donor1[2]
      don4 <- donor2[2]
      donprops1 <- paste(recip,arun,2,sep=".")
      donprops1 <- admixturesources5[donprops1,21:80]
      donprops1 <- donprops1[as.numeric(donprops1)>0]
      #donprops1 <- donprops1/sum(donprops1)
      donprops2 <- paste(recip,arun,2,sep=".")
      donprops2 <- admixturesources5[donprops2,21:80]
      donprops2 <- donprops2[as.numeric(donprops2)>0]
      #donprops2 <- donprops2/sum(donprops2)
      donprops3 <- paste(recip,arun,3,sep=".")
      donprops3 <- admixturesources5[donprops3,21:80]
      donprops3 <- donprops3[as.numeric(donprops3)>0]
      #donprops3 <- donprops3/sum(donprops3)
      donprops4 <- paste(recip,arun,2,sep=".")
      donprops4 <- admixturesources5[donprops4,21:80]
      donprops4 <- donprops4[as.numeric(donprops4)>0]
      #donprops4 <- donprops4/sum(donprops4)
      date <- makeDate(as.numeric(dates[1]))
      if(gsub("[0-9]","",date) != "")
      {
        date <- -as.numeric(gsub("B","",date))
      } else
      {
        date <- as.numeric(date)
      }
      
      prop1 <- round((as.numeric(props[1])*100)/20)
      if(prop1 > 5) prop1 <- 5
      prop1 <- linewds[prop1]
      
      prop2 <- round(((1-as.numeric(props[1]))*100)/20)
      if(prop2 > 5) prop2 <- 5
      prop2 <- linewds[prop2]
      
      prop3 <- round((as.numeric(props[2])*100)/20)
      if(prop3 > 5) prop3 <- 5
      prop3 <- linewds[prop3]
      
      prop4 <- round(((1-as.numeric(props[2]))*100)/20)
      if(prop4 > 5) prop4 <- 5
      prop4 <- linewds[prop4]
      
      x1 <- recipcolumn
      y1 <- (1:length(pops3))[pops3==recip]
      
      bothdons1 <- c(don2)
      bothdonprops1 <- c(as.numeric(donprops2)*(1-as.numeric(props[1])))
      names(bothdonprops1) <- c(names(donprops2))
      bothdons2 <- c(don1)
      bothdonprops2 <- c(as.numeric(donprops1)*(as.numeric(props[1])))
      names(bothdonprops2) <- c(names(donprops1))
      bothdons3 <- c(don4)
      bothdonprops3 <- c(as.numeric(donprops4)*(1-as.numeric(props[2])))
      names(bothdonprops3) <- c(names(donprops4))
      bothdons4 <- c(don3)
      bothdonprops4 <- c(as.numeric(donprops3)*(as.numeric(props[2])))
      names(bothdonprops4) <- c(names(donprops3))
      
      plot_mat1[i,names(bothdonprops1)] <- bothdonprops1
      plot_mat2[i,names(bothdonprops2)] <- bothdonprops2
      
      date_mat <- rbind(date_mat,c(i,event,date,0))
  } 
  if(event == "2D")
  {
      don1 <- donor1[3]
      don2 <- donor2[3]
      don3 <- donor1[4]
      don4 <- donor2[4]
      donprops1 <- paste(recip,arun,5,sep=".")
      donprops1 <- admixturesources5[donprops1,21:80]
      donprops1 <- donprops1[as.numeric(donprops1)>0]
      #donprops1 <- donprops1/sum(donprops1)
      donprops2 <- paste(recip,arun,6,sep=".")
      donprops2 <- admixturesources5[donprops2,21:80]
      donprops2 <- donprops2[as.numeric(donprops2)>0]
      #donprops2 <- donprops2/sum(donprops2)
      donprops3 <- paste(recip,arun,7,sep=".")
      donprops3 <- admixturesources5[donprops3,21:80]
      donprops3 <- donprops3[as.numeric(donprops3)>0]
     #donprops3 <- donprops3/sum(donprops3)
      donprops4 <- paste(recip,arun,8,sep=".")
      donprops4 <- admixturesources5[donprops4,21:80]
      donprops4 <- donprops4[as.numeric(donprops4)>0]
      #donprops4 <- donprops4/sum(donprops4)
      
      date <- sapply(as.numeric(dates[c(2,3)]),makeDate)
      for(j in 1:2)
      {
        if(gsub("[0-9]","",date[j]) != "")
        {
          date[j] <- -as.numeric(gsub("B","",date[j]))
        } else
        {
          date[j] <- as.numeric(date[j])
        }
      }
      prop1 <- round((as.numeric(props[3])*100)/20)
      if(prop1 > 5) prop1 <- 5
      prop1 <- linewds[prop1]
      
      prop2 <- round(((1-as.numeric(props[3]))*100)/20)
      if(prop2 > 5) prop2 <- 5
      prop2 <- linewds[prop2]
      
      prop3 <- round((as.numeric(props[4])*100)/20)
      if(prop3 > 5) prop3 <- 5
      prop3 <- linewds[prop3]
      
      prop4 <- round(((1-as.numeric(props[4]))*100)/20)
      if(prop4 > 5) prop4 <- 5
      prop4 <- linewds[prop4]
      
      x1 <- recipcolumn
      y1 <- (1:length(pops3))[pops3==recip]
      
     bothdons1 <- c(don2)
     bothdonprops1 <- c(as.numeric(donprops2)*(1-as.numeric(props[3])))
     names(bothdonprops1) <- c(names(donprops2))
     bothdons2 <- c(don1)
     bothdonprops2 <- c(as.numeric(donprops1)*(as.numeric(props[3])))
     names(bothdonprops2) <- c(names(donprops1))
     bothdons3 <- c(don4)
     bothdonprops3 <- c(as.numeric(donprops4)*(1-as.numeric(props[4])))
     names(bothdonprops3) <- c(names(donprops4))
     bothdons4 <- c(don3)
     bothdonprops4 <- c(as.numeric(donprops3)*(as.numeric(props[4])))
     names(bothdonprops4) <- c(names(donprops3))
     
     plot_mat1[i,names(bothdonprops1)] <- bothdonprops1
     plot_mat2[i,names(bothdonprops2)] <- bothdonprops2 
     date_mat <- rbind(date_mat,c(i,event,date[1],date[2]))
  }
 }
  
  
### DEFINE COLOURS ##
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

pop_list2 <- pop_list
pop_list2$V2 <- factor(pop_list2$V2,levels=ancreg_list)
pop_list2 <- pop_list2[order(pop_list2$V2),]
final_pop_order <- as.character(pop_list2$V1[1:48])
final_pop_order1 <- as.character(pop_list2$V1)


plot_mat1 <- plot_mat1[,final_pop_order1]
plot_mat2 <- plot_mat2[,final_pop_order1]

srcreg <- sapply(colnames(plot_mat1),function(x) {
  x <- popkey$RegionM[popkey$Ethnic_Group==x]
  return(x)
})
srcreg <- as.vector(srcreg)
regions <- as.character(ancreg_list)
srcreg <- factor(srcreg,levels=regions)

srccols <- sapply(as.character(srcreg),function(x)
{
  x <- pcolshex[regions==x]
  return(x)
})

srccols <- as.vector(srccols)
srccols[srccols=="black"] <- "grey"

poprowkey <- as.vector(sapply(rownames(plot_mat1),function(x){strsplit(x,split="\\.")[[1]][1]}))
poprunkey <- as.vector(sapply(rownames(plot_mat1),function(x){strsplit(x,split="\\.")[[1]][2]}))
altruns <- c("nolocal","nolocalnoeast","nolocalnosouth","nolocalnowest","nolocalnomalawi")
altruns2 <- c("nolocal","nolocalnosouth","nolocalnowest","nolocalnomalawi")
plot_mat_new <- t(cbind(plot_mat1,0.1,plot_mat2))
#plotpops <- cls[1:48]

plot_labels <- c("Full Analysis",
                 "No local\nregion",
                 "No local,\neast or south",
                 "No local\n or west",
                 "No local\n or Malawi")


pdf(paste("figure/GLOBETROTTERnolocalSourcesFullColour.pdf",sep=""),height=10,width=8)
nplots <- length(altruns2) + 2
layout(matrix(c(1:nplots),1,nplots))
par(mar=c(0,0,2,1))

## plot labels
tmp_mat <- matrix(0,nrow=nrow(plot_mat_new),
                  ncol=length(plotpops))
colnames(tmp_mat) <- rev(plotpops)
bp <- barplot(tmp_mat,
        yaxt="n",xlab="",beside=F,
        cex.main=0.75,border=NA,horiz=T,
        xaxt="n",main="")
axis(2,at = bp, labels=tidyNames(colnames(tmp_mat),fula=T),
     pos=1,las=2,lwd=0)

for(i in c("main",altruns2))
{
 tmp_mat <- matrix(0,nrow=nrow(plot_mat_new),
                   ncol=length(plotpops))
 colnames(tmp_mat) <- rev(plotpops)
 for(j in plotpops) 
 {
   tmpname <- paste(j,i,sep=".")
   if(!tmpname%in%colnames(plot_mat_new) & i == "nolocalnosouth")
   {
     tmpname <- paste(j,"nolocalnoeast",sep=".")
   }
   if(tmpname%in%colnames(plot_mat_new))
   {
    tmp_mat[,j] <- plot_mat_new[,tmpname]  
   }
 }
   
  barplot(tmp_mat,
          col=c(srccols,"white",srccols),
          yaxt="n",xlab="",beside=F,
          cex.main=0.75,border=NA,horiz=T,
          xaxt="n",main="")
 mtext(3,text=plot_labels[c("main",altruns2) == i],line=-1)
}
dev.off()


srccols <- sapply(as.character(srcreg),function(x)
{
  if(length(grep("Niger-Congo",x)>0))
  {
    x <- pcolshex[regions==x]
  } else
  {
    x <- "grey"
  }
  return(x)
})


srccols[colnames(plot_mat1)== "MALAWI"] <- "black"
srccols[colnames(plot_mat1)== "BANTU"] <- "red"
srccols[colnames(plot_mat1)== "SEMI.BANTU"] <- "red"

# srccols <- sapply(as.character(srcreg),function(x)
# {
#     x <- pcolshex[regions==x]
#   return(x)
# })



srccols <- as.vector(srccols)

pdf(paste("figure/GLOBETROTTERnolocalSourcesMalawiCameroon.pdf",sep=""),height=10,width=8)
nplots <- length(altruns2) + 2
layout(matrix(c(1:nplots),1,nplots))
par(mar=c(0,0,2,1))

## plot labels
tmp_mat <- matrix(0,nrow=nrow(plot_mat_new),
                  ncol=length(plotpops))
colnames(tmp_mat) <- rev(plotpops)
bp <- barplot(tmp_mat,
              yaxt="n",xlab="",beside=F,
              cex.main=0.75,border=NA,horiz=T,
              xaxt="n",main="")
axis(2,at = bp, labels=tidyNames(colnames(tmp_mat),fula=T),
     pos=1,las=2,lwd=0)

for(i in c("main",altruns2))
{
  tmp_mat <- matrix(0,nrow=nrow(plot_mat_new),
                    ncol=length(plotpops))
  colnames(tmp_mat) <- rev(plotpops)
  for(j in plotpops) 
  {
    tmpname <- paste(j,i,sep=".")
    if(!tmpname%in%colnames(plot_mat_new) & i == "nolocalnosouth")
    {
      tmpname <- paste(j,"nolocalnoeast",sep=".")
    }
    if(tmpname%in%colnames(plot_mat_new))
    {
      tmp_mat[,j] <- plot_mat_new[,tmpname]  
    }
  }
  
  barplot(tmp_mat,
          col=c(srccols,"white",srccols),
          yaxt="n",xlab="",beside=F,
          cex.main=0.75,border=NA,horiz=T,
          xaxt="n",main="")
  mtext(3,text=plot_labels[c("main",altruns2) == i],line=-1)
}
dev.off()
