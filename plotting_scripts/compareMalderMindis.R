### SCRIPT TO GET MINUIMUM INFERRED GENETIC DISTANCE OVER WHICH TO INFER
### ADMIXTURE USING MALDER/ALDER
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
library(xtable)
setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
malder_out <- paste0("~/R/Copy/Rprojects/AfricaPOPGEN/manuscript/f3tables/AllPopsMalderFinalALLmindis")

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
## LOAD POPORDER FILE ##
popplot <- scan("/mnt/kwiat/home/popgen/scripts/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplotorder <- popplot
## CHOOSE ONLY THE AFRICAN POPS FROM POPPLOT ORDER
popplot <- popplot[c(1:16,18:length(popplot))]
popplotorder <- popplotorder[c(1:16,18:49)]
poplist <- factor(popplotorder,levels=popplotorder)
poplist<- as.matrix(poplist)
n_pops <- nrow(poplist)

### colours and regions
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
ancreg_printlist <- c("West African Niger-Congo",
                      "Central West African Niger-Congo",
                      "East African Niger-Congo",
                      "East African Nilo-Saharan",
                      "East African Afroasiatic",
                      "South African Niger-Congo",
                      "KhoeSan",
                      "Eurasia")

## PULL IN THE COMPUTED MIN DISTANCES
mindist <- read.table("/mnt/kwiat/data/bayes/users/george/popgen/analysis3/alder/MalariaGENallPopsAlderMinDistances.txt", header=F, as.is=T)
colnames(mindist) <- c("target","ref","mindist")

## WORK OUT THE MINIMUM DISTANCES
dists<- mindist$mindist[!is.infinite(mindist$mindist)]
mdist<- max(dists)


pdf("figures/AllPopsALDERmindistComparisons.pdf",height=5,width=10)
popsx <- popplot
popsx <- popsx[popsx!="SEMI.BANTU"]
par(mar=c(8,4,2,1))
plot(0,0, xlim=c(1,length(popsx)),ylim=c(0,2),
     type="n",ylab="minimum distance inferred from the data (cM)",xlab="",
     axes=F)
axis(2,las=2)
abline(h=0.5,lty=2,col="red")
allpopcomps <- c()
for(pop in poplist)
{
    tmp <- mindist[mindist$target==pop&!is.infinite(mindist$mindist),]
    tmppnts <- cbind(pop,popsx,NA)    
    for(j in 1:nrow(tmp))
    {
        tmppnts[tmppnts[,2]==tmp$ref[j],3] <- tmp$mindist[j]
    }
    allpopcomps <- rbind(allpopcomps,tmppnts)
}    
   
allcols <- c()
meds <- c()
popmat <-matrix(0,nc=n_pops,nr=length(popsx))
colnames(popmat) <- poplist
rownames(popmat) <- popsx
for(pop in popsx)
{
    tmp <- as.numeric(allpopcomps[allpopcomps[,2]==pop,3])
    popmat[pop,allpopcomps[allpopcomps[,2]==pop,1]] <- tmp
    tmp[is.na(tmp)] <- 2
    boxreg <- getPopRegion(tidyNames(pop,fula=T),popkey)
    if(pop == "KARRETJIE") boxreg <- "South_Africa_KhoeSan"
    boxcol <- pcolshex[ancreg_list==boxreg]
    allcols <- c(allcols,boxcol)
    meds <- rbind(meds,c(pop,median(tmp)))
    boxplot(tmp,add=T,at=(1:length(popsx))[popsx==pop],axes=F,col=boxcol,pch=20)
}
for(i in 1:length(popsx)) axis(1,at=i,las=2,lwd=0,
                               labels=tidyNames(popsx[i],fula=T,khoesan=T),col.axis=allcols[i])
dev.off()

##############################################################################
## LOOK AT WHETHER PAIRWISE CHUNKLENGTHS CORRELATE WITH LENGTH OF 
## SHARED LD: THE IDEA HERE IS TO SEE IF GROUPS THAT SHARE MORE
## ANCESTRY TEND TO HAVE LONGER SHARED LD
## correlation of popmat v copyingvectors
mat_file <- paste0("data/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalALL2.chunklengths.out")
tmpmat <- read.table(mat_file,header=T,row.names=1)

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

colnames(tmpmat) <- gsub("\\.","\\-",colnames(tmpmat))
rownames(tmpmat) <- gsub("\\.","\\-",rownames(tmpmat))

comat <- t(rowsAsMapClusts(final_clusts2,t(tmpmat),sum))
comat <- rowsAsMapClusts(final_clusts2,comat,mean)
comat <- comat/rowSums(comat)

cvmat <- comat[rownames(popmat),colnames(popmat)]
popmat[is.na(popmat)] <- 2

cormat <- cor(cvmat,popmat)

diag(popmat) <- NA
diag(cvmat) <- NA

pdf("figures/AllPopsALDERmindistComparisonsVesusCP.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
plot(c(cvmat),c(popmat),pch=20,axes=F,
     xlab="pairwise copying proportions (CHROMOPAINTER)",
     ylab="pairwise LD correlation (cM)")
axis(1)
axis(2,las=2)
box()
tt <- cor.test(c(cvmat),c(popmat))
legend("bottomright",legend=expression(italic(R)^2 == 0.53 ~ italic(P)<"0.00001"),bty="n")
dev.off()

#image(cormat[,rev(colnames(cormat))],na.rm=T,axes=F)

####################################################################
## REPEAT THE ABOVE BOXPLOTS, BUT SPLIT BY REGION
pdf("figures/AllPopsALDERmindistComparisonsSplitByReg.pdf",height=16,width=10)
layout(matrix(c(1:8)))
popsx <- popplot
popsx <- popsx[popsx!="SEMI.BANTU"]

for(reg in ancreg_list[ancreg_list!="Eurasia"])
{
    par(mar=c(0,4,0,1))
    plot(0,0, xlim=c(1,length(popsx)),ylim=c(0,2),
         type="n",ylab="min dist (cM)",xlab="",
         axes=F)
    axis(2,las=2)
    abline(h=0.5,lty=2,col="red")
    allpopcomps <- c()
    for(pop in poplist)
    {
        tmp <- mindist[mindist$target==pop&!is.infinite(mindist$mindist),]
        tmppnts <- cbind(pop,popsx,NA)    
        for(j in 1:nrow(tmp))
        {
            tmppnts[tmppnts[,2]==tmp$ref[j],3] <- tmp$mindist[j]
        }
        allpopcomps <- rbind(allpopcomps,tmppnts)
    }    
    
    popmat <-matrix(0,nc=n_pops,nr=length(popsx))
    colnames(popmat) <- poplist
    rownames(popmat) <- popsx
    for(pop in popsx)
    {
        tmp <- as.numeric(allpopcomps[allpopcomps[,2]==pop,3])
        popmat[pop,allpopcomps[allpopcomps[,2]==pop,1]] <- tmp
        tmppops <- allpopcomps[allpopcomps[,2]==pop,1]
        tmppops <- tmppops[!is.na(tmp)]
        tmp[is.na(tmp)] <- 2
        tmpreg <- sapply(tmppops, function(x)getPopRegion(tidyNames(x,fula=T),popkey))
        tmpreg[tmppops=="KARRETJIE"] <- "South_Africa_KhoeSan"
        #for(reg in unique(tmpreg))
        #{
            boxcol <- pcolshex[ancreg_list==reg] 
            bat <- (1:length(popsx))[popsx==pop]
            boxplot(tmp[tmpreg==reg],add=T,at=bat,axes=F,col=boxcol,pch=20)    
        #}    
    }
    legend("topright",legend= ancreg_printlist[ancreg_list==reg],bty="n")
}

par(mar=c(12,4,0,1))
plot(0,0, xlim=c(1,length(popsx)),ylim=c(0,0),
     type="n",ylab="",xlab="",
     axes=F)
for(i in 1:length(popsx)) axis(1,at=i,las=2,lwd=0,labels=tidyNames(popsx[i], fula=T, khoesan=T),
                               col.axis=allcols[i],line=-2)
dev.off()

#######################################################################
## AND THE OPPOSITE -- IE, WHAT IS THE AVERAGE LD FOR A GIVEN POP
pdf("figures/AllPopsALDERmindistComparisonsSplitByPop.pdf",height=5,width=8)
popsx <- popplot
popsx <- popsx[popsx!="SEMI.BANTU"]
par(mar=c(8,4,2,1))
plot(0,0, xlim=c(1,length(poplist)),ylim=c(0,2),
     type="n",ylab="minimum distance inferred from the data (cM)",xlab="",
     axes=F)
axis(2,las=2)
abline(h=0.5,lty=2,col="red")

allcols <- c()
meds <- c()
popmat <-matrix(0,nc=n_pops,nr=length(popsx))
colnames(popmat) <- poplist
rownames(popmat) <- popsx
for(pop in poplist)
{
    tmp <- as.numeric(allpopcomps[allpopcomps[,1]==pop,3])
    popmat[pop,allpopcomps[allpopcomps[,1]==pop,1]] <- tmp
    tmp[is.na(tmp)] <- 2
    boxreg <- getPopRegion(tidyNames(pop,fula=T),popkey)
    if(pop == "KARRETJIE") boxreg <- "South_Africa_KhoeSan"
    boxcol <- pcolshex[ancreg_list==boxreg]
    allcols <- c(allcols,boxcol)
    meds <- rbind(meds,c(pop,median(tmp),mean(tmp),sd(tmp),se(tmp)))
    bat <- (1:length(poplist))[poplist==pop]
    boxplot(tmp,add=T,at=bat,axes=F,col=boxcol,pch=20,horizontal=F)
}
for(i in 1:length(poplist))
{
    axis(1,at=i,las=2,lwd=0,
         labels=tidyNames(poplist[i], fula=T,khoesan=T),
         col.axis=allcols[i])
}
dev.off()

############################################################################
## I THINK THAT THE MEDIANS ARE A BETTER CHOICE HERE
## plot means and 2*se
pdf("figures/AllPopsALDERmindistComparisonsSplitByPopMeans.pdf",height=10,width=5)
popsx <- popplot
popsx <- popsx[popsx!="SEMI.BANTU"]
par(mar=c(4,8,2,1))
    plot(0,0, xlim=c(0,2),ylim=c(0,length(poplist)),
         type="n",xlab="minimum distance inferred from the data (cM)",ylab="",
         axes=F)
    axis(1,las=1)
    abline(v=0.5,lty=2,col="red")
    
    for(i in length(poplist):1)
    {
        axis(2,at=i,las=2,lwd=0,
             labels=tidyNames(rev(poplist)[i], fula=T,khoesan=T),
             col.axis=rev(allcols)[i])
    }
    
    for(i in 1:nrow(meds))
    {
        x0 <- as.numeric(meds[i,3])-2*as.numeric(meds[i,5])
        x1 <- as.numeric(meds[i,3])+2*as.numeric(meds[i,5])
        y <- rev(1:nrow(meds))[i]
        arrows(x0,y,x1,y,angle=90, code=3,length=0.05)
    }
    points(meds[,3],n_pops:1,pch=21,bg=allcols)
dev.off()




