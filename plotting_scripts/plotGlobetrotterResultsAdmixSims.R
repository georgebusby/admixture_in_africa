### SCRIPT TO PLOT GLOBETROTTER RESULTS FOR PAPER ###
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
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$Ethnic_Group[popkey$Ethnic_Group=="SEMI-BANTU"] <- "SEMI.BANTU"

### DEFINE SOME PLOTTING VARIABLES ##
dateLabelCex=1
datelines=c(-1000,-500,0,500,1000,1500)
yAxisLim=c(-1000,2000)
pcolshex <- c("#0000CD", "#03B4CC", "#A65628", "#FF7F00", "#984EA3", "#4DAF4A", "#CCCC00")[c(1,2,4,5,3,6,7)]
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
popplot <- scan("/mnt/kwiat/home/popgen/scripts/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplot <- popplot[popplot!="SEMI-BANTU"]
popplotorder <- popplot
popplot <- factor(popplot,levels=popplotorder)

## PULL IN DATA
admixturesources2 <- read.table("data/MalariaGenGlobetrotterSimsAdmixtureSources3.txt",header=T,row.names=1,as.is=T)
final.res2plot <- read.table("data/MalariaGenGlobetrotterSims2plot.txt",header=T,row.names=1,as.is=T)
admixturesources2 <- admixturesources2[,levels(popplot)]
dateboots <- read.table("data/MalariaGenGlobetrotterSimsOneDateBootstraps.txt",header=T,row.names=1,as.is=T)
date2boots <- read.table("data/MalariaGenGlobetrotterSimsTwoDateBootstraps.txt",header=T,row.names=1,as.is=T)

## DEFINE THE REGION AND COLOURS FOR EACH OF THE DIFFERENT SOURCES
srcreg <- getPopRegion(tidyNames(colnames(admixturesources2),fula=T),popkey)
srcreg[tidyNames(colnames(admixturesources2),fula=T)=="KARRETJIE"] <- "South_Africa_KhoeSan"
regions <- as.character(ancreg_list)
srcreg <- factor(srcreg,levels=regions)
srccols <- pcolshex[srcreg]
srccols <- as.vector(srccols)
srccols[srccols=="black"] <- "grey"

###########################################################
## SORT OUT RESULTS TO PLOT
pltable <- final.res2plot
pltable <- pltable[!pltable$Result%in%c("U","NA"),]
pltable$Result[pltable$Result=="1D(2D)"] <- "1D"
pltable$Result[pltable$Result=="1MW(2D)"] <- "1MW"

nores <- popplot[!popplot%in%pltable$Cluster]
noprop <- FALSE ## USE ACTUAL PROPORTIONS IN SOURCE BAR CHARTS?

## TABLE OF SIM DETAILS
src1list <- list(c("JOLA","CEU"),c("JOLA","JUHOAN"),c("JOLA","MALAWI"),
                 c("JUHOAN","CEU"),c("JUHOAN","JOLA"),c("JUHOAN","MALAWI"),
                 c("MALAWI","CEU"),c("MALAWI","JOLA"),c("MALAWI","JUHOAN"))
names(src1list) <- c("JOLACEU","JOLAJUHOAN","JOLAMALAWI",
                     "JUHOANCEU","JUHOANJOLA","JUHOANMALAWI",
                     "MALAWICEU","MALAWIJOLA","MALAWIJUHOAN")

sim_pops <- pops <- as.character(pltable$Cluster)
sim_srcs <- sapply(sim_pops,function(x)strsplit(x,split="[0-9]")[[1]][1])
sim_gens <- sapply(sim_pops,function(x)as.numeric(gsub("[A-Z]","",strsplit(x,split="g")[[1]][1])))
sim_props <- sapply(sim_pops,function(x)as.numeric(gsub("[a-z]","",strsplit(x,split="00")[[1]][2])))
sim_src1 <- sapply(sim_srcs,function(x) as.character(unlist(src1list[x])[1]))
sim_src2 <- sapply(sim_srcs,function(x) as.character(unlist(src1list[x])[2]))

sim_tab <- data.frame(cbind(sim_pops,sim_srcs,sim_src1,sim_src2,sim_gens,sim_props))
sim_tab$sim_src1 <- factor(sim_tab$sim_src1,levels=c("JOLA","MALAWI","JUHOAN"))
sim_tab$sim_src2 <- factor(sim_tab$sim_src2,levels=c("CEU","JOLA","MALAWI","JUHOAN"))
sim_tab <- sim_tab[!duplicated(sim_tab),]
sim_tab <- sim_tab[order(sim_tab$sim_src2,sim_tab$sim_src2),]
sim_pops <- as.character(sim_tab$sim_pops)
## PLOT JUST THE CEU AND JUHOAN MINOR SOURCES SIMS
#sim_pops <- sim_pops[sim_tab$sim_src2%in%c("CEU","JUHOAN")]
plgens <- as.numeric(sapply(pltable$Cluster,function(x)strsplit(gsub("[A-Z]","",x),split="g")[[1]][1]))
over100gens <- pltable$Cluster%in%sim_pops&pltable$Analysis=="ancient" & plgens>100
under100gens <- pltable$Cluster%in%sim_pops&pltable$Analysis=="main" & plgens<=100     
pltable <- rbind(pltable[over100gens,],pltable[under100gens,])

pltable$Result <- "1D"

## I'VE ADDED A FUNCTION TO PULL THESE THINGS OUT
#pltable$Result <- "1D"
tempresults <- addGTresults(pltable,rev_pops="")
all_dates <- tempresults[[1]]
all_plot_mat <- tempresults[[2]]
all_src_mat <- tempresults[[3]]
dboots <- tempresults[[4]]
allsources <- tempresults[[5]]

##################################################################################
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
##################################################################################
## MAKE A MATRIX OF SIMULATED SOURCE COPPYING VECTORS
## FOR THIS WE NEED THE ORIGINAL COPYING VECTORS FOR EACH OF THE
## SIM DONOR POPS
if(!exists("orig_predmat"))
{
    orig_predmat <- read.table("data/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalALL2.chunklengths.out",
                               header=T,row.names=1,as.is=T)
    colnames(orig_predmat) <- gsub("\\.","\\-",colnames(orig_predmat))
    rownames(orig_predmat) <- gsub("\\.","\\-",rownames(orig_predmat))
    
    nlength1 <- t(rowsAsMapClusts(final_clusts2,t(orig_predmat),stat=sum))
    nlength1 <- rowsAsMapClusts(final_clusts2,nlength1,stat=mean)
    predmat <- nlength1/rowSums(nlength1)
    ## get columns of predmat in correct order
    srcorder <- colnames(all_src_mat)
    predmat <- predmat[,srcorder]
    
    sim_plot_mat <- matrix(0,nrow=nrow(all_plot_mat),ncol=length(sim_pops))
    colnames(sim_plot_mat) <- sim_pops
    for(i in sim_pops)
    {
        row <- sim_tab$sim_pops==i
        majorsource <- as.character(sim_tab$sim_src1[row])
        minorsource <- as.character(sim_tab$sim_src2[row])
        srcprop <- as.numeric(as.character(sim_tab$sim_props[row]))/100
        simsourcevec <- c((predmat[minorsource,]*(1-srcprop)),0.1,
                          (predmat[majorsource,]*srcprop))
        sim_plot_mat[,i] <- simsourcevec
    }
}
##################################################################################

## ORDER GROUPS BY THE NORMAL POPKEY ORDER
# roword <- rownames(all_dates)
# roword <- as.character(sapply(roword,function(x){gsub("\\_a","",x)}))
# roword <- factor(roword,levels=popkey$Ethnic_Group)
# roword <- order(roword)
# all_dates <- all_dates[roword,]

## MAKE SOME NICE REGION LABELS
regions2 <- c("Western NC","Central West NC","Southern NC",
              "Eastern NC", "Eastern NS/AA", "Southern KS",
              "Eurasia")
regions2 <- toupper(regions2)
##################################################################################
## GENERATE A MATRIX OF SRC COPYING VECTORS, FOR ALL EVENTS,    
## WE'LL SELECT WHICH TO PLOT LATER 
all_src_mat <- matrix(0,ncol=length(popkey$Ethnic_Group),nrow=0)
all_predmat <- matrix(0,ncol=length(popkey$Ethnic_Group),nrow=length(popkey$Ethnic_Group))
rownames(all_predmat) <- colnames(all_predmat) <- colnames(all_src_mat) <- colnames(admixturesources2)
for(i in 1:ncol(all_plot_mat))
{
    pop <- colnames(all_plot_mat)[i]
    pop1 <- as.character(sim_tab$sim_src1[sim_tab$sim_pops==gsub("\\_a","",pop)])
    reg <- as.character(popkey$RegionM[popkey$Ethnic_Group==gsub("\\.","\\-",pop1)])
    predmat <- read.table(paste0("/mnt/kwiat/data/bayes/users/george/popgen/analysis3/globetrotter/input/",
                                 "MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinal",reg,".chunklengths.out"),
                          header=T,row.names=1)
    predmat <- rowsAsMapClusts(final_clusts2,predmat,mean)
    predmat <- predmat/rowSums(predmat)
    rownames(predmat) <- gsub("\\-","\\.",rownames(predmat))
    s1 <- unlist(all_plot_mat[1:60,i])
    s2 <- unlist(all_plot_mat[62:121,i])
    names(s1) <- names(s2) <- colnames(admixturesources2)
    cv1 <- apply(predmat[names(s1),]*s1,2,sum)
    cv2 <- apply(predmat[names(s2),]*s2,2,sum)
    add <- colnames(all_src_mat)[!colnames(all_src_mat)%in%names(cv1)]
    cv1 <- c(cv1,rep(0,length(add)))
    names(cv1)[names(cv1)==""] <- add
    cv1 <- cv1[colnames(all_src_mat)]
    all_src_mat <- rbind(all_src_mat,cv1)
    rownames(all_src_mat)[nrow(all_src_mat)] <- paste(pop,"minor",sep="_")
    cv2 <- c(cv2,rep(0,length(add)))
    names(cv2)[names(cv2)==""] <- add
    cv2 <- cv2[colnames(all_src_mat)]
    all_src_mat <- rbind(all_src_mat,cv2)
    rownames(all_src_mat)[nrow(all_src_mat)] <- paste(pop,"major",sep="_")
    if(length(strsplit(pop,split="\\_")[[1]]) == 1)
    {
        if(sum(all_predmat[pop1,])==0 )
        {
            tp <- predmat[pop1,]
            add <- colnames(all_predmat)[!colnames(all_predmat)%in%names(tp)]
            tp <- c(tp,rep(0,length(add)))
            names(tp)[names(tp)==""] <- add
            all_predmat[pop1,] <- tp[colnames(all_predmat)]
        }
    }
}

##################################################################################
## PULL IN THE RESULTS OF THE MIXTURE MODEL (RUN ELSEWHERE)
## HAVEN'T DONE THIS FOR THE SIMS
mix_dir <- "/mnt/kwiat/well/human/george/admix_sims/globetrotter/props/"
# GENERATE A MIXTURE MATRIX, IF NONE ALREADY
# mixmat <- matrix(0,ncol=length(popplot),nrow=nrow(sim_tab))
# rownames(mixmat) <- levels(sim_tab$sim_pops)
# colnames(mixmat) <- popplot
# for(i in levels(sim_tab$sim_pops))
# {
#     mixture <- read.table(paste0(mix_dir,i,"nolocal.globetrotter.main.nnls.txt"),skip=1,header=T)
#     mixmat[i,match(colnames(mixture),colnames(mixmat))] <- as.numeric(mixture)
# }
# write.table(mixmat,"data/AdmixtureSimsNNLS.txt",col.names=T,row.names=T,quote=F)
mixmat <- read.table("data/AdmixtureSimsNNLS.txt",header=T,row.names=1,as.is=T)

##################################################################################
## PLOT
rev_pops <- ""
pdf("figures/GLOBETROTTERSimresultsAncientMainCombined1D.pdf",height=15,width=10)
    layout(matrix(c(1,1,1,
                    2,2,2,
                    3,3,3,
                    4,4,4,
                    5,5,5,
                    6,7,8),3,6),
           widths=c(4.25,0.75,0.75,0.75,0.75,2.75),heights=c(3,3,3))
    par(mar=c(4,15,4,0.5))
    n_pops <- length(sim_pops)
    d_pops <- gsub("\\_a","",rownames(all_dates))
    
    poplabpos <- 6000
    ev1pos1 <- 5600
    ev1pos2 <- 4600
    ev2pos1 <- 3600
    ev2pos2 <- 2600

    ## EMPTY PLOT FOR DATES
    x_labs3 <- c(2000,0,-2000,-5000,-10000,-15000)
    plot(0,0,xlim=rev(range(x_labs3)),
         ylim=c(n_pops,1),type="n",axes=F,xlab="",ylab="")
    text(poplabpos,-1,labels=LETTERS[1],adj=0,las=1,cex=2,lwd=3,xpd=T)
    x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE","15000\nBCE")
    x_max <- min(x_labs3)
    axis(1,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5)
    for(j in x_labs3) abline(v=j,lty=2)
    mtext("Date of Admixture",1,line=3)
    ## ADD DATE IN GENERATIONS?
    #x_labs4 <- c(100,150,300,400)
    #x_at4 <- sapply(x_labs4,makeDate1)
    #for(i in 1:length(x_labs4))
    #  text(y=-1,x=x_at4[i],labels=paste(x_labs4[i],"gens"),xpd=T,srt=35,adj=0,cex=1)
    #for(j in sapply(x_labs4,makeDate1)) abline(v=j,lty=3)
    
    ## DO WE WANT TO CHANGE THE ORDERING? 
    ## EDIT THIS IN SOME WAY TO CHANGE THE ORDER 
    dpopsorder <- order(factor(d_pops,levels(sim_tab$sim_pops)))
    #dpops <- d_pops[dpopsorder]
    #all_dates <- all_dates[dpopsorder,]
    all_plot_mat <- all_plot_mat[,dpopsorder]
    ###################################################################
    ## PLOT DATES AND EVENT SOURCE ANCESTRY
    ## RUN THE CODE TWICE, THE FIRST TIME GETS THE ORDERING AND THE 
    ## SECOND TIME ACTUALLY PLOTS THE DATA
    for(run in 1:2)
    {
        if(run == 1)
        {
            poporder <- as.character(sort(factor(as.character(sim_pops),levels(sim_tab$sim_pops))))
            popordertab <- c()
        }
        if(run == 2)
        {
            poporder <- as.character(sort(factor(as.character(sim_pops),levels(sim_tab$sim_pops))))
            #poporder <- popordertab[,1]
            #missing_pops <- popplot[1:48]
            #missing_pops <- as.character(missing_pops[!missing_pops%in%poporder])
            #poporder <- c(poporder,missing_pops)
        }
        
        for(i in 1:nrow(all_dates))
        {
            pop <- rownames(all_dates)[i]
            pop1 <- gsub("\\_a","",pop)
            poppos <- (1:length(poporder))[poporder%in%pop1]
            ## plot simulated dates
            popgen <- makeDate(as.numeric(as.character(sim_tab$sim_gens[sim_tab$sim_pops==pop1])),add_BCE=F)
            points(popgen,poppos,pch=15,cex=2,col="red")
            
            res <- all_dates[pop,5]
            d <- as.numeric(all_dates[pop,1])
            dh <- as.numeric(all_dates[pop,3])
            dl <- as.numeric(all_dates[pop,2])
            prop <- as.numeric(all_dates[pop,4])
            if(pop == "ARI") d <- -d ## not sure why this wasn't reversed before
            if(run == 2) points(d,poppos,pch=20, cex=2)
            #if(run == 2) lines(x=c(dl,dh),y=c(poppos,poppos))
            if(run == 2) arrows(dl,poppos,dh,poppos,code=3,length=0.025,angle=90)
           
            
            if(res %in% c("1D","1MW"))
            {
                pop2 <- gsub("\\_a","",pop)
                anc1 <- pltable[pltable$Cluster==pop2,"best.source1"]
                anc2 <- pltable[pltable$Cluster==pop2,"best.source2"]
                anc3 <- pltable[pltable$Cluster==pop2,"best.source1.ev2"]
                anc4 <- pltable[pltable$Cluster==pop2,"best.source2.ev2"]
                pcol1 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc1])]
                pcol2 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc2])]
                pcol3 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc3])]
                pcol4 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc4])]
                if(run == 2) points(ev1pos1,poppos,pch=15,col=pcol1,xpd=T,cex=2)
                if(run == 2) points(ev1pos2,poppos,pch=15,col=pcol2,xpd=T,cex=2)
                if(run == 2) points(ev1pos2,poppos,pch=22,xpd=T,cex=2,lwd=1.5)

                
                if(res == "1MW")
                {
                    if(run == 2) points(ev2pos1,poppos,pch=15,col=pcol3,xpd=T,cex=2)
                    if(run == 2) points(ev2pos2,poppos,pch=15,col=pcol4,xpd=T,cex=2)
                    if(run == 2) points(ev2pos2,poppos,pch=22,xpd=T,cex=2,lwd=1.5)
                } 
                if(pop2 %in% rev_pops)
                {
                    addancs <- c(pop,anc2,anc1,d)
                } else
                {
                    addancs <- c(pop,anc1,anc2,d)
                }
                if(run == 1) popordertab <- rbind(popordertab,addancs)
            }
            if(res %in% c("2D"))
            {
                pop2 <- gsub("\\_a","",pop)
                anc1 <- pltable[pltable$Cluster==pop2,"best.source1.date1"]
                anc2 <- pltable[pltable$Cluster==pop2,"best.source2.date1"]
                anc3 <- pltable[pltable$Cluster==pop2,"best.source1.date2"]
                anc4 <- pltable[pltable$Cluster==pop2,"best.source2.date2"]
                pcol1 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc1])]
                pcol2 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc2])]
                pcol3 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc3])]
                pcol4 <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==anc4])]
                if(run == 2) points(ev1pos1,poppos,pch=15,col=pcol1,xpd=T,cex=2)
                if(run == 2) points(ev1pos2,poppos,pch=15,col=pcol2,xpd=T,cex=2)
                if(run == 2) points(ev1pos2,poppos,pch=22,xpd=T,cex=2,lwd=1.5)
                if(run == 2) points(ev2pos1,poppos,pch=15,col=pcol3,xpd=T,cex=2)
                if(run == 2) points(ev2pos2,poppos,pch=15,col=pcol4,xpd=T,cex=2)
                if(run == 2) points(ev2pos2,poppos,pch=22,xpd=T,cex=2,lwd=1.5)
                if(pop2 %in% rev_pops)
                {
                    addancs <- c(pop,anc2,anc1,d)
                } else
                {
                    addancs <- c(pop,anc1,anc2,d)
                }
                if(run == 1) popordertab <- rbind(popordertab,addancs)
            }
        }
        if(run == 1)
        {
            popordertab <- popordertab[grep("_a",popordertab[,1],invert = T),]
            popordertab <- cbind(popordertab,
                                 as.character(sapply(popordertab[,2],function(x){
                                     popkey$RegionM[popkey$Ethnic_Group==x]})))
            popordertab <- cbind(popordertab,
                                 as.character(sapply(popordertab[,3],function(x){
                                     popkey$RegionM[popkey$Ethnic_Group==x]})))
            popordertab[,5] <- factor(popordertab[,5],levels=ancreg_list)
            popordertab[,6] <- factor(popordertab[,6],levels=ancreg_list)
            popordertab <- popordertab[order(popordertab[,6],popordertab[,5],
                                             round((1950-as.numeric(popordertab[,4]))/29)),]
        }
    }
    
    ###################################################################
    ## PLOT Y-AXIS NAMES AND COLOURS
    y_ax_cols <- c()
    for(i in poporder) 
    {
        if(i=="SEMI-BANTU") i <- "SEMI.BANTU"
        i <- as.character(sim_tab$sim_src1[sim_tab == i])
        ii <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
        y_ax_cols <- c(y_ax_cols,ii)
    }
    for(i in 1:n_pops)
    {
        poplab <- paste0(as.character(sim_tab$sim_src1[sim_tab$sim_pops==poporder[i]]),"v",
                        as.character(sim_tab$sim_src2[sim_tab$sim_pops==poporder[i]]),": ",
                        as.character(sim_tab$sim_props[sim_tab$sim_pops==poporder[i]]))
                                     
        axis(2,pos=poplabpos,at=i,labels=poplab,col.axis=y_ax_cols[i],las=2,tck=0,lwd=0,line=-0.5)
    }
    for(i in seq(0.5,(n_pops+1),1)) abline(h=i,lty=3,lwd=0.5)
    lines(x=c(median(c(ev1pos2,ev2pos1)),median(c(ev1pos2,ev2pos1))),y=c(0,n_pops+1),lwd=1,xpd=T)
    
    text(y=-1,x=median(c(ev1pos1,ev1pos2)),labels="1ST EVENT",xpd=T,srt=35,adj=0)
    text(y=-1,x=median(c(ev2pos1,ev2pos2)),labels="2ND EVENT",xpd=T,srt=35,adj=0)
    
    ###################################################################
    ## MAKE A DUMMY PLOT TO GET SOME POSITIONS FOR THE BARPLOTS
    x1 <- mixmat[poporder,]
    bp <- barplot(t(as.matrix(x1)),
                  yaxt="n",xlab="",beside=F,
                  main="",cex.main=0.75,border=NA,horiz=T,
                  col=srccols,xaxt="n",add=F,plot=F)

    ###################################################################
    ## PLOT TRUTH 
    plot_mat <- matrix(0,nrow=nrow(all_plot_mat),ncol=length(poporder))
    colnames(plot_mat) <- as.character(poporder)
    for(i in 1:ncol(plot_mat))
    {
        pop2 <- as.character(poporder[i])
        if(pop2 %in% colnames(all_plot_mat)) plot_mat[,pop2] <- unlist(sim_plot_mat[,pop2])
    }
    
    par(mar=c(4,0,4,1))
    plot(0,0,xlim=c(0,1),
         ylim=c(bp[length(bp)],bp[1]),type="n",axes=F,xlab="",ylab="")
    barplot(plot_mat,col=c(srccols,"white",srccols),
            yaxt="n",xlab="",beside=F,
            main="",cex.main=0.75,border=NA,horiz=T,
            xaxt="n",add=T)
    text(y=-1,x=0.2,labels="SIMULATED\nADMIXTURE",xpd=T,srt=35,adj=0)
    
    ###################################################################
    ## PLOT NNLS MIXTURE MODEL
    par(mar=c(4,0,4,0.5))
    plot(0,0,xlim=c(0,1),
         ylim=c(bp[length(bp)],bp[1]),type="n",axes=F,xlab="",ylab="")
    bp <- barplot(t(as.matrix(x1)),
                  yaxt="n",xlab="",beside=F,
                  main="",cex.main=0.75,border=NA,horiz=T,
                  col=srccols,xaxt="n",add=T,plot=T)
    #text(0,-1,labels=LETTERS[2],adj=0,las=1,cex=2,lwd=3,xpd=T)
    text(y=-1,x=0.2,labels="MIXTURE\nMODEL",xpd=T,srt=35,adj=0)
    ##################################################################
    ## PLOT ADMIXTURE SOURCES
    ## FIRST EVENTS
    plot_mat <- matrix(0,nrow=nrow(all_plot_mat),ncol=length(poporder))
    colnames(plot_mat) <- as.character(poporder)
    for(i in 1:ncol(plot_mat))
    {
        pop2 <- as.character(poporder[i])
        if(pop2 %in% colnames(all_plot_mat)) plot_mat[,pop2] <- unlist(all_plot_mat[,pop2])
    }
    par(mar=c(4,0,4,1))
    plot(0,0,xlim=c(0,1),
         ylim=c(bp[length(bp)],bp[1]),type="n",axes=F,xlab="",ylab="")
    barplot(plot_mat,col=c(srccols,"white",srccols),
            yaxt="n",xlab="",beside=F,
            main="",cex.main=0.75,border=NA,horiz=T,
            xaxt="n",add=T)
    text(y=-1,x=0.2,labels="1ST EVENT\nSOURCES",xpd=T,srt=35,adj=0)
    ## SECOND EVENTS
    plot_mat <- matrix(0,nrow=nrow(all_plot_mat),ncol=length(poporder))
    colnames(plot_mat) <- poporder
    for(i in 1:ncol(plot_mat))
    {
        pop2 <- paste0(poporder[i],"_a")
        if(pop2 %in% colnames(all_plot_mat)) plot_mat[,i] <- unlist(all_plot_mat[,pop2])
    }
    par(mar=c(4,0,4,1))
    plot(0,0,xlim=c(0,1),
         ylim=c(bp[length(bp)],bp[1]),type="n",axes=F,xlab="",ylab="")
    barplot(plot_mat,col=c(srccols,"white",srccols),
            yaxt="n",xlab="",beside=F,
            main="",cex.main=0.75,border=NA,horiz=T,
            xaxt="n",add=T)
    text(y=-1,x=0.2,labels="2ND EVENT\nSOURCES",xpd=T,srt=35,adj=0)

    ###################################################################
    ## LEGEND
    par(mar=c(1,1.5,4,0))
    plot(0,0,axes=F,xlab="",ylab="",type="n")
    legend_text <- c("West African Niger-Congo",
                     "Central West African Niger-Congo",
                     "East African Niger-Congo",
                     "South African Niger-Congo",
                     "Nilo-Saharan / Afro-Asiatic",
                     "KhoeSan",
                     "Eurasia",
                     "main event ancestry")
    l <- legend("top",legend=legend_text,pch=c(rep(15,7),22),
                col=c(pcolshex[c(1:3,5,4,6,7)],"black"),bty="n",
                ncol=1,xpd=T,pt.cex=3,x.intersp=1,y.intersp=1.25,
                pt.lwd=2,cex=1.25, title="Ancestry Region")
    par(mar=c(5,5,2,1))

    ###################################################################
    ### MALDER AND GLOBETROTTER COMPARISONS
#     maldata <- read.table("data/AllPopsMalderFinalHAPMAP.txt",header=T,as.is=T)
#     gtdata <- read.csv("data/AfricaGlobetrotterFinalResults.csv",header=T)
#     gtdata <- gtdata[grep("null",gtdata$Analysis,invert=T),]
#     gtdata$Result[is.na(gtdata$Result)] <- "1D"
#     ## MAKE A TABLE OF GLOBETROTTER DATES
#     gtres <- c()
#     for(i in 1:nrow(gtdata))
#     {
#         pop <- as.character(gtdata$Cluster[i])
#         res <- gtdata$Result[i]
#         if(res %in% c("1D","1MW","U"))
#         {
#             date <- gtdata$date.1D[i]
#             gtres <- rbind(gtres,c(pop,round(date)))
#         } else 
#         {
#             date1 <- gtdata$date.2D.1[i]
#             date2 <- gtdata$date.2D.2[i]
#             gtres <- rbind(gtres,c(pop,round(date1)),c(pop,round(date2)))
#         }
#     }
# 
#     gtres <- gtres[order(gtres[,1],as.numeric(gtres[,2])),]
#     colnames(gtres) <- c("pop","gtdate")
#     gtres[,2] <- as.numeric(gtres[,2])
# 
#     ## COMBINE WITH INFERENCE FROM MALDER
#     datecomps <- c()
#     for(i in unique(gtres[,"pop"]))
#     {
#         maldate <- sort(maldata$Date.Gens[maldata$EthnicGroup==i])
#         
#         if(length(maldate)==1)
#         {
#             gtdate <- gtdata$date.1D[gtdata$Cluster==i]
#             gtdatel <- gtdata$date.1D.L[gtdata$Cluster==i]
#             gtdateh <- gtdata$date.1D.H[gtdata$Cluster==i]  
#             maldateci <- sort(maldata$Date.Gens.CI[maldata$EthnicGroup==i]) 
#             maldatel <- maldate - maldateci
#             maldateh <- maldate + maldateci
#             datecomps <- rbind(datecomps,c(i,maldate,maldatel,maldateh,
#                                            round(gtdate),round(gtdatel),round(gtdateh)))
#         }
#         
#         if(length(maldate)==2)
#         {
#             gtdate <- sort(c(gtdata$date.2D.1[gtdata$Cluster==i],
#                              gtdata$date.2D.2[gtdata$Cluster==i]))
#             gtdatel <- sort(c(gtdata$date.2D.1.L[gtdata$Cluster==i],
#                               gtdata$date.2D.2.L[gtdata$Cluster==i]))
#             gtdateh <- sort(c(gtdata$date.2D.1.H[gtdata$Cluster==i],
#                               gtdata$date.2D.2.H[gtdata$Cluster==i]))
#             maldateci <- sort(maldata$Date.Gens.CI[maldata$EthnicGroup==i]) 
#             maldatel <- c(maldate[1] - maldateci[1],maldate[2] - maldateci[2])
#             maldateh <- c(maldate[1] + maldateci[1],maldate[2] + maldateci[2])
#             
#             
#             datecomps <- rbind(datecomps,c(i,maldate[1],maldatel[1],maldateh[1],
#                                            round(gtdate[1]),round(gtdatel[1]),round(gtdateh[1])))
#             datecomps <- rbind(datecomps,c(i,maldate[2],maldatel[2],maldateh[2],
#                                            round(gtdate[2]),round(gtdatel[2]),round(gtdateh[2])))
#             
#         }
#     }
# 
#     colnames(datecomps) <- c("pop","maldate","maldatel","maldateh","gtdate","gtdatel","gtdateh")
#     datecomps <- data.frame(datecomps)
#     
#     pnts <- getPopSymbols(tidyNames(datecomps$pop,fula=F),leginfo)
#     pointpch <-  pnts$pch2plot
#     pointcol <- pnts$col2plot
#     pointrim <- pnts$rim2plot
#     axmax <- max(sapply(apply(datecomps[,2:3],1,as.character),as.numeric))
#     plot(sapply(as.numeric(as.character(datecomps$gtdate)),makeDate,add_BCE=F),
#          sapply(as.numeric(as.character(datecomps$maldate)),makeDate,add_BCE=F),
#          xlim=sapply(c(0,240),makeDate,add_BCE=F),
#          ylim=sapply(c(400,0),makeDate,add_BCE=F),
#          ylab="Date inferred by MALDER",xlab="Date inferred by GLOBETROTTER",
#          type="n",axes=F,main="")
#     abline(a=0,b=1,lwd=2,col="red")
#     x_labs3 <- c(2000,0,-2000,-5000,-10000)
#     x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE")
#     axis(1,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5)
#     axis(2,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5,las=2)
#     for(j in x_labs3) abline(v=j,lty=2)
#     for(j in x_labs3) abline(h=j,lty=2)
#     box()
#     mtext(3,line=0,text=LETTERS[2],adj=0,las=1,cex=1.5,lwd=3,xpd=T)
#     mtext(3,line=0,text="Number of events inferred\nby MALDER",adj=0.5,las=1,cex=0.5,lwd=1,xpd=T)
# 
#     ## ADD CONFIDENCE INTERVALS
#     for(i in 1:nrow(datecomps))
#     {
#         xs <- c(as.numeric(as.character(datecomps[i,"gtdatel"])),
#                 as.numeric(as.character(datecomps[i,"gtdateh"])))
#         ys <- c(as.numeric(as.character(datecomps[i,"maldate"])),
#                 as.numeric(as.character(datecomps[i,"maldate"])))
#         lines(sapply(xs,makeDate,add_BCE=F),
#               sapply(ys,makeDate,add_BCE=F))
#         xs <- c(as.numeric(as.character(datecomps[i,"gtdate"])),
#                 as.numeric(as.character(datecomps[i,"gtdate"])))
#         ys <- c(as.numeric(as.character(datecomps[i,"maldatel"])),
#                 as.numeric(as.character(datecomps[i,"maldateh"])))
#         lines(sapply(xs,makeDate,add_BCE=F),
#               sapply(ys,makeDate,add_BCE=F))
#         
#     }
#     ## ADD POINTS
#     points(sapply(as.numeric(as.character(datecomps$gtdate)),makeDate,add_BCE=F),
#            sapply(as.numeric(as.character(datecomps$maldate)),makeDate,add_BCE=F),
#            pch=as.numeric(pointpch),
#            col=pointrim,bg=pointcol,cex=1,lwd=1)
#     ###################################################################
#     ## COMPARISONS OF DATES WITH GLOBETROTTER INFERENCE ##
# 
#     #??????best_ald ??????
#     datecomps2 <- c()
#     for(i in unique(gtres[,"pop"]))
#     {
#         maldate <- sort(maldata$Date.Gens[maldata$EthnicGroup==i])
#         numevents <- sum(gtres[,"pop"]==i)
#         
#         if(numevents == 2 & length(maldate) == 2)
#         {
#             gtdate <- sort(c(gtdata$date.2D.1[gtdata$Cluster==i],
#                              gtdata$date.2D.2[gtdata$Cluster==i]))
#             gtdatel <- sort(c(gtdata$date.2D.1.L[gtdata$Cluster==i],
#                               gtdata$date.2D.2.L[gtdata$Cluster==i]))
#             gtdateh <- sort(c(gtdata$date.2D.1.H[gtdata$Cluster==i],
#                               gtdata$date.2D.2.H[gtdata$Cluster==i]))
#             
#             maldateci <- sort(maldata$Date.Gens.CI[maldata$EthnicGroup==i]) 
#             maldatel <- maldate - maldateci
#             maldateh <- maldate + maldateci
#             
#             datecomps2 <- rbind(datecomps2,c(i,maldate[1],maldatel[1],maldateh[1],
#                                              round(gtdate[1]),round(gtdatel[1]),round(gtdateh[1])))
#             datecomps2 <- rbind(datecomps2,c(i,maldate[2],maldatel[2],maldateh[2],
#                                              round(gtdate[2]),round(gtdatel[2]),round(gtdateh[2])))
#         }
#         if(numevents == 2 & length(maldate) == 1)
#         {
#             gtdate <- c(gtdata$date.1D[gtdata$Cluster==i])
#             gtdatel <- c(gtdata$date.1D.L[gtdata$Cluster==i])
#             gtdateh <- c(gtdata$date.1D.H[gtdata$Cluster==i])
#             maldateci <- sort(maldata$Date.Gens.CI[maldata$EthnicGroup==i]) 
#             maldatel <- maldate - maldateci
#             maldateh <- maldate + maldateci
#             datecomps2 <- rbind(datecomps2,c(i,maldate,maldatel,maldateh,
#                                              round(gtdate[1]),round(gtdatel[1]),round(gtdateh[1])))
#         }
#         if(numevents == 1 & length(maldate) == 2)
#         {
#             gtdate <- c(gtdata$date.1D[gtdata$Cluster==i])
#             gtdatel <- c(gtdata$date.1D.L[gtdata$Cluster==i])
#             gtdateh <- c(gtdata$date.1D.H[gtdata$Cluster==i])
#             aldate <- as.character(best_ald$Date.Gens[best_ald$adm.pop==tidyNames(i)])
#             if(aldate != "")
#             {
#                 aldate1 <- round(as.numeric(strsplit(aldate,split=" ")[[1]][1]))
#                 aldateci <- round(as.numeric(strsplit(aldate,split=" ")[[1]][3]))
#                 aldate <- aldate1
#                 aldatel <- aldate - aldateci
#                 aldateh <- aldate + aldateci
#                 datecomps2 <- rbind(datecomps2,c(i,aldate,aldatel,aldateh,
#                                                  round(gtdate),round(gtdatel),round(gtdateh)))
#             }  
#         }
#         if(numevents == 1 & length(maldate) == 1)
#         {
#             gtdate <- c(gtdata$date.1D[gtdata$Cluster==i])
#             gtdatel <- c(gtdata$date.1D.L[gtdata$Cluster==i])
#             gtdateh <- c(gtdata$date.1D.H[gtdata$Cluster==i])
#             maldateci <- sort(maldata$Date.Gens.CI[maldata$EthnicGroup==i]) 
#             maldatel <- maldate - maldateci
#             maldateh <- maldate + maldateci
#             datecomps2 <- rbind(datecomps2,c(i,maldate,maldatel,maldateh,
#                                              round(gtdate),round(gtdatel),round(gtdateh)))
#         }  
#     }
#     colnames(datecomps2) <- c("pop","maldate","maldatel","maldateh","gtdate","gtdatel","gtdateh")
#     datecomps2 <- data.frame(datecomps2)
# 
# 
#     pnts <- getPopSymbols(tidyNames(datecomps2$pop,fula=F),leginfo)
#     pointpch <-  pnts$pch2plot
#     pointcol <- pnts$col2plot
#     pointrim <- pnts$rim2plot
# 
#     axmax <- max(sapply(apply(datecomps2[,2:3],1,as.character),as.numeric))
#     plot(sapply(as.numeric(as.character(datecomps2$gtdate)),makeDate,add_BCE=F),
#          sapply(as.numeric(as.character(datecomps2$maldate)),makeDate,add_BCE=F),
#          xlim=sapply(c(0,240),makeDate,add_BCE=F),
#          ylim=sapply(c(400,0),makeDate,add_BCE=F),
#          ylab="Date inferred by MALDER",xlab="Date inferred by GLOBETROTTER",
#          type="n",axes=F,main="")
#     abline(a=0,b=1,lwd=2,col="red")
#     x_labs3 <- c(2000,0,-2000,-5000,-10000)
#     x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE")
#     axis(1,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5)
#     axis(2,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5,las=2)
#     for(j in x_labs3) abline(v=j,lty=2)
#     for(j in x_labs3) abline(h=j,lty=2)
# 
#     box()
#     mtext(3,line=0,text=LETTERS[3],adj=0,las=1,cex=1.5,lwd=3,xpd=T)
#     mtext(3,line=0,text="Number of events inferred\nby GLOBETROTTER",adj=0.5,las=1,cex=0.5,lwd=1,xpd=T)
# 
#     for(i in 1:nrow(datecomps2))
#     {
#         xs <- c(as.numeric(as.character(datecomps2[i,"gtdatel"])),
#                 as.numeric(as.character(datecomps2[i,"gtdateh"])))
#         ys <- c(as.numeric(as.character(datecomps2[i,"maldate"])),
#                 as.numeric(as.character(datecomps2[i,"maldate"])))
#         lines(sapply(xs,makeDate,add_BCE=F),
#               sapply(ys,makeDate,add_BCE=F))
#         xs <- c(as.numeric(as.character(datecomps2[i,"gtdate"])),
#                 as.numeric(as.character(datecomps2[i,"gtdate"])))
#         ys <- c(as.numeric(as.character(datecomps2[i,"maldatel"])),
#                 as.numeric(as.character(datecomps2[i,"maldateh"])))
#         lines(sapply(xs,makeDate,add_BCE=F),
#               sapply(ys,makeDate,add_BCE=F))
#     }
# 
# 
#     points(sapply(as.numeric(as.character(datecomps2$gtdate)),makeDate,add_BCE=F),
#            sapply(as.numeric(as.character(datecomps2$maldate)),makeDate,add_BCE=F),
#            pch=as.numeric(pointpch),col=pointrim,bg=pointcol,cex=1,lwd=1)
# 

dev.off()





