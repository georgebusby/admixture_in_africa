### SCRIPT TO PLOT ALDER HEATMAP ###
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

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
## LOAD LEGEND INFO
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
leginfo <- read.table(leginfo_file, header = T, comment.char = "")

## LOAD POPORDER FILE ##
popplot <- scan("/mnt/kwiat/home/popgen/scripts/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplotorder <- popplot
popplotorder <- popplotorder[c(1:16,18:61)]

pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
######################################
## PULL IN DATA FROM SERVER

poplist <- factor(popplotorder,levels=popplotorder)
poplist<- as.matrix(poplist)
n_pops <- nrow(poplist)
all_res <- c()
for(i in 50:60) #(1:nrow(poplist)))
{
    pop1 <- as.character(poplist[i,1])
    print(paste("getting alder results for pop:", pop1))
    
    for(j in (1:nrow(poplist)))
    {
        if(i != j)
        {
            popA <- as.character(poplist[j,1])
            file_root <- paste0("/mnt/kwiat/data/bayes/users/george/popgen/analysis3/alder/output/",pop1,"/",pop1,popA)
            file_root2 <- paste0("/mnt/kwiat/data/bayes/users/george/popgen/analysis3/alder/output/",pop1,popA,"_alderonerefF")
            ######################################
            ## FILE READ IN
            if(file.info(paste0(file_root2,".ampexp"))$size > 0)
            {
                alder_amp <- read.table(paste0(file_root2,".ampexp"),header=F,sep="\t")
                alder_amp <- strsplit(as.character(alder_amp[nrow(alder_amp),]),split="\\ ")[[1]][c(8,10)]
                alder_decay <- read.table(paste0(file_root2,".decay"),header=F,sep="\t")
                alder_decay <- strsplit(as.character(alder_decay[nrow(alder_decay),]),split="\\ ")[[1]][c(15,17)]
                alder_mixture <- read.table(paste0(file_root2,".mixture"),header=F,sep="\t")
                alder_mixture <- strsplit(as.character(alder_mixture[nrow(alder_mixture),]),split="\\ ")[[1]][c(2,4)]
                all_res <- rbind(all_res,c(pop1,popA,as.numeric(alder_amp),as.numeric(alder_decay),as.numeric(alder_mixture)))
            } else 
            {
                all_res <- rbind(all_res,c(pop1,popA,rep(NA,6)))
            }
        }
    }
}

all_res_bak <- all_res
all_res <- data.frame(all_res)
colnames(all_res) <- c("pop1","ref.pop","amp","ampCI","decay","decayCI","prop","propCI")

pdf("figures/MalariaGENAlderAmplitudeComparisonsSampledPops.pdf",height=10,width=10)
par(mar=c(1,5,2,1))
let_cnt <- 1

pops2plot <- c("JOLA","FULAI","AKANS","MOSSI","YRI","MALAWI","CHONYI","MZIGUA","JUHOAN")
layout(matrix(c(1:length(pops2plot))))
for(p in pops2plot)
{
  tt <- all_res[which(all_res$pop1==p),]
  tt<-tt[order(as.numeric(as.character(tt$amp)),decreasing=T),]
  
  pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
  
  plot_things <- getPopSymbols(tidyNames(tt$ref.pop,fula=F),leginfo)
  plot_col <- plot_things$col2plot
  plot_pch <- plot_things$pch2plot
  plot_rims <- plot_things$rim2plot
  y_max <- as.numeric(as.character(tt$amp))+as.numeric(as.character(tt$ampCI))
  y_max <- y_max[is.finite(y_max)]
  y_max <- max(y_max,na.rm=T)
  plot(1:nrow(tt),ylim=c(0,y_max),
       axes=F,ylab=expression(paste("amplitude x", 10^-4)),
       xlab="",main=tidyNames(p,fula=T,khoesan=T))
  box()
  #mtext(3,text=LETTERS[let_cnt],adj=0)
  for(j in 1:nrow(tt))
  {
    a<-as.numeric(as.character(tt$amp[j]))
    a1<-as.numeric(as.character(tt$ampCI[j]))
    lines(c(j,j),c(a-a1,a+a1),lwd=2,col=plot_col[j])
  }
  points(1:nrow(tt),
         as.numeric(as.character(tt$amp)),
         pch=as.numeric(plot_pch),
         bg=plot_col,col=plot_rims,cex=2)
  
  ys <-seq(0,y_max,length.out=3)
  axis(2,at=ys,labels=round(ys*10000,1),las=2)
  legend("topright",bty="n",text.col=plot_col[1:3],adj=c(0,0),cex=1,ncol=3,
         legend=paste0(1:3,": ",tidyNames(tt$ref.pop,fula=T,khoesan=T)[1:3]))
  let_cnt <- let_cnt+1
}
dev.off()


########################################
### PLOT A HEAT MAP OF THE RANK MATRIX OF THE REF AMPLITUDES ##
dummy_dat <- matrix(0,ncol=length(popplotorder),nrow=length(popplotorder))
colnames(dummy_dat) <- rownames(dummy_dat) <- popplotorder

for(i in levels(all_res$pop1))
{
    for(j in levels(all_res$ref.pop[all_res$pop1==i]))
    {
        k <- as.numeric(as.character(all_res$amp[all_res$pop1==i&all_res$ref.pop==j]))
        if(length(k)==0) k <- 0
        dummy_dat[i,j] <- k
    }
}


########################################
## GET RANK OF AMPLITUDES FOR EACH POP
library(RColorBrewer)
for(i in 1:nrow(dummy_dat)) dummy_dat[i,] <- rank(-as.numeric(dummy_dat[i,]))
## CAP AT 15
dummy_dat[dummy_dat>15] <- 15
dummy_dat <- dummy_dat[rev(popplotorder[1:48]),]

topcols <- rev(brewer.pal(9,"Greens"))
othercols <- c(rep(topcols[length(topcols)],5),"grey")

y_ax_cols <- c()
for(i in colnames(dummy_dat)) 
{
    ii <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
    y_ax_cols <- c(y_ax_cols,ii)
}

pdf("figures/MalariaGENAlderAmplitudeComparisonsRankMatrix.pdf",height=10,width=11)
layout(matrix(c(1,2),1,2),widths=c(10,1))
par(mar=c(0.5,6,6,0.5))
image(1:ncol(dummy_dat),
      1:nrow(dummy_dat),
      t(dummy_dat),col=c(topcols,othercols),
      ylab="",xlab="",axes=F)

for(i in 1:ncol(dummy_dat))
{
    axis(3,at=i,las=2,
         labels=tidyNames(colnames(dummy_dat)[i],fula=T,khoesan=T),
         cex.axis=0.75,lwd=0, col.axis=y_ax_cols[i])
}
mtext(1,line=8,text="reference population amplitude (rank)")
mtext(2,line=8,text="test population")
for(i in 1:nrow(dummy_dat))
{
    axis(2,at=i,las=2,
         labels=tidyNames(rownames(dummy_dat)[i],fula=T, khoesan=T),
         tck=0,lwd=0,cex.axis=0.75,col.axis=rev(y_ax_cols[1:48])[i])
}
abline(h=1:nrow(dummy_dat)-0.5,col="grey",lwd=0.5)
abline(v=1:ncol(dummy_dat)-0.5,col="grey",lwd=0.5)
## ADD SCALE BAR
sc_bar <- c(11:1)
par(mar=c(8,0.5,20,3.5))
image(1,1:length(sc_bar),
      t(sc_bar),col=c(topcols,"grey"),axes=F,ylab="",xlab="",
      main="")
axis(4,at=c(12,sc_bar),
     labels=c("RANK",1:9,"10-14","15+"),tck=0,lwd=0,las=2,xpd=T,adj=0.5)
abline(v=1:length(sc_bar)-0.5,col="grey",lwd=0.5)

dev.off()


# ######################################
# 
# ## why do Jola and Malawi choose Juhoan first?
# ## note that most African groups have similar
# ## eurasian signals to Pickrell but some don't
# ## what does this mean?
# 
# ### plot amplitudes with Yoruba v TSI
# 
# pop_x <- "JUHOAN"
# pop_y <- "YRI"
# afr_pops <- poplist[1:48,1]
# afr_pops <- afr_pops[afr_pops!=pop_y]
# afr_pops <- afr_pops[afr_pops!=pop_x]
# 
# pdf(paste0("alder/AfricanPopsAmplitudeComparisonsF_",pop_x,"_",pop_y,".pdf"),
#     height=5,width=5)
# par(mar=c(4,5,2,1))
# plot1 <- as.numeric(as.character(all_res$amp[all_res$ref.pop==pop_x&as.character(all_res$pop1)%in%afr_pops]))
# plot1ci <- as.numeric(as.character(all_res$ampCI[all_res$ref.pop==pop_x&as.character(all_res$pop1)%in%afr_pops]))#/sqrt(22)
# plot2 <- as.numeric(as.character(all_res$amp[all_res$ref.pop==pop_y&as.character(all_res$pop1)%in%afr_pops]))
# plot2ci <- as.numeric(as.character(all_res$ampCI[all_res$ref.pop==pop_y&as.character(all_res$pop1)%in%afr_pops]))#/sqrt(22)
# plot_pops <- as.character(afr_pops)
# plot_col <- sapply(as.character(plot_pops),function(x){as.character(info3$popcols.info.Country.[info3$EthnicGroup==tidyNames(x)])})
# plot_pch <- sapply(as.character(plot_pops),function(x){as.character(info3$poppch[info3$EthnicGroup==tidyNames(x)])})
# plot_rim <- sapply(as.character(plot_pops),function(x){as.character(info3$rim[info3$EthnicGroup==tidyNames(x)])})
# plot_rims <- rep("black",length(plot_rim))
# plot_rims[plot_rim==1] <- plot_col[plot_rim==1]
# 
# max_lim <- max(max(plot1[!is.na(plot2)],na.rm=T),max(plot2[!is.na(plot2)],na.rm=T))
# plot(plot1,plot2,pch=as.numeric(plot_pch),bg=plot_cols,
#      xlab=bquote(.(paste0("amplitude (", tidyNames(pop_x,fula=F)," reference) x ")) ~ 10^-4),
#      ylab=bquote(.(paste0("amplitude (", tidyNames(pop_y,fula=F)," reference) x ")) ~ 10^-4),
#      xlim=c(0,max_lim),ylim=c(0,max_lim), axes=F,type="n")
# box()
# axis(1,at=seq(0,max_lim,length=5),labels=round(seq(0,max_lim,length=5)*1000,2))
# axis(2,at=seq(0,max_lim,length=5),labels=round(seq(0,max_lim,length=5)*1000,2),las=2)
# abline(a=0,b=1)
# for(i in 1:length(plot1ci)) lines(c(plot1[i]+plot1ci[i],plot1[i]-plot1ci[i]),c(plot2[i],plot2[i]),col=plot_col[i])
# for(i in 1:length(plot1ci)) lines(c(plot1[i],plot1[i]),c(plot2[i]+plot2ci[i],plot2[i]-plot2ci[i]),col=plot_col[i])
# points(plot1,plot2,pch=as.numeric(plot_pch),bg=plot_col,col=plot_rims)
# 
# dev.off()
# 
# 
# 
# 
# ########################################
# plot1 <- as.numeric(as.character(all_res$decay[all_res$ref.pop=="TSI"&as.character(all_res$pop1)%in%afr_pops]))
# plot2 <- as.numeric(as.character(all_res$decay[all_res$ref.pop=="JUHOAN"&as.character(all_res$pop1)%in%afr_pops]))
# plot_pops <- as.character(afr_pops)
# plot_col <- sapply(as.character(plot_pops),function(x){as.character(info3$popcols.info.Country.[info3$EthnicGroup==tidyNames(x)])})
# plot_pch <- sapply(as.character(plot_pops),function(x){as.character(info3$poppch[info3$EthnicGroup==tidyNames(x)])})    
# plot(plot1,plot2,pch=as.numeric(plot_pch),bg=plot_cols,xlab="decay (TSI reference)",
#      ylab="decay (JUHOANSI reference)")
# abline(a=0,b=1)
# 
