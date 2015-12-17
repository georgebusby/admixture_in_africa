## PLOTS FULL SUPPLEMENTARY MAP
## RUN THE BEGINNING OF FIGURE1.R

setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
fsanalyname <- 'MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalForce.A'
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
mapinfo_file <- "data/MalariaGenAdmixturePopulationKey2.txt"
latlong_file <- "data/MalariaGenAdmixturePopulationKeyLatLongs.txt"
poppos_file <- "data/MalariaGenAdmixturePopulationKeyMapPositions.txt"
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"
popnums_file <- "data/MalariaGenPopulationKeyCPanalysisPopNumbers.txt"
pca_file <- "data/Africa300KPCS.txt"
tree_file <- paste0("data/",fsanalyname,".mcmc.tree.xml")
mat_file <- paste0("data/",fsanalyname,"CoAncestry.txt")

## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)

##########################################################
##########################################################
pdf("figures/MalariaGenMAPofPops.pdf",height=9,width=9)
layout(matrix(c(1,1,5,1,1,5,2,3,4),3,3,byrow=T),heights=c(3,3,3))
par(mar=c(0.5,0.5,0.5,0.5))
plot_legends<- FALSE
plot_letter <- FALSE
pt_cex <- 2
pt_lwd <- 0.75
source("plotting_scripts/africamap.R")

##########################################################
## PLOT LEGENDS
afr_cnts <- list(c("Gambia","Mali","BurkinaFaso","Ghana","Nigeria","Cameroon"),
                 c("Sudan","Ethiopia","Somalia","Tanzania","Kenya"),
                 c("Angola","Namibia","Botswana","SouthAfrica","Malawi"))
names(afr_cnts) <- c("WEST/CENTRAL","EASTERN","SOUTHERN")

for(i in 1:length(afr_cnts))
{
    leginfo <-read.table(leginfo_file,header=T,comment.char="")
    legnumbers <- read.table(pop_nums_file,header=F)
    legnumbers[,1] <- tidyNames(legnumbers[,1])
    pop_nums <- c()
    pop_orig <- c()
    for(j in as.character(leginfo$EthnicGroup))
    {
        pn <- legnumbers[legnumbers$V1==j,2]
        po <- legnumbers[legnumbers$V1==j,3]
        if(length(pn)>0)
        {
            pop_nums <- c(pop_nums,pn)    
            pop_orig <- c(pop_orig,as.character(po))
        }  else
        {
            pop_nums <- c(pop_nums,"X")
            pop_orig <- c(pop_orig,"X")
        }
    }
    leginfo <- cbind(leginfo,pop_nums,pop_orig)
    leginfo <- leginfo[leginfo$Country%in%unlist(afr_cnts[i]),]
    leginfo$Country <- factor(leginfo$Country,levels=unlist(afr_cnts[i]))
    leginfo <- leginfo[order(leginfo$Country,leginfo$Colour,leginfo$EthnicGroup),]
    if(i==3)
    {
        leginfo$Country <- factor(leginfo$Country,levels=c("Tanzania","Malawi","SouthAfrica","Namibia","Botswana","Angola"))
        leginfo <- leginfo[order(leginfo$Country,decreasing=F),]
    }
    legrim <- rep("#000000",nrow(leginfo))
    legrim[leginfo$rim==1]  <- as.character(leginfo$Colour)[leginfo$rim==1]
    leginfo <- leginfo[order(leginfo$pop_orig),]
    par(mar=c(0,0,0,0))
    plot(0,0,type="n",axes=F,xlab="",ylab="")
    legend_text <- tidyNames(as.character(leginfo$EthnicGroup),khoesan=T)
    legend_text[leginfo$pop_orig == "MALGEN"] <- paste0(legend_text[leginfo$pop_orig == "MALGEN"],"\t(MG)")
    legend_text[leginfo$pop_orig == "SCHLEBUSCH"] <- paste0(legend_text[leginfo$pop_orig == "SCHLEBUSCH"],"\t(SC)")
    legend_text[leginfo$pop_orig == "PAGANI"] <- paste0(legend_text[leginfo$pop_orig == "PAGANI"],"\t(PA)")
    legend_text[leginfo$pop_orig == "PETERSON"] <- paste0(legend_text[leginfo$pop_orig == "PETERSON"],"\t(PE)")
    legend_text[leginfo$pop_orig == "PETERSON/SCHLEBUSCH"] <- paste0(legend_text[leginfo$pop_orig == "PETERSON/SCHLEBUSCH"],"\t(PS)")
    legend_text[leginfo$pop_orig == "1KGP"] <- paste0(legend_text[leginfo$pop_orig == "1KGP"],"\t(1K)")
    legend_text <- c(gsub("\\_","\n",paste0(legend_text," [",as.character(leginfo$pop_nums),"]")),"\n")
    legend("top",ncol=2,
           legend=legend_text,
           pch=as.numeric(leginfo$poppch),
           pt.bg=as.character(leginfo$Colour),
           pt.cex=c(rep(pt_cex,length(leginfo$poppch)),0),
           pt.lwd=0.5,y.intersp = 0.9,
           col=legrim,cex=1,
           text.col=as.character(leginfo$Colour),xpd=T,
           bty="n",title=names(afr_cnts)[i],
           title.col="black",title.adj=0)
}
##########################################################
## PLOT EURASIA AND KEY
par(mar=c(0,0,4,0))
plot(0,0,type="n",axes=F,xlab="",ylab="")
leginfo <-read.table(leginfo_file,header=T,comment.char="")
legnumbers <- read.table(pop_nums_file,header=F)
legnumbers[,1] <- tidyNames(legnumbers[,1])
pop_nums <- c()
pop_orig <- c()
for(j in as.character(leginfo$EthnicGroup))
{
    pn <- legnumbers[legnumbers$V1==j,2]
    po <- legnumbers[legnumbers$V1==j,3]
    if(length(pn)>0)
    {
        pop_nums <- c(pop_nums,pn)    
        pop_orig <- c(pop_orig,as.character(po))
    } else 
    {
        pop_nums <- c(pop_nums,"X")
        pop_orig <- c(pop_orig,"X")
    }
}
leginfo <- cbind(leginfo,pop_nums,pop_orig)
leginfo <- leginfo[leginfo$Region%in%c("AMERICAS","ASIA","EUROPE"),]
leginfo <- leginfo[leginfo$pop_orig!="X",]
leginfo <- leginfo[order(leginfo$EthnicGroup),]
leginfo <- leginfo[order(leginfo$Region,decreasing=T),]
#leginfo$Country <- factor(leginfo$Country,levels=unlist(afr_cnts[i]))

legend_text <-as.character(leginfo$EthnicGroup)
legend_text[leginfo$pop_orig == "1KGP"] <- paste0(legend_text[leginfo$pop_orig == "1KGP"],"\t(1K)")
legend_text <- c(gsub("\\_","\n",paste0(legend_text," [",as.character(leginfo$pop_nums),"]")),"\n")
# legend("bottomleft",ncol=1,
#        legend=legend_text[1:12],
#        #pch = rep(20,12),
#        #pt.bg=rep("white",12),
#        #pt.cex=rep(0,12),
#        #pt.lwd=0,
#        col="white",
#        text.col="black",
#        y.intersp = 0.9,
#        cex=1,
#        xpd=T,
#        bty="n",title="EURASIA (not shown)",
#        title.col="black",title.adj=0)

     
legend("topleft",legend=c("EURASIA (not shown)",
                          legend_text,
                          "KEY TO POPULATION ORIGINS",
                          "1K = ONE THOUSAND GENOME PROJECT",
                          "MG = MALARIAGEN POPULATIONS",
                          "PA = PAGANI ET AL 2012 AJHG ",
                          "PE = PETERSON ET AL 2014 PLOS GENETICS",
                          "SC = SCHLEBUSCH ET AL 2014 SCIENCE",
                          "PS = POPULATIONS FROM PE AND SC"),    
       bty="n",
       title.col="black",title.adj=0)
  
dev.off()
