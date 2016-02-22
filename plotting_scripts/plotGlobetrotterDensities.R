### SCRIPT TO PLOT GLOBETROTTER RESULTS FOR PAPER ###
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
library(xtable)
library(grid)
library(gridExtra)
library(ggplot2)

setwd(paste0(main_dir,"popgen/"))
###########################################################
## DEFINE DATAFILES
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$Ethnic_Group[popkey$Ethnic_Group=="SEMI-BANTU"] <- "SEMI.BANTU"
## LOAD LEGEND INFO
leginfo <- read.table(leginfo_file, header = T, comment.char = "")

### DEFINE SOME PLOTTING VARIABLES ##
dateLabelCex=1
datelines=c(-1000,-500,0,500,1000,1500)
yAxisLim=c(-1000,2000)
## OLD VERSION OF THE COLOURS
#pcolshex <- c("#0000CD", "#03B4CC", "#A65628", "#FF7F00", "#984EA3", "#4DAF4A", "#CCCC00")[c(1,2,4,5,3,6,7)]
## NEW COLOURS THAT DIFFERENTIATE AA AND NS
pcolshex <- c("#0000CD","#03B4CC","#FF7F00","#984EA3","#FF69B4","#A65628","#4DAF4A","#CCCC00")
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
popplot <- scan("/mnt/kwiat/home/popgen/scripts/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplot <- popplot[popplot!="SEMI-BANTU"]
popplotorder <- popplot
popplot <- factor(popplot,levels=popplotorder)

## PULL IN DATA
admixturesources2 <- read.table("data/MalariaGenGlobetrotterAdmixtureSources3.txt",header=T,row.names=1,as.is=T)
#final.res2plot2 <- read.table("data/MalariaGenGlobetrotter2plot.txt",header=T,row.names=1,as.is=T)
final.res2plot <- read.table("data/MalariaGenGlobetrotter2plotFinal.txt",header=T,row.names=1,as.is=T)
admixturesources2 <- admixturesources2[,levels(popplot)]
dateboots <- read.table("data/MalariaGenGlobetrotterOneDateBootstraps.txt",header=T,row.names=1,as.is=T)
date2boots <- read.table("data/MalariaGenGlobetrotterTwoDateBootstraps.txt",header=T,row.names=1,as.is=T)
best_ald <- read.table("data/MalariaGenBestAlder.txt",header=T,as.is=T,comment.char="",fill=T)


## DEFINE THE REGION AND COLOURS FOR EACH OF THE DIFFERENT SOURCES
srcreg <- getPopRegion(tidyNames(colnames(admixturesources2),fula=T),popkey)
regions <- as.character(ancreg_list)
srcreg[tidyNames(colnames(admixturesources2),fula=T) == "KARRETJIE"] <- "South_Africa_KhoeSan"
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

## TO MAKE SIMILAR SOURCES ALIGN IN THE FIGURE BELOW, WE MIGHT SOMETIMES
## WANT TO DEFINE THE MAJOR SOURCE AS THE ONE THAT CONTRIBUTES THE 
## LEAST TO AN EVENT. rev_pops ALLOWS US TO REVERSE WHICH IS THE MAJOR
## AND WHICH IS THE MINOR SOURCE IN THE OUTPUT TABLE
rev_pops <- c("KHWE","AMAXHOSA","SEBANTU","MALAWI")
pops <- as.character(pltable$Cluster)
## I'VE ADDED A FUNCTION TO PULL THESE THINGS OUT
tempresults <- addGTresults(pltable,rev_pops)
all_dates <- tempresults[[1]]
all_plot_mat <- tempresults[[2]]
all_src_mat <- tempresults[[3]]
dboots <- tempresults[[4]]
allsources <- tempresults[[5]]


###########################################################
## NEED TO MAKE SURE THAT WE'RE USING THE SAME DATE BOOTS AS
## THOSE REPORTED IN FIGURE 4

dboots2get <- final.res2plot[,c(1,2,4)]
dboots2get$boots <- "dateboots"
dboots2get$boots[dboots2get$Result=="2D"] <- "date2boots"
dboots2get$pick <- 1
dboots2get$pick[dboots2get$Result=="2D"] <- "both"
test <- dboots2get$Cluster%in%popkey$Ethnic_Group[popkey$RegionM=="Western_Africa_Niger-Congo"]
dboots2get$boots[test] <- "date2boots"
dboots2get$pick[test] <- "2"

## NOW UPDATE dboots
for(i in 1:nrow(dboots2get))
{
    run <- dboots2get$Analysis[i]
    tmp <- get(dboots2get$boots[i])
    tmp <- tmp[tmp$pop==dboots2get$Cluster[i]&tmp$X.main.==run,]
    if(dboots2get$boots[i] == "dateboots")
    {
        bts <- tmp$date1.est.boot
        if(length(bts) < 100) bts <- c(bts,rep(mean(bts),100-length(bts)))
            
    } else if(dboots2get$boots[i] == "date2boots" & dboots2get$pick[i] == 2)
    {
        bts <- tmp$date2.est.boot
        if(length(bts) < 100) bts <- c(bts,rep(mean(bts),100-length(bts)))
    } else
    {
        bts <- tmp$date1.est.boot
        if(length(bts) < 100) bts <- c(bts,rep(mean(bts),100-length(bts)))
        bts2 <- tmp$date2.est.boot
        if(length(bts2) < 100) bts2 <- c(bts2,rep(mean(bts2),100-length(bts2)))
        bts <- rbind(bts,bts2)
    }
    dboots[dboots[,1]==dboots2get$Cluster[i],6:ncol(dboots)] <- bts
}





###########################################################
## PLOT DATE DENSITIES
## RESHAPE DATEBOOTS FOR GGPLOT
dboots2 <- c()
for(i in 1:nrow(dboots))
{
    tboots <- cbind(dboots[i,1],dboots[i,4],dboots[i,5],as.numeric(dboots[i,6:ncol(dboots)]))
    colnames(tboots) <- c("pop","don1","don2","dates")
    dboots2 <- rbind(dboots2,tboots)
}
dboots2 <- data.frame(dboots2,stringsAsFactors=F)
dboots2$dates <- 1950-round((sapply(dboots2$dates,as.numeric) * (29+1))-65)
for(i in c("pop","don1","don2"))
{
    newregs <- dboots2[,i]
    newregs <- as.character(sapply(newregs,function(x){
        popkey$RegionM[popkey$Ethnic_Group==x]
    }))
    dboots2 <- cbind(dboots2,newregs)
    colnames(dboots2)[ncol(dboots2)] <- paste0(i,".reg")
}
################################################################################
## WHAT ABOUT SPLITTING UP SOURCES INTO ACUTAL MIXTURES?
## RESHAPE DATEBOOTS FOR GGPLOT
## GET THE ANCESTRY FOR EACH SOURCE
event_srcs <- apply(allsources[,6:ncol(allsources)],2,as.numeric)
## ROUND TO THE NEAREST 0.01
event_srcs[event_srcs < 0.025] <- 0
event_srcs <- event_srcs/rowSums(event_srcs)
event_srcs <- event_srcs*100
event_srcs <- apply(event_srcs,2,round)
event_srcs <- event_srcs*as.numeric(allsources[,4])

## NOW GET SOURCE CONRITBUTIONS IN PROPORTION TO EFFECT ON EVENTS
#for(i in 1:nrow(event_srcs)) event_srcs[i,] <- event_srcs[i,] * as.numeric(allsources[i,4])



## NOW FOR EACH EVENT, MULTIPLY THE DATES BY EACH OF THE COLUMNS...
dboots3 <- c()
event_pop <- unique(dboots[,1])

for(pop in event_pop)
{
    res <- unique(dboots[dboots[,1]==pop,3])
    don1 <- dboots[dboots[,1]==pop,4]
    don2 <- dboots[dboots[,1]==pop,5]
        
    event_row <- allsources[,1] == pop & allsources[,2] == res
    facts <- event_srcs[event_row,]
    date_row <- dboots[,1] == pop
    if(sum(date_row) == 2)   dates <- apply(apply(dboots[date_row,6:ncol(dboots)],2,as.numeric),2,round)
    if(sum(date_row) == 1)   dates <- sapply(sapply(dboots[date_row,6:ncol(dboots)],as.numeric),round)
    for(j in 1:nrow(facts))
    {
        tfacts <- facts[j,facts[j,]>0]
        for(k in 1:length(tfacts))
        {
            if( j %in% c(1,2) ) date_row <- 1
            if( j %in% c(3,4) ) date_row <- 2
            if( j %in% c(1,3) ) src <- "min"
            if( j %in% c(2,4) ) src <- "maj"
            comp_prop <- as.numeric(round(tfacts[k]))
            if(comp_prop<1) comp_prop <- 1
            if(!is.null(nrow(dates))) dateprops <- rep(dates[date_row,],comp_prop)
            if(is.null(nrow(dates)))  dateprops <- rep(dates,comp_prop)
            pop_reg <- getPopRegion(tidyNames(pop,fula=T),popkey)
            if(pop == "KARRETJIE") pop_reg <- "South_Africa_KhoeSan"
            don_comp <- names(tfacts[k])
            don_comp_reg <- getPopRegion(tidyNames(don_comp,fula=T),popkey)
            if(don_comp == "KARRETJIE") don_comp_reg <- "South_Africa_KhoeSan"
            tboots <- cbind(pop,pop_reg,
                            res,don1[date_row],don2[date_row],
                            don_comp,don_comp_reg,
                            src,dateprops)
            colnames(tboots) <- c("pop","pop.reg","res","don1","don2","don1.component",
                                  "don1.component.reg","src","dates")
            dboots3 <- rbind(dboots3,tboots)
        }
    }
}

dboots3 <- data.frame(dboots3,stringsAsFactors=F)
dboots3$dates <- sapply(sapply(dboots3$dates,as.numeric),makeDate,add_BCE=F)
plotXlimggplot <- c(1950,-2000)

###################################################
### AND PLOT
dboots4 <- dboots3
dboots4$don1.component.cnt <- dboots4$don1.component.reg
#dboots4$don1.component <- factor(dboots4$don1.component,levels=popplot,labels=popplot)
reg_labs <- c("West Africa NC","Central West Africa NC","East Africa NC",
              "South Africa NC","Nilo-Saharan","Afroasiatic",
              "Khoesan","Eurasia")
dboots4$pop.reg <- factor(dboots4$pop.reg,levels=ancreg_list[c(1:3,6,4,5,7,8)],labels=reg_labs)
dens_info <- cbind(reg_labs,pcolshex[c(1:3,6,4,5,7,8)],ancreg_list[c(1:3,6,4,5,7,8)])
###################################################
## REPLACE OLD WITH NEW
for(i in c(1:8))
{
    old_reg <- dens_info[i,3]
    new_reg <- dens_info[i,1]
    dboots4$don1.component.cnt[dboots4$don1.component.reg==old_reg] <- new_reg
}


## ADD A FEW EXTRA LEVELS TO don1.component.reg
## ie SPLIT UP EURASIA
## HIGHLIGHT FULAI AND SEMI-BANTU
## NORTH EUROPE
pops <- c("CEU","GBR","FIN")
new_pops <- "North Europe"
new_col <- "#eeee13"
dboots4$don1.component.cnt[dboots4$don1.component%in%pops] <- new_pops
dens_info <- rbind(dens_info,c(new_pops,new_col))

## SOUTH EUROPE
pops <- c("IBS","TSI")
new_pops <- "South Europe"
new_col <- "#CCCC00"
dboots4$don1.component.cnt[dboots4$don1.component%in%pops] <- new_pops
dens_info <- rbind(dens_info,c(new_pops,new_col))

## SOUTH ASIA
pops <- c("KHV","GIH")
new_pops <- "South Asia"
new_col <- "#787833"
dboots4$don1.component.cnt[dboots4$don1.component%in%pops] <- new_pops
dens_info <- rbind(dens_info,c(new_pops,new_col))

## EAST ASIA
pops <-c("CHB","CHS","CHD","CDX","JPT","PELII")
new_pops <- "East Asia"
new_col <- "#8B8878"
dboots4$don1.component.cnt[dboots4$don1.component%in%pops] <- new_pops
dens_info <- rbind(dens_info,c(new_pops,new_col))

## SEMI-BANTU
# pops <- c("BANTU","SEMI.BANTU")
# new_pops <- "Bantu/Semi-Bantu"
# new_col <- "#FF001E"
# dboots4$don1.component.cnt[dboots4$don1.component%in%pops] <- new_pops
# dens_info <- rbind(dens_info,c(new_pops,new_col))

## GET RID OF THE EURASIA LEVEL
dens_info <- dens_info[dens_info[,"reg_labs"]!="Eurasia",]


## MAKE INTO A FACTOR
dboots4$don1.component.cnt <- factor(dboots4$don1.component.cnt,
                                     levels=dens_info[,"reg_labs"],
                                     labels=dens_info[,"reg_labs"])
## sort out colours
dens_cols <- dens_info[,2]

##################################################
## SPLIT MAJOR AND MINOR SOURCES
dbootsMin <- dboots4[dboots4[,"src"]=="min",]
dbootsMaj <- dboots4[dboots4[,"src"]=="maj",]

## MAKE PRETTY NAMES FOR MINOR SOURCE REGIONS
dbootsMin$don1.component.reg <- sapply(dbootsMin$don1.component.reg,as.character)
for(i in unique(dbootsMin$don1.component.reg)[grep("Niger-Congo",unique(dbootsMin$don1.component.reg))])
{
    dbootsMin$don1.component.reg <- gsub(i,"All Niger Congo",dbootsMin$don1.component.reg)
}
new_regs <- unique(dbootsMin$don1.component.reg)
new_regs <- gsub("\\_"," ",new_regs)
new_regs <- gsub("South Africa KhoeSan","Khoesan",new_regs)
new_regs <- gsub("East Africa Nilo-Saharan","Nilo-Saharan",new_regs)
new_regs <- gsub("East Africa Afro-Asiatic","Afroasiatic",new_regs)
new_regstab <- cbind(unique(dbootsMin$don1.component.reg),new_regs)

for(i in 2:nrow(new_regstab))
{
    dbootsMin$don1.component.reg <- gsub(new_regstab[i,1],new_regstab[i,2],dbootsMin$don1.component.reg)
}

dbootsMin$don1.component.reg <- factor(dbootsMin$don1.component.reg,levels=new_regs[c(1,4,3,5,2)])

## NOW FOR MAJOR SOURCES
dbootsMaj$don1.component.reg <- sapply(dbootsMaj$don1.component.reg,as.character)
for(i in unique(dbootsMaj$don1.component.reg)[grep("Niger-Congo",unique(dbootsMaj$don1.component.reg))])
{
    dbootsMaj$don1.component.reg <- gsub(i,"All Niger Congo",dbootsMaj$don1.component.reg)
}

for(i in 2:nrow(new_regstab))
{
    dbootsMaj$don1.component.reg <- gsub(new_regstab[i,1],new_regstab[i,2],dbootsMaj$don1.component.reg)
}
dbootsMaj$don1.component.reg <- factor(dbootsMaj$don1.component.reg,levels=new_regs[c(1,4,3,5,2)])    


p7 <- ggplot(dbootsMin, aes(x=dates)) +
    geom_density(aes(fill=don1.component.cnt,colour=don1.component.cnt,y=..count..), position="stack") +
    geom_density(data=dbootsMaj,
                 aes(fill=don1.component.cnt,colour=don1.component.cnt,y=-..count..), position="stack") +
    xlab("Date of Admixture") + ylab("Donor Region") +
    ggtitle("Recipient Region") +
    scale_fill_manual(values=as.character(dens_cols),name="") +
    scale_colour_manual(values=as.character(dens_cols),name="") +
    theme_bw() +
    guides(linetype=FALSE,
           colour = guide_legend(nrow = 3),
           fill = guide_legend(nrow = 3)) +
    coord_cartesian(xlim = plotXlimggplot) +
    facet_grid(don1.component.reg~pop.reg) + #don1.component.reg
    theme(strip.text.y = element_text(size= 10 ),
          strip.text.x = element_text(size= 10 ),
          strip.background = element_rect(colour="black", fill="white"),
          axis.text.x = element_text(size = 10 ),
          axis.title.x = element_text(size = 16 ),
          axis.text.y = element_text(size = 10 ),
          axis.title.y = element_text(size = 16 ),
          plot.title = element_text(size = 16),
          panel.margin.x = unit(0.5, "lines")) +
    scale_x_reverse(breaks=c(1000,0,-1000),
                    labels=c("1000\nCE","0","1000\nBCE")) +
    theme(axis.ticks.length = unit(0.001, "mm")) +
    theme(plot.margin = unit(c(0,1,0.5,0), "cm")) +
    scale_y_continuous(expand = c(0,0)) +
    theme(axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_vline(xintercept = c(1950,1000,0,-1000,-2000), linetype = "longdash",size=0.1) +
    theme(legend.position='bottom',
          legend.justification = 'right',
          legend.text = element_text(size=16))


pdf(paste("figures/PropsAndDatesDensitiesAllRegBothSourcesLanguage.pdf",sep=""),width=12,height=9)
print(p7)
dev.off()

