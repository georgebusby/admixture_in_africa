### SCRIPT TO PLOT MALDER RESULTS FOR PAPER -- CAN USE RESULTS FROM DIFFERENT GENETIC MAPS ###
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
## LOAD POPORDER FILE ##
popplot <- scan("/mnt/kwiat/home/popgen/scripts/poplists/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplotorder <- popplot
sim_pops <- scan("/mnt/kwiat/well/human/george/admix_sims/chromopainter/admix_sims.pops",what="char")


## CHOOSE ONLY THE AFRICAN POPS FROM POPPLOT ORDER
popplotorder <- popplotorder[c(1:16,18:49)]
poplist <- factor(popplotorder,levels=popplotorder)
poplist<- as.matrix(poplist)
n_pops <- nrow(poplist)
analys <- c("HAPMAP") #,"CeuMap","AfrMap")
malder_file <- "~/R/Copy/Rprojects/AfricaPOPGEN/manuscript/f3tables/AllPopsMalderFinalAdmixSimsALL.txt"
malder_plot <- read.table(malder_file,header=T, as.is=T)
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
leginfo <- read.table(leginfo_file,header=T,comment.char="", as.is = T)
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)

### DEFINE COLOURS ##
pcolshex <- c("#0000CD", "#03B4CC", "#A65628", "#FF7F00", "#984EA3", "#4DAF4A", "#CCCC00")[c(1,2,4,5,3,6,7)]
ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

src1list <- list(c("JOLA","CEU"),c("JOLA","JUHOAN"),c("JOLA","MALAWI"),
                 c("JUHOAN","CEU"),c("JUHOAN","JOLA"),c("JUHOAN","MALAWI"),
                 c("MALAWI","CEU"),c("MALAWI","JOLA"),c("MALAWI","JUHOAN"))
names(src1list) <- c("JOLACEU","JOLAJUHOAN","JOLAMALAWI",
                     "JUHOANCEU","JUHOANJOLA","JUHOANMALAWI",
                     "MALAWICEU","MALAWIJOLA","MALAWIJUHOAN")

###################################################################
## DEFINE A PLOTTING ORDER FOR THE SIMS
sim_srcs <- sapply(sim_pops,function(x)strsplit(x,split="[0-9]")[[1]][1])
sim_gens <- sapply(sim_pops,function(x)as.numeric(gsub("[A-Z]","",strsplit(x,split="g")[[1]][1])))
sim_props <- sapply(sim_pops,function(x)as.numeric(gsub("[a-z]","",strsplit(x,split="00")[[1]][2])))
sim_src1 <- sapply(sim_srcs,function(x) as.character(unlist(src1list[x])[1]))
sim_src2 <- sapply(sim_srcs,function(x) as.character(unlist(src1list[x])[2]))

sim_tab <- data.frame(cbind(sim_srcs,sim_src1,sim_src2,sim_gens,sim_props))
sim_tab$sim_src1 <- factor(sim_tab$sim_src1,levels=c("JOLA","MALAWI","JUHOAN"))
sim_tab$sim_src2 <- factor(sim_tab$sim_src2,levels=c("CEU","JOLA","MALAWI","JUHOAN"))
sim_pops <- sim_pops[order(sim_tab$sim_src2,sim_tab$sim_src2)]

## DO THE DIFFERENT MAPS INFER THE SAME NUMBER OF EVENTS? ##
pop_vec <- sim_pops #unique(malder_plot$EthnicGroup)
n_pops <- length(pop_vec)

### FOR A GIVEN ANALYSIS, SUBSET THE MAIN DATABASE AND PLOT
for(analy in analys)
{
    ## ORDER pop_vec BY THE DATES AND ANCESTRIES IN THE HAPMAP ##
    malder <- malder_plot[malder_plot$Rec.Map==analy,]
    sim_dates <- c()
    ##################################################################
    ## PLOT DATES AND COMPARISONS
    pdf(paste0("figures/AllPopsMalderTimesMapComps",analy,"AdmixSims.pdf"),height=12,width=9)
        layout(matrix(c(1,1,1,2,3,4),3,2),widths=c(6,3)) #,widths=c(6,3))
        par(mar=c(4,6,3,1))
        x_labs3 <- c(2000,0,-2000,-5000,-10000,-15000)
        x_labs3char <- c("2000\nCE",0,"2000\nBCE","5000\nBCE","10000\nBCE","15000\nBCE")
        x_max <- min(x_labs3)
        plot(0,0,xlim=rev(range(x_labs3)),
             ylim=c(n_pops,1),type="n",axes=F,xlab="",ylab="")
        axis(1,at=x_labs3,labels=x_labs3char,cex.axis=1,tick=F,padj=0.5)
        for(j in x_labs3) abline(v=j,lty=2)
        mtext("Date of Admixture",1,line=3)
        
        ## COLOURS FOR EACH GROUP
        y_ax_cols <- rep("black",length(pop_vec))
#         for(i in pop_vec) 
#         {
#             if(i=="SEMI-BANTU") i <- "SEMI-BANTU"
#             ii <- pcolshex[ancreg_list==as.character(popkey$RegionM[popkey$Ethnic_Group==i])]
#             y_ax_cols <- c(y_ax_cols,ii)
#         }
        
        ## SET THE SCENE
#         for(i in 1:n_pops) axis(2,pos=4000,at=i,
#                                 labels=tidyNames(pop_vec[i],fula=T),
#                                 col.axis=y_ax_cols[i],las=2,tck=0,lwd=0,line=-0.5)
         
        for(i in seq(0.5,(n_pops+1),1)) abline(h=i,lty=3,lwd=0.5)
        ## put thick line between different groups
        pop_vec2 <- as.character(sapply(pop_vec,function(x){strsplit(x,split="[0-9]")[[1]][1]}))
        for(i in unique(pop_vec2)) abline(h=max((1:length(pop_vec))[pop_vec2==i])+0.5,lty=1,lwd=1,xpd=T)
        

        lines(x=c(3200,3200),y=c(0,n_pops+1),lwd=1,xpd=T)
        lines(x=c(4000,4000),y=c(0,n_pops+1),lwd=2,xpd=T)
        
        text(y=-1,x=5000,labels="SIM PROPS",xpd=T,srt=35,adj=0,cex=1)
        text(y=-1,x=4400,labels="SIM SRCS",xpd=T,srt=35,adj=0,cex=1)
        text(y=-1,x=3600,labels="1ST EVENT",xpd=T,srt=35,adj=0,cex=1)
        text(y=-1,x=2800,labels="2ND EVENT",xpd=T,srt=35,adj=0,cex=1)
        text(y=-5,x=5000,labels=LETTERS[1],adj=0,las=1,cex=2,lwd=3,xpd=T)
        ## NOW PLOT DATES AND ANCESTRIES ##
        
        trupropa <- 5200
        trupropb <- 5050
        trupropc <- 4900
        truposa <- 4600
        truposb <- 4200
        ev1posa <- 3800
        ev1posb <- 3400
        ev2posa <- 3000
        ev2posb <- 2600
        evcex <- 1.4

        for(i in (1:length(pop_vec)))
        {
            pop1 <- as.character(pop_vec)[i]
            srcs <- strsplit(pop1,split="[0-9]")[[1]][1]
            sgens <- as.numeric(gsub("[A-Z]","",strsplit(pop1,split="g")[[1]][1]))
            sprops <- as.numeric(gsub("[a-z]","",strsplit(pop1,split="00")[[1]][2]))
            src1 <- as.character(unlist(src1list[srcs])[1])
            src2 <- as.character(unlist(src1list[srcs])[2])
            srcanc1 <- getPopRegion(tidyNames(src1,fula=T),popkey)
            srcanc2 <- getPopRegion(tidyNames(src2,fula=T),popkey)
            
            ## PLOT TRUE ANCESTRIES
            pcol1 <- pcolshex[ancreg_list==srcanc1]
            pcol2 <- pcolshex[ancreg_list==srcanc2]
            
            points(truposa,i,pch=15,col=pcol1,xpd=T,cex=evcex)
            points(truposb,i,pch=15,col=pcol2,xpd=T,cex=evcex)
            
            ## PLOT TRUE PROPS
            points(trupropa,i,pch=15,col="black",xpd=T,cex=1)
            if(sprops %in% c(90,95)) points(trupropb,i,pch=15,col="black",xpd=T,cex=1)
            if(sprops %in% c(95)) points(trupropc,i,pch=15,col="black",xpd=T,cex=1)
            
            if(pop1%in%malder$EthnicGroup)
            {
                
                tt <- malder[malder$EthnicGroup==pop1,]
                tt <- tt[order(tt$Date.Gens),]
                res_tmp <- nrow(tt)
                if(nrow(tt)>0)
                {
                    for(j in 1:nrow(tt))
                    {
                        date <- tt$Date.Gens[j]
                        dateci <- tt$Date.Gens.CI[j]
                        donpop1 <- strsplit(tt$Test.Pops[j],split="\\;")[[1]][1]
                        donanc1 <- getPopRegion(tidyNames(donpop1,fula=T),popkey)
                        donpop2 <- strsplit(tt$Test.Pops[j],split="\\;")[[1]][2]
                        donanc2 <- getPopRegion(tidyNames(donpop2,fula=T),popkey)
                        mainanc <- tt$Main.Anc[j]
                        mainancp <- tt$Main.Anc.p[j]
                    
                        ## switch to match the order of the true ancestries
                        if(donanc1 == srcanc2 | donanc2 == srcanc1)
                        {
                            donanc1 <- getPopRegion(tidyNames(donpop2,fula=T),popkey)
                            donanc2 <- getPopRegion(tidyNames(donpop1,fula=T),popkey)
                        }
                        
                        ## PLOT DATE
                        ## plot true date
                        points(makeDate(sgens,add_BCE=F),i,pch=21,bg="red", col="red", cex=1)
                        
                        arrows(makeDate(date-dateci,add_BCE=F),i,
                               makeDate(date+dateci,add_BCE=F),i,
                               angle=90,length=0.025,code=3)
                        pcol1 <- pcolshex[ancreg_list==donanc1]
                        pcol2 <- pcolshex[ancreg_list==donanc2]
                        pnt_bg <- "black"
                        if(mainancp>0.001) pnt_bg <- "grey"
                        if(mainancp>0.05) pnt_bg <- "white"
                        points(makeDate(date,add_BCE=F),i,pch=21,bg=pnt_bg, cex=1.2)
                        
                        sim_dates <- rbind(sim_dates,
                                           c(src1,src2,sgens,sprops,
                                             date,dateci,donpop1,donpop2,mainancp,j))
                        
                                                
                        ## PLOT EVENT ANCESTRIES
                        if(j == 1)
                        {
                            points(ev1posa,i,pch=15,col=pcol1,xpd=T,cex=evcex)
                            points(ev1posb,i,pch=15,col=pcol2,xpd=T,cex=evcex)
    #                         if(donanc1 == mainanc) points(ev1posa,i,pch=22,xpd=T,cex=1,lwd=0.5)
    #                         if(donanc2 == mainanc) points(ev1posb,i,pch=22,xpd=T,cex=1,lwd=0.5)
                        }
                        
                        if(j == 2)
                        {    points(ev2posa,i,pch=15,col=pcol1,xpd=T,cex=evcex)
                             points(ev2posb,i,pch=15,col=pcol2,xpd=T,cex=evcex)
    #                          if(donanc1 == mainanc) points(ev2posa,i,pch=22,xpd=T,cex=1,lwd=0.5)
    #                          if(donanc2 == mainanc) points(ev2posb,i,pch=22,xpd=T,cex=1,lwd=0.5)
                        }     
                    }
                }
            }
        }
    
        ## PLOT LEGEND 
        par(mar=c(1,1.5,3,0))
        plot(0,0,axes=F,xlab="",ylab="",type="n")
        legend_text <- c("West African Niger-Congo",
                         "Central West African Niger-Congo",
                         "East African Niger-Congo",
                         "South African Niger-Congo",
                         "Nilo-Saharan / Afro-Asiatic",
                         "KhoeSan",
                         "Eurasia",
                         "high confidence date (P<0.001)",
                         "low confidence date (P<0.05)",
                         "non-significant date (P>0.05)")
        l <- legend("top",legend=legend_text,pch=c(rep(15,7),21,21,21),
                    col=c(pcolshex[c(1:3,5,4,6,7)],"black","black","black"),bty="n",
                    pt.bg=c(pcolshex[c(1:3,5,4,6,7)],"black","grey","white"),
                    ncol=1,xpd=T,pt.cex=3,x.intersp=1,y.intersp=1.25,
                    pt.lwd=2,cex=1.25, title="Ancestry Region")
    
        ## PLOT DATE COMPARISONS
        colnames(sim_dates) <- c("src1","src2","true.date","prop","sim.date","sim.date.ci","sim.src1","sim.src2","p","event")
        sim_dates <- data.frame(sim_dates)
                
        ## FOR POPS WHERE TWO DATES ARE INFERRED ONLY USE THE ONE THAT
        ## IS CLOSEST TO THE SIM DATE
        
        dups <- sim_dates[duplicated(sim_dates[,1:4]),]
        for(i in 1:nrow(dups))
        {
            test <- sim_dates$src1==dups$src1[i]&sim_dates$src2==dups$src2[i] &sim_dates$true.date==dups$true.date[i]&sim_dates$prop==dups$prop[i]
            tmp <- sim_dates[test,]
            tgen <- as.numeric(as.character(unique(sim_dates$true.date[test])))
            tgens <- as.numeric(as.character(sim_dates$sim.date[test]))
            ## remove whichever date is furthest away from gens
            to_keep <- which.min(abs(tgens-tgen))
            sim_dates[test,] <- tmp[rep(to_keep,2),]
        }
        sim_dates <- sim_dates[!duplicated(sim_dates),]
    
        ## PLOT 
        xs <- sapply(as.numeric(as.character(sim_dates$true.date)),makeDate,add_BCE=F) 
        ys <- sapply(as.numeric(as.character(sim_dates$sim.date)),makeDate,add_BCE=F)

        par(mar=c(6,5,2,1))
        plot(xs,ys,xlab="",ylab="",
             pch=20,xlim=c(2000,-12000),ylim=c(-12000,2000),
             axes=F,type="n")
        x_labs <- abs(round(unique(xs),-1))
#         for(i in x_labs)
#         {
#             if(nchar(i)>3)
#             {
#                 ii <- paste(substr(i,1,1),substr(i,2,nchar(i)),sep=",")
#             x_labs[x_labs == i] <- ii
#             }
#         }
        x_labs <- paste(x_labs,"BCE")
        x_labs <- paste(x_labs,c("[100 gens]","[200 gens]","[300 gens]","[400 gens]"),sep="\n")
        
        for(i in 1:4)
        {
            text(y=-13000,x=unique(xs)[i],labels=x_labs[i],xpd=T,srt=315,adj=c(0,1),cex=1)
        }
        mtext(1,text="Simulated Date",line=5)
        mtext(2,text="MALDER Date Inference",line=2.5)
        y_at <- c(2000,0,-2000,-5000,-10000)
        y_at_labs <- c("2000\nCE","0","2000\nBCE","5000\nBCE","10000\nBCE")
        axis(2,line=0.5,at=y_at,labels=y_at_labs,las=2,hadj=0.5,padj=0.5,lwd=0)
        abline(h=y_at,lty=2)
        abline(a=0,b=1,col="red",lwd=2)
        box()
        
        box_width <- 200
        for(box_date in unique(xs))
        {
            bp<- boxplot(ys[xs==box_date]~xs[xs==box_date],plot=F,
                         boxwex=0.5,add=T,at=box_date,axes=F,col="blue")
            ## add box
            rect(xleft=box_date+box_width,
                 ybottom=bp$stats[2],
                 xright=box_date-box_width,
                 ytop=bp$stats[4],lwd=1.5)
            ## add median
            lines(x=c(box_date+box_width,box_date-box_width),
                  y=c(bp$stats[3],bp$stats[3]),lwd=1.5,col="blue",lend=2)
            ## add whiskers
            lines(x=c(box_date,box_date),
                  y=c(bp$stats[5],bp$stats[4]))
            lines(x=c(box_date,box_date),
                  y=c(bp$stats[1],bp$stats[2]))
        }
        
        points(xs,ys,pch=20)
        mtext(3,line=0,text=LETTERS[2],adj=0,las=1,cex=1.5,lwd=3,xpd=T)


    dev.off()                
}


