## SCRIPT TO TAKE THE MALDER TEXT FILE AND GENERATE A TEX FORMAT FILE
### SCRIPT TO PLOT GLOBETROTTER RESULTS ON MAP ###
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
library(xtable)

setwd(paste0(main_dir,"popgen/"))

popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$Ethnic_Group[popkey$Ethnic_Group=="SEMI-BANTU"] <- "SEMI.BANTU"

mald <- read.table("~/R/Copy/Rprojects/AfricaPOPGEN/manuscript/f3tables/AllPopsMalderFinalALLmindis.txt",
                   header=T,as.is=T)
plot_colmns <- c("EthnicGroup","amplitude","Date.Gens","Date.Gens.CI",
                 "pop1",grep("tpop1",colnames(mald),value=T),
                 "pop2",grep("tpop2",colnames(mald),value=T),"N.evs")

ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )
reg_labs <- c("West Africa NC","Central West Africa NC","East Africa NC",
              "Nilo-Saharan","Afroasiatic","South Africa NC",
              "Khoesan","Eurasia")
regtab <- cbind(ancreg_list,reg_labs)

for(analy in c("HAPMAP","AfrMap","CeuMap"))
{
    mtmp <- mald[mald$Rec.Map==analy,plot_colmns]
    mtmp[mtmp$N.evs==0,2:ncol(mtmp)] <- NA
    mtmp$EthnicGroup <- tidyNames(mtmp$EthnicGroup,fula=T,khoesan=T)
    mtmp$pop1 <- tidyNames(mtmp$pop1,fula=T,khoesan=T)
    mtmp$pop2 <- tidyNames(mtmp$pop2,fula=T,khoesan=T)
    mtmp$Date.Gens <- sapply(round(mtmp$Date.Gens),as.character)
    mtmp$Date.Gens.CI <- sapply(round(mtmp$Date.Gens.CI),as.character)
    mtmp$amplitude <- sapply(signif(mtmp$amplitude,2),as.character)
    newreg <- c()
    for(i in mtmp$tpop1)
    {
        if(is.na(i))
        {
            newreg <- c(newreg,NA)
        } else
        {
            newreg <- c(newreg,as.character(regtab[regtab[,1]==i,2]))
        }
    }    
    mtmp$tpop1 <- newreg
    newreg <- c()
    for(i in mtmp$tpop2)
    {
        if(is.na(i))
        {
            newreg <- c(newreg,NA)
        } else
        {
            newreg <- c(newreg,as.character(regtab[regtab[,1]==i,2]))
        }
    }   
    mtmp$tpop2 <- newreg
    
    ctmp <- colnames(mtmp)
    ctmp <- gsub("\\.","\\-",ctmp)
    ctmp <- gsub("tpop1","Pop1Anc",ctmp)
    ctmp <- gsub("tpop2","Pop2Anc",ctmp)
    ctmp <- gsub("pop1","Pop1",ctmp)
    ctmp <- gsub("pop2","Pop2",ctmp)
    for(j in 1:nrow(regtab)) ctmp <- gsub(regtab[j,1],paste(regtab[j,2],"(Z)",sep=""),ctmp)
    colnames(mtmp) <- ctmp
    
    mtmp <- mtmp[,1:(ncol(mtmp)-1)]
    
    mtmp <- mtmp[!duplicated(mtmp),]
    
    mtmp <- xtable(mtmp,digits=2)
    
    ## PRINT TO TABLE
    newlines <- c()
    for(i in (2:nrow(mtmp))) if(mtmp$EthnicGroup[i]!=mtmp$EthnicGroup[(i-1)]) newlines <- c(newlines,i)
    #newlines <- sort(c(newlines,newlines - 1))
    #newlines <- newlines[newlines != 0]
    
    addtorow          <- list()
    addtorow$pos      <- list()
    addtorow$pos[[1]] <- c(-1)
    addtorow$pos[[2]] <- c(newlines-1)
    
    addtorow$command[[1]]  <- "all results"
    addtorow$command[[2]] <- "\\hline \n"
    print(mtmp, floating=FALSE,
          tabular.environment="longtable", 
          comment=FALSE,
          include.colnames=TRUE,
          include.rownames=FALSE,
          caption.placement="top",
          file=paste0("data/AllPopsMalderMinDis",analy,".tex"),
          booktabs=TRUE,add.to.row=addtorow)
}

