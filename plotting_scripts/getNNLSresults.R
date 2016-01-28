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
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
leginfo_file <- "data/MalariaGenAdmixturePopulationKey.txt"
## LOAD POPKEY FILE ##
popkey <- read.table(popkey_file,header=T,as.is=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$Ethnic_Group[popkey$Ethnic_Group=="SEMI-BANTU"] <- "SEMI.BANTU"


### COMBINE NNLS RESULTS, INFERRED USING GLOBETROTTER
pops <- popkey$Ethnic_Group
mgpops <- pops[1:48]

mixmat <- matrix(0,nc=length(pops),nr=length(mgpops))
rownames(mixmat) <- mgpops
colnames(mixmat) <- pops

for(pop in mgpops)
{
    infile <- paste0("/mnt/kwiat/well/human/george/globetrotter/props/",pop,".globetrotter.main.nnls.txt")
    mc <- read.table(infile,skip=1,header=T)
    mixmat[pop,colnames(mc)] <- as.numeric(mc)
}

mixmat[mixmat<0.01] <- 0
mixmat <- mixmat/rowSums(mixmat)

write.table(mixmat,"data/Malgen23EthnicGroups1KGNoAmericaFinalAnalysisNonLocalNNLS.txt",
            col.names=T,row.names=T)


