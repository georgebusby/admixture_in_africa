#########################################################################
### SCRIPT TO PULL OUT GLOBETROTTER RESULTS FROM MULTIPLE OUTPUTS     ###
setwd("~/repos/popgen/")

############################################################
## SOURCE SOME USEFUL FUNCTIONS FROM copyselection PACKAGE ##
## ~~~~~~~~~~~       !!! IN DEVELOPMENT !!!     ~~~~~~~~~~ ##
############################################################
main_dir <- "~/repos/" ## where the code is kept
source(paste0(main_dir,"popgen/packages_dev/functionWriter.R"))
source(paste0(main_dir,"popgen/packages_ext/FinestructureLibrary.R"))
source(paste0(main_dir,"popgen/packages_ext/FinestructureLibrary_GB.R"))
library("xtable")
       
###########################################################
## DEFINE VARIABLES AND FILES
setwd(paste0(main_dir,"popgen/"))
popkey_file <- "data/MalariaGenAdmixturePopulationOverviewNSAA.txt"
popkey <- read.table(popkey_file,header=T)
popkey$Ethnic_Group <- toupper(popkey$Ethnic_Group)
popkey$Ethnic_Group[popkey$Ethnic_Group=="SEMI-BANTU"] <- "SEMI.BANTU"
popplot <- scan("data/MalariaGen23EthnicGroups1KGSouthAfricaFinalAnalysisPopsSSAonlyOrder.txt",what="char")
popplot <- popplot[popplot!="SEMI.BANTU"]

ancreg_list <- c("Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                 "East_Africa_Niger-Congo","East_Africa_Nilo-Saharan","East_Africa_Afro-Asiatic",
                 "South_Africa_Niger-Congo","South_Africa_KhoeSan","Eurasia" )

in_dir <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/globetrotter/input/"
res_dir <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/globetrotter/props/"
dates_dir <- "/mnt/kwiat/well/human/george/globetrotter/dates/"
param_dir <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/globetrotter/paramfiles/"
analysis <- "MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinal"
c_suffix <- ".globetrotter.main_curves.txt"
p_suffix <- ".globetrotter.main.txt"

##########################################################
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
##########################################################
## DEFINE VARIABLES TO STORE MAIN RESULTS
resmat <- resmatnull <- resmatalt <- matrix(NA,nrow=0,ncol=22)
dateboots <- matrix(NA,nrow=0,ncol=4)
date2boots <- matrix(NA,nrow=0,ncol=4)
datebootsnull <- matrix(NA,nrow=0,ncol=4)
## DEFINE MATRICES TO STORE SOURCE COPYING VECTORS
admixturesources2 <- matrix(nrow=0,ncol=length(popplot))
colnames(admixturesources2) <- popplot
admixturesources3 <- admixturesources4 <- matrix(nrow=0,ncol=length(popplot)+20)

## DEFINE ALTERNATIVE RUNS -- THIS IS WHERE I HAVE SUCCESSIVELY REMOVED DONORS
altruns <- c("nolocal","nolocalnomalawi","nolocalnowest","nolocalnosouth","nolocalnoeast")
##########################################################
## GET RESULTS
for(pop in popplot[1:48])
{
    print(paste("getting GLOBETROTTER reults for:",pop))
    pop <- gsub("-",".",pop)
    pop1 <- pop
    ## NEED TO USE SEMI.BANTUII, WHERE I DISALLOWED BANTU FROM COPYING
    if(pop1 == "SEMI.BANTU" ) pop1 <- "SEMI.BANTUII"
    reg <- as.character(popkey$RegionM[popkey$Ethnic_Group==pop])
    if(reg == "East_Africa_Afro-Asiatic") reg <- "East_Africa_Nilo-Saharan"
    nlength <- read.table(paste0(in_dir,analysis,reg,".chunklengths.out"),header=T,row.names=1)
    nlength <- rowsAsMapClusts(final_clusts2,nlength,stat=mean)
    predmat <- nlength/rowSums(nlength)
    donors <- gsub("-",".",popplot)
    a <- "main"
    infile <- paste0(res_dir,pop1,".globetrotter.",a,".txt")
    infile_null <- paste0(res_dir,pop1,".globetrotter.null.txt")
    infile_dates <- paste0(dates_dir,pop1,".globetrotter.",a,".dates.bt")
    infile_2dates <- paste0(dates_dir,pop1,".globetrotter.",a,".2dates.bt")
    infile_dates_null <- paste0(dates_dir,pop1,".globetrotter.null.dates.bt")
    infile_2dates_null <- paste0(dates_dir,pop1,".globetrotter.null.2dates.bt")
    if(file.exists(infile))
    {
        ## GET MAIN DATA
        ll <- getGlobetrotter(infile,donors)
        ## GET ONE DATE BOOTSTRAPS
        ld <- getGlobetrotterDates(infile_dates)
        ## GET TWO DATE BOOTSTRAPS
        ld2 <- getGlobetrotterDates(infile_2dates)
        
        if(nrow(ld)>0) dateboots <- rbind(dateboots,cbind(pop,"main",ld))
        if(nrow(ld2)>0) date2boots <- rbind(date2boots,cbind(pop,"main",ld2))
        ## ADD TO THE MAIN RESULTS TABLE
        resmat <- rbind(resmat,c(pop,ll[[1]],"main"))
        
        ## GET THE RESULTS OF THE NULL RUN
        if(file.exists(infile_null))
        {
            lln <- getGlobetrotter(infile=infile_null,donors)
            resmat <- rbind(resmat,c(pop,lln[[1]],"main.null"))
            ld <- getGlobetrotterDates(infile_dates_null)
            if(nrow(ld)>0)
            {
                ld <- cbind(pop,"main.null",ld)
                colnames(ld) <- colnames(dateboots)
                dateboots <- rbind(dateboots,ld)
            }
            ld2 <- getGlobetrotterDates(infile_2dates_null)
            if(nrow(ld2)>0)
            {
                ld2 <- cbind(pop,"main.null",ld2)
                colnames(ld2) <- colnames(date2boots)
                date2boots <- rbind(date2boots,ld2)
            }
            
        }
        ## PULL OUT THE SOURCES AND GENERATE COPYING VECTORS USING PREDMAT
        if(!is.null(ll[[2]]))
        {
            finn <- ll[[2]]
            ## REORDER finn TO MATCH PREDMAT; THERE ARE 8 SOURCES THAT TWO FOR EACH TYPE OF EVENT
            cvs <- finn[,2:ncol(finn)]
            rownames(cvs) <- paste(paste(pop,a,sep="."),1:8,sep=".")
            ## ADD TO THE MAIN MATRIX
            while(ncol(cvs)!=length(popplot)) cvs <- cbind(cvs,0)
            colns <- gsub("-",".",popplot)
            colnames(cvs)[colnames(cvs) == ""] <- colns[!colns%in%colnames(cvs)]
            admixturesources2 <- rbind(admixturesources2,cvs[,colns])
            lln_tmp <- c()
            for(k in 1:nrow(cvs))
            {
                lln_tmp <- rbind(lln_tmp,ll[[1]])  
            }
            admixturesources3 <- rbind(admixturesources3,cbind(lln_tmp,cvs[,colns]))
        }
        ## ADD NULL SOURCES
        if(!is.null(lln[[2]]))
        {
            finn <- lln[[2]]
            ## REORDER finn TO MATCH PREDMAT; THERE ARE 8 SOURCES THAT TWO FOR EACH TYPE OF EVENT
            cvs <- finn[,2:ncol(finn)]
            rownames(cvs) <- paste(paste(pop,"null",sep="."),1:8,sep=".")
            ## ADD TO THE MAIN MATRIX
            while(ncol(cvs)!=length(popplot)) cvs <- cbind(cvs,0)
            colns <- gsub("-",".",popplot)
            colnames(cvs)[colnames(cvs) == ""] <- colns[!colns%in%colnames(cvs)]
            admixturesources2 <- rbind(admixturesources2,cvs[,colns])
            lln_tmp <- c()
            for(k in 1:nrow(cvs))
            {
                lln_tmp <- rbind(lln_tmp,lln[[1]])  
            }
            admixturesources3 <- rbind(admixturesources3,cbind(lln_tmp,cvs[,colns]))
        }
        
        
        ## GET RESULTS FOR ALTERNATIVE RUNS
        for(j in altruns)
        {
            if(j=="nolocalnoeast" & pop == "MALAWI")
            {
                infile2 <- paste0(paste0(res_dir,pop1,j,".globetrotter.null.txt"))
            } else
            {
                infile2 <- paste0(paste0(res_dir,pop1,j,".globetrotter.main.txt"))
            }
            infile_dates2 <- paste0(dates_dir,pop1,j,".globetrotter.dates.bt")
            if(file.exists(infile2))
            {
                lln <- getGlobetrotter(infile2,donors)
                resmatalt <- rbind(resmatalt,c(pop,lln[[1]],j))
                if(!is.null(lln[[2]]))
                {
                    finn <- lln[[2]]
                    ## reorder finn to have same cols as predmat
                    cvs <- finn[,2:ncol(finn)]
                    rownames(cvs) <- paste(paste(pop,j,sep="."),1:8,sep=".")
                    ## add to admixture sources results matrix
                    colns <- gsub("-",".",popplot)
                    while(ncol(cvs)!=length(colns)) cvs <- cbind(cvs,0)
                    colnames(cvs)[colnames(cvs) == ""] <- colns[!colns%in%colnames(cvs)]
                    lln_tmp <- c()
                    for(k in 1:nrow(cvs))
                    {
                        lln_tmp <- rbind(lln_tmp,lln[[1]])  
                    }
                    admixturesources4 <- rbind(admixturesources4,cbind(lln_tmp,cvs[,colns]))
                    
                }
            }
        }
    }
}

resmat <- data.frame(resmat,stringsAsFactors = FALSE)
resmat$result <- toupper(resmat$result)

#############################################################
## NOW WE HAVE OUR RESULTS TABLE AND WE CAN WORK WITH IT TO 
## GENERATE OUR FINAL RESULTS
#############################################################
## GEN ADMIXTURE P-VALUE BASED ON NULL BOOTSTRAPS
resmat$pval <- 1
resmat$nump <- 0
nullruns <- grep("null",paste(dateboots$pop,dateboots[,2]))
nullruns <- dateboots[nullruns,]
for(i in as.character(unique(nullruns$pop)))
{
    bb <- nullruns$date1.est.boot[nullruns$pop==i]
    rowindex <- paste(resmat$V1,resmat$V22,sep=".") == paste(i,"main",sep=".")
    resmat$pval[rowindex] <- sum(bb>400 | bb < 1) / length(bb)
    resmat$nump[rowindex] <- length(bb)
    rowindex <- paste(resmat$V1,resmat$V22,sep=".") == paste(i,"main.null",sep=".")
    resmat$pval[rowindex] <- sum(bb>400 | bb < 1) / length(bb)
}

#############################################################
## GENERATE A FINAL TABLE OF RESULTS
restypes <- c("ONE-DATE","ONE-DATE-MULTIWAY","MULTIPLE-DATES","UNCERTAIN" ,"NO-ADMIXTURE" ,"UNKNOWN")
numericcols <- c("maxR2fit.1date","fit.quality.1event","fit.quality.2events",
                 "maxScore.2events","proportion.date1.source1","proportion.date2.source1")
finaltable <- matrix(NA,nrow=0,ncol=ncol(resmat))
for(i in restypes)
{
    ii <- resmat[resmat$result==i,]
    if(nrow(ii) > 0)
    {
        ii[,numericcols] <- apply(signif(apply(ii[,numericcols],2,as.numeric),3),2,as.character)
        ii <- convert.factors.to.strings.in.dataframe(ii)
        ii <- ii[order(ii$bestmatch.event1.source1,ii$gen.1date,decreasing=T),]
    }
    finaltable <- rbind(finaltable,ii)
}

finaltable <- cbind(finaltable[,c(1,22)],finaltable[,c(2:21,23:(ncol(finaltable)))])
colnames(finaltable) <- c("Cluster","Analysis","Result","Date","alpha",
                          "max(R1)","FQ1","FQ2","best.source1","best.source2",
                          "alpha2","best.source1.ev2","best.source2.ev2",
                          "Date.2a","Date2b","maxScore2ev","alpha2.date1",
                          "best.source1.date1","best.source2.date1",
                          "alpha2.date2","best.source1.date2","best.source2.date2","pval","nump")
#############################################################
## LOOK AT DATE BOOTSTRAPS AND GENERATE CONFIDENCE INTERVALS
dateboots2 <- dateboots3 <- c()
for(i in unique(dateboots$pop))
{
    for(analy in c("main","main.null"))
    {
        ds <- round(dateboots$date1.est.boot[dateboots$pop==i & dateboots[,2] == analy])
        ii <- sapply(quantile(ds,c(0.975,0.025)),makeDate)
        dateboots2 <- rbind(dateboots2,c(i,analy,ii))
        ii <- round(quantile(ds,c(0.025,0.975)))
        dateboots3 <- rbind(dateboots3,c(i,analy,ii))
    }
}
#############################################################
## LOOK AT TWO DATE BOOTSTRAPS AND GENERATE CONFIDENCE INTERVALS
date2boots2 <- date2boots3 <- c()
for(i in unique(date2boots$pop))
{
    for(analy in c("main","main.null"))
    {
        ds <- round(date2boots$date1.est.boot[date2boots$pop==i & date2boots[,2] == analy])
        ds1 <- round(date2boots$date1.est.boot[date2boots$pop==i& date2boots[,2] == analy])
        ds2 <- round(date2boots$date2.est.boot[date2boots$pop==i& date2boots[,2] == analy])
        ii <- sapply(quantile(ds1,c(0.975,0.025)),makeDate)
        iii <- sapply(quantile(ds2,c(0.975,0.025)),makeDate)
        date2boots2 <- rbind(date2boots2,c(i,analy,ii,iii))
        ii <- round(quantile(ds1,c(0.025,0.975)))
        iii <- round(quantile(ds2,c(0.025,0.975)))
        date2boots3 <- rbind(date2boots3,c(i,analy,ii,iii))
    }
}

finaltable <- convert.factors.to.strings.in.dataframe(finaltable)
finaltable$Date.CI <- finaltable$dateL <- finaltable$dateH <- NA

#############################################################
## FILL IN DATES AND CIS IN THE MAIN RESULTS TABLE

rowindex <- match(paste(dateboots2[,1],dateboots2[,2],sep="."),
                      paste(finaltable$Cluster,finaltable$Analysis,sep="."),nomatch=0)
finaltable$Date.CI[rowindex] <- paste0("(",dateboots2[,3],"-",dateboots2[,4],")")
finaltable$dateL[rowindex] <- dateboots3[,3]
finaltable$dateH[rowindex] <- dateboots3[,4]
    
finaltable$Date2a.CI <- finaltable$date2aL <- finaltable$date2aH <- finaltable$Date2b.CI <- finaltable$date2bL <- finaltable$date2bH <- NA
rowindex <- match(paste(date2boots2[,1],date2boots2[,2],sep="."),
                      paste(finaltable$Cluster,finaltable$Analysis,sep="."),nomatch=0)
    
finaltable$Date2a.CI[rowindex] <- paste0("(",date2boots2[,3],"-",date2boots2[,4],")")
finaltable$date2aL[rowindex] <- date2boots3[,3]
finaltable$date2aH[rowindex] <- date2boots3[,4]
finaltable$Date2b.CI[rowindex] <- paste0("(",date2boots2[,5],"-",date2boots2[,6],")")
finaltable$date2bL[rowindex] <- date2boots3[,5]
finaltable$date2bH[rowindex] <- date2boots3[,6]
#############################################################
## ALSO ROUND DATE COLUMN
finaltable$Date <- round(as.numeric(finaltable$Date))
finaltable$Date.2a <- round(as.numeric(finaltable$Date.2a))
finaltable$Date2b <- round(as.numeric(finaltable$Date2b))
#############################################################
### MAKE SOME PRETTY RESULTS TABLES AND STORE FOR PLOTTING
allcols <- c("Cluster","Analysis","pval","maxScore2ev","max.R1.","FQ1","FQ2","Result")
onedatecols <- c(allcols,"Date","dateL","dateH","alpha","best.source1","best.source2")
onemultcols <- c(onedatecols,"alpha2","best.source1.ev2","best.source2.ev2") 
twodatecols <- c(allcols,"Date.2a","date2aH","date2aL","best.source1.date1","best.source2.date1","alpha2.date1","Date2b","date2bH","date2bL","best.source1.date2","best.source2.date2","alpha2.date2")

#############################################################
#############################################################
#############################################################
## NOW ALSO DEFINE NO-ADMIXTURE WHEN THE REDUCTION IN THE NULL
## INFERENCE R^2 IS GREATER THAN 1/3
res.tabA <- convert.factors.to.strings.in.dataframe(finaltable) #[,twodatecols]
## DEFINE NO ADMIXTURE FROM P-VALUES
res.tabA$Result[res.tabA$pval>0.05] <- "NO-ADMIXTURE"
nullruns <- grep("null",res.tabA$Analysis)
mainruns <- grep("null",res.tabA$Analysis,invert=T)
mall <- rep(NA,nrow(res.tabA))
for(i in mainruns)
{
    pop <- res.tabA$Cluster[i]
    analy <- res.tabA$Analysis[i]
    res <- res.tabA$Result[i]
    r1 <- as.numeric(res.tabA[,"max(R1)"][res.tabA$Cluster==pop&res.tabA$Analysis==analy])
    r1n <- as.numeric(res.tabA[,"max(R1)"][res.tabA$Cluster==pop&res.tabA$Analysis==paste(analy,"null",sep=".")])
    if(length(r1n)>0 & res !="UNKNOWN")
    {
        m <- r1n/r1 # tests again whether multiple dates is correct inference
        mall[i] <- unique(m) #if(m<(2/3)) 
    }
}
## NOW WE SWITCH TO UNCERTAIN ANY EVENT WHERE MAX(R1) < 0.175 OR 
## NULL R^2 IS REDUCED BY MORE THAN 2/3
res.tabA$Result[which(as.numeric(res.tabA[,"max(R1)"])<0.175 | mall<(1/3))] <- "UNCERTAIN"
#############################################################
#############################################################
## ORDER RESULTS BY EVENT TYPE AND DATE
res.tabAll <- res.tabA
rt1 <- res.tabA[!res.tabA$Result%in%c("UNCERTAIN","NO-ADMIXTURE"),]
rt2 <- res.tabA[res.tabA$Result%in%c("UNCERTAIN","NO-ADMIXTURE"),]
res.ord <- rt1$best.source1.date1
rt1 <- rt1[order(res.ord,as.numeric(rt1$Date.2a)),]
res.ord <- c(res.ord,rep("XXX",nrow(rt2)))
res.tabA <- rbind(rt1,rt2)
res.res <- res.tabA$Result
res.res <- gsub("ONE-DATE-MULTIWAY","1MW",res.res)
res.res <- gsub("ONE-DATE","1D",res.res)
res.res <- gsub("MULTIPLE-DATES","2D",res.res)
res.res <- gsub("UNCERTAIN","U",res.res)
res.res <- gsub("NO-ADMIXTURE","NA",res.res)
res.res <- gsub("UNKNOWN","NA",res.res)
res.tabA$Result <- res.res
res.tabA[,"maxScore2ev"] <- round(as.numeric(as.character(res.tabA[,"maxScore2ev"])),2)
res.tabA[,"pval"] <- round(as.numeric(as.character(res.tabA[,"pval"])),2)


#############################################################
## SWITCH RESULT IF MULTIPLE DATE ARE UNREASONABLE ie CI IS LESS THAN 2
res.tabA$FinalResult <- res.tabA$Result
res.tabA$FinalAnaly <- "main"
min_gens <- 0
## for multiple dates, if either the point estimate or CI include 3
## don't bother with nulls
## then switch result

## IF TWO DATES:
## CHECK IF FIRST EVENT HAS SMALL CI --
##   IF THIS <= 5, THEN USE SECOND EVENT ONLY?
## IF NULL NULL AND MAIN RUN GIVE DIFFERENT RESULTS:
## FOR KHOESAN : USE NULL RUN AN FORCE SINGLE DATE
##      JUSTIFICATION IS THAT THE 2D HS CI INCLUDING 1
pops <- popkey$Ethnic_Group[popkey$RegionM=="South_Africa_KhoeSan"]
test1mw <- res.tabA$Cluster%in%pops & res.tabA$Result == "2D" & res.tabA$FQ1 <= 0.975
res.tabA$FinalResult[test1mw] <- "1MW"
test1mw <- res.tabA$Cluster%in%pops & res.tabA$Result == "2D" & res.tabA$FQ1 > 0.975
res.tabA$FinalResult[test1mw] <- "1D"
res.tabA$FinalAnaly[res.tabA$Cluster%in%pops] <- "main.null"

## FOR SOUTH NC: USE NULL RUN AND FORCE SINGLE DAT
##       SEBANTU AND AMAXHOSA -- REMOVE OLD DATE
pops <- popkey$Ethnic_Group[popkey$RegionM=="South_Africa_Niger-Congo"]
test1mw <- res.tabA$Cluster%in%pops & res.tabA$Result == "2D" & res.tabA$FQ1 <= 0.975
res.tabA$FinalResult[test1mw] <- "1MW"
test1mw <- res.tabA$Cluster%in%pops & res.tabA$Result == "2D" & res.tabA$FQ1 > 0.975
res.tabA$FinalResult[test1mw] <- "1D"
res.tabA$FinalResult[res.tabA$Cluster=="SWBANTU"] <- "2D"
res.tabA$FinalAnaly[res.tabA$Cluster%in%pops] <- "main.null"

## FROR AFROASIATIC SPEAKERS, USE NULL
pops <- popkey$Ethnic_Group[popkey$RegionM=="East_Africa_Afro-Asiatic"]
res.tabA$FinalAnaly[res.tabA$Cluster%in%pops] <- "main.null"

## FOR GAMBIA/MALI : FORCE TWO DATE, BUT ONLY REPORT SECOND DATE
pops <- popkey$Ethnic_Group[popkey$RegionM=="Western_Africa_Niger-Congo"]
res.tabA$FinalResult[res.tabA$Cluster%in%pops] <- "1D"

## SWAP SECOND DATE COLUMNS FOR FIRST DATE
temp <- res.tabA[res.tabA$Cluster%in%pops,]
date1cols <- c("Date","alpha","best.source1","best.source2","dateH","dateL","Date.CI")
date2cols <- c("Date2b","alpha2.date2","best.source1.date2","best.source2.date2","date2bH","date2bL","Date2b.CI")
res.tabA[res.tabA$Cluster%in%pops,date1cols] <- temp[,date2cols]
res.tabA[res.tabA$Cluster%in%pops,date2cols] <- temp[,date1cols]

###################################################
# matching <- c()
# for(i in unique(res.tabA$Cluster))
# {
#     if(length(unique(res.tabA$Result[res.tabA$Cluster==i])) != 1) matching <- c(matching,i)
# }
# 
# res.tabA[res.tabA$Cluste%in%matching,1:10]
###################################################
## THINK ABOUT COMPARING TWO DATE INF TO ONE DATE
## WHAT ARE THE TWO DATES?
## HOW DOES THE SECOND DATE RELATE TO FIRST
## WHAT ARE THE SOURCES? ARE THEY THE SAME FOR BOTH EVENTS?

test <- res.tabA$Result=="2D"
tmp.res <- res.tabA[test==T,]
tmp.res <- tmp.res[grep("null",tmp.res$Analysis,invert=T),]
test <- (as.numeric(tmp.res$date2aL)>min_gens)
test[is.na(test)] <- F
tmp.res <- tmp.res[test==F,]
if(nrow(tmp.res)>0)
{
    tmp.res$FinalResult <- "1D"
    tmp.res$FinalResult[tmp.res$FQ1 <= 0.975] <- "1MW"
    res.tabA[rownames(res.tabA)%in%rownames(tmp.res),] <- tmp.res
}

## then generate a final table will all results
## plus the final results table 
colnames(res.tabA)[colnames(res.tabA)=="Date"] <- "date.1D"
colnames(res.tabA)[colnames(res.tabA)=="dateL"] <- "date.1D.L"
colnames(res.tabA)[colnames(res.tabA)=="dateH"] <- "date.1D.H"
colnames(res.tabA)[colnames(res.tabA)=="Date.2a"] <- "date.2D.1"
colnames(res.tabA)[colnames(res.tabA)=="date2aH"] <- "date.2D.1.H"
colnames(res.tabA)[colnames(res.tabA)=="date2aL"] <- "date.2D.1.L"
colnames(res.tabA)[colnames(res.tabA)=="Date2b"] <- "date.2D.2"
colnames(res.tabA)[colnames(res.tabA)=="date2bH"]  <- "date.2D.2.H"
colnames(res.tabA)[colnames(res.tabA)=="date2bL"] <- "date.2D.2.L"

res.tabA <- res.tabA[order(res.tabA$FinalResult,res.tabA$best.source1),]
write.csv(res.tabA,file=paste0(main_dir,"popgen/figures/AfricaGlobetrotterFinalResults.csv"),quote=F,row.names=F)

#############################################################
## NOW WORK ON GETTING MAKING NICE TEX TABLES
res.tabB <- res.tabA
res.tabB$Result <- res.tabB$FinalResult
rt1 <- res.tabB[!res.tabB$FinalResult%in%c("U","NA"),]
rt2 <- res.tabB[res.tabB$FinalResult%in%c("U","NA"),]
res1.ord <- rt1[,"best.source1"]
res1.ord[rt1$Result=="2D"] <- rt1$best.source1.date1[rt1$Result=="2D"]
res.ord.levs <- popplot
res1.ord <- factor(res1.ord,levels=res.ord.levs)
res2.ord <- rt1[,"best.source2"]
res2.ord[rt1$Result=="2D"] <- rt1$best.source2.date1[rt1$Result=="2D"]
res2.ord <- factor(res2.ord,levels=res.ord.levs)

rt1 <- rt1[order(res2.ord,res1.ord,as.numeric(rt1[,"date.1D"])),]
#res.ord <- c(res.ord,rep("XXX",nrow(rt2)))
res.tabB <- rbind(rt1,rt2)
## PULL OUT THE FINAL RESULTS
#rows2keep <- grep("main.null",res.tabB[,2],invert=F) ## OLD VERSION
res2keep <- res.tabB[,c("Cluster","FinalAnaly")]
res2keep <- res2keep[duplicated(res2keep),]
rows2keep <- c()
for(i in 1:nrow(res2keep))
{
    keeper <- (1:nrow(res.tabB))[res.tabB$Cluster==res2keep[i,1] & res.tabB$Analysis==res2keep[i,2]]
    rows2keep <- c(rows2keep,keeper)
}

res.tabB <- res.tabB[rows2keep,]

#############################################################
## THE FINAL MAIN RESULTS TABLE
write.csv(res.tabB,file=paste0(main_dir,"popgen/figures/AfricaGTBESTtable.csv"),quote=F)
#############################################################
res.tabcols <- c("Cluster","Analysis","pval","Result","FinalResult","max.R1.","FQ1","FQ2","maxScore2ev" ,
                 "date.1D","Date.CI", "alpha","best.source1","best.source2",
                 "alpha2","best.source1.ev2","best.source2.ev2",
                 "date.2D.1","Date2a.CI","alpha2.date1","best.source1.date1","best.source2.date1",
                 "date.2D.2","Date2b.CI","alpha2.date2","best.source1.date2","best.source2.date2")
colnames(res.tabB)[colnames(res.tabB)=="max(R1)"] <- "max.R1."
res.tabB <- res.tabB[,res.tabcols]

test <- res.tabB$Result!=res.tabB$FinalResult
#res.tabB$Result[test] <- paste(res.tabB$FinalResult[test],"(",res.tabB$Result[test],")",sep="")
res.tabB$Result[test] <- res.tabB$FinalResult[test]
res.tabB <- res.tabB[,colnames(res.tabB)!="FinalResult"]

for(i in c("date.1D","date.2D.1","date.2D.2"))
{
    res.tabB[,i] <- sapply(res.tabB[,i],function(x)makeDate(round(as.numeric(x)),add_BCE=F))
    res.tabB[,i][which(res.tabB[,i]<0)] <- paste0(-res.tabB[,i][which(res.tabB[,i]<0)],"B")
    res.tabB[,i] <- as.character(res.tabB[,i])
}

#############################################################
## SOME FURTHER COSMETIC CHANGES:: change p = 0 to p<0.01
res.tabB$pval[res.tabB$pval==0] <- "<0.01"
## only include best guess info
onedatecols <- c("date.1D","Date.CI","alpha","best.source1","best.source2" )
onemwcols <- c( "alpha2","best.source1.ev2","best.source2.ev2")
twodatecols <- c("date.2D.1","Date2a.CI","alpha2.date1","best.source1.date1","best.source2.date1",
                 "date.2D.2","Date2b.CI","alpha2.date2","best.source1.date2","best.source2.date2")

res.tabB[res.tabB$Result%in%c("1D","1D(2D)"),c(onemwcols,twodatecols)] <- NA
res.tabB[res.tabB$Result%in%c("1MW","1MW(2D)"),twodatecols] <- NA
res.tabB[res.tabB$Result=="2D",c(onedatecols,onemwcols)] <- NA

final.res2plot <- res.tabB
## DO I NEED TO OUTPUT THIS SOMEWHERE?
write.table(final.res2plot,file=paste0(main_dir,"popgen/data/MalariaGenGlobetrotter2plotFinal.txt"))
write.table(admixturesources2,paste0(main_dir,"popgen/data/MalariaGenGlobetrotterAdmixtureSources2.txt"))
write.table(admixturesources3,paste0(main_dir,"popgen/data/MalariaGenGlobetrotterAdmixtureSources3.txt"))
write.table(admixturesources4,paste0(main_dir,"popgen/data/MalariaGenGlobetrotterAdmixtureSources4.txt"))
write.table(dateboots,paste0(main_dir,"popgen/data/MalariaGenGlobetrotterOneDateBootstraps.txt"))
write.table(date2boots,paste0(main_dir,"popgen/data/MalariaGenGlobetrotterTwoDateBootstraps.txt"))

## combine the date CI colummns
res.tabB <- final.res2plot
res.tabB$date.1D <- gsub("NAnewlineNA",NA,paste(res.tabB$date.1D,res.tabB$Date.CI,sep="newline"))
res.tabB$date.2D.1 <- gsub("NAnewlineNA",NA,paste(res.tabB$date.2D.1,res.tabB$Date2a.CI,sep="newline"))
res.tabB$date.2D.2 <- gsub("NAnewlineNA",NA,paste(res.tabB$date.2D.2,res.tabB$Date2b.CI,sep="newline"))
res.tabB <- res.tabB[,!colnames(res.tabB)%in%c("Date.CI","Date2a.CI","Date2b.CI")]

## put in proper names
namecols <- c("Cluster",grep("source",colnames(res.tabB),value=T))
for(i in namecols)
{
    tmp <- tidyNames(res.tabB[,i], fula=T, khoesan=T)
    res.tabB[,i] <- tmp
}
#############################################################
res.tabB <- xtable(res.tabB,align="|r|r|cccccccccccccccccccccc|")
newlines <- c()
newlineord <- as.character(res.tabB$best.source1)
newlineord[is.na(newlineord)] <- res.tabB$best.source1.ev2[is.na(newlineord)]
newlineord[is.na(newlineord)] <- res.tabB$best.source1.date1[is.na(newlineord)]
for(i in (2:length(newlineord))) if(newlineord[i]!=newlineord[(i-1)]) newlines <- c(newlines,i)

addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(-1)
addtorow$pos[[2]] <- c(newlines-1)

addtorow$command[[1]]  <- "main results table"
addtorow$command[[2]] <- "\\hline \n"
print(res.tabB, floating=FALSE,
      tabular.environment="longtable", 
      comment=FALSE,
      include.colnames=FALSE,
      include.rownames=FALSE,
      caption.placement="top",file=paste0(main_dir,"popgen/figures/AfricaGTBESTtable.tex"),
      booktabs=TRUE,add.to.row=addtorow)#,

#system(paste0("sed -i 's|newline|\newline|g' ",main_dir,"popgen/figures/AfricaGTBESTtable.tex"))

#############################################################
## FINALLY RECORD ALL RESULTS
res.tabA <- res.tabA[order(res.tabA$Cluster,res.tabA$Analysis,res.tabA$best.source1),]
rescls <- res.tabA$Cluster
colnames(res.tabA)[colnames(res.tabA)=="max(R1)"] <- "max.R1."
res.tabA <- res.tabA[,res.tabcols]

test <- res.tabA$Result!=res.tabA$FinalResult
res.tabA$Result[test] <- paste(res.tabA$FinalResult[test],"(",res.tabA$Result[test],")",sep="")
res.tabA <- res.tabA[,colnames(res.tabA)!="FinalResult"]

for(i in c("date.1D","date.2D.1","date.2D.2"))
{
    res.tabA[,i] <- sapply(res.tabA[,i],function(x)makeDate(round(as.numeric(x)),add_BCE=F))
    res.tabA[,i][which(res.tabA[,i]<0)] <- paste0(-res.tabA[,i][which(res.tabA[,i]<0)],"B")
    res.tabA[,i] <- as.character(res.tabA[,i])
}
res.tabA$pval[res.tabA$pval==0] <- "<0.01"

tmp <- res.tabA
## combine the date CI colummns
res.tabA$date.1D <- gsub("NAnewlineNA",NA,paste(res.tabA$date.1D,res.tabA$Date.CI,sep="newline"))
res.tabA$date.2D.1 <- gsub("NAnewlineNA",NA,paste(res.tabA$date.2D.1,res.tabA$Date2a.CI,sep="newline"))
res.tabA$date.2D.2 <- gsub("NAnewlineNA",NA,paste(res.tabA$date.2D.2,res.tabA$Date2b.CI,sep="newline"))
res.tabA <- res.tabA[,!colnames(res.tabA)%in%c("Date.CI","Date2a.CI","Date2b.CI")]

## put in proper names
namecols <- c("Cluster",grep("source",colnames(res.tabA),value=T))
for(i in namecols)
{
    tmp <- tidyNames(res.tabA[,i], fula=T, khoesan=T)
    res.tabA[,i] <- tmp
}


res.tabA <- xtable(res.tabA,align="|r|r|cccccccccccccccccccccc|")

### print all events
newlines <- c()
for(i in (2:length(rescls))) if(rescls[i]!=rescls[(i-1)]) newlines <- c(newlines,i)
#newlines <- sort(c(newlines,newlines - 1))
#newlines <- newlines[newlines != 0]

addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(-1)
addtorow$pos[[2]] <- c(newlines-1)

addtorow$command[[1]]  <- "all results"
addtorow$command[[2]] <- "\\hline \n"
print(res.tabA, floating=FALSE,
      tabular.environment="longtable", 
      comment=FALSE,
      include.colnames=FALSE,
      include.rownames=FALSE,
      caption.placement="top",file=paste0(main_dir,"popgen/figures/AfricaGTALLtable.tex"),
      booktabs=TRUE,add.to.row=addtorow)#,




