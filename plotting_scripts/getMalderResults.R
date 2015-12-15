### SCRIPT TO GET MALDER RESULTS FOR PAPER -- CAN GET RESULTS FROM DIFFERENT GENETIC MAPS ###
### ONCE RUN WE CAN USE THE OUTPUT TO PLOT SEPARATELY IN 
### plotMalderResults.R
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
popplotorder <- popplotorder[c(1:16,18:49)]
poplist <- factor(popplotorder,levels=popplotorder)
poplist<- as.matrix(poplist)
n_pops <- nrow(poplist)
in_dir <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/alder/outputfinal/"

analys <- c("","CeuMap","AfrMap")

## GET MALDER FOR DIFFERENT MAP ANALYSES ##
malder_dates <- c()
malder_results <- c()
for(analy in analys)
{
    #data_suff <- paste0("_malder",analy,".data")
    data_suff <- paste0("_malderMinDis",analy,".data")
    ## GET RESULTS FOR ALL POPS IN POPLIST ##
    all_res2 <- c()
    for(i in 1:nrow(poplist))
    {
        pop1 <- as.character(poplist[i,1])
        file_name <- paste0(in_dir,pop1,data_suff)
        if(file.exists(file_name))
        {
            x <- getMalder(file_name,pop1)
            all_res2 <- rbind(all_res2,x)
        } else 
        {
            print(paste("no result for", pop1))
        }
    }
    all_res2 <- data.frame(all_res2)
    colnames(all_res2)[1] <- "pop"
    ## GET RESULTS FOR COMPARISON PLOTS ##
    malder_d <- matrix(0,nr=nrow(poplist),nc=4)
    rownames(malder_d) <- poplist
    colnames(malder_d) <- c("date1","date1.ci","date2","date2.ci")
    for(i in (1:nrow(poplist)))
    {
        pop1 <- as.character(poplist[i,1])
        if(pop1%in%as.character(unique(all_res2$pop)))
        {
            tt <- all_res2[all_res2[,"pop"]==pop1,]
            res_tmp <- unique(as.character(tt$result))
            if(res_tmp == "0")
            {
                d <- 0
                di <- 0
            }
            if(res_tmp == "1")
            {
                d <- unique(as.numeric(as.character(tt$time0)))[1]
                di <- unique(as.numeric(as.character(tt$time0.ci)))[1]
            }
            if(res_tmp == "2")
            {
                d <- unique(as.numeric(as.character(tt$time0)))[1]
                di <- unique(as.numeric(as.character(tt$time0.ci)))[1]
                d1 <- unique(as.numeric(as.character(tt$time1)))[1]
                di1 <- unique(as.numeric(as.character(tt$time1.ci)))[1]
                malder_d[i,c(3,4)] <- c(d1,di1)    
            }
            malder_d[i,c(1,2)] <- c(d,di)
            
        }
    }
    ### SAVE IF LOOKING AT NORMAL MAP RESULTS ###
    if(analy == "") analy <- "HAPMAP"
    malder_dates <- cbind(malder_dates,analy,malder_d)
    malder_results <- rbind(malder_results,cbind(all_res2[,1],analy,all_res2[,2:ncol(all_res2)]))
}

colnames(malder_results)[1] <- "pop"
## find the pop and their major region for each comparison
## this takes a bit of time, so try to only run once
if(is.null(malder_results$pop1[1]))
{
    malder_results$pop1 <- tidyNames(sapply(strsplit(as.character(malder_results$test.pops),
                                                     split="\\;"), "[[",1), fula=T)
    malder_results$tpop1 <- sapply(malder_results$pop1,function(x) getPopRegion(x,popkey))
    malder_results$pop2 <- tidyNames(sapply(strsplit(as.character(malder_results$test.pops),
                                                     split="\\;"), "[[",2), fula=T)
    malder_results$tpop2 <- sapply(malder_results$pop2,function(x) getPopRegion(x,popkey))
}


### COMPUTE Z-SCORES  ##
all_z <- c()
for(i in (1:nrow(poplist)))
{
    for(analy in c("HAPMAP",analys[2:3]))
    {
        pop1 <- as.character(poplist[i,1])
        print(paste("getting top comps for", pop1, "for the", analy, "analysis"))
        tt <- malder_results[malder_results$pop==pop1 & malder_results$analy==analy,]
        tt <- tt[order(tt$tpop1,tt$tpop2,as.numeric(as.character(tt$amp0)),decreasing=T),]
        if(nrow(tt)>0)
        {
            ## find largest amplitude for all pairs of major regions
            unique_comp <- unique(tt[,c("tpop1","tpop2")])
            for(j in 1:nrow(unique_comp))
            {
                test <- tt$tpop1==unique_comp[j,1] & tt$tpop2==unique_comp[j,2]
                ttt <- tt[test,]
                #ttt <- ttt[order(as.numeric(as.character(ttt$amp0)),decreasing=T),][1,]
                ttt <- cbind("ev1",analy,ttt[1,])
                all_z <- rbind(all_z,ttt)
                if(as.character(ttt$result)=="2")
                {
                    test <- tt$tpop1==unique_comp[j,1] & tt$tpop2==unique_comp[j,2]
                    ttt <- tt[test,]
                    ttt <- ttt[order(as.numeric(as.character(ttt$amp1)),decreasing=T),][1,]
                    ttt <- cbind("ev2",analy,ttt)
                    colnames(ttt) <- colnames(all_z)
                    all_z <- rbind(all_z,ttt)
                }
            }
        }
    }
}

colnames(all_z)[1] <- "event"
## now i have a table with the top amplitudes for each comparison for each population
## calculate z-scores

getZ <- function(tt,group="Eurasia",amp=0)
{
    ev_amp <- paste0("amp",amp)
    tt <- tt[order(as.numeric(as.character(tt[,ev_amp])),decreasing=T),]
    ## GET MAX CURVE OF CONTAINING AT LEAST ONE POP FROM GROUP OF INTEREST
    etest <- tt$tpop1%in%group | tt$tpop2%in%group
    cmax_self <- tt[etest,][1,]
    cmax_selfamp <- as.numeric(as.character(cmax_self$amp0))
    cmax_selfampci <- as.numeric(as.character(cmax_self$amp0.ci))

    ## GET MAX CURVE OF CONTAINING NO POPS FROM GROUP OF INTEREST
    ntest <- !tt$tpop1%in%group & !tt$tpop2%in%group
    cmax_notself <- tt[ntest,][1,]
    cmax_notselfamp <- as.numeric(as.character(cmax_notself$amp0))
    cmax_notselfampci <- as.numeric(as.character(cmax_notself$amp0.ci))
    
    ## COMPUTE Z-SCORE
    num <- cmax_selfamp - cmax_notselfamp
    den <- sqrt((cmax_selfampci^2) + (cmax_notselfampci^2))
    ez <- num/den 
    return(rbind(cbind(ez,cmax_self),cbind(ez,cmax_notself)))
}

compZ <- function(cmax1,cmax2,amp=0){
    if(amp == 0)
    {
        cmax1amp <- as.numeric(as.character(cmax1$amp0))
        cmax1ampci <- as.numeric(as.character(cmax1$amp0.ci))
        cmax2amp <- as.numeric(as.character(cmax2$amp0))
        cmax2ampci <- as.numeric(as.character(cmax2$amp0.ci))
    }
    
    if(amp == 1)
    {
        cmax1amp <- as.numeric(as.character(cmax1$amp1))
        cmax1ampci <- as.numeric(as.character(cmax1$amp1.ci))
        cmax2amp <- as.numeric(as.character(cmax2$amp1))
        cmax2ampci <- as.numeric(as.character(cmax2$amp1.ci))
    }
    num <- cmax1amp - cmax2amp
    den <- sqrt((cmax1ampci^2) + (cmax2ampci^2))
    z <- num/den 
    return(z)
}
    


## COMPARE THE TOP CURVE RESULT TO CURVES CONTAINING POPS FROM OTHER ANCESTRIES
#z_res <- matrix(nr=0,nc=4+length(ancreg_list))
#colnames(z_res) <- c("pop","analy","event","anc",ancreg_list)
z_res <- c()
for(i in (1:nrow(poplist)))
{
    for(analy in c("HAPMAP",analys[2:3]))
    {
        pop1 <- as.character(poplist[i,1])
        tres <- all_z[all_z$pop==pop1 & all_z$analy == analy,]
        n_evs <- unique(as.character(tres$result))
        for(k in 1:n_evs)
        {
            event <- paste0("ev",k)
            eventsortcol <- paste0("amp",k-1)
            tt <- tres[tres$event==event,]
            if(k == 1) best_anc_cols <- c(1:6,7:12,19:22)
            if(k == 2) best_anc_cols <- c(1:6,13:18,19:22)
            
            if(nrow(tt)>0)
            {
                tt <- tt[order(as.numeric(as.character(tt[,eventsortcol])),decreasing=T),]
                ## WE ARE INTERESTED IN COMPARING THE TOP DONORS
                best_anc <- tt[1,]
                ## THE TOP CURVE HAS TWO ANCESTRIES
                ## WE WANT TO KNOW IF THIS CURVE IS SIGNIFICANTLY DIFFERENT
                ## FROM OTHER CURVES INVOLVING THESE TWO ANCESTRIES, TO DO
                ## THIS WE TAKE EACH ANCESTRY IN TURN AND LOOK AT ALL OTHER CURVES
                ## THAT CONTAIN THIS ANCESTRY
                tmp_mat <- matrix(NA,nr=1,nc=length(ancreg_list)*2)
                colnames(tmp_mat) <- c(paste("tpop1",ancreg_list,sep="."),paste("tpop2",ancreg_list,sep="."))
                for(anc in 1:2)
                {
                    if(anc == 1 )
                    {
                        best_ancs <- best_anc$tpop1
                        anc_cols <- "tpop2"   
                    }
                    if(anc == 2 )
                    {
                        best_ancs <- best_anc$tpop2
                        anc_cols <- "tpop1"
                    }
                    newtt <- tt[tt$tpop1 == best_ancs | tt$tpop2 == best_ancs ,]
                    newtt <- newtt[newtt$tpop1 != newtt$tpop2,]
                    
                    for(j in 1:nrow(newtt))
                    {
                        z <- compZ(newtt[1,],newtt[j,])
                        test_anc <- as.character(newtt[j,c("tpop1","tpop2")])
                        test_anc <- test_anc[!test_anc%in%best_ancs]
                        tmp_mat[1,paste(anc_cols,test_anc,sep=".")] <- z
                    }
                }
                ## FOR EACH OF THE TWO MAIN ANCESTRIES, GET THE Z-SCORE FOR DIFFERENCE BETWEEN
                ## BEST MAIN ANC CURVE AND BEST CURVE NOT INCLUDING MAIN ANC
                tpop1Z <- getZ(tt,group=best_anc$tpop1,amp=(k-1))$ez[1]
                tpop2Z <- getZ(tt,group=best_anc$tpop2,amp=(k-1))$ez[1]
            }
            
            z_res <- rbind(z_res,c(t(apply(best_anc[1,best_anc_cols],1,as.character)),tmp_mat,tpop1Z,tpop2Z))
        }
    }
}

colnames(z_res) <- c(colnames(best_anc)[c(1:6,7:12,19:22)],colnames(tmp_mat),c("tpop1.Z","tpop2.Z"))
z_res <- data.frame(z_res)
## z_res DESCRIBES THE COMPARISONS OF THE TOP ANCESTRY CURVES WITH CURVES CONTAINING OTHER
## ANCESTRIES, IF THE VALUE IN A COLUMN IS 0 < VALUE < 2, THEN THE CURVES FROM THAT ANCESTRY
## ARE NOT DIFFERENT FROM THE TOP ANCESTRY
main_anc <- other_anc1 <- other_anc2 <- c()
for(i in 1:nrow(z_res))
{
    tmp <- as.numeric(apply(z_res[i,c("tpop1.Z","tpop2.Z")],2,as.character))
    tmp <- apply(z_res[i,c("tpop1","tpop2")],2,as.character)[tmp>2]
    if(length(tmp)==0) main_anc <- c(main_anc,"none")
    if(length(tmp)==1) main_anc <- c(main_anc,tmp)
    if(length(tmp)==2) main_anc <- c(main_anc,paste(tmp,collapse=";"))
    
    for(tpop in c("tpop1","tpop2"))
    {
        tmp <- z_res[i,grep(tpop,colnames(z_res))]
        tmp <- apply(tmp[2:(length(tmp)-1)],2,as.numeric)
        tmp <- tmp[!is.na(tmp)]
        tmp <- tmp[tmp<2]
        other_ances <- gsub("\\.","\\-",paste(gsub(paste0(tpop,"\\."),"",names(tmp)),collapse=";"))
        if(tpop == "tpop1") other_anc1 <- c(other_anc1,other_ances)
        if(tpop == "tpop2") other_anc2 <- c(other_anc2,other_ances)
    }
}

z_res <- cbind(z_res,main_anc,other_anc1,other_anc2)

colnames(z_res)[c(3,10,11,6,7,35,2,5)] <- c("EthnicGroup","Date.Gens","Date.Gens.CI","Test.Pops","amplitude",
                                            "Main.Anc","Rec.Map", "N.evs")
write.table(z_res,paste0(malder_out,".txt"),quote=F,col.names=T,row.names=F)


############################################################
## FINALLY, GENERATE A TABLE FOR THE SOM
print_cols <- colnames(t_res)[c(1,2,3,16,15,18,19,17,20)]
tmp.tab <- t_res[,print_cols]
tmp.tab$main_anc_p <- signif(tmp.tab$main_anc_p,5)
colnames(tmp.tab) <- c("EthnicGroup","Date.Gens","Date.Gens.CI","Test.Pops","amplitude",
                       "Main.Anc","Main.Anc.p","Rec.Map", "N.evs")

tmp.tab$Date.Gens <- as.character(round(as.numeric(as.character(tmp.tab$Date.Gens))))
tmp.tab$Date.Gens.CI <- as.character(round(as.numeric(as.character(tmp.tab$Date.Gens.CI))))
tmp.tab$amplitude <- as.character(signif(as.numeric(as.character(tmp.tab$amplitude)),digits=3))


tmp.tab$Main.Anc.p[tmp.tab$Main.Anc.p<10e-5] <- "<10^{-5}"

## NEED TO ZERO OUT THE RELEVANT LINES FOR THE O EVENTS

tmp.tab <- xtable(tmp.tab)
align(tmp.tab) <- "l|l|c|c|l|c|l|c|c|c|"

print(tmp.tab,include.rownames=F,floating=F,include.colnames=T,
      tabular.environment="longtable",
      file=paste0(malder_out,".tex"),
       hline.after=c(-1,0,nrow(tmp.tab)),size="tiny",
      sanitize.text.function = function(x){x}) 

system(paste0("sed -i 's#_#-#g' ",malder_out,".tex"))
