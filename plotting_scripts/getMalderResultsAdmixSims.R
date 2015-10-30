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
popkey_file <- "data/MalariaGenAdmixturePopulationOverview.txt"
malder_out <- paste0("~/R/Copy/Rprojects/AfricaPOPGEN/manuscript/f3tables/AllPopsMalderFinalAdmixSimsALL")
sim_pops <- scan("/mnt/kwiat/well/human/george/admix_sims/chromopainter/admix_sims.pops",what="char")

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
in_dir <- "/mnt/kwiat/well/human/george/admix_sims/malder/output/"
analys <- c("") #,"CeuMap","AfrMap")

## GET MALDER FOR DIFFERENT MAP ANALYSES ##
malder_dates <- c()
malder_results <- c()
for(analy in analys)
{
    data_suff <- paste0("_alderallrefsmalder",analy,".data")
    ## GET RESULTS FOR ALL POPS IN POPLIST ##
    all_res2 <- c()
    for(i in sim_pops)
    {
        pop1 <- i
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
    malder_d <- matrix(0,nr=length(sim_pops),nc=4)
    rownames(malder_d) <- sim_pops
    colnames(malder_d) <- c("date1","date1.ci","date2","date2.ci")
    for(i in sim_pops)
    {
        pop1 <- i
        if(pop1%in%as.character(unique(all_res2$pop1)))
        {
            tt <- all_res2[all_res2[,"pop"]==pop1,]
            res_tmp <- unique(as.character(tt$result))
            d <- unique(as.numeric(as.character(tt$time0)))[1]
            di <- unique(as.numeric(as.character(tt$time0.ci)))[1]
            malder_d[i,c(1,2)] <- c(d,di)
            if(res_tmp == "multi")
            {
                d <- unique(as.numeric(as.character(tt$time1)))[1]
                di <- unique(as.numeric(as.character(tt$time1.ci)))[1]
                malder_d[i,c(3,4)] <- c(d,di)
            }
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
for(i in sim_pops)
{
    for(analy in c("HAPMAP")) #,analys[2:3]))
    {
        pop1 <- i
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
                if(as.character(ttt$result)=="multi")
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


## now i have a table with the top amplitudes for each comparison for each population
## calculate z-scores

getZ <- function(tt,group=eurasia,ev_amp="amp0")
{
    tt <- tt[order(as.numeric(as.character(tt[,ev_amp])),decreasing=T),]
    etest <- tt$tpop1%in%group | tt$tpop2%in%group
    e <- tt[etest,][1,]
    ntest <- !tt$tpop1%in%group & !tt$tpop2%in%group
    n <- tt[ntest,][1,]
    ez <- as.numeric(as.character(n[,ev_amp])) - as.numeric(as.character(e[,ev_amp])) / sqrt((as.numeric(as.character(n[,paste0(ev_amp,".ci")])))^2 + (as.numeric(as.character(e[,paste0(ev_amp,".ci")])))^2 )
    return(rbind(cbind(ez,e),cbind(ez,n)))
}



t_res <- c()
for(i in sim_pops)
{
    for(analy in c("HAPMAP"))#,analys[2:3]))
    {
        pop1 <- i
        tt <- all_z[all_z$pop==pop1 & all_z$analy == analy,]
        if(nrow(tt)>0)
        {
            n_evs <- unique(as.character(tt$result))
            ee <- getZ(tt,group="Eurasia")
            eez <- as.numeric(as.character(ee$ez[1]))
            if(is.na(eez)) eez <- 0
            nc <- getZ(tt,group="Western_Africa_Niger-Congo")
            ncz <- as.numeric(as.character(nc$ez[1]))
            nc1 <- getZ(tt,group="Central_West_Africa_Niger-Congo")
            nc1z <- as.numeric(as.character(nc1$ez[1]))
            nc2 <- getZ(tt,group="South_Africa_Niger-Congo")
            nc2z <- as.numeric(as.character(nc2$ez[1]))
            nc3 <- getZ(tt,group="East_Africa_Niger-Congo")
            nc3z <- as.numeric(as.character(nc2$ez[1]))
            nc4 <- getZ(tt,group="East_Africa_Nilo-Saharan")
            nc4z <- as.numeric(as.character(nc3$ez[1]))
            nc5 <- getZ(tt,group="South_Africa_KhoeSan")
            nc5z <- as.numeric(as.character(nc4$ez[1]))
            best_anc <- tt[order(as.numeric(as.character(tt$amp0)),decreasing=T),][1,]
            
            f_res <- c(pop1,
                       as.character(nc$time0[1]),
                       as.character(nc$time0.ci[1]),
                       eez,ncz,nc1z,nc2z,nc3z,nc4z,nc5z,
                       as.character(best_anc$pop1),
                       as.character(best_anc$tpop1),
                       as.character(best_anc$pop2),
                       as.character(best_anc$tpop2),
                       as.character(best_anc$amp0),
                       as.character(best_anc$test.pops),
                       analy)
            t_res <- rbind(t_res,f_res)
            if(n_evs == "multi")
            {
                tt <- all_z[all_z$pop==pop1 & all_z$analy == analy,]
                ee <- getZ(tt,group="Eurasia",ev_amp="amp1")
                eez <- as.numeric(as.character(ee$ez[1]))
                nc <- getZ(tt,group="Western_Africa_Niger-Congo",ev_amp="amp1")
                ncz <- as.numeric(as.character(nc$ez[1]))
                nc1 <- getZ(tt,group="Central_West_Africa_Niger-Congo",ev_amp="amp1")
                nc1z <- as.numeric(as.character(nc1$ez[1]))
                nc2 <- getZ(tt,group="East_Africa_Niger-Congo",ev_amp="amp1")
                nc2z <- as.numeric(as.character(nc2$ez[1]))
                nc3 <- getZ(tt,group="South_Africa_Niger-Congo",ev_amp="amp1")
                nc3z <- as.numeric(as.character(nc2$ez[1]))
                nc4 <- getZ(tt,group="East_Africa_Nilo-Saharan",ev_amp="amp1")
                nc4z <- as.numeric(as.character(nc3$ez[1]))
                nc5 <- getZ(tt,group="South_Africa_KhoeSan",ev_amp="amp1")
                nc5z <- as.numeric(as.character(nc4$ez[1]))
                best_anc <- tt[order(as.numeric(as.character(tt$amp1)),decreasing=T),][1,]
                f_res <- c(pop1,
                           as.character(ee$time1[1]),
                           as.character(ee$time1.ci[1]),
                           eez,ncz,nc1z,nc2z,nc3z,nc4z,nc5z,
                           as.character(best_anc$pop1),
                           as.character(best_anc$tpop1),
                           as.character(best_anc$pop2),
                           as.character(best_anc$tpop2),
                           as.character(best_anc$amp1),
                           as.character(best_anc$test.pops),
                           analy)
                t_res <- rbind(t_res,f_res)
            }
        }
    }
}

t_res <- data.frame(t_res)
colnames(t_res) <- c("pop1","date","dateCI",
                     "eurasia.Z","wanc.Z","cwanc.Z","sanc.Z","eanc.Z",
                     "eans.Z","saks.Z","donpop1","ancestry1",
                     "donpop2","ancestry2","amp","test.pops", "analy")
####
## find col with most neg z-score for the two inferred ancestries
z_cols <- c("eurasia.Z","wanc.Z","cwanc.Z","sanc.Z","eanc.Z","eans.Z","saks.Z")
anc_table <- data.frame(z_cols,
                        c("Eurasia","Western_Africa_Niger-Congo","Central_West_Africa_Niger-Congo",
                          "South_Africa_Niger-Congo","East_Africa_Niger-Congo","East_Africa_Nilo-Saharan",
                          "South_Africa_KhoeSan"))
colnames(anc_table) <- c("Z","reg")
main_anc <- c()
main_anc_p <- c()
for(i in 1:nrow(t_res))
{
    ii <- t_res[i,z_cols]
    ianc <- c(as.character(t_res$ancestry1[i]),as.character(t_res$ancestry2[i]))
    ianc2 <- as.character(sapply(ianc,function(x){anc_table$Z[x==as.character(anc_table$reg)]}))
    ii <- names(sort(ii[ianc2],decreasing=T)[1])
    ip <- z2p(as.numeric(as.character(t_res[i,ii])))
    ii <- as.character(anc_table$reg[ii==as.character(anc_table$Z)])
    main_anc <- c(main_anc,ii)
    main_anc_p <- c(main_anc_p,ip)
}

t_res <- cbind(t_res,main_anc,main_anc_p)

## FINALLY, GENERATE A TABLE FOR THE SOM
print_cols <- colnames(t_res)[c(1,2,3,16,15,18,19,17)]
tmp.tab <- t_res[,print_cols]
tmp.tab$main_anc_p <- signif(tmp.tab$main_anc_p,5)
colnames(tmp.tab) <- c("EthnicGroup","Date.Gens","Date.Gens.CI","Test.Pops","amplitude",
                       "Main.Anc","Main.Anc.p","Rec.Map")

tmp.tab$Date.Gens <- as.character(round(as.numeric(as.character(tmp.tab$Date.Gens))))
tmp.tab$Date.Gens.CI <- as.character(round(as.numeric(as.character(tmp.tab$Date.Gens.CI))))
tmp.tab$amplitude <- as.character(signif(as.numeric(as.character(tmp.tab$amplitude)),digits=3))
write.table(tmp.tab,paste0(malder_out,".txt"),quote=F,col.names=T,row.names=F)


tmp.tab$Main.Anc.p[tmp.tab$Main.Anc.p<10e-5] <- "<10^{-5}"

tmp.tab <- xtable(tmp.tab)
align(tmp.tab) <- "l|l|c|c|l|c|l|c|c|"

print(tmp.tab,include.rownames=F,floating=F,include.colnames=T,
      tabular.environment="longtable",
      file=paste0(malder_out,".tex"),
      hline.after=c(-1,0,nrow(tmp.tab)),size="tiny",
      sanitize.text.function = function(x){x}) 

system(paste0("sed -i 's#_#-#g' ",malder_out,".tex"))
