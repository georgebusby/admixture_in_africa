## READS IN REAL DATA propsanddates/ RESULTS AND RANKS BASED ON EVIDENCE FOR 0,1, OR 2 EVENTS
## FOR GEORGE POPS!

## AvRelMSEIter4/AvRelMSEIter5:  Yoruba (ONLY pop with NaN in "Scores.txt")
## AvRelSMEIter5DateCIii: SanNamibia/Yoruba failed (replaced SanN with SanK; replaced Yoruba with AvRelMSEIter5DateCI Yoruba)

#########################################
## INPUT:

out.file="/ugi/home/ucbtgrh/WTCHG/George/gavinmaps/GeorgePopsFinal1EventsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIviSummaryWebInfoAugust2012SPLITPvalueFinal.txt"

                      ## for selecting analysis (and getting everything but propoportion, dates, intercepts, avrelMSE, and side-of-event-info):
admixture.conclusion.file="/ugi/home/ucbtgrh/WTCHG/George/summaryresults/AdmixtureHistoryPaper1EventFinalResultsSummaryPCNoShrinkVaryGridsAugust2012SPLITPvalue.txt"
ignore.nonsplit=1	## preferentially take split pops?
                 ## print curves from donor pops involved in second event, even if they do not contribute to PC1 assuming 1 event?
print.curves.secondevent='yes'
#curves.cutoff.secondevent=0.01  ## minimum weight put on second event pops for consideration -- any second event pops w/ mixing-prop > side.of.event.cutoff/(curves.cutoff.secondevent*admixprop.prop) is included (currently set so that any second event pop with admixprop*mixing-prop > 10% is included)
max.pops.multievent=10  ## no more than this many curves per side of event, in 2-date and complex cases

#################################
## FIXED INPUT:

                ## 1-EVENT:
fst.filesPOST="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txtDatesIntercepts"
fst.filesPOSTii="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txt"
fst.filesPOSTiii="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txtScores"
fst.filesPOSTiv="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txtResults"

                ## MULTIPLE-EVENTS, SINGLE-DATE:
fst.filesPOST.multievents.1="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txt"
fst.filesPOST.multievents.2="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txt1DatePC2"

                ## 2-EVENTS:
fst.filesPOST.2dates.1="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txt2DatesPC1"
fst.filesPOST.2dates.2="1EventfstplotGeorgePopsChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txt2DatesPC2"
#fst.filesPOST.1st.init="2x1stEventfstplotChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemove.txt"
#fst.filesPOST.2nd.init="2x2ndEventfstplotChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemove.txt"

                ## RAW PROPS: (NOW NEEDED FOR MULTIPLE OR 2 EVENTS)
rawprops.pathwayHGDP="/ugi/home/ucbtgrh/WTCHG/George/ChromoPainterOUTPUTAugust2012SPLIT/LeaveOneOutProps/"
rawprops.suffixHGDP="GeorgePopsFinalAllChromLeaveOneOutPropsFixedNeMut.chunklengths.out"
rawprops.pathwayHGDP.sep="/ugi/home/ucbtgrh/WTCHG/George/ChromoPainterOUTPUTAugust2012SPLITSeparateAnalysis/LeaveOneOutProps/"
rawprops.pathway.sep="/ugi/home/ucbtgrh/WTCHG/George/ChromoPainterOUTPUTAugust2012SPLITSeparateAnalysis/PaintingSamples/"
rawprops.suffixHGDP.sep="PopsAllChromLeaveOneOutPropsFixedNeMut.chunklengths.out"
rawprops.suffix.sep="PopsAllChromPaintingSamples.chunklengths.out"

#################################
##################################
## FIXED INPUT:

date.summary.files.POST.vec="AdmixCurves1EventOneDateBootstrapChromoChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveBetterMixingSummary.txt"
date.summary.files.POST.vec.auxfiles="AdmixCurves1EventOneDateBootstrapChromoChunkLengthsNoSelfCopyWEIGHTINGGenDistCutoffNoChunkRemoveBetterMixingSummary.txt"
date.summary.files.POST.vec.2event="AdmixCurves2EventOneDateBootstrapChromoChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveBetterMixingSummary.txt"
date.summary.files.POST.vec.2event.twodategrid="AdmixCurves2EventAncientRecentAdmixtureOneDateBootstrapChromoChunkLengthsWEIGHTINGGenDistCutoffNoChunkRemoveBetterMixingSummary.txt"
date.summary.files.POST.vec.auxfiles.2event="AdmixCurves2EventOneDateBootstrapChromoChunkLengthsNoSelfCopyWEIGHTINGGenDistCutoffNoChunkRemoveBetterMixingSummary.txt"
date.summary.files.POST.vec.auxfiles.2event.twodategrid="AdmixCurves2EventAncientRecentAdmixtureOneDateBootstrapChromoChunkLengthsNoSelfCopyWEIGHTINGGenDistCutoffNoChunkRemoveBetterMixingSummary.txt"
fst.filesPRE.full="/ugi/scratch/ucbtgrh/WTCHG/George/ChromoPainterOUTPUTAugust2012SPLIT/plotINFO/"
fst.filesPRE.sep="/ugi/scratch/ucbtgrh/WTCHG/George/ChromoPainterOUTPUTAugust2012SPLITSeparateAnalysis/plotINFO/"
pathway.full="/ugi/scratch/ucbtgrh/WTCHG/George/ChromoPainterOUTPUTAugust2012SPLIT/propsanddates/"
pathway.sep="/ugi/scratch/ucbtgrh/WTCHG/George/ChromoPainterOUTPUTAugust2012SPLITSeparateAnalysis/propsanddates/"
pop.infile="/ugi/home/ucbtgrh/WTCHG/George/inputfiles/AllGeorgePopsFinalAugust2012SPLIT.pops.list"
geo.info.file="/ugi/home/ucbtgrh/WTCHG/George/datainfo/MedPeopleAdmixtureLatLongGAVINGeorge.txt"

                     ## sim info (for getting "p-values" for 2-date and multiway admixture):
pvalue.infileSIMS="/ugi/home/ucbtgrh/WTCHG/George/summaryresults/SIMSMarch2013.txt"
excluded.pops=c("Yoruba7GenProp20BrahuiHan30GenProp20FakeSwitchIndNOOVERLAP2Event","Yoruba7GenProp20BrahuiHan30GenProp50FakeSwitchIndNOOVERLAP2Event","Yoruba7GenProp20BrahuiHan150GenProp20FakeSwitchIndNOOVERLAP2Event","Yoruba7GenProp20BrahuiHan150GenProp50FakeSwitchIndNOOVERLAP2Event","BasqueFrenchFakeSwitchIndNOOVERLAP1000GenProp50","FrenchFrenchFakeSwitchIndNOOVERLAP1000GenProp50","BrahuiBrahuiFakeSwitchIndNOOVERLAP1000GenProp50","HanJapaneseFakeSwitchIndNOOVERLAP1000GenProp50","SindhiPathanFakeSwitchIndNOOVERLAP1000GenProp50")   ## don't include these pops when calculating thresholds for finding 2-events or multi-way admixture
dates.filesPRESIMS="/ugi/scratch/ucbtgrh/WTCHG/George/sims/ChromoPainterOUTPUTAugust2012SPLIT/propsanddates/"
dates.filesPOSTSIMS="GeorgePopsPHASEAdmixCurves1EventOneDateBootstrapChromoChunkLengthsNoSelfCopyWEIGHTINGGenDistCutoffNoChunkRemoveBetterMixingSummary.txt"
fit.qualityPRESIMS="/ugi/scratch/ucbtgrh/WTCHG/George/sims/ChromoPainterOUTPUTAugust2012SPLIT/plotINFO/"
fit.qualityPOSTSIMS="1EventfstplotGeorgePopsPHASEChunkLengthsNoSelfCopyWEIGHTINGGenDistCutoffNoChunkRemoveAvRelMSEIter5DateCIvi.txtResults"
bootstrap.filesPRESIMS="/ugi/scratch/ucbtgrh/WTCHG/George/sims/ChromoPainterOUTPUTMarch2013/bootstraps/"
bootstrap.filesPOSTSIMS="EverythingBootstrapSummary.txt"
fit.qualityBOTH.thres=0.985
pvalue.thres=0
r2.1event.thres.2eventtest=0.975       ## only consider curves whose R^2 for one-event is less than this when calculating "two-date" score (should be 0.975, but I blew it for this submission)
r2.cutoff=0.0		## only consider curves whose R^2 for two-events is at least this when calculating "two-date" score
two.date.null.real.r2.fac=1/3		 ## null two-event-score divided by "real" two-event-score must be >= than this to declare two-dates

gen.time=28
prop.cutoff=0.001
side.of.event.cutoff=0    # 0.001
max.donors=3         ## max number of donors listed per side of event
round.r2=3
sd.fac=1/20         ## for proportion CI
round.date=0
round.val.int=10
round.fit=2
round.prop=5
round.coeff=3
round.avrel=5
round.pval=3

################################
## PREVIOUS JUNK:

# #sim.files=c("/home.oak/hellenth/George/inputfiles/AllGeorgePopsFinal.pops.listII","/home.oak/hellenth/George/inputfiles/SplitRedo.pops.list")
# #aux.pop.files=c("/home.oak/hellenth/George/inputfiles/AllGeorgePopsFinal.pops.listII","/home.oak/hellenth/George/inputfiles/SplitRedo.pops.list","/home.oak/hellenth/George/inputfiles/ArabSeparateAnalysisPops.txt","/home.oak/hellenth/George/inputfiles/SplitArabSeparateAnalysis.pops.list","/home.oak/hellenth/George/inputfiles/EthiopianSeparateAnalysisPops.txt","/home.oak/hellenth/George/inputfiles/MediterraneanSeparateAnalysisPops.txt","/home.oak/hellenth/George/inputfiles/PakistanSeparateAnalysisPops.txt","/home.oak/hellenth/George/inputfiles/SplitPakistanSeparateAnalysis.pops.list","/home.oak/hellenth/George/inputfiles/SanSeparateAnalysisPops.txt","/home.oak/hellenth/George/inputfiles/SlavicSeparateAnalysisPopsII.txt","/home.oak/hellenth/George/inputfiles/SlavicSeparateAnalysisPopsIIINoLith.txt","/home.oak/hellenth/George/inputfiles/SplitSlavicSeparateAnalysis.pops.list","/home.oak/hellenth/George/inputfiles/SplitSlavicSeparateAnalysis.pops.list","/home.oak/hellenth/George/inputfiles/
# AllGeorgePopsFinal2Dates.pops.list","/home.oak/hellenth/George/inputfiles/SplitRedo2Dates.pops.list","/home.oak/hellenth/George/inputfiles/ArabSeparateAnalysisPops2Dates.txt","/home.oak/hellenth/George/inputfiles/SplitArabSeparateAnalysis2Dates.pops.list","/home.oak/hellenth/George/inputfiles/MediterraneanSeparateAnalysisPops2Dates.txt","/home.oak/hellenth/George/inputfiles/PakistanSeparateAnalysisPops2Dates.txt","/home.oak/hellenth/George/inputfiles/SplitPakistanSeparateAnalysis.pops.list","/home.oak/hellenth/George/inputfiles/SanSeparateAnalysisPops2Dates.txt","/home.oak/hellenth/George/inputfiles/AllGeorgePopsFinalMultiOnly.pops.list","/home.oak/hellenth/George/inputfiles/SplitRedoMultiOnly.pops.list","/home.oak/hellenth/George/inputfiles/ArabSeparateAnalysisPopsMultiOnly.txt","/home.oak/hellenth/George/inputfiles/SplitArabSeparateAnalysis.pops.list","/home.oak/hellenth/George/inputfiles/MediterraneanSeparateAnalysisPopsMultiOnly.txt","/home.oak/hellenth/George/inputfiles/SanSeparateAnalysisPopsMultiOnly.
# txt","/home.oak/hellenth/George/inputfiles/SlavicIVSeparateAnalysisPopsMultiOnly.txt","/home.oak/hellenth/George/inputfiles/SlavicVSeparateAnalysisPopsMultiOnly.txt","/home.oak/hellenth/George/inputfiles/SplitSlavicIVSeparateAnalysisMultiOnly.pops.list","/home.oak/hellenth/George/inputfiles/SplitSlavicVSeparateAnalysisMultiOnly.pops.list")
# #aux.pop.files.addon=c("FullAnalysis","FullAnalysis","ArabSeparateAnalysis","ArabSeparateAnalysis","EthiopianSeparateAnalysis","MediterraneanSeparateAnalysis","PakistanSeparateAnalysisII","PakistanSeparateAnalysisII","SanSeparateAnalysis","SlavicSeparateAnalysisIV","SlavicSeparateAnalysisV","SlavicSeparateAnalysisIV","SlavicSeparateAnalysisV","FullAnalysis","FullAnalysis","ArabSeparateAnalysis","ArabSeparateAnalysis","MediterraneanSeparateAnalysis","PakistanSeparateAnalysisII","PakistanSeparateAnalysisII","SanSeparateAnalysis","FullAnalysis","FullAnalysis","ArabSeparateAnalysis","ArabSeparateAnalysis","MediterraneanSeparateAnalysis","SanSeparateAnalysis","SlavicSeparateAnalysisIV","SlavicSeparateAnalysisV","SlavicSeparateAnalysisIV","SlavicSeparateAnalysisV")
# #num.events=c(1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3)     ## "3"=multi

###############################
## PROGRAM:

options(scipen=999)

full.pop.list=as.character(read.table(pop.infile,header=TRUE)$pop)
full.pop.nhaps=as.integer(read.table(pop.infile,header=TRUE)$nhaps)

             ## get sim names/haps:
admixture.conclusion.table=read.csv(admixture.conclusion.file,sep=",",header=TRUE)
pop.vec=as.character(admixture.conclusion.table$pop)
pop.haps=as.integer(admixture.conclusion.table$nhaps)
pop.conclusion=as.character(admixture.conclusion.table$admixture)
pvalue.vec=as.double(admixture.conclusion.table$pvalue)
rsquared.1event.max=as.double(admixture.conclusion.table$rsquared.1event.max)
rsquared.2event.max=as.double(admixture.conclusion.table$rsquared.2event.max)
r2.2event.minus.1event.rescaled.max=as.double(admixture.conclusion.table$r2.2event.minus.1event.rescaled.max)
r2.2event.minus.1event.rescaled.maxNULL=as.double(admixture.conclusion.table$r2.2event.minus.1event.rescaled.maxNULL)
fit.quality1=as.double(admixture.conclusion.table$fit.quality1)
fit.quality2.vec=as.double(admixture.conclusion.table$fit.quality2)
fit.quality2.minus1.rescaled=as.double(admixture.conclusion.table$fit.quality2.minus1.rescaled)
fit.qualityBOTH.vec=as.double(admixture.conclusion.table$fit.qualityBOTH)
fit.qualityBOTH.vecNULL=as.double(admixture.conclusion.table$fit.qualityBOTHNULL)
grid.type=as.character(admixture.conclusion.table$grid.type)
first.source.best.copyvec=as.character(admixture.conclusion.table$first.source.best.copyvec)
second.source.best.copyvec=as.character(admixture.conclusion.table$second.source.best.copyvec)
first.source.best.copyvec.PC1=as.character(admixture.conclusion.table$first.source.best.copyvec.PC1)
second.source.best.copyvec.PC1=as.character(admixture.conclusion.table$second.source.best.copyvec.PC1)
first.source.best.copyvec.PC2=as.character(admixture.conclusion.table$first.source.best.copyvec.PC2)
second.source.best.copyvec.PC2=as.character(admixture.conclusion.table$second.source.best.copyvec.PC2)
pop.analysis=as.character(admixture.conclusion.table$analysis)
pop.analysis[pop.analysis=="Final"]="FullAnalysis"
pop.analysis[pop.analysis=="Split"]="FullAnalysisSplit"

            ## sort:
pop.haps=pop.haps[sort(pop.vec,index.return=TRUE)$ix]
pop.conclusion=pop.conclusion[sort(pop.vec,index.return=TRUE)$ix]
pvalue.vec=pvalue.vec[sort(pop.vec,index.return=TRUE)$ix]
rsquared.1event.max=rsquared.1event.max[sort(pop.vec,index.return=TRUE)$ix]
rsquared.2event.max=rsquared.2event.max[sort(pop.vec,index.return=TRUE)$ix]
r2.2event.minus.1event.rescaled.max=r2.2event.minus.1event.rescaled.max[sort(pop.vec,index.return=TRUE)$ix]
r2.2event.minus.1event.rescaled.maxNULL=r2.2event.minus.1event.rescaled.maxNULL[sort(pop.vec,index.return=TRUE)$ix]
fit.quality1=fit.quality1[sort(pop.vec,index.return=TRUE)$ix]
fit.quality2.vec=fit.quality2.vec[sort(pop.vec,index.return=TRUE)$ix]
fit.quality2.minus1.rescaled=fit.quality2.minus1.rescaled[sort(pop.vec,index.return=TRUE)$ix]
fit.qualityBOTH.vec=fit.qualityBOTH.vec[sort(pop.vec,index.return=TRUE)$ix]
fit.qualityBOTH.vecNULL=fit.qualityBOTH.vecNULL[sort(pop.vec,index.return=TRUE)$ix]
grid.type=grid.type[sort(pop.vec,index.return=TRUE)$ix]
first.source.best.copyvec=first.source.best.copyvec[sort(pop.vec,index.return=TRUE)$ix]
second.source.best.copyvec=second.source.best.copyvec[sort(pop.vec,index.return=TRUE)$ix]
first.source.best.copyvec.PC1=first.source.best.copyvec.PC1[sort(pop.vec,index.return=TRUE)$ix]
second.source.best.copyvec.PC1=second.source.best.copyvec.PC1[sort(pop.vec,index.return=TRUE)$ix]
first.source.best.copyvec.PC2=first.source.best.copyvec.PC2[sort(pop.vec,index.return=TRUE)$ix]
second.source.best.copyvec.PC2=second.source.best.copyvec.PC2[sort(pop.vec,index.return=TRUE)$ix]
pop.analysis=pop.analysis[sort(pop.vec,index.return=TRUE)$ix]
pop.vec=sort(pop.vec)
#pop.vec=c("Finnish","Russian","Kalash","Adygei")
#pop.haps=(1:length(pop.vec))*2

         ## get geo.info:
geo.info=read.table(geo.info.file,header=TRUE)

             ## GET PVALUES FOR 2-date, multiway:
                         ## get sim information:
pvalue.readSIMS=read.csv(pvalue.infileSIMS,sep=',',header=TRUE)
pop.vecSIMS=unique(as.character(pvalue.readSIMS$sim))
fit.quality.1event.allSIMS=rsquared.2event.rescaled.allSIMS=NULL
for (i in 1:length(pop.vecSIMS))
{
	if (as.double(pvalue.readSIMS$p.value)[as.character(pvalue.readSIMS$sim)==pop.vecSIMS[i] & as.character(pvalue.readSIMS$analysis.type)=="one-date-recent-grid"]<=pvalue.thres && sum(is.element(pop.vecSIMS[i],excluded.pops))==0)
	{
		boot.resultsSIMS=read.table(paste(bootstrap.filesPRESIMS,pop.vecSIMS[i],bootstrap.filesPOSTSIMS,sep=''),header=TRUE)
		load(paste(fit.qualityPRESIMS,pop.vecSIMS[i],fit.qualityPOSTSIMS,sep=''))
		date.readSIMS=read.table(paste(dates.filesPRESIMS,pop.vecSIMS[i],dates.filesPOSTSIMS,sep=''),skip=1)
		rsquared.2event.rescaledSIMS=max(((as.double(date.readSIMS[,4])-as.double(date.readSIMS[,3]))/(1.0-as.double(date.readSIMS[,3])))[as.double(date.readSIMS[,4])>=r2.cutoff & as.double(date.readSIMS[,3])<r2.1event.thres.2eventtest])
		fit.quality.both.vec=c(fit.quality+fit.quality2,boot.resultsSIMS$fit.quality.1+boot.resultsSIMS$fit.quality.2)
		fit.quality.vec=c(fit.quality,boot.resultsSIMS$fit.quality.1)[fit.quality.both.vec>fit.qualityBOTH.thres]
		rsquared.2event.rescaled.vec=c(rsquared.2event.rescaledSIMS,boot.resultsSIMS$max.rsquared.2event.rescaled)[fit.quality.both.vec>fit.qualityBOTH.thres]
		fit.quality.1event.allSIMS=c(fit.quality.1event.allSIMS,fit.quality.vec)
		rsquared.2event.rescaled.allSIMS=c(rsquared.2event.rescaled.allSIMS,rsquared.2event.rescaled.vec)
	}
}
                         ## calculate 2-date and multi-way pvalues:
min.fit.qualityBOTH.vec=twodate.pvalue.vec=multiway.pvalue.vec=rep(NA,length(pop.vec))
for (i in 1:length(pop.vec))
{
				## check if any other pops are "uncertain" (based on NULL sim fit-quality):
	min.fit.qualityBOTH.vec[i]=min(c(fit.qualityBOTH.vecNULL[i],fit.qualityBOTH.vec[i]))
                                ## get two-date pvalue:
	twodate.pvalue.vec[i]=length(rsquared.2event.rescaled.allSIMS[rsquared.2event.rescaled.allSIMS>=r2.2event.minus.1event.rescaled.max[i]])/length(rsquared.2event.rescaled.allSIMS)
	if (twodate.pvalue.vec[i]<0.05 && (r2.2event.minus.1event.rescaled.maxNULL[i]/r2.2event.minus.1event.rescaled.max[i])<two.date.null.real.r2.fac) twodate.pvalue.vec[i]=NA
                                ## get multi-way pvalue:
	multiway.pvalue.vec[i]=length(fit.quality.1event.allSIMS[fit.quality.1event.allSIMS<=fit.quality1[i]])/length(fit.quality.1event.allSIMS)
}

         ## read in data and extract 1.event-0.event "likelihood diff"; 1-2 event "R-squared"; and 1-date/2-date intercepts:
full.info.mat=num.events.vec=avrelmse=proportion.est.mat=fit.quality.mat=fit.quality.rescaled=rsquared.1event.vec=rsquared.2event.vec=rsquared.1event.vec.max=rsquared.2event.vec.max=rsquared.2event.vec.max.rescaled=date.est.mat=date.est.ci.mat=date.est.ci.matII=proportion.est.mat.ci=grid.type.vec=best.copyvec.matrix=pvalue.vec.final=twodate.pvalue.vec.final=multiway.pvalue.vec.final=min.fit.quality.vec.final=NULL
count=1
for (i in 1:length(pop.vec))
  {
    split.pop="no"
    if (length(strsplit(pop.vec[i],split="Split")[[1]])>1) split.pop="yes"

    print(noquote(c(i,pop.vec[i],length(pop.vec))))

                     ## GET APPROPRIATE FILE NAMES:
    file.addonFST=""
    fst.filesPRE=fst.filesPRE.full
    rawprops.donorsPRE=rawprops.pathwayHGDP
    rawprops.donorsPOST=rawprops.suffixHGDP
    rawprops.recipientPRE=rawprops.pathwayHGDP
    rawprops.recipientPOST=rawprops.suffixHGDP
    if (pop.analysis[i]!="FullAnalysis" && pop.analysis[i]!="FullAnalysisSplit")
    {
 	file.addonFST=strsplit(pop.analysis[i],split="Split")[[1]][1]
	fst.filesPRE=fst.filesPRE.sep
	rawprops.donorsPRE=rawprops.pathwayHGDP.sep
    	rawprops.donorsPOST=paste("George",file.addonFST,rawprops.suffixHGDP.sep,sep='')
    	rawprops.recipientPRE=rawprops.pathway.sep
    	rawprops.recipientPOST=paste("George",file.addonFST,rawprops.suffix.sep,sep='')
    }

                        ## GET DATE, PROPORTION, AND SIDE-OF-EVENT-INFO:
    if (pop.conclusion[i]=="one-event" || pop.conclusion[i]=="no-admixture" || pop.conclusion[i]=="uncertain")
      {
        if (pop.analysis[i]=="FullAnalysis" || pop.analysis[i]=="FullAnalysisSplit")
          {
            file.addonFST=""
            file.addon=paste("PopsFinal",date.summary.files.POST.vec,sep='')
            ci.filein=paste(pathway.full,pop.vec[i],"George",file.addon,sep='')
          }
        if (pop.analysis[i]!="FullAnalysis" && pop.analysis[i]!="FullAnalysisSplit")
          {
            file.addon=paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"Pops",date.summary.files.POST.vec.auxfiles,sep='')
            file.addonFST=strsplit(pop.analysis[i],split="Split")[[1]][1]
            ci.filein=paste(pathway.sep,pop.vec[i],"George",file.addon,sep='')
          }

        file.inI=paste(fst.filesPRE,pop.vec[i],file.addonFST,fst.filesPOST,sep='')
        file.inII=paste(fst.filesPRE,pop.vec[i],file.addonFST,fst.filesPOSTii,sep='')
        file.inIII=paste(fst.filesPRE,pop.vec[i],file.addonFST,fst.filesPOSTiii,sep='')
        file.inIV=paste(fst.filesPRE,pop.vec[i],file.addonFST,fst.filesPOSTiv,sep='')

               ## side-of-event info:
        #load(file.inIV)
        #fit.quality.mat=rbind(fit.quality.mat,c(fit.quality,fit.quality2))
        proportion.est.ci=NULL
        x=read.table(file.inII,nrows=1,skip=1)
        avrelmse=rbind(avrelmse,c(as.double(x[2]),NA))
	prop.val.est=as.double(x[1])
        proportion.est.mat=rbind(proportion.est.mat,c(prop.val.est,NA))
        x=read.table(file.inII,nrows=1,skip=3,as.is=TRUE)
        proportion.est.ci=c(as.double(strsplit(strsplit(as.character(x[1]),split="-")[[1]][1],split="(",fixed=TRUE)[[1]][2]),as.double(strsplit(strsplit(as.character(x[1]),split="-")[[1]][2],split=")",fixed=TRUE)[[1]][1]))

        x=read.table(file.inII,nrows=2,as.is=TRUE)
        side.of.event.1=as.double(x[2,][3:length(x[2,])])
        event.side.1.all=as.character(x[1,][3:length(x[1,])])
        x=read.table(file.inII,skip=2,nrows=2,as.is=TRUE)
        side.of.event.2=as.double(x[2,][3:length(x[2,])])
        event.side.2.all=as.character(x[1,][3:length(x[1,])])
    	donor.pops=unique(c(event.side.1.all,event.side.2.all))
    	ourmix.pop1=ourmix.pop2=rep(0,length(donor.pops))
    	ourmix.pop1[match(event.side.1.all,donor.pops)]=side.of.event.1
    	ourmix.pop2[match(event.side.2.all,donor.pops)]=side.of.event.2
    	ourmix=prop.val.est*ourmix.pop1+(1-prop.val.est)*ourmix.pop2
    	donor.pops.withcurves=donor.pops[ourmix>prop.cutoff]
        event.side.1.all=event.side.1.all[side.of.event.1>prop.cutoff]
        side.of.event.1=side.of.event.1[side.of.event.1>prop.cutoff]
        event.side.2.all=event.side.2.all[side.of.event.2>prop.cutoff]
        side.of.event.2=side.of.event.2[side.of.event.2>prop.cutoff]
	side.of.event.1=side.of.event.1/sum(side.of.event.1)
	side.of.event.2=side.of.event.2/sum(side.of.event.2)
	curves.ind.side1=rep("no",length(event.side.1.all))
	curves.ind.side1[match(intersect(event.side.1.all,donor.pops.withcurves),event.side.1.all)]="yes"
	curves.ind.side2=rep("no",length(event.side.2.all))
	curves.ind.side2[match(intersect(event.side.2.all,donor.pops.withcurves),event.side.2.all)]="yes"
        full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,strsplit(pop.analysis[i],split="Split")[[1]][1],"1",as.character(event.side.1.all)[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],proportion.est.mat[count,1]*side.of.event.1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,strsplit(pop.analysis[i],split="Split")[[1]][1],"2",as.character(event.side.2.all)[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],(1.0-proportion.est.mat[count,1])*side.of.event.2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
        event.side.1=event.side.1.all[side.of.event.1>prop.cutoff]
        event.side.2=event.side.2.all[side.of.event.2>prop.cutoff]
        if (length(event.side.1)==0) event.side.1=event.side.1.all[1]
        if (length(event.side.2)==0) event.side.2=event.side.2.all[1]
        event.side.1[event.side.1=="BantuSouthAfrica"]="BantuSA"
        event.side.1[event.side.1=="BantuKenya"]="BantuK"
        event.side.1[event.side.1=="Han.NChina"]="Han-NChina"
        event.side.2[event.side.2=="BantuSouthAfrica"]="BantuSA"
        event.side.2[event.side.2=="BantuKenya"]="BantuK"
        event.side.2[event.side.2=="Han.NChina"]="Han-NChina"
        event.side.1.final=event.side.1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix][1]
        if (length(event.side.1)>1 && max.donors>1)
          {
            for (k in 2:min(length(event.side.1),max.donors)) event.side.1.final=paste(event.side.1.final," ",event.side.1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix][k],sep='')
          }
        if (length(event.side.1)>max.donors) event.side.1.final=paste(event.side.1.final," +",sep='')
        event.side.2.final=event.side.2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix][1]
        if (length(event.side.2)>1 && max.donors>1)
          {
            for (k in 2:min(length(event.side.2),max.donors)) event.side.2.final=paste(event.side.2.final," ",event.side.2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix][k],sep='')
          }
        if (length(event.side.2)>max.donors) event.side.2.final=paste(event.side.2.final," +",sep='')
        #side.of.event.mat[i,(((k-1)*2+1):((k-1)*2+2))]=c(event.side.1.final,event.side.2.final)            
        proportion.est.mat.ci=rbind(proportion.est.mat.ci,c(proportion.est.ci,NA,NA))
          
                       ## 1-2 event "R-squared", final dates and intercepts:    
        x=read.table(file.inI,header=TRUE)
        intercepts.all=as.double(x$intercept.est)
        intercepts.all.fitted=as.double(x$intercept.est.fitted)
        rsquared.1event.vec.start=as.double(x$rsquared.1event)
        rsquared.2event.vec.start=as.double(x$r.squared.est.2event)
        rsquared.1event.vec=c(rsquared.1event.vec,round(mean(rsquared.1event.vec.start),round.r2))
        rsquared.2event.vec=c(rsquared.2event.vec,round(mean(rsquared.2event.vec.start),round.r2))
        #rsquared.1event.vec.max=c(rsquared.1event.vec.max,round(max(rsquared.1event.vec.start),round.r2))
        #rsquared.2event.vec.max=c(rsquared.2event.vec.max,round(max(rsquared.2event.vec.start),round.r2))
        #rsquared.2event.vec.max.rescaled=c(rsquared.2event.vec.max.rescaled,round(max((rsquared.2event.vec.start-rsquared.1event.vec.start)/(1.0-rsquared.1event.vec.start)),round.r2))

	rsquared.1event.vec.max=c(rsquared.1event.vec.max,round(rsquared.1event.max[i],round.r2))
	rsquared.2event.vec.max=c(rsquared.2event.vec.max,round(rsquared.2event.max[i],round.r2))
	rsquared.2event.vec.max.rescaled=c(rsquared.2event.vec.max.rescaled,round(r2.2event.minus.1event.rescaled.max[i],round.r2))
	fit.quality.mat=rbind(fit.quality.mat,c(fit.quality1[i],fit.quality2.vec[i]))
	fit.quality.rescaled=c(fit.quality.rescaled,fit.quality2.minus1.rescaled[i])
	grid.type.vec=c(grid.type.vec,grid.type[i])
	pvalue.vec.final=c(pvalue.vec.final,round(pvalue.vec[i],round.pval))
	twodate.pvalue.vec.final=c(twodate.pvalue.vec.final,round(twodate.pvalue.vec[i],round.pval))
	multiway.pvalue.vec.final=c(multiway.pvalue.vec.final,round(multiway.pvalue.vec[i],round.pval))
	min.fit.quality.vec.final=c(min.fit.quality.vec.final,round(min(c(fit.qualityBOTH.vec[i],fit.qualityBOTH.vecNULL[i])),round.fit))

               ## "best copy-vec" info:
	best.copyvec.matrix=rbind(best.copyvec.matrix,c(first.source.best.copyvec[i],second.source.best.copyvec[i]))
        
                        ## 1-date/2-date confidence intervals:
        x=read.table(ci.filein,skip=1,nrows=1,as.is=TRUE)
        date.est=as.double(x[5])
        date.est.ci=as.character(x[6])
        date.est.ciII=c(floor(as.double(strsplit(as.character(strsplit(as.character(x[6]),split="-")[[1]][1]),"(",fixed=TRUE)[[1]][2])),ceiling(as.double(strsplit(as.character(strsplit(as.character(x[6]),split="-")[[1]][2]),")",fixed=TRUE)[[1]][1])))
        date.est.mat=rbind(date.est.mat,c(date.est,NA))
        date.est.ci.mat=rbind(date.est.ci.mat,date.est.ci)
        date.est.ci.matII=rbind(date.est.ci.matII,c(date.est.ciII,NA,NA))
        
        count=count+1
        num.events.vec=c(num.events.vec,1)
      }

    if (pop.conclusion[i]=="two-dates" || pop.conclusion[i]=="two-dates-multiple-events" || pop.conclusion[i]=="two-dates-uncertain" || pop.conclusion[i]=="two-dates-multiple-events-uncertain")
      {
        if (pop.analysis[i]=="FullAnalysis" || pop.analysis[i]=="FullAnalysisSplit")
          {
            if (grid.type[i]!="two-date") dates.suffix=paste("GeorgePopsFinal",date.summary.files.POST.vec.2event,sep='')
            if (grid.type[i]=="two-date") dates.suffix=paste("GeorgePopsFinal",date.summary.files.POST.vec.2event.twodategrid,sep='')
            fst.filesPOST.1st=fst.filesPOST.2dates.1
            fst.filesPOST.2nd=fst.filesPOST.2dates.2
          }
        if (pop.analysis[i]!="FullAnalysis" && pop.analysis[i]!="FullAnalysisSplit")
          {
            if (grid.type[i]!="two-date") dates.suffix=paste("George",strsplit(pop.analysis[i],split="Split")[[1]][1],"Pops",date.summary.files.POST.vec.auxfiles.2event,sep='')
            if (grid.type[i]=="two-date") dates.suffix=paste("George",strsplit(pop.analysis[i],split="Split")[[1]][1],"Pops",date.summary.files.POST.vec.auxfiles.2event.twodategrid,sep='')
            
            fst.filesPOST.1st=paste(strsplit(pop.analysis[i],split="Split")[[1]][1],fst.filesPOST.2dates.1,sep='')
            fst.filesPOST.2nd=paste(strsplit(pop.analysis[i],split="Split")[[1]][1],fst.filesPOST.2dates.2,sep='')
          }

                       ## proportions:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=1,skip=1)
        y=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),nrows=1,skip=1)
        avrelmse=rbind(avrelmse,c(as.double(x[2]),as.double(y[2])),c(as.double(x[2]),as.double(y[2])),c(as.double(x[2]),as.double(y[2])))
        proportion.est.mat=rbind(proportion.est.mat,c(as.double(x[1]),as.double(y[1])),c(as.double(x[1]),as.double(y[1])),c(as.double(x[1]),as.double(y[1])))
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=1,skip=3,as.is=TRUE)
        proportion.est.ci.1st=c(as.double(strsplit(strsplit(as.character(x[1]),split="-")[[1]][1],split="(",fixed=TRUE)[[1]][2]),as.double(strsplit(strsplit(as.character(x[1]),split="-")[[1]][2],split=")",fixed=TRUE)[[1]][1]))
        y=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),nrows=1,skip=3,as.is=TRUE)
        proportion.est.ci.2nd=c(as.double(strsplit(strsplit(as.character(y[1]),split="-")[[1]][1],split="(",fixed=TRUE)[[1]][2]),as.double(strsplit(strsplit(as.character(y[1]),split="-")[[1]][2],split=")",fixed=TRUE)[[1]][1]))
        proportion.est.mat.ci=rbind(proportion.est.mat.ci,c(proportion.est.ci.1st,proportion.est.ci.2nd),c(proportion.est.ci.1st,proportion.est.ci.2nd),c(proportion.est.ci.1st,proportion.est.ci.2nd))

                       ## get raw copy-props (recipient):
	x=read.table(paste(rawprops.recipientPRE,pop.vec[i],rawprops.recipientPOST,sep=''),header=TRUE)
	x.names=read.table(paste(rawprops.recipientPRE,pop.vec[i],rawprops.recipientPOST,sep=''),as.is=TRUE,nrow=1)[-1]
	donor.pops.all=x.names[x.names!="Self"]
	raw.copyprops.mat=matrix(as.matrix(x[,2:dim(x)[2],]),ncol=dim(x)[2]-1)
	raw.copyprops.mat=raw.copyprops.mat[,match(donor.pops.all,x.names)]
	raw.copyprops.mat=raw.copyprops.mat/apply(raw.copyprops.mat,1,sum)
	raw.copyprops=apply(raw.copyprops.mat,2,mean)
	donor.pops.nhaps=full.pop.nhaps[match(donor.pops.all,full.pop.list)]

                       ## get raw copy-props (donors):
	predmat=NULL
	for (k in 1:length(donor.pops.all))
  	{
		x=read.table(paste(rawprops.donorsPRE,donor.pops.all[k],rawprops.donorsPOST,sep=''),header=TRUE)
    		x.names=read.table(paste(rawprops.donorsPRE,donor.pops.all[k],rawprops.donorsPOST,sep=''),as.is=TRUE,nrow=1)[-1]
    		x.names[x.names=="Self"]=donor.pops.all[k]
    		raw.copypropsHGDP.i=apply(matrix(as.matrix(x[,2:dim(x)[2],]),ncol=dim(x)[2]-1),2,mean)
    		predmat=rbind(predmat,raw.copypropsHGDP.i[match(donor.pops.all,x.names)])
  	}
	rownames(predmat)=donor.pops.all
	colnames(predmat)=donor.pops.all
	predmat.unscaled=predmat
	for (k in 1:nrow(predmat)) predmat[k,]=predmat[k,]/sum(predmat[k,])

                       ## side-of-event info, event 1:
                                              ## regular mixing props:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=2,as.is=TRUE)
        side.of.event.1=as.double(x[2,][3:length(x[2,])])
	prop.val.est=as.double(x[2,1])
        event.side.1.all=as.character(x[1,][3:length(x[1,])])
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),skip=2,nrows=2,as.is=TRUE)
        side.of.event.2=as.double(x[2,][3:length(x[2,])])
        event.side.2.all=as.character(x[1,][3:length(x[1,])])
    	donor.pops=unique(c(event.side.1.all,event.side.2.all))
    	ourmix.pop1=ourmix.pop2=rep(0,length(donor.pops))
    	ourmix.pop1[match(event.side.1.all,donor.pops)]=side.of.event.1
    	ourmix.pop2[match(event.side.2.all,donor.pops)]=side.of.event.2
    	ourmix=prop.val.est*ourmix.pop1+(1-prop.val.est)*ourmix.pop2
    	donor.pops.withcurves=donor.pops[ourmix>prop.cutoff]
        event.side.1.all=event.side.1.all[side.of.event.1>prop.cutoff]
        side.of.event.1=side.of.event.1[side.of.event.1>prop.cutoff]
        event.side.2.all=event.side.2.all[side.of.event.2>prop.cutoff]
        side.of.event.2=side.of.event.2[side.of.event.2>prop.cutoff]
	side.of.event.1=side.of.event.1/sum(side.of.event.1)
	side.of.event.2=side.of.event.2/sum(side.of.event.2)
	curves.ind.side1=rep("no",length(event.side.1.all))
	curves.ind.side1[match(intersect(event.side.1.all,donor.pops.withcurves),event.side.1.all)]="yes"
	curves.ind.side2=rep("no",length(event.side.2.all))
	curves.ind.side2[match(intersect(event.side.2.all,donor.pops.withcurves),event.side.2.all)]="yes"
        full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"1",as.character(event.side.1.all)[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],proportion.est.mat[count,1]*side.of.event.1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"2",as.character(event.side.2.all)[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],(1.0-proportion.est.mat[count,1])*side.of.event.2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
                                             ## copy-vec diff mixing props:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=2,as.is=TRUE)
        side.of.event.1.1st=as.double(x[2,][3:length(x[2,])])
	alpha.PC1=as.double(x[2,1])
        event.side.1.all.1st=as.character(x[1,][3:length(x[1,])])
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),skip=2,nrows=2,as.is=TRUE)
        side.of.event.2.1st=as.double(x[2,][3:length(x[2,])])
        event.side.2.all.1st=as.character(x[1,][3:length(x[1,])])
	copyvec.source1.PC1=copyvec.source2.PC1=rep(0,dim(predmat)[2])
	for (k in 1:length(event.side.1.all.1st))
  	    copyvec.source1.PC1=copyvec.source1.PC1+side.of.event.1.1st[k]*predmat[donor.pops.all==event.side.1.all.1st[k],]
	for (k in 1:length(event.side.2.all.1st))
  	    copyvec.source2.PC1=copyvec.source2.PC1+side.of.event.2.1st[k]*predmat[donor.pops.all==event.side.2.all.1st[k],]
	#copyvec.source1.PC1=copyvec.source1.PC1/(donor.pops.nhaps-2)
	#copyvec.source2.PC1=copyvec.source2.PC1/(donor.pops.nhaps-2)
	copyvec.diff.PC1=copyvec.source1.PC1-copyvec.source2.PC1
	copyvec.diff.PC1=copyvec.diff.PC1/sum(abs(copyvec.diff.PC1))
	copyvec.diff.PC1=copyvec.diff.PC1/(donor.pops.nhaps-2)
	#copyvec.diff.PC1=(copyvec.source1.PC1/sum(copyvec.source1.PC1))-(copyvec.source2.PC1/sum(copyvec.source2.PC1))
	#copyvec.diff.PC1=sqrt(alpha.PC1*(1-alpha.PC1))*copyvec.diff.PC1
	#donor.pops.withcurves.twodates.event1=donor.pops.all[abs(copyvec.diff.PC1)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC1)>side.of.event.cutoff)>0) donor.pops.withcurves.twodates.event1=donor.pops.all[abs(copyvec.diff.PC1)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC1)>side.of.event.cutoff)==0) donor.pops.withcurves.twodates.event1=c(donor.pops.all[copyvec.diff.PC1==min(copyvec.diff.PC1)],donor.pops.all[copyvec.diff.PC1==max(copyvec.diff.PC1)])
        event.side.1.all.1st=donor.pops.all[copyvec.diff.PC1>side.of.event.cutoff]
	side.of.event.1.1st=copyvec.diff.PC1[copyvec.diff.PC1>side.of.event.cutoff]
	if (length(side.of.event.1.1st)>max.pops.multievent)
	{
		event.side.1.all.1st=event.side.1.all.1st[side.of.event.1.1st>=sort(side.of.event.1.1st,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.1.1st=side.of.event.1.1st[side.of.event.1.1st>=sort(side.of.event.1.1st,decreasing=TRUE)[max.pops.multievent]]
	}
        event.side.2.all.1st=donor.pops.all[copyvec.diff.PC1 < -1*side.of.event.cutoff]
	side.of.event.2.1st=abs(copyvec.diff.PC1[copyvec.diff.PC1 < -1*side.of.event.cutoff])
	if (length(side.of.event.2.1st)>max.pops.multievent)
	{
		event.side.2.all.1st=event.side.2.all.1st[side.of.event.2.1st>=sort(side.of.event.2.1st,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.2.1st=side.of.event.2.1st[side.of.event.2.1st>=sort(side.of.event.2.1st,decreasing=TRUE)[max.pops.multievent]]
	}
	side.of.event.1.1st=side.of.event.1.1st/sum(side.of.event.1.1st);side.of.event.2.1st=side.of.event.2.1st/sum(side.of.event.2.1st)
	side.of.event.1.1st=side.of.event.1.1st[side.of.event.1.1st>prop.cutoff]/sum(side.of.event.1.1st[side.of.event.1.1st>prop.cutoff]);side.of.event.2.1st=side.of.event.2.1st[side.of.event.2.1st>prop.cutoff]/sum(side.of.event.2.1st[side.of.event.2.1st>prop.cutoff]);
	donor.pops.withcurves.twodates.event1=unique(c(event.side.1.all.1st,event.side.2.all.1st))

                       ## side-of-event info, event 2:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),nrows=2,as.is=TRUE)
        side.of.event.1.2nd=as.double(x[2,][3:length(x[2,])])
	alpha.PC2=as.double(x[2,1])
        event.side.1.all.2nd=as.character(x[1,][3:length(x[1,])])
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),skip=2,nrows=2,as.is=TRUE)
        side.of.event.2.2nd=as.double(x[2,][3:length(x[2,])])
        event.side.2.all.2nd=as.character(x[1,][3:length(x[1,])])
	copyvec.source1.PC2=copyvec.source2.PC2=rep(0,dim(predmat)[2])
	for (k in 1:length(event.side.1.all.2nd))
  	    copyvec.source1.PC2=copyvec.source1.PC2+side.of.event.1.2nd[k]*predmat[donor.pops.all==event.side.1.all.2nd[k],]
	for (k in 1:length(event.side.2.all.2nd))
  	    copyvec.source2.PC2=copyvec.source2.PC2+side.of.event.2.2nd[k]*predmat[donor.pops.all==event.side.2.all.2nd[k],]
	#copyvec.source1.PC2=copyvec.source1.PC2/(donor.pops.nhaps-2)
	#copyvec.source2.PC2=copyvec.source2.PC2/(donor.pops.nhaps-2)
	#copyvec.diff.PC2=(copyvec.source1.PC2/sum(copyvec.source1.PC2))-(copyvec.source2.PC2/sum(copyvec.source2.PC2))
	copyvec.diff.PC2=copyvec.source1.PC2-copyvec.source2.PC2
	copyvec.diff.PC2=copyvec.diff.PC2/sum(abs(copyvec.diff.PC2))
	copyvec.diff.PC2=copyvec.diff.PC2/(donor.pops.nhaps-2)
	#copyvec.diff.PC2=sqrt(alpha.PC2*(1-alpha.PC2))*copyvec.diff.PC2
	#donor.pops.withcurves.twodates.event2=donor.pops.all[abs(copyvec.diff.PC2)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC2)>side.of.event.cutoff)>0) donor.pops.withcurves.twodates.event2=donor.pops.all[abs(copyvec.diff.PC2)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC2)>side.of.event.cutoff)==0) donor.pops.withcurves.twodates.event2=c(donor.pops.all[copyvec.diff.PC2==min(copyvec.diff.PC2)],donor.pops.all[copyvec.diff.PC2==max(copyvec.diff.PC2)])
        event.side.1.all.2nd=donor.pops.all[copyvec.diff.PC2>side.of.event.cutoff]
	side.of.event.1.2nd=copyvec.diff.PC2[copyvec.diff.PC2>side.of.event.cutoff]
	if (length(side.of.event.1.2nd)>max.pops.multievent)
	{
		event.side.1.all.2nd=event.side.1.all.2nd[side.of.event.1.2nd>=sort(side.of.event.1.2nd,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.1.2nd=side.of.event.1.2nd[side.of.event.1.2nd>=sort(side.of.event.1.2nd,decreasing=TRUE)[max.pops.multievent]]
	}
        event.side.2.all.2nd=donor.pops.all[copyvec.diff.PC2 < -1*side.of.event.cutoff]
	side.of.event.2.2nd=abs(copyvec.diff.PC2[copyvec.diff.PC2 < -1*side.of.event.cutoff])
	if (length(side.of.event.2.2nd)>max.pops.multievent)
	{
		event.side.2.all.2nd=event.side.2.all.2nd[side.of.event.2.2nd>=sort(side.of.event.2.2nd,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.2.2nd=side.of.event.2.2nd[side.of.event.2.2nd>=sort(side.of.event.2.2nd,decreasing=TRUE)[max.pops.multievent]]
	}
	side.of.event.1.2nd=side.of.event.1.2nd/sum(side.of.event.1.2nd);side.of.event.2.2nd=side.of.event.2.2nd/sum(side.of.event.2.2nd)
	side.of.event.1.2nd=side.of.event.1.2nd[side.of.event.1.2nd>prop.cutoff]/sum(side.of.event.1.2nd[side.of.event.1.2nd>prop.cutoff]);side.of.event.2.2nd=side.of.event.2.2nd[side.of.event.2.2nd>prop.cutoff]/sum(side.of.event.2.2nd[side.of.event.2.2nd>prop.cutoff]);
	donor.pops.withcurves.twodates.event2=unique(c(event.side.1.all.2nd,event.side.2.all.2nd))

                       ## final curve/proportion info:
	curves.ind.side1.1st=rep("no",length(event.side.1.all.1st))
	if (print.curves.secondevent!='yes') curves.ind.side1.1st[match(intersect(event.side.1.all.1st,donor.pops.withcurves),event.side.1.all.1st)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side1.1st[match(intersect(event.side.1.all.1st,donor.pops.withcurves.twodates.event1),event.side.1.all.1st)]="yes"
	curves.ind.side2.1st=rep("no",length(event.side.2.all.1st))
	if (print.curves.secondevent!='yes') curves.ind.side2.1st[match(intersect(event.side.2.all.1st,donor.pops.withcurves),event.side.2.all.1st)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side2.1st[match(intersect(event.side.2.all.1st,donor.pops.withcurves.twodates.event1),event.side.2.all.1st)]="yes"
	curves.ind.side1.2nd=rep("no",length(event.side.1.all.2nd))
	if (print.curves.secondevent!='yes') curves.ind.side1.2nd[match(intersect(event.side.1.all.2nd,donor.pops.withcurves),event.side.1.all.2md)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side1.2nd[match(intersect(event.side.1.all.2nd,donor.pops.withcurves.twodates.event2),event.side.1.all.2nd)]="yes"
	curves.ind.side2.2nd=rep("no",length(event.side.2.all.2nd))
	if (print.curves.secondevent!='yes') curves.ind.side2.2nd[match(intersect(event.side.2.all.2nd,donor.pops.withcurves),event.side.2.all.2nd)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side2.2nd[match(intersect(event.side.2.all.2nd,donor.pops.withcurves.twodates.event2),event.side.2.all.2nd)]="yes"
        full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"3",as.character(event.side.1.all.1st)[sort(side.of.event.1.1st,index.return=TRUE,decreasing=TRUE)$ix],proportion.est.mat[count,1]*side.of.event.1.1st[sort(side.of.event.1.1st,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1.1st[sort(side.of.event.1.1st,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"4",as.character(event.side.2.all.1st)[sort(side.of.event.2.1st,index.return=TRUE,decreasing=TRUE)$ix],(1.0-proportion.est.mat[count,1])*side.of.event.2.1st[sort(side.of.event.2.1st,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2.1st[sort(side.of.event.2.1st,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
        full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"5",as.character(event.side.1.all.2nd)[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],proportion.est.mat[count,2]*side.of.event.1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"6",as.character(event.side.2.all.2nd)[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],(1.0-proportion.est.mat[count,2])*side.of.event.2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
	#full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"3",as.character(event.side.1.all.2nd)[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],0.5*side.of.event.1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"2Event",sep=''),"4",as.character(event.side.2.all.2nd)[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],0.5*side.of.event.2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
      
                       ## 1-2 event "R-squared", final dates and intercepts:
        if (pop.analysis[i]=="FullAnalysis" || pop.analysis[i]=="FullAnalysisSplit") x=read.table(paste(pathway.full,pop.vec[i],dates.suffix,sep=''),skip=1,as.is=TRUE)
        if (pop.analysis[i]!="FullAnalysis" && pop.analysis[i]!="FullAnalysisSplit") x=read.table(paste(pathway.sep,pop.vec[i],dates.suffix,sep=''),skip=1,as.is=TRUE)
        date.est=c(x[1,5],x[1,7])
        date.est.ci=as.character(x[1,6])
        date.est.ciII=c(floor(as.double(strsplit(as.character(strsplit(as.character(x[1,6]),split="-")[[1]][1]),"(",fixed=TRUE)[[1]][2])),ceiling(as.double(strsplit(as.character(strsplit(as.character(x[1,6]),split="-")[[1]][2]),")",fixed=TRUE)[[1]][1])),floor(as.double(strsplit(as.character(strsplit(as.character(x[1,8]),split="-")[[1]][1]),"(",fixed=TRUE)[[1]][2])),ceiling(as.double(strsplit(as.character(strsplit(as.character(x[1,8]),split="-")[[1]][2]),")",fixed=TRUE)[[1]][1])))
        date.est.mat=rbind(date.est.mat,date.est,date.est,date.est)
        date.est.ci.mat=rbind(date.est.ci.mat,date.est.ci,date.est.ci,date.est.ci)
        date.est.ci.matII=rbind(date.est.ci.matII,date.est.ciII,date.est.ciII,date.est.ciII)
        #rsquared.1event.vec.start=as.double(x[,3])
        #rsquared.2event.vec.start=as.double(x[,4])
        #rsquared.1event.vec=c(rsquared.1event.vec,round(mean(rsquared.1event.vec.start),round.r2),round(mean(rsquared.1event.vec.start),round.r2))
        #rsquared.2event.vec=c(rsquared.2event.vec,round(mean(rsquared.2event.vec.start),round.r2),round(mean(rsquared.2event.vec.start),round.r2))
        #rsquared.1event.vec.max=c(rsquared.1event.vec.max,round(max(rsquared.1event.vec.start),round.r2),round(max(rsquared.1event.vec.start),round.r2))
        #rsquared.2event.vec.max=c(rsquared.2event.vec.max,round(max(rsquared.2event.vec.start),round.r2),round(max(rsquared.2event.vec.start),round.r2))
        #rsquared.2event.vec.max.rescaled=c(rsquared.2event.vec.max.rescaled,round(max((rsquared.2event.vec.start-rsquared.1event.vec.start)/(1.0-rsquared.1event.vec.start)),round.r2),round(max((rsquared.2event.vec.start-rsquared.1event.vec.start)/(1.0-rsquared.1event.vec.start)),round.r2))

	rsquared.1event.vec.max=c(rsquared.1event.vec.max,round(rsquared.1event.max[i],round.r2),round(rsquared.1event.max[i],round.r2),round(rsquared.1event.max[i],round.r2))
	rsquared.2event.vec.max=c(rsquared.2event.vec.max,round(rsquared.2event.max[i],round.r2),round(rsquared.2event.max[i],round.r2),round(rsquared.2event.max[i],round.r2))
	rsquared.2event.vec.max.rescaled=c(rsquared.2event.vec.max.rescaled,round(r2.2event.minus.1event.rescaled.max[i],round.r2),round(r2.2event.minus.1event.rescaled.max[i],round.r2),round(r2.2event.minus.1event.rescaled.max[i],round.r2))
	fit.quality.mat=rbind(fit.quality.mat,c(fit.quality1[i],fit.quality2.vec[i]),c(fit.quality1[i],fit.quality2.vec[i]),c(fit.quality1[i],fit.quality2.vec[i]))
	fit.quality.rescaled=c(fit.quality.rescaled,fit.quality2.minus1.rescaled[i],fit.quality2.minus1.rescaled[i],fit.quality2.minus1.rescaled[i])
	grid.type.vec=c(grid.type.vec,grid.type[i],grid.type[i],grid.type[i])
	pvalue.vec.final=c(pvalue.vec.final,round(pvalue.vec[i],round.pval),round(pvalue.vec[i],round.pval),round(pvalue.vec[i],round.pval))
	twodate.pvalue.vec.final=c(twodate.pvalue.vec.final,round(twodate.pvalue.vec[i],round.pval),round(twodate.pvalue.vec[i],round.pval),round(twodate.pvalue.vec[i],round.pval))
	multiway.pvalue.vec.final=c(multiway.pvalue.vec.final,round(multiway.pvalue.vec[i],round.pval),round(multiway.pvalue.vec[i],round.pval),round(multiway.pvalue.vec[i],round.pval))
	min.fit.quality.vec.final=c(min.fit.quality.vec.final,round(min(c(fit.qualityBOTH.vec[i],fit.qualityBOTH.vecNULL[i])),round.fit),round(min(c(fit.qualityBOTH.vec[i],fit.qualityBOTH.vecNULL[i])),round.fit),round(min(c(fit.qualityBOTH.vec[i],fit.qualityBOTH.vecNULL[i])),round.fit))

               ## "best copy-vec" info:
	#best.copyvec.matrix=rbind(best.copyvec.matrix,c(paste(first.source.best.copyvec.PC1[i],"??",sep=''),paste(second.source.best.copyvec.PC1[i],"??",sep='')),c(paste(first.source.best.copyvec.PC2[i],"??",sep=''),paste(second.source.best.copyvec.PC2[i],"??",sep='')))
	best.copyvec.matrix=rbind(best.copyvec.matrix,c(strsplit(first.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1],strsplit(second.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1]),c(strsplit(first.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1],strsplit(second.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1]),c("??","??"))

        count=count+3
        num.events.vec=c(num.events.vec,2,2,2)
      }

    if (pop.conclusion[i]=="multiple-events" || pop.conclusion[i]=="complex?")
      {
        if (pop.analysis[i]=="FullAnalysis" || pop.analysis[i]=="FullAnalysisSplit")
          {
            file.addonFST=""
            file.addon=paste("PopsFinal",date.summary.files.POST.vec,sep='')
            fst.filesPOST.1st=fst.filesPOST.multievents.1
            fst.filesPOST.2nd=fst.filesPOST.multievents.2
            ci.filein=paste(pathway.full,pop.vec[i],"George",file.addon,sep='')
          }
        if (pop.analysis[i]!="FullAnalysis" && pop.analysis[i]!="FullAnalysisSplit")
          {
         
            file.addon=paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"Pops",date.summary.files.POST.vec.auxfiles,sep='')
            file.addonFST=strsplit(pop.analysis[i],split="Split")[[1]][1]
            fst.filesPOST.1st=paste(strsplit(pop.analysis[i],split="Split")[[1]][1],fst.filesPOST.multievents.1,sep='')
            fst.filesPOST.2nd=paste(strsplit(pop.analysis[i],split="Split")[[1]][1],fst.filesPOST.multievents.2,sep='')
            ci.filein=paste(pathway.sep,pop.vec[i],"George",file.addon,sep='')
          }
        file.inI=paste(fst.filesPRE,pop.vec[i],file.addonFST,fst.filesPOST,sep='')

                       ## proportions:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=1,skip=1)
        y=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),nrows=1,skip=1)
        avrelmse=rbind(avrelmse,c(as.double(x[2]),as.double(y[2])),c(as.double(x[2]),as.double(y[2])),c(as.double(x[2]),as.double(y[2])))
        proportion.est.mat=rbind(proportion.est.mat,c(as.double(x[1]),as.double(y[1])),c(as.double(x[1]),as.double(y[1])),c(as.double(x[1]),as.double(y[1])))
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=1,skip=3,as.is=TRUE)
        proportion.est.ci.1st=c(as.double(strsplit(strsplit(as.character(x[1]),split="-")[[1]][1],split="(",fixed=TRUE)[[1]][2]),as.double(strsplit(strsplit(as.character(x[1]),split="-")[[1]][2],split=")",fixed=TRUE)[[1]][1]))
        y=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),nrows=1,skip=3,as.is=TRUE)
        proportion.est.ci.2nd=c(as.double(strsplit(strsplit(as.character(y[1]),split="-")[[1]][1],split="(",fixed=TRUE)[[1]][2]),as.double(strsplit(strsplit(as.character(y[1]),split="-")[[1]][2],split=")",fixed=TRUE)[[1]][1]))
        proportion.est.mat.ci=rbind(proportion.est.mat.ci,c(proportion.est.ci.1st,proportion.est.ci.2nd),c(proportion.est.ci.1st,proportion.est.ci.2nd),c(proportion.est.ci.1st,proportion.est.ci.2nd))

                       ## get raw copy-props (recipient):
	x=read.table(paste(rawprops.recipientPRE,pop.vec[i],rawprops.recipientPOST,sep=''),header=TRUE)
	x.names=read.table(paste(rawprops.recipientPRE,pop.vec[i],rawprops.recipientPOST,sep=''),as.is=TRUE,nrow=1)[-1]
	donor.pops.all=x.names[x.names!="Self"]
	raw.copyprops.mat=matrix(as.matrix(x[,2:dim(x)[2],]),ncol=dim(x)[2]-1)
	raw.copyprops.mat=raw.copyprops.mat[,match(donor.pops.all,x.names)]
	raw.copyprops.mat=raw.copyprops.mat/apply(raw.copyprops.mat,1,sum)
	raw.copyprops=apply(raw.copyprops.mat,2,mean)
	donor.pops.nhaps=full.pop.nhaps[match(donor.pops.all,full.pop.list)]

                       ## get raw copy-props (donors):
	predmat=NULL
	for (k in 1:length(donor.pops.all))
  	{
		x=read.table(paste(rawprops.donorsPRE,donor.pops.all[k],rawprops.donorsPOST,sep=''),header=TRUE)
    		x.names=read.table(paste(rawprops.donorsPRE,donor.pops.all[k],rawprops.donorsPOST,sep=''),as.is=TRUE,nrow=1)[-1]
    		x.names[x.names=="Self"]=donor.pops.all[k]
    		raw.copypropsHGDP.i=apply(matrix(as.matrix(x[,2:dim(x)[2],]),ncol=dim(x)[2]-1),2,mean)
    		predmat=rbind(predmat,raw.copypropsHGDP.i[match(donor.pops.all,x.names)])
  	}
	rownames(predmat)=donor.pops.all
	colnames(predmat)=donor.pops.all
	predmat.unscaled=predmat
	for (k in 1:nrow(predmat)) predmat[k,]=predmat[k,]/sum(predmat[k,])

                       ## side-of-event info, event 1:
                                              ## regular mixing props:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=2,as.is=TRUE)
        side.of.event.1=as.double(x[2,][3:length(x[2,])])
	prop.val.est=as.double(x[2,1])
        event.side.1.all=as.character(x[1,][3:length(x[1,])])
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),skip=2,nrows=2,as.is=TRUE)
        side.of.event.2=as.double(x[2,][3:length(x[2,])])
        event.side.2.all=as.character(x[1,][3:length(x[1,])])
    	donor.pops=unique(c(event.side.1.all,event.side.2.all))
    	ourmix.pop1=ourmix.pop2=rep(0,length(donor.pops))
    	ourmix.pop1[match(event.side.1.all,donor.pops)]=side.of.event.1
    	ourmix.pop2[match(event.side.2.all,donor.pops)]=side.of.event.2
    	ourmix=prop.val.est*ourmix.pop1+(1-prop.val.est)*ourmix.pop2
    	donor.pops.withcurves=donor.pops[ourmix>prop.cutoff]
        event.side.1.all=event.side.1.all[side.of.event.1>prop.cutoff]
        side.of.event.1=side.of.event.1[side.of.event.1>prop.cutoff]
        event.side.2.all=event.side.2.all[side.of.event.2>prop.cutoff]
        side.of.event.2=side.of.event.2[side.of.event.2>prop.cutoff]
	side.of.event.1=side.of.event.1/sum(side.of.event.1)
	side.of.event.2=side.of.event.2/sum(side.of.event.2)
	curves.ind.side1=rep("no",length(event.side.1.all))
	curves.ind.side1[match(intersect(event.side.1.all,donor.pops.withcurves),event.side.1.all)]="yes"
	curves.ind.side2=rep("no",length(event.side.2.all))
	curves.ind.side2[match(intersect(event.side.2.all,donor.pops.withcurves),event.side.2.all)]="yes"
        full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"1",as.character(event.side.1.all)[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],proportion.est.mat[count,1]*side.of.event.1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1[sort(side.of.event.1,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"2",as.character(event.side.2.all)[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],(1.0-proportion.est.mat[count,1])*side.of.event.2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2[sort(side.of.event.2,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
                                             ## copy-vec diff mixing props:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),nrows=2,as.is=TRUE)
        side.of.event.1.1st=as.double(x[2,][3:length(x[2,])])
	alpha.PC1=as.double(x[2,1])
        event.side.1.all.1st=as.character(x[1,][3:length(x[1,])])
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.1st,sep=''),skip=2,nrows=2,as.is=TRUE)
        side.of.event.2.1st=as.double(x[2,][3:length(x[2,])])
        event.side.2.all.1st=as.character(x[1,][3:length(x[1,])])
	copyvec.source1.PC1=copyvec.source2.PC1=rep(0,dim(predmat)[2])
	for (k in 1:length(event.side.1.all.1st))
  	    copyvec.source1.PC1=copyvec.source1.PC1+side.of.event.1.1st[k]*predmat[donor.pops.all==event.side.1.all.1st[k],]
	for (k in 1:length(event.side.2.all.1st))
  	    copyvec.source2.PC1=copyvec.source2.PC1+side.of.event.2.1st[k]*predmat[donor.pops.all==event.side.2.all.1st[k],]
	#copyvec.source1.PC1=copyvec.source1.PC1/(donor.pops.nhaps-2)
	#copyvec.source2.PC1=copyvec.source2.PC1/(donor.pops.nhaps-2)
	#copyvec.diff.PC1=(copyvec.source1.PC1/sum(copyvec.source1.PC1))-(copyvec.source2.PC1/sum(copyvec.source2.PC1))
	copyvec.diff.PC1=copyvec.source1.PC1-copyvec.source2.PC1
	copyvec.diff.PC1=copyvec.diff.PC1/sum(abs(copyvec.diff.PC1))
	copyvec.diff.PC1=copyvec.diff.PC1/(donor.pops.nhaps-2)
	#copyvec.diff.PC1=sqrt(alpha.PC1*(1-alpha.PC1))*copyvec.diff.PC1
	#donor.pops.withcurves.multievents.event1=donor.pops.all[abs(copyvec.diff.PC1)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC1)>side.of.event.cutoff)>0) donor.pops.withcurves.multievents.event1=donor.pops.all[abs(copyvec.diff.PC1)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC1)>side.of.event.cutoff)==0) donor.pops.withcurves.multievents.event1=c(donor.pops.all[copyvec.diff.PC1==min(copyvec.diff.PC1)],donor.pops.all[copyvec.diff.PC1==max(copyvec.diff.PC1)])
        event.side.1.all.1st=donor.pops.all[copyvec.diff.PC1>side.of.event.cutoff]
	side.of.event.1.1st=copyvec.diff.PC1[copyvec.diff.PC1>side.of.event.cutoff]
	if (length(side.of.event.1.1st)>max.pops.multievent)
	{
		event.side.1.all.1st=event.side.1.all.1st[side.of.event.1.1st>=sort(side.of.event.1.1st,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.1.1st=side.of.event.1.1st[side.of.event.1.1st>=sort(side.of.event.1.1st,decreasing=TRUE)[max.pops.multievent]]
	}
        event.side.2.all.1st=donor.pops.all[copyvec.diff.PC1 < -1*side.of.event.cutoff]
	side.of.event.2.1st=abs(copyvec.diff.PC1[copyvec.diff.PC1 < -1*side.of.event.cutoff])
	if (length(side.of.event.2.1st)>max.pops.multievent)
	{
		event.side.2.all.1st=event.side.2.all.1st[side.of.event.2.1st>=sort(side.of.event.2.1st,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.2.1st=side.of.event.2.1st[side.of.event.2.1st>=sort(side.of.event.2.1st,decreasing=TRUE)[max.pops.multievent]]
	}
	side.of.event.1.1st=side.of.event.1.1st/sum(side.of.event.1.1st);side.of.event.2.1st=side.of.event.2.1st/sum(side.of.event.2.1st)
	side.of.event.1.1st=side.of.event.1.1st[side.of.event.1.1st>prop.cutoff]/sum(side.of.event.1.1st[side.of.event.1.1st>prop.cutoff]);side.of.event.2.1st=side.of.event.2.1st[side.of.event.2.1st>prop.cutoff]/sum(side.of.event.2.1st[side.of.event.2.1st>prop.cutoff]);
	donor.pops.withcurves.multievents.event1=unique(c(event.side.1.all.1st,event.side.2.all.1st))

                       ## side-of-event info, event 2:
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),nrows=2,as.is=TRUE)
        side.of.event.1.2nd=as.double(x[2,][3:length(x[2,])])
	alpha.PC2=as.double(x[2,1])
        event.side.1.all.2nd=as.character(x[1,][3:length(x[1,])])
        x=read.table(paste(fst.filesPRE,pop.vec[i],fst.filesPOST.2nd,sep=''),skip=2,nrows=2,as.is=TRUE)
        side.of.event.2.2nd=as.double(x[2,][3:length(x[2,])])
        event.side.2.all.2nd=as.character(x[1,][3:length(x[1,])])
	copyvec.source1.PC2=copyvec.source2.PC2=rep(0,dim(predmat)[2])
	for (k in 1:length(event.side.1.all.2nd))
  	    copyvec.source1.PC2=copyvec.source1.PC2+side.of.event.1.2nd[k]*predmat[donor.pops.all==event.side.1.all.2nd[k],]
	for (k in 1:length(event.side.2.all.2nd))
  	    copyvec.source2.PC2=copyvec.source2.PC2+side.of.event.2.2nd[k]*predmat[donor.pops.all==event.side.2.all.2nd[k],]
	#copyvec.source1.PC2=copyvec.source1.PC2/(donor.pops.nhaps-2)
	#copyvec.source2.PC2=copyvec.source2.PC2/(donor.pops.nhaps-2)
	#copyvec.diff.PC2=(copyvec.source1.PC2/sum(copyvec.source1.PC2))-(copyvec.source2.PC2/sum(copyvec.source2.PC2))
	copyvec.diff.PC2=copyvec.source1.PC2-copyvec.source2.PC2
	copyvec.diff.PC2=copyvec.diff.PC2/sum(abs(copyvec.diff.PC2))
	copyvec.diff.PC2=copyvec.diff.PC2/(donor.pops.nhaps-2)
	#copyvec.diff.PC2=sqrt(alpha.PC2*(1-alpha.PC2))*copyvec.diff.PC2
	#donor.pops.withcurves.multievents.event2=donor.pops.all[abs(copyvec.diff.PC2)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC2)>side.of.event.cutoff)>0) donor.pops.withcurves.multievents.event2=donor.pops.all[abs(copyvec.diff.PC2)>side.of.event.cutoff]
    	#if (sum(abs(copyvec.diff.PC2)>side.of.event.cutoff)==0) donor.pops.withcurves.multievents.event2=c(donor.pops.all[copyvec.diff.PC2==min(copyvec.diff.PC2)],donor.pops.all[copyvec.diff.PC2==max(copyvec.diff.PC2)])
        event.side.1.all.2nd=donor.pops.all[copyvec.diff.PC2>side.of.event.cutoff]
	side.of.event.1.2nd=copyvec.diff.PC2[copyvec.diff.PC2>side.of.event.cutoff]
	if (length(side.of.event.1.2nd)>max.pops.multievent)
	{
		event.side.1.all.2nd=event.side.1.all.2nd[side.of.event.1.2nd>=sort(side.of.event.1.2nd,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.1.2nd=side.of.event.1.2nd[side.of.event.1.2nd>=sort(side.of.event.1.2nd,decreasing=TRUE)[max.pops.multievent]]
	}
        event.side.2.all.2nd=donor.pops.all[copyvec.diff.PC2 < -1*side.of.event.cutoff]
	side.of.event.2.2nd=abs(copyvec.diff.PC2[copyvec.diff.PC2 < -1*side.of.event.cutoff])
	if (length(side.of.event.2.2nd)>max.pops.multievent)
	{
		event.side.2.all.2nd=event.side.2.all.2nd[side.of.event.2.2nd>=sort(side.of.event.2.2nd,decreasing=TRUE)[max.pops.multievent]]
		side.of.event.2.2nd=side.of.event.2.2nd[side.of.event.2.2nd>=sort(side.of.event.2.2nd,decreasing=TRUE)[max.pops.multievent]]
	}
	side.of.event.1.2nd=side.of.event.1.2nd/sum(side.of.event.1.2nd);side.of.event.2.2nd=side.of.event.2.2nd/sum(side.of.event.2.2nd)
	side.of.event.1.2nd=side.of.event.1.2nd[side.of.event.1.2nd>prop.cutoff]/sum(side.of.event.1.2nd[side.of.event.1.2nd>prop.cutoff]);side.of.event.2.2nd=side.of.event.2.2nd[side.of.event.2.2nd>prop.cutoff]/sum(side.of.event.2.2nd[side.of.event.2.2nd>prop.cutoff]);
	donor.pops.withcurves.multievents.event2=unique(c(event.side.1.all.2nd,event.side.2.all.2nd))

                       ## final curve/proportion info:
	curves.ind.side1.1st=rep("no",length(event.side.1.all.1st))
	if (print.curves.secondevent!='yes') curves.ind.side1.1st[match(intersect(event.side.1.all.1st,donor.pops.withcurves),event.side.1.all.1st)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side1.1st[match(intersect(event.side.1.all.1st,donor.pops.withcurves.multievents.event1),event.side.1.all.1st)]="yes"
	curves.ind.side2.1st=rep("no",length(event.side.2.all.1st))
	if (print.curves.secondevent!='yes') curves.ind.side2.1st[match(intersect(event.side.2.all.1st,donor.pops.withcurves),event.side.2.all.1st)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side2.1st[match(intersect(event.side.2.all.1st,donor.pops.withcurves.multievents.event1),event.side.2.all.1st)]="yes"
	curves.ind.side1.2nd=rep("no",length(event.side.1.all.2nd))
	if (print.curves.secondevent!='yes') curves.ind.side1.2nd[match(intersect(event.side.1.all.2nd,donor.pops.withcurves),event.side.1.all.2md)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side1.2nd[match(intersect(event.side.1.all.2nd,donor.pops.withcurves.multievents.event2),event.side.1.all.2nd)]="yes"
	curves.ind.side2.2nd=rep("no",length(event.side.2.all.2nd))
	if (print.curves.secondevent!='yes') curves.ind.side2.2nd[match(intersect(event.side.2.all.2nd,donor.pops.withcurves),event.side.2.all.2nd)]="yes"
	if (print.curves.secondevent=='yes') curves.ind.side2.2nd[match(intersect(event.side.2.all.2nd,donor.pops.withcurves.multievents.event2),event.side.2.all.2nd)]="yes"
         full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"3",as.character(event.side.1.all.1st)[sort(side.of.event.1.1st,index.return=TRUE,decreasing=TRUE)$ix],proportion.est.mat[count,1]*side.of.event.1.1st[sort(side.of.event.1.1st,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1.1st[sort(side.of.event.1.1st,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"4",as.character(event.side.2.all.1st)[sort(side.of.event.2.1st,index.return=TRUE,decreasing=TRUE)$ix],(1.0-proportion.est.mat[count,1])*side.of.event.2.1st[sort(side.of.event.2.1st,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2.1st[sort(side.of.event.2.1st,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
        full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"5",as.character(event.side.1.all.2nd)[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],proportion.est.mat[count,2]*side.of.event.1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"6",as.character(event.side.2.all.2nd)[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],(1.0-proportion.est.mat[count,2])*side.of.event.2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
	#full.info.mat=rbind(full.info.mat,cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"3",as.character(event.side.1.all.2nd)[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],0.5*side.of.event.1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side1.2nd[sort(side.of.event.1.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]),cbind(as.character(pop.vec[i]),split.pop,paste(strsplit(pop.analysis[i],split="Split")[[1]][1],"MultiEvent",sep=''),"4",as.character(event.side.2.all.2nd)[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],0.5*side.of.event.2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],curves.ind.side2.2nd[sort(side.of.event.2.2nd,index.return=TRUE,decreasing=TRUE)$ix],pop.conclusion[i]))
     
                       ## 1-2 event "R-squared", final dates and intercepts:    
        x=read.table(file.inI,header=TRUE)
        intercepts.all=as.double(x$intercept.est)
        intercepts.all.fitted=as.double(x$intercept.est.fitted)
        rsquared.1event.vec.start=as.double(x$rsquared.1event)
        rsquared.2event.vec.start=as.double(x$r.squared.est.2event)
        rsquared.1event.vec=c(rsquared.1event.vec,round(mean(rsquared.1event.vec.start),round.r2),round(mean(rsquared.1event.vec.start),round.r2),round(mean(rsquared.1event.vec.start),round.r2))
        rsquared.2event.vec=c(rsquared.2event.vec,round(mean(rsquared.2event.vec.start),round.r2),round(mean(rsquared.2event.vec.start),round.r2),round(mean(rsquared.2event.vec.start),round.r2))
        #rsquared.1event.vec.max=c(rsquared.1event.vec.max,round(max(rsquared.1event.vec.start),round.r2),round(max(rsquared.1event.vec.start),round.r2))
        #rsquared.2event.vec.max=c(rsquared.2event.vec.max,round(max(rsquared.2event.vec.start),round.r2),round(max(rsquared.2event.vec.start),round.r2))
        #rsquared.2event.vec.max.rescaled=c(rsquared.2event.vec.max.rescaled,round(max((rsquared.2event.vec.start-rsquared.1event.vec.start)/(1.0-rsquared.1event.vec.start)),round.r2),round(max((rsquared.2event.vec.start-rsquared.1event.vec.start)/(1.0-rsquared.1event.vec.start)),round.r2))

	rsquared.1event.vec.max=c(rsquared.1event.vec.max,round(rsquared.1event.max[i],round.r2),round(rsquared.1event.max[i],round.r2),round(rsquared.1event.max[i],round.r2))
	rsquared.2event.vec.max=c(rsquared.2event.vec.max,round(rsquared.2event.max[i],round.r2),round(rsquared.2event.max[i],round.r2),round(rsquared.2event.max[i],round.r2))
	rsquared.2event.vec.max.rescaled=c(rsquared.2event.vec.max.rescaled,round(r2.2event.minus.1event.rescaled.max[i],round.r2),round(r2.2event.minus.1event.rescaled.max[i],round.r2),round(r2.2event.minus.1event.rescaled.max[i],round.r2))
	fit.quality.mat=rbind(fit.quality.mat,c(fit.quality1[i],fit.quality2.vec[i]),c(fit.quality1[i],fit.quality2.vec[i]),c(fit.quality1[i],fit.quality2.vec[i]))
	fit.quality.rescaled=c(fit.quality.rescaled,fit.quality2.minus1.rescaled[i],fit.quality2.minus1.rescaled[i],fit.quality2.minus1.rescaled[i])
	grid.type.vec=c(grid.type.vec,grid.type[i],grid.type[i],grid.type[i])
	pvalue.vec.final=c(pvalue.vec.final,round(pvalue.vec[i],round.pval),round(pvalue.vec[i],round.pval),round(pvalue.vec[i],round.pval))
	twodate.pvalue.vec.final=c(twodate.pvalue.vec.final,round(twodate.pvalue.vec[i],round.pval),round(twodate.pvalue.vec[i],round.pval),round(twodate.pvalue.vec[i],round.pval))
	multiway.pvalue.vec.final=c(multiway.pvalue.vec.final,round(multiway.pvalue.vec[i],round.pval),round(multiway.pvalue.vec[i],round.pval),round(multiway.pvalue.vec[i],round.pval))
	min.fit.quality.vec.final=c(min.fit.quality.vec.final,round(min(c(fit.qualityBOTH.vec[i],fit.qualityBOTH.vecNULL[i])),round.fit),round(min(c(fit.qualityBOTH.vec[i],fit.qualityBOTH.vecNULL[i])),round.fit),round(min(c(fit.qualityBOTH.vec[i],fit.qualityBOTH.vecNULL[i])),round.fit))

               ## "best copy-vec" info:
	#best.copyvec.matrix=rbind(best.copyvec.matrix,c(paste(first.source.best.copyvec.PC1[i],"??",sep=''),paste(second.source.best.copyvec.PC1[i],"??",sep='')),c(paste(first.source.best.copyvec.PC2[i],"??",sep=''),paste(second.source.best.copyvec.PC2[i],"??",sep='')))
	best.copyvec.matrix=rbind(best.copyvec.matrix,c(strsplit(first.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1],strsplit(second.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1]),c(strsplit(first.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1],strsplit(second.source.best.copyvec.PC1[i],split=" [",fixed=TRUE)[[1]][1]),c("??","??"))
  
                        ## 1-date/2-date confidence intervals:
        x=read.table(ci.filein,skip=1,nrows=1,as.is=TRUE)
        date.est=as.character(x[5])
        date.est.ci=as.character(x[6])
        date.est.ciII=c(floor(as.double(strsplit(as.character(strsplit(as.character(x[6]),split="-")[[1]][1]),"(",fixed=TRUE)[[1]][2])),ceiling(as.double(strsplit(as.character(strsplit(as.character(x[6]),split="-")[[1]][2]),")",fixed=TRUE)[[1]][1])))
        date.est.mat=rbind(date.est.mat,c(date.est,NA),c(date.est,NA),c(date.est,NA))
        date.est.ci.mat=rbind(date.est.ci.mat,date.est.ci,date.est.ci,date.est.ci)
        date.est.ci.matII=rbind(date.est.ci.matII,c(date.est.ciII,NA,NA),c(date.est.ciII,NA,NA),c(date.est.ciII,NA,NA))
        
        count=count+3
        num.events.vec=c(num.events.vec,3,3,3)
      }
  }
warnings()
proportion.est.mat.ci[proportion.est.mat.ci<0]=0
full.info.mat[,5][full.info.mat[,5]=="EthiopianJewII"]="EthiopianJew"
date.est.mat[!is.na(date.est.mat)]=round(1950-gen.time*(round(as.double(date.est.mat[!is.na(date.est.mat)]),round.date)+1),round.date)
date.est.ci.matII[!is.na(date.est.ci.matII)]=round(1950-gen.time*(as.double(date.est.ci.matII[!is.na(date.est.ci.matII)])+1),round.date)
date.est.mat.pre=date.est.mat
date.est.ci.matII.pre=date.est.ci.matII
date.est.mat[!is.na(date.est.mat.pre) & as.double(date.est.mat.pre)<0]=paste(abs(as.double(date.est.mat.pre[!is.na(date.est.mat.pre) & as.double(date.est.mat.pre)<0])),"BCE",sep='')
date.est.mat[!is.na(date.est.mat.pre) & as.double(date.est.mat.pre)>0]=paste(abs(as.double(date.est.mat.pre[!is.na(date.est.mat.pre) & as.double(date.est.mat.pre)>0])),"CE",sep='')
date.est.ci.matII[!is.na(date.est.ci.matII.pre) & as.double(date.est.ci.matII.pre)<0]=paste(abs(as.double(date.est.ci.matII.pre[!is.na(date.est.ci.matII.pre) & as.double(date.est.ci.matII.pre)<0])),"BCE",sep='')
date.est.ci.matII[!is.na(date.est.ci.matII.pre) & as.double(date.est.ci.matII.pre)>0]=paste(abs(as.double(date.est.ci.matII.pre[!is.na(date.est.ci.matII.pre) & as.double(date.est.ci.matII.pre)>0])),"CE",sep='')
date.est.ci.matII=date.est.ci.matII[,c(2,1,4,3)]
pvalue.vec.final[pvalue.vec.final==0]="<0.01"
twodate.pvalue.vec.final[twodate.pvalue.vec.final==0]="<0.01"
twodate.pvalue.vec.final[is.na(twodate.pvalue.vec.final)]="--"
multiway.pvalue.vec.final[multiway.pvalue.vec.final==0]="<0.01"
grid.type.vec[grid.type.vec=="two-date"]="multiple-date"

                ## PRINT OUT IN GAVIN'S PREFERRED FORMAT:
pop.vecII=pop.vec
pop.vecIII=pop.vec
pop.vecII[pop.vecII=="IndianSplitII"]="Myanmar"
pop.vecIII[pop.vecIII=="IndianSplitII"]="Myanmar"
pop.vecII[pop.vecII=="EthiopianJewII"]="EthiopianJew"
for (i in 1:length(pop.vecII))
{
  pop.vecII[i]=strsplit(pop.vecII[i],split="Split")[[1]][1]
  pop.vecIII[i]=strsplit(pop.vecIII[i],split="Split")[[1]][1]
}
write("{",file=out.file,ncolumns=1)
write(paste("\t","\"file_version\": 2,",sep=''),file=out.file,ncolumns=1,append=TRUE)
write(paste("\t","\"geography\": [",sep=''),file=out.file,ncolumns=1,append=TRUE)
unique.pop.vecII=sort(unique(pop.vecII))
unique.pop.vecIII=sort(unique(pop.vecIII))
for (i in 1:length(unique.pop.vecII))
  {
    if (i<length(unique.pop.vecII)) to.print.val=paste("\t\t","{ \"name\": \"",unique.pop.vecII[i],"\", \"latitude\": ",round(geo.info$lat[as.character(geo.info$pop)==unique.pop.vecIII[i]],1),", \"longitude\": ",round(geo.info$long[as.character(geo.info$pop)==unique.pop.vecIII[i]],1)," },",sep='')
    if (i==length(unique.pop.vecII)) to.print.val=paste("\t\t","{ \"name\": \"",unique.pop.vecII[i],"\", \"latitude\": ",round(geo.info$lat[as.character(geo.info$pop)==unique.pop.vecIII[i]],1),", \"longitude\": ",round(geo.info$long[as.character(geo.info$pop)==unique.pop.vecIII[i]],1)," }",sep='')
        write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
  }
write(paste("\t","],",sep=''),file=out.file,ncolumns=1,append=TRUE)
write(paste("\t","\"admixture\": [",sep=''),file=out.file,ncolumns=1,append=TRUE)
count=1
unique.pop.vec=unique(pop.vec)
for (i in 1:length(unique.pop.vec))
{
   print(i)
               ## if there is a "split version of pop, ignore this"
   split.ind=1
   if (ignore.nonsplit==1 && sum(grep(paste(unique.pop.vec[i],"Split",sep=''),unique.pop.vec))>0) split.ind=0
   if (ignore.nonsplit==1 && split.ind==0)
     {
       to.print.finalII.i=full.info.mat[as.character(full.info.mat[,1])==unique.pop.vec[i],]
       unique.events=unique(as.character(to.print.finalII.i[,3]))
       count=count+length(unique.events)+length(grep("2Event",unique.events,value=TRUE))+length(grep("MultiEvent",unique.events,value=TRUE))
     }
   #print(c(i,unique.pop.vec[i],count,date.est.mat[count,],split.ind))
   if (split.ind==1)
     {
	       ## select "winning" analysis:
       admixture.conclusion.pop=admixture.conclusion.table[as.character(admixture.conclusion.table$pop)==unique.pop.vec[i],]
       names(admixture.conclusion.pop)=names(admixture.conclusion.table)
       analysis.choices=as.character(admixture.conclusion.pop$analysis)
                   ## (i) first limit yourself to "separate" (or "full") analyses: 
       good.analysis.ind=rep(0,length(analysis.choices))
       for (j in 1:length(analysis.choices))
         {
           #if (length(strsplit(analysis.choices[j],split="Separate")[[1]])>1) good.analysis.ind[j]=1
           if (length(strsplit(analysis.choices[j],split="Separate")[[1]])==1) good.analysis.ind[j]=1
         }
       if (sum(good.analysis.ind)>0) admixture.conclusion.pop=admixture.conclusion.pop[good.analysis.ind==1,]
                     ## (ii) then pick from among separate analyses with highest joint R^2 and fit-quality (only applicable to Slavic pops):
       sort.value=as.double(admixture.conclusion.pop$rsquared.1event.max)+as.double(admixture.conclusion.pop$fit.qualityBOTH)
       if (length(sort.value)>1)   ##  (or just make sure SlavicIV wins)
         {
           if ((strsplit(strsplit(as.character(admixture.conclusion.pop$analysis[1]),split="SeparateAnalysis")[[1]][2],split="Split")[[1]][1]=="IV") && (strsplit(strsplit(as.character(admixture.conclusion.pop$analysis[2]),split="SeparateAnalysis")[[1]][2],split="Split")[[1]][1]=="V")) sort.value=c(2,1)
           if ((strsplit(strsplit(as.character(admixture.conclusion.pop$analysis[1]),split="SeparateAnalysis")[[1]][2],split="Split")[[1]][1]=="V") && (strsplit(strsplit(as.character(admixture.conclusion.pop$analysis[2]),split="SeparateAnalysis")[[1]][2],split="Split")[[1]][1]=="IV")) sort.value=c(1,2)
         }
       analysis.type=strsplit(as.character(admixture.conclusion.pop$analysis)[sort.value==max(sort.value)][1],split="Split")[[1]]
       admixture.conclusion=as.character(admixture.conclusion.pop$admixture)[sort.value==max(sort.value)][1]
       if (analysis.type=="Final" || analysis.type=="") analysis.type="FullAnalysis"
       winning.analysis=analysis.type
       if (admixture.conclusion=="two-dates" || admixture.conclusion=="two-dates-multiple-events" || admixture.conclusion=="two-dates-uncertain" || admixture.conclusion=="two-dates-multiple-events-uncertain") winning.analysis=paste(winning.analysis,"2Event",sep='')
       #if (admixture.conclusion=="multiple-events" || admixture.conclusion=="complex?") winning.analysis=paste(winning.analysis,"MultiEvent",sep='')
       admixture.conclusion.toprint=admixture.conclusion
       if (admixture.conclusion.toprint=="no-admixture") admixture.conclusion.toprint="no admixture"
       if (admixture.conclusion.toprint=="uncertain") admixture.conclusion.toprint="uncertain"
       if (admixture.conclusion.toprint=="one-event") admixture.conclusion.toprint="one date"
       if (admixture.conclusion.toprint=="multiple-events" || admixture.conclusion.toprint=="complex?") admixture.conclusion.toprint="one date, multiway"
       if (admixture.conclusion.toprint=="two-dates" || admixture.conclusion.toprint=="two-dates-multiple-events" || admixture.conclusion.toprint=="two-dates-uncertain" || admixture.conclusion.toprint=="two-dates-multiple-events-uncertain") admixture.conclusion.toprint="multiple dates"

                ## print:
       to.print.finalII.i=full.info.mat[as.character(full.info.mat[,1])==unique.pop.vec[i],]
       write(paste("\t\t","{",sep=''),file=out.file,append=TRUE,ncolumns=1)
       pop.toprint=as.character(unique.pop.vec[i])
       if (pop.toprint=="IndianSplitII") pop.toprint="Myanmar"
       if (pop.toprint=="EthiopianJewII") pop.toprint="EthiopianJew"
       if (ignore.nonsplit==1) pop.toprint=strsplit(pop.toprint,split="Split")[[1]][1]
       write(paste("\t\t\t","\"target\": \"",pop.toprint,"\",",sep=''),file=out.file,append=TRUE,ncolumns=1)
       write(paste("\t\t\t","\"num_ind\": ",as.integer(pop.haps[match(unique.pop.vec[i],pop.vec)])/2,",",sep=''),file=out.file,append=TRUE,ncolumns=1)
       write(paste("\t\t\t","\"main_event\": \"",winning.analysis,"\",",sep=''),file=out.file,append=TRUE,ncolumns=1)
       write(paste("\t\t\t","\"main_conclusion\": \"",admixture.conclusion.toprint,"\",",sep=''),file=out.file,append=TRUE,ncolumns=1)
       write(paste("\t\t\t","\"split\": \"",unique(to.print.finalII.i[,2])[1],"\",",sep=''),file=out.file,append=TRUE,ncolumns=1)
       write(paste("\t\t\t","\"analysis\": {",sep=''),file=out.file,append=TRUE,ncolumns=1)
       files.to.make=unique(as.character(to.print.finalII.i[,3]))
       for (j in 1:length(files.to.make))
         {
           to.print.finalII.i.j=to.print.finalII.i[as.character(to.print.finalII.i[,3])==files.to.make[j],]

	   conclusion.toprint=unique(to.print.finalII.i.j[,dim(to.print.finalII.i.j)[2]])[1]
       	   if (conclusion.toprint=="no-admixture") conclusion.toprint="no admixture"
       	   if (conclusion.toprint=="uncertain") conclusion.toprint="uncertain"
       	   if (conclusion.toprint=="one-event") conclusion.toprint="one date"
       	   if (conclusion.toprint=="multiple-events" || conclusion.toprint=="complex?") conclusion.toprint="one date, multiway"
           if (conclusion.toprint=="two-dates" || conclusion.toprint=="two-dates-multiple-events" || conclusion.toprint=="two-dates-uncertain" || conclusion.toprint=="two-dates-multiple-events-uncertain") conclusion.toprint="multiple dates"
           write(paste("\t\t\t\t","\"",strsplit(files.to.make[j],split="MultiEvent")[[1]][1],"\": {",sep=''),file=out.file,append=TRUE,ncolumns=1)
           to.print.val=paste("\t\t\t\t\t","\"details\": {\"conclusion\": \"",conclusion.toprint,"\", \"grid_type\": \"",grid.type.vec[count],"\"",sep='')
           if (num.events.vec[count]==1)
	     {
		to.print.val=paste(to.print.val,", \"admixture_pvalue\": \"",pvalue.vec.final[count],"\", \"two_date_pvalue\": \"",twodate.pvalue.vec.final[count],"\", \"multiway_pvalue\": \"",multiway.pvalue.vec.final[count],"\", \"min_fit_qualityBOTH\": ",min.fit.quality.vec.final[count],", \"rsquared_max_1\": ",as.double(rsquared.1event.vec.max[count]),", \"rsquared_max_2\": ",as.double(rsquared.2event.vec.max[count]),", \"rsquared_2eventrescaled_max\": ",rsquared.2event.vec.max.rescaled[count],", \"av_rel_mse\": ",round(as.double(avrelmse[count,1]),round.avrel),", \"fit_quality_1\": ",round(as.double(fit.quality.mat[count,1]),round.fit),", \"fit_quality_2\": ",round(as.double(fit.quality.mat[count,2]),round.fit),", \"fit_quality_2eventrescaled\": ",round(fit.quality.rescaled[count],round.fit),", \"date_est\": \"",date.est.mat[count,1],"\", \"date_est_lower\": \"",date.est.ci.matII[count,1],"\", \"date_est_upper\": \"",date.est.ci.matII[count,2],"\", \"proportion_est\": ",as.double(proportion.est.mat[count,1]),", \"
proportion_est_lower\": ",as.double(proportion.est.mat.ci[count,1]),", \"proportion_est_upper\": ",as.double(proportion.est.mat.ci[count,2]),"},",sep='')
		count.add=1
	     }
           if (num.events.vec[count]==2) 
             {
               to.print.val=paste(to.print.val,", \"admixture_pvalue\": \"",pvalue.vec.final[count],"\", \"two_date_pvalue\": \"",twodate.pvalue.vec.final[count],"\", \"multiway_pvalue\": \"",multiway.pvalue.vec.final[count],"\", \"min_fit_qualityBOTH\": ",min.fit.quality.vec.final[count],", \"rsquared_max_1\": ",as.double(rsquared.1event.vec.max[count]),", \"rsquared_max_2\": ",as.double(rsquared.2event.vec.max[count]),", \"rsquared_2eventrescaled_max\": ",rsquared.2event.vec.max.rescaled[count],", \"fit_quality_1\": ",round(as.double(fit.quality.mat[count,1]),round.fit),", \"fit_quality_2\": ",round(as.double(fit.quality.mat[count,2]),round.fit),", \"fit_quality_2eventrescaled\": ",round(fit.quality.rescaled[count],round.fit),sep='')
               for (k in 1:num.events.vec[count]) to.print.val=paste(to.print.val,", \"date_est_",k,"\": \"",date.est.mat[count,k],"\", \"date_est_lower_",k,"\": \"",date.est.ci.matII[count,(2*(k-1)+1)],"\", \"date_est_upper_",k,"\": \"",date.est.ci.matII[count,(2*(k-1)+2)],"\", \"proportion_est_",k,"\": ",as.double(proportion.est.mat[count,k]),", \"proportion_est_lower_",k,"\": ",as.double(proportion.est.mat.ci[count,(2*(k-1)+1)]),", \"proportion_est_upper_",k,"\": ",as.double(proportion.est.mat.ci[count,(2*(k-1)+2)]),", \"av_rel_mse_",k,"\": ",round(as.double(avrelmse[count,k]),round.avrel),sep='')
               to.print.val=paste(to.print.val,"},",sep='')
	       count.add=3
             }
           if (num.events.vec[count]==3) 
             {	
               to.print.val=paste(to.print.val,", \"admixture_pvalue\": \"",pvalue.vec.final[count],"\", \"two_date_pvalue\": \"",twodate.pvalue.vec.final[count],"\", \"multiway_pvalue\": \"",multiway.pvalue.vec.final[count],"\", \"min_fit_qualityBOTH\": ",min.fit.quality.vec.final[count],", \"rsquared_max_1\": ",as.double(rsquared.1event.vec.max[count]),", \"rsquared_max_2\": ",as.double(rsquared.2event.vec.max[count]),", \"rsquared_2eventrescaled_max\": ",rsquared.2event.vec.max.rescaled[count],", \"av_rel_mse\": ",round(as.double(avrelmse[count,1]),round.fit),", \"fit_quality_1\": ",round(as.double(fit.quality.mat[count,1]),round.fit),", \"fit_quality_2\": ",round(as.double(fit.quality.mat[count,2]),round.fit),", \"fit_quality_2eventrescaled\": ",round(fit.quality.rescaled[count],round.fit),", \"date_est\": \"",date.est.mat[count,1],"\", \"date_est_lower\": \"",date.est.ci.matII[count,1],"\", \"date_est_upper\": \"",date.est.ci.matII[count,2],"\"",sep='')
               for (k in 1:2) to.print.val=paste(to.print.val,", \"proportion_est_",k,"\": ",as.double(proportion.est.mat[count,k]),", \"proportion_est_lower_",k,"\": ",as.double(proportion.est.mat.ci[count,(2*(k-1)+1)]),", \"proportion_est_upper_",k,"\": ",as.double(proportion.est.mat.ci[count,(2*(k-1)+2)]),", \"av_rel_mse_",k,"\": ",round(as.double(avrelmse[count,k]),round.avrel),sep='')
               to.print.val=paste(to.print.val,"},",sep='')
	       count.add=3
             }

           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   to.print.val=paste("\t\t\t\t\t","\"events\":[",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   to.print.val=paste("\t\t\t\t\t\t","{ \"name\": \"EventOne\",",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           to.print.val=paste("\t\t\t\t\t\t\t","\"best_matches\": { \"side1\": \"",best.copyvec.matrix[count,1],"\", \"side2\": \"",best.copyvec.matrix[count,2],"\" },",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   to.print.val=paste("\t\t\t\t\t\t\t","\"views\":{",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   to.print.val=paste("\t\t\t\t\t\t\t\t","\"sourceView\":{",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   to.print.val=paste("\t\t\t\t\t\t\t\t\t","\"type\": \"source\",",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   to.print.val=paste("\t\t\t\t\t\t\t\t\t","\"sides\":{",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           to.print.val=paste("\t\t\t\t\t\t\t\t\t\t","\"side1\":[",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           to.print.finalII.i.j.side1=matrix(to.print.finalII.i.j[as.integer(to.print.finalII.i.j[,4])==1,],ncol=dim(to.print.finalII.i.j)[2])
           for (k in 1:dim(to.print.finalII.i.j.side1)[1])
             {
               to.print.val=paste("\t\t\t\t\t\t\t\t\t\t\t","{ \"source\": \"",as.character(to.print.finalII.i.j.side1[k,5]),"\", \"proportion\": ",round(as.double(to.print.finalII.i.j.side1[k,6]),round.prop),", \"curves\": \"",as.character(to.print.finalII.i.j.side1[k,7]),"\"",sep='')
               if (k < dim(to.print.finalII.i.j.side1)[1]) to.print.val=paste(to.print.val," },",sep='')  
               if (k == dim(to.print.finalII.i.j.side1)[1]) to.print.val=paste(to.print.val," }",sep='')  
               write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
             }
           write(paste("\t\t\t\t\t\t\t\t\t\t","],",sep=''),file=out.file,append=TRUE,ncolumns=1)
           to.print.val=paste("\t\t\t\t\t\t\t\t\t\t","\"side2\":[",sep='')
           write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           to.print.finalII.i.j.side2=matrix(to.print.finalII.i.j[as.integer(to.print.finalII.i.j[,4])==2,],ncol=dim(to.print.finalII.i.j)[2])
           for (k in 1:dim(to.print.finalII.i.j.side2)[1])
             {
               to.print.val=paste("\t\t\t\t\t\t\t\t\t\t\t","{ \"source\": \"",as.character(to.print.finalII.i.j.side2[k,5]),"\", \"proportion\": ",round(as.double(to.print.finalII.i.j.side2[k,6]),round.prop),", \"curves\": \"",as.character(to.print.finalII.i.j.side2[k,7]),"\"",sep='')
               if (k < dim(to.print.finalII.i.j.side2)[1]) to.print.val=paste(to.print.val," },",sep='')  
               if (k == dim(to.print.finalII.i.j.side2)[1]) to.print.val=paste(to.print.val," }",sep='')  
               write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
             }
           write(paste("\t\t\t\t\t\t\t\t\t\t","]",sep=''),file=out.file,append=TRUE,ncolumns=1)
           write(paste("\t\t\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
           if (num.events.vec[count]==1)
	   {
		write(paste("\t\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
	   }
           if (num.events.vec[count]==2 || num.events.vec[count]==3)
	   {
		write(paste("\t\t\t\t\t\t\t\t","},",sep=''),file=out.file,append=TRUE,ncolumns=1)
	   	#to.print.val=paste("\t\t\t\t\t\t\t\t","{ \"name\": \"EventOne\",",sep='')
	   	to.print.val=paste("\t\t\t\t\t\t\t\t","\"contrastView\":{",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   	to.print.val=paste("\t\t\t\t\t\t\t\t\t","\"type\": \"contrast\",",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	#to.print.val=paste("\t\t\t\t\t\t\t","\"best_matches\": { \"side1\": \"",best.copyvec.matrix[count+1,1],"\", \"side2\": \"",best.copyvec.matrix[count+1,2],"\" },",sep='')
           	#write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   	to.print.val=paste("\t\t\t\t\t\t\t\t\t","\"sides\":{",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	to.print.val=paste("\t\t\t\t\t\t\t\t\t\t","\"side1\":[",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	to.print.finalII.i.j.side1=matrix(to.print.finalII.i.j[as.integer(to.print.finalII.i.j[,4])==3,],ncol=dim(to.print.finalII.i.j)[2])
           	for (k in 1:dim(to.print.finalII.i.j.side1)[1])
             	{
			to.print.val=paste("\t\t\t\t\t\t\t\t\t\t\t","{ \"source\": \"",as.character(to.print.finalII.i.j.side1[k,5]),"\", \"proportion\": ",round(as.double(to.print.finalII.i.j.side1[k,6]),round.prop),", \"curves\": \"",as.character(to.print.finalII.i.j.side1[k,7]),"\"",sep='')
               		if (k < dim(to.print.finalII.i.j.side1)[1]) to.print.val=paste(to.print.val," },",sep='')  
               		if (k == dim(to.print.finalII.i.j.side1)[1]) to.print.val=paste(to.print.val," }",sep='')  
               		write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
             	}
           	write(paste("\t\t\t\t\t\t\t\t\t\t","],",sep=''),file=out.file,append=TRUE,ncolumns=1)
           	to.print.val=paste("\t\t\t\t\t\t\t\t\t\t","\"side2\":[",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	to.print.finalII.i.j.side2=matrix(to.print.finalII.i.j[as.integer(to.print.finalII.i.j[,4])==4,],ncol=dim(to.print.finalII.i.j)[2])
           	for (k in 1:dim(to.print.finalII.i.j.side2)[1])
             	{
			to.print.val=paste("\t\t\t\t\t\t\t\t\t\t\t","{ \"source\": \"",as.character(to.print.finalII.i.j.side2[k,5]),"\", \"proportion\": ",round(as.double(to.print.finalII.i.j.side2[k,6]),round.prop),", \"curves\": \"",as.character(to.print.finalII.i.j.side2[k,7]),"\"",sep='')
               		if (k < dim(to.print.finalII.i.j.side2)[1]) to.print.val=paste(to.print.val," },",sep='')  
               		if (k == dim(to.print.finalII.i.j.side2)[1]) to.print.val=paste(to.print.val," }",sep='')  
               		write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
             	}
           	write(paste("\t\t\t\t\t\t\t\t\t\t","]",sep=''),file=out.file,append=TRUE,ncolumns=1)
           	write(paste("\t\t\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t","},",sep=''),file=out.file,append=TRUE,ncolumns=1)
	   	to.print.val=paste("\t\t\t\t\t\t","{ \"name\": \"EventTwo\",",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	to.print.val=paste("\t\t\t\t\t\t\t","\"best_matches\": { \"side1\": \"",best.copyvec.matrix[count+2,1],"\", \"side2\": \"",best.copyvec.matrix[count+2,2],"\" },",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   	to.print.val=paste("\t\t\t\t\t\t\t","\"views\":{",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   	to.print.val=paste("\t\t\t\t\t\t\t\t","\"contrastView\":{",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   	to.print.val=paste("\t\t\t\t\t\t\t\t\t","\"type\": \"contrast\",",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
	   	to.print.val=paste("\t\t\t\t\t\t\t\t\t","\"sides\":{",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	to.print.val=paste("\t\t\t\t\t\t\t\t\t\t","\"side1\":[",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	to.print.finalII.i.j.side1=matrix(to.print.finalII.i.j[as.integer(to.print.finalII.i.j[,4])==5,],ncol=dim(to.print.finalII.i.j)[2])
           	for (k in 1:dim(to.print.finalII.i.j.side1)[1])
             	{
			to.print.val=paste("\t\t\t\t\t\t\t\t\t\t\t","{ \"source\": \"",as.character(to.print.finalII.i.j.side1[k,5]),"\", \"proportion\": ",round(as.double(to.print.finalII.i.j.side1[k,6]),round.prop),", \"curves\": \"",as.character(to.print.finalII.i.j.side1[k,7]),"\"",sep='')
               		if (k < dim(to.print.finalII.i.j.side1)[1]) to.print.val=paste(to.print.val," },",sep='')  
               		if (k == dim(to.print.finalII.i.j.side1)[1]) to.print.val=paste(to.print.val," }",sep='')  
               		write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
             	}
           	write(paste("\t\t\t\t\t\t\t\t\t\t","],",sep=''),file=out.file,append=TRUE,ncolumns=1)
           	to.print.val=paste("\t\t\t\t\t\t\t\t\t\t","\"side2\":[",sep='')
           	write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
           	to.print.finalII.i.j.side2=matrix(to.print.finalII.i.j[as.integer(to.print.finalII.i.j[,4])==6,],ncol=dim(to.print.finalII.i.j)[2])
           	for (k in 1:dim(to.print.finalII.i.j.side2)[1])
             	{
			to.print.val=paste("\t\t\t\t\t\t\t\t\t\t\t","{ \"source\": \"",as.character(to.print.finalII.i.j.side2[k,5]),"\", \"proportion\": ",round(as.double(to.print.finalII.i.j.side2[k,6]),round.prop),", \"curves\": \"",as.character(to.print.finalII.i.j.side2[k,7]),"\"",sep='')
               		if (k < dim(to.print.finalII.i.j.side2)[1]) to.print.val=paste(to.print.val," },",sep='')  
               		if (k == dim(to.print.finalII.i.j.side2)[1]) to.print.val=paste(to.print.val," }",sep='')  
               		write(to.print.val,file=out.file,append=TRUE,ncolumns=1)
             	}
           	write(paste("\t\t\t\t\t\t\t\t\t\t","]",sep=''),file=out.file,append=TRUE,ncolumns=1)
           	write(paste("\t\t\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
		write(paste("\t\t\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
	   }
	   write(paste("\t\t\t\t\t","]",sep=''),file=out.file,append=TRUE,ncolumns=1)
           if (j<length(files.to.make)) write(paste("\t\t\t\t","},",sep=''),file=out.file,append=TRUE,ncolumns=1)
           if (j==length(files.to.make)) write(paste("\t\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
           count=count+count.add
         }
       write(paste("\t\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)

       if (i < length(unique.pop.vec)) write(paste("\t\t","},",sep=''),file=out.file,append=TRUE,ncolumns=1)
       if (i == length(unique.pop.vec)) write(paste("\t\t","}",sep=''),file=out.file,append=TRUE,ncolumns=1)
     }
 }
write(paste("\t","]",sep=''),file=out.file,append=TRUE,ncolumns=1)
write("}",file=out.file,append=TRUE,ncolumns=1)
warnings()

print(quantile(fit.quality.1event.allSIMS,c(0.01,0.025,0.05)))
print(quantile(rsquared.2event.rescaled.allSIMS,c(0.95,0.975,0.99)))

q(save='no')

