###################################################################
###################################################################
##################
##################
###################################################################
###################################################################

## ALL LOCAL PAINTINGS ARE STORED IN AN hdf5 FILE
hdf5file <- '/mnt/kwiat/well/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
## ALL LENGTH INFO IS STRORED IN ANOTHER hdf5 FILE

library("rhdf5")
library("bigmemory")


## SWITCHES V RECOMBINATION RATE


## for each snp,
##   for each haplotype
##     a) see who they're copying
##     b) see who copier is copying
##     c) get lengths

## FUNCTIONS
hap2sampleindex <- function(hap,nsamps=10){
    ## finds the first sample index for a haplotype
    sample <- (hap*nsamps)-(nsamps-1)
    return(sample)
}


## GET SAMPLE INFORMATION
inds <- data.frame(t(h5read(hdf5file,"/paintings/samples/individuals")))
colnames(inds) <- c("ind","pop","X")
regions <-data.frame(t(h5read(hdf5file,"/paintings/samples/regions")))
colnames(regions) <- c("pop","country","region")
H5close()

## GENERATE SOME INDIVIDUAL, POPULATION, AND REGION INDEXES
happops <- c()
for(i in 1:nrow(inds)) happops <- c(happops,rep(as.character(inds$pop[i]),2))
hapregs <- c()
for(i in happops) hapregs <- c(hapregs,as.character(regions$region[regions$pop==i]))
hapregs <- as.factor(hapregs)
## GET SNP INFO
snps <- data.frame(t(h5read(hdf5file,paste("/paintings/chrom",chrom,"/snps",sep=""))))
colnames(snps) <- c("chrom","rsid","pos","a0","a1")
map <- data.frame(t(h5read(hdf5file,paste("/paintings/chrom",chrom,"/map",sep=""))))
colnames(map) <- c("pos","cm.Mb")

## DEFINE WHERE TO FIND THE PAINTINGS + INFO IN THE HDF5 FILE
nonlocalpnts <- paste("/paintings/chrom",chrom,"/nonlocal",sep="")
localpnts <- paste("/paintings/chrom",chrom,"/local",sep="")

n_samps <- 10 
n_haps <- length(happops)
n_snps <- nrow(snps)
n_paintedhaps <- length(happops) * n_samps
hap <- 1


paintings <- h5read(hdf5file,nonlocalpnts,index=list(hap,1:n_snps))
regpaintings <- hapregs[paintings]
paintingsindex <- hap2sampleindex(paintings)

paintingsindexlist <- list(paintingsindex,1:n_snps)
whocopies <- diag(h5read(hdf5file,localpnts,index=paintingsindexlist))
H5close()


## the haplotype that people copy
test <- lines4[s,1:10]
whocopied <- donor_hap_vec[test]
## now i need to find the index for each haplotype in test
## let's look at the first sample for the time being ...
test <- (test*10)-9
whocopiedcopies1 <- h5read(hdf5file,paste("/paintings/chrom",chrom,"/nonlocal",sep=""),
                           index=list(test,s))
H5close()

whocopied ## this is the region where the actual samples come from
donor_hap_vec[whocopiedcopies] ## this is the region that the copied samples copy from 




###################################################################
###################################################################
## SOME IMPORTANT VARIABLES
options(scipen=999,digits=20)
## INFORMATION ON SNPS - CHROMOSOME,POSITION,ALLELES ETC
## DIRECTORY WITH FILE WITH SNP INFO IN THEM
snp_dir <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
snp_file_pre <- "AllPops330KChrom"
snp_file_pos <- "phased.legend.gz"
## NUMBER OF SAMPLES TO USE TO GENERATE LIKELIHOODS
n_samps <- 10 # IF 1 THEN JUST USES THE FIRST
useallsamps <- TRUE ## USE ALL SAMPLES?
samp2use <- 1 ## IF useallsamps == FALSE, THEN USE THIS SAMPLE
###################################################################
###################################################################
###################################################################
## INFORMATION ON THE SAMPLES
## POPFILE MUST HAVE TWO COLUMNS: "Ethnic_Group" AND "Region"
## THIS FILE WILL BE READ TO GROUP POPULATIONS INTO REGIONS
pop_file <- "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/analysislists/populationOverviewCopyProbs.txt"
popkey <- read.table(pop_file,header=T)
###################################################################
###################################################################
###################################################################
### FUNCTIONS ###
## THE MAIN LIKELIHOOD FUNTION:
## THIS TAKES:
##      v: COPYING PROBS PER HAPLOTYPE; THESE ARE ESTIMATED FROM GENOME-WIDE SAMPLES
##      data: A VECTOR OF HAPLOTYPES AT A GIVEN SNP
##      nsamps: THE NUMBER OF SAMPLES OF EACH HAPLOTYPE AT THIS SNP
##      useall: SHOLD ALL PAINTED HAPLOTPYES BE USED TO GENERATE LIKELIHOOD?
##      slike: IF useall = FALSE ; THEN USE THIS SAMPLE TO GENERATE LIKELIHOOD
##      NOTE: THE COLUMNS OF V SHOULD CORRESPOND TO THE COPYING PROBS
##          OF COPYING FROM THE REGION WITH THE SAME INDEX
##          THAT IS, COLUMN 1 OF v SHOULD DESCRIBE THE PROB OF COPYING
##          FROM THE REGION IDENTIFIED BY 1 IN data
##      NOTE: THAT DATA SHOULD BE n_haps * n_samps long

loglik <- function(v,data,nsamps=n_samps,useall=useallsamps,slike=samp2use)
{
    ## GENERATE AN IDENTITY MATRIX WHERE 1ST COL IS 1:nhaps
    ## AND NEXT nsamps COLS ARE THE IDENTITY OF THE DONOR THAT
    ## EACH HAPLOTYPE COPIES AT THIS SNP
    nhaps <- length(data)/nsamps
    ind <- matrix(cbind(1:nhaps,matrix(data,ncol=nsamps,byrow=T)),ncol=nsamps+1)
    ## MAKE AN EMPTY MATRIX AND FILL WITH THE PROBABILITIES
    ## OF COPYING FROM THAT REGION AT THAT SNP ie BY USING v
    prob <- matrix(NA, nrow=nhaps,ncol=nsamps)
    for(k in 2:(nsamps+1)) prob[,(k-1)] <- v[ind[,c(1,k)]]
    
    ## NOW GENERATE LOG LIKELIHOOD ACROSS ALL HAPLOTYPES
    ## CURRENT IMPLEMENTATION USES A SINGLE PAINTING: slike DEFINES THE SAMPLE TO USE
    if(useall == FALSE ) prob <- sum(log(prob[,slike]))
    #if(useall == TRUE ) prob <- sum(log(rowSums(prob)/ncol(prob)))
    ## THIS IS GBJB'S PREVIOUS VERSION, WHICH LOOKS LIKE IT WORKS
    ## BUT IS PERHAPS INCORRECT THING TO DO
    if(useall == TRUE) prob <- mean(apply(prob,2,function(x)sum(log(x))))
    return(prob)
}

## THIS FUNCTION OPTIMISES THE VALUE OF LAMBDA, BASED ON MLE USING OPTIM
##  RECALL THAT v MUST HAVE n_regions COLUMNS, WITH EACH COLUMN DESCRIBING
##  THE COPYING PROBS FROM A GIVEN REGION INDEX. THESE ARE THE INDICES IN data
##      lambda: IS THE VALUE WE WANT TO OPTIMISE (THE PROPORIONAL INCREASE IN 
##          COPYING FROM A GIVEN DONOR REGION)
##      colindex: THE COLUMN OF v THAT WE WANT TO ADJUST TO GET OPTIMUM VALUE
##          OF LAMBDA - IE THE REGION OF INTEREST

par.loglik <- function(v,data,nsamps=n_samps,useall=useallsamps,slike=samp2use,lambda,colindex)
{
    ## WE NEED TO CHOOSE A REGION TO PIVOT THE LAMBDA VALUE AROUND
    ## WE WILL USE COLUMN 1, UNLESS THIS IS THE REGION THAT THE POP
    ## OF INTEREST COMES FROM, IN WHICH CASE THE COPYING PROPS WILL
    ## BE 0 (AS WE DIS-ALLOW) SELF-COPYING. IN THIS CASE WE USE COL 2
    ## AS THE "base_col"
    base_col <- 1
    if(sum(v[,base_col]) == 0) base_col <- 2
    ## INTIATE A MATRIX TO STORE OUR ADJUSTED v VALUES
    x <- adjv <- matrix(0,ncol=ncol(v),nrow(v))
    ## x IS THE PROBS DIVIDED BY THE BASE COLUMN, AND LOGGED
    x <- log(v/v[,base_col])
    
    ## NOW ADJUST COPYING PROBS BY LAMBDA
    x[,colindex] <- x[,colindex] + lambda
    ## GET OUT OF LOG SPACE
    adjv <- apply(x,2,function(i) exp(i))
    ## RENORMALISE
    adjv <- adjv / rowSums(adjv)
    ## ESTIMATE THE LOG-LIKELIHOOD
    llik <- -(loglik(v=adjv,data=data,nsamps=nsamps,useall=useall,slike=slike) + dnorm(lambda,0,10,log=TRUE))
    return(llik)
}
###################################################################
###################################################################
###################################################################
### PROGRAM ###
## 00 LOAD PAINTING SAMPLES
## LOAD SAMPLES: NOTE THESE ARE GENERATED USING THE PY-PROG
## AND SHOULD BE SNPS AS COLUMNS AND HAPS AS COLUMNS
#lines3 <- read.big.matrix(in_file,type="char",sep=" ")
lines3 <- as.big.matrix(read.table(in_file,colClasses="integer"))

## 01 GET POPULATION INFO PLUS INFO ON REGIONS, NUMBERS OF HAPS ETC.
ids <- read.table(id_file)    
pops <- as.character(ids[,2])
regions <- sapply(pops,function(x){as.character(popkey$Region[popkey$Ethnic_Group==x])})
ids <- cbind(ids,regions)
n_regs <- length(unique(regions))
region_ids <- levels(popkey$Region)
n_haps <- sum(ids$V2==pop)*2
## FOR EACH DONOR HAPLOTYPE WE WANT TO KNOW THE IDENTITY OF THE RGION THAT IS
## COMES FROM. THIS ALIGNS THE DONORS IN THE *samples.out FILES WITH THEIR
## REGIONAL IDENTITY
donor_hap_vec <- c()
for(i in 1:nrow(ids)) donor_hap_vec <- c(donor_hap_vec,ids$regions[i],ids$regions[i])
pop_hap_vec <- c()
for(i in 1:nrow(ids)) pop_hap_vec <- c(pop_hap_vec,as.character(ids$V2[i]),as.character(ids$V2[i]))
id_hap_vec <- c()
for(i in 1:nrow(ids)) id_hap_vec <- c(id_hap_vec,as.character(ids$V1[i]),as.character(ids$V1[i]))

## 02 LOAD SNP INFO
snps <- c()
for(chrom in 1:22)
{
    if(chrom < 10) chrom  <- paste0("0",chrom)
    snp_file <- paste(snp_dir,snp_file_pre,chrom,snp_file_pos,sep="")
    snp <- read.table(snp_file,header=F)
    snps <- rbind(snps,snp)
}    
colnames(snps) <- c("chrom","rsid","pos","a0","a1")

## 03 LEAVE-ONE-OUT COPYING PROBS
## COMPUTE GENOME-WIDE COPYING-PROBS BASED ON NUMBER OF 
## SNPS COPIED FROM EACH REGION ACROSS ALL CHROMSOMES EXCEPT
## THE ONE THAT WE'RE ANALYSING
## WE AVERAGE THESE ACROSS ALL SAMPLES ????

chrom <- mainchrom
rows <- snps$chrom!=chrom
#lines4 <- lines3[rows,]
print(paste0("estimating likelihoods for chromosome: ", chrom, " in ", pop))
ind_copy_probs <- matrix(0,nrow=n_haps*n_samps,ncol=n_regs)
colnames(ind_copy_probs) <- 1:n_regs
ind_vec <- seq(1,ncol(lines3),n_samps)
for(i in 1:n_haps)
{
    print(paste("esimating probs for hap:", i))
    cp <- matrix(0,nrow=n_samps,ncol=n_regs)
    colnames(cp) <- 1:n_regs
    for(j in 1:n_samps)
    {
        tmp_cp <- table(donor_hap_vec[lines3[rows,(ind_vec[i]:(ind_vec[i]+9))[j]]])
        tmp_cp <- tmp_cp/sum(tmp_cp)
        cp[j,names(tmp_cp)] <- tmp_cp
    }
    ind_copy_probs[((i-1)*10+1):((i-1)*10+10),] <- cp
}

## THESE ARE THE MAIN COPYING PROPORTIONS
ind_copy_probs <- ind_copy_probs/rowSums(ind_copy_probs)

## REMOVE REGION THAT IS ARE NOT COPIED FROM A REGION ID VECTOR
region_ids2 <- region_ids[colSums(ind_copy_probs)!=0]
n_regs2 <- length(region_ids2)
self_reg <- region_ids[!region_ids%in%region_ids2]

## 04 GET SNPS ON CURRENT CHROMOSOME AND ESTIMATE NULL LIKELIHOODS
snps2 <- snps[snps$chrom==chrom,]
n_snps <- nrow(snps2)
lines4 <- as.big.matrix(lines3[snps$chrom==chrom,])

## DEFINE WHETHER WE AVERAGE COPY PROBS ACROSS ALL SAMPLES OR JUST USE ONE SAMPLE
v <- c()
if(useallsamps == TRUE) for(i in seq(1,nrow(ind_copy_probs),by=10))   v <- rbind(v,apply(ind_copy_probs[i:(i+9),],2,mean))
if(useallsamps == FALSE) v <- ind_copy_probs[seq(samp2use,nrow(ind_copy_probs),by=10),]

## NULL SNP LIKELIHOODS
snp_liks <- apply(lines4,1,function(x){loglik(v,donor_hap_vec[x],nsamps=n_samps,useall=useallsamps,slike=samp2use)})

## WE CAN CHECK IF A CHUNK IS 'GENUINE' BY LOOKING AT A DONOR AT A SNP
## AND SEEING WHO THEY COPY AT THE SAME SNP
## NEED TO CHECK DIMENSIONS OF THE FINAL FILE -- THEY MAY BE SWAPPED ROUND ...

## let's look at th first snp on chrom 22
s <- 1000
## the haplotype that people copy
test <- lines4[s,1:10]
whocopied <- donor_hap_vec[test]
## now i need to find the index for each haplotype in test
## let's look at the first sample for the time being ...
test <- (test*10)-9
whocopiedcopies1 <- h5read(hdf5file,paste("/paintings/chrom",chrom,"/nonlocal",sep=""),
                          index=list(test,s))
H5close()

whocopied ## this is the region where the actual samples come from
donor_hap_vec[whocopiedcopies] ## this is the region that the copied samples copy from 


## NOW FIND MLE OF LAMBDA
mle <- matrix(0,nrow=n_snps,ncol=2*n_regs2)
cnames <- c()
for(i in region_ids2) cnames <- c(cnames,i,i)
colnames(mle) <- paste(cnames,rep(c("likelihood","lambda"),n_regs2),sep=".")
for(i in 1:n_snps)
{
    pcdone <- signif((i/n_snps)*100,2)
    if(pcdone%%10 == 10) print(paste(pcdone," % through snps"))
    for(reg_index in 1:length(region_ids))
    {
        reg_id <- region_ids[reg_index]
    	if(reg_id != self_reg)
    	{
            	lambda <- 0
            	opt1 <- optim(lambda,par.loglik,
            	              data=donor_hap_vec[lines4[i,]],
            	              nsamps=n_samps,useall=useallsamps,
            	              slike=samp2use,colindex=reg_index,v=v,
            	              method="Nelder-Mead")
                mle[i,grep(reg_id,colnames(mle))]  <- c(-opt1$value,opt1$par)
    	}
    }
}

mle <- cbind(snps2,snp_liks,mle)
colnames(mle)[colnames(mle)=="snp_liks"] <- "null.likelihood"
write.table(mle,file=out_file,quote=F,col.names=T,row.names=F)
