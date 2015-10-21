###################################################################
###################################################################
##################
## SCRIPT TO GENERATE MULTI-VARIATE NORM TEST FOR INDIVIDUAL SNPS
## COMPARED TO GENOME(CHROMOSOME)-WIDE PAINTINGS
##################
###################################################################
###################################################################
## NB THIS IS AN UPDATE OF CHRIS'S Copy_Deviant_v5.R

## EXAMPLE R  < /well/malariagen/malariagen/human/george/copy_selection/MVN_CopySelection.R FULAI --no-save > /well/malariagen/malariagen/human/george/copy_selection/Copy_Selection_v5_FULAI.out 

options(stringsAsFactors=FALSE);
temp <- commandArgs()
pop <- temp[2]

## NB: TO DO - EDIT SO THAT INPUT IS THE SAME AS LRT_CopySelection.R ##
#in_file <- temp[4]
#id_file <- temp[5]
#out_file<- temp[6]

#Get the population for analysis

in_dir <- "/well/malariagen/malariagen/human/george/copy_selection/"
samp_dir <- "/data/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/"
snp_dir <- "/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/"
pop_file <- paste(indir,"populationOverviewCopyProbs.txt",sep="")
popkey <- read.table(pop_file,header=T)


out_file <- paste(in_dir,"Copy_Selection_deviant_v5_",pop,".csv",sep="")
header=-1

if(!file.exists(out_file))
{
    tmp <- read.table(paste(samp_dir,pop,"nolocalAllChromsPP.samples.out.gz",sep=""),header=FALSE,as.is=TRUE,nrow=header);
    tmp <- as.matrix(tmp);

#Just one of the samples
    subset <- seq(1,ncol(tmp),by=1);
    tmp <- t(as.matrix(tmp[,subset]));
## x11(type="Xlib");
## rho <- cor(haps)
## image(t(cor(haps)));

#Read the population and individual information
    id_file <- paste(in_dir,pop,"nolocal.idfile.txt",sep="");
    ids <- read.table(id_file)
    ids <- ids[rep(1:nrow(ids),each=2),];
    own_reg <- as.character(popkey$Region[popkey$Ethnic_Group==pop])
    snp_file <- paste(snp_dir,"AllPops330Kphased.legend.gz",sep="");
    snps <- read.table(snp_file,header=F);
    colnames(snps) <- c("chr","rsid","pos","a0","a1")
    haps <- array(NA,dim(tmp));

##Reassign copying into broad regions
    regions <- unique(as.character(popkey[,2]));
    for(i in 1:length(regions))
    {
        haps[matrix(tmp %in% which(ids[,2] %in% popkey[which(popkey[,2] == regions[i]),1]),
                    ncol=ncol(haps),byrow=FALSE)] = i;
        print(regions[i]);
    }

    #table(regions[c(haps)])
    groups <- regions;

# #Reassign copying into groups
# groups <- unique(as.character(popkey[,1]));
# for(i in 1:length(groups)){
#      haps[matrix(tmp %in% which(ids[,2] %in% groups[i]),
#                ncol=ncol(haps),byrow=FALSE)] = i;
#      print(groups[i]);
# }
# table(groups[c(haps)]);	

#Set things up
    nsamp = 10;
    nind = nrow(haps) / nsamp;
    nsnps = ncol(haps);
    first = seq(1,nrow(haps),by=nsamp);
#Get individuals averages
    avg <- array(NA, c(nrow(haps),length(groups)));
    for(i in 1:length(groups)) avg[,i] = rowSums(haps == i);
    avg = avg / ncol(haps);
# #Get the counts per site
# counts = array(NA,c(nsamp,nsnps,length(groups)));
# for(j in 1:nsamp){
#     index = first + (j-1);
#     for(i in 1:length(groups))
#         counts[j,,i] = colSums(haps[index,] == i);
#     print(j);
# }
# props = counts / nind;

#Residual
    res <- haps;
    for(i in 1:nrow(res))
    {
        for(j in 1:length(groups))
        {
            res[i,haps[i,] == j] = (1 - avg[i,j])
        }
    print(i);
    }

#Get the total residual deviation
    deviant = array(0,c(nsamp,nsnps,length(groups)));
    for(j in 1:nsamp)
    {
      index = first + (j-1);
      for(i in 1:length(groups))
      {
          tmp = res[index,];
          tmp[haps[index,] != i] = 0;
          deviant[j,,i] = colSums(tmp);
      }
      print(j);
    }

    #x11(type="Xlib",width=15);
    #barplot(t(counts[1,,]),beside=FALSE,border=NA,space=0,col=2:7);
    random <- sample(1:nsnps,500);
    tmp <- deviant[,,which(colSums(avg) != 0)];
    x <- array(0,dim(tmp)[-1])
    for(i in 1:nsamp) x <- x + tmp[i,,] 
    x = x / nsamp;
    mu <- colSums(x[random,]) / nrow(x[random,]);
    sigma <- cov(x[random,]);
# anc = 6;
# plot(-log10(pnorm(x[,anc],mu[anc],sd=sqrt(sigma[anc,anc]))))
    prior.sigma = array(0,dim(sigma))
    diag(prior.sigma) = 0.0001;
    sigma <- sigma + prior.sigma;
    output <- array(NA,c(nsnps,2));
    for(i in 1:nsnps)
          output[i,1] <- t(x[i,] - mu) %*% solve(sigma) %*% (x[i,] - mu);
    output[,2] <- pchisq(output[,1],ncol(sigma),lower=FALSE)

    marg.pval = array(NA,dim(x));
    for(i in 1:ncol(marg.pval))
          marg.pval[,i] = pchisq((x[,i]-mu[i])^2/diag(sigma)[i],1,lower=FALSE);

#if(DISPLAY) plot(-log10(output[,2]))
    colnames(output) = c("chisq","pval");
    colnames(x) <- regions[regions!=own_reg]
    colnames(marg.pval) <- paste(regions[regions!=own_reg],"_pval",sep="")
    res <- data.frame(snps[1:nsnps,],output,x,marg.pval);
    #res <- res[order(res[,7]),]

#if(DISPLAY) qqplot(rchisq(nrow(output),ncol(sigma)),output[,1])
#if(DISPLAY) abline(0,1)
    options(digits = 5);
    write.csv(res,outfile);
}

