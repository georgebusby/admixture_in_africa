



inds <- read.table("/mnt/kwiat/data/bayes/users/george/phased/inds.withnewIDS",header=F,as.is=T)
indsold <- scan("/mnt/kwiat/data/bayes/users/george/phased/inds",what="char")
indsold2 <- read.table("~/R/Copy/Rprojects/AfricaPOPGEN/data/MalariaGen23EthnicGroups1KGSouthAfricaFinalIndsAndPops.txt",as.is=T)

indsold3 <- read.table("~/Data/1000G/Omni25_genotypes_2141_samples.b37.sample",as.is=T,header=1)

indsnew <-c()
for(i in 1:length(indsold))
{
    oldind <- indsold[i]
    tmp <- inds[inds[,1]==oldind,]
    if(nrow(tmp) == 0)
    {
        tmp <- sapply(indsold3[indsold3[,1]==oldind,c(1,6,7)],as.character)
        if(tmp[1]!="character(0)")
        {
            tmp <- c(oldind,oldind,tmp[2],tmp[2],tmp[3])
        } else
        {
            tmp <- sapply(indsold2[indsold2[,1]==oldind,c(1,2)],as.character)
            tmp <- c(oldind,oldind,tmp[2],tmp[2],tmp[2])
        }
    }
    indsnew <- rbind(indsnew,tmp)
}

indsnew[grep("ARI",indsnew[,3]),3:5] <- "ARI"
indsnew[grep("ASHANTI",indsnew[,4]),4] <- "AKANS"

indsnew[grep("EUR",indsnew[,4]),4] <- indsnew[grep("EUR",indsnew[,4]),3]
indsnew[grep("AMR",indsnew[,4]),4] <- indsnew[grep("AMR",indsnew[,4]),3]
indsnew[grep("ASN",indsnew[,4]),4] <- indsnew[grep("ASN",indsnew[,4]),3]
indsnew[grep("AFR",indsnew[,4]),4] <- indsnew[grep("AFR",indsnew[,4]),3]

indsnew[is.na(indsnew[,5]),5] <- 0
indsnew[indsnew[,5]=="M",5] <- 1
indsnew[indsnew[,5]=="F",5] <- 2
indsnew[(indsnew[,5]!=1 & indsnew[,5] != 2 & indsnew[,5]!=0),5] <- 0

colnames(indsnew) <- c("oldID","newID","Country","EthnicGp","Sex")

popkey <- read.table("~/repos/popgen/data/MalariaGenAdmixturePopulationKey.txt",header=T,as.is=T,comment.char="")

for(i in unique(indsnew$Country))
{
    if(i %in% popkey$EthnicGroup)
    {
        ii <- popkey$Country[popkey$EthnicGroup == i]
        indsnew$Country[indsnew$Country == i] <- ii
    } else if (i %in% c("YRI","LWK","KHOMANI","GUIGHANAKGAL","JUHOAN","SWBANTU"))
    {
        if(i == "YRI") ii <- "Nigeria"
        if(i == "LWK") ii <- "Kenya"
        if(i == "KHOMANI") ii <- "SouthAfrica"
        if(i == "GUIGHANAKGAL") ii <- "Botswana"
        if(i == "JUHOAN") ii <- "Namibia"
        if(i == "SWBANTU") ii <- "Namibia"
        indsnew$Country[indsnew$Country == i] <- ii
    }
}



oldsamp <- read.table("/mnt/kwiat/data/bayes/users/george/phased/haps/MalariaGen23EthnicGroupsSouthAfricaFinalChrom22phased.sample.gz",header=F,skip=2,fill=T,as.is=T)

iinds <- oldsamp[,1:3]
fams <- oldsamp[,11:14]
fams[is.na(fams)]<- 0
colnames(fams) <- c("Family","Relationship","Father","Mother")

indsnewfinal <- cbind(indsnew[,2],indsnew[,2],iinds[,3],indsnew,fams)
colnames(indsnewfinal)[1:3] <- c("ID_1","ID_2","missing")

indsnewfinal <- rbind(c(0,0,0,rep("D",ncol(indsnewfinal)-3)),
                      indsnewfinal)
write.table(indsnewfinal,file="data/MalariaGen23EthnicGroups1KGSouthAfricaFinalPHASINGsamples.txt",col.names=T,row.names=F,quote=F)



