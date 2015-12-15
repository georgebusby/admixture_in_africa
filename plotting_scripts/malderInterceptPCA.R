## try getting more info out of the malder results

m <- all_res2[all_res2[,"pop1"]=="JOLA",]

pops <- m[,"test.pops"]
pop1 <- sapply(pops,function(x)strsplit(x,split="\\;")[[1]][1])
pop2 <- sapply(pops,function(x)strsplit(x,split="\\;")[[1]][2])

icepts <- cbind(pop1,pop2,m[,"amp0"])
colnames(icepts) <- c("pop1","pop2","amp0")
icepts[,"amp0"] <- as.numeric(icepts[,"amp0"])
allpops <- unique(c(pop1,pop2))
imat <- matrix(0,nr=length(allpops),nc=length(allpops))
rownames(imat)<- allpops
colnames(imat) <- allpops
for(i in 1:nrow(icepts))
{
    imat[icepts[i,"pop1"],icepts[i,"pop2"]] <- as.numeric(icepts[i,"amp0"])
}


test=imat
for(i in 1:ncol(test)){
    for(j in 1:ncol(test)){
        imat[i,j]=test[i,j]-colMeans(test)[j]-rowMeans(test)[i]+mean(test);
    }
}

pcnum <- 1
qfitted=eigen(imat)$vectors[,pcnum]*sqrt(eigen(imat)$values[pcnum])
quality=1-sum((qfitted %*% t(qfitted)-imat)^2)/sum(imat^2)

q2 <- Im(qfitted)/sum(Im(qfitted))
pdf("figures/trialInterceptPCAJola.pdf",width=10,height=3)
bp=barplot(q2)
axis(1,at=bp,labels=allpops,las=2,lwd=0)
dev.off()