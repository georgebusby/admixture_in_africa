##########################################################################
#### SCRIPT TO GENERATE CHUNKLENGTH COPYING VECTOR FILES FOR ADMIXTURE ###
#### SIMULATIONS;
#### SIM CHUNKLLENGTH FILE ARE BROUGHT IN AND A SINGLE COPYING VECTOR IS 
#### GENERATED FOR EACH SIM, WHICH IS THE ADDED TO THE MASTER CV FILE FOR
#### EACH POP1/POP2 PAIR

#### NOTE THAT BECAUSE THE SIMS WERE PAINTED WITHOUT INCLUDING POP2
#### THIS POPULATION IS REMOVED FROM THE ORIGINAL CV PAINTING FILE
#### AND THEN RE-NORMALISED

pops <- c("MALAWI","JOLA","JUHOAN","CEU")
gens <- c(100,200,300,400)
props <- c(95,90,80)

length_dir <- "/mnt/kwiat/well/human/george/admix_sims/chromopainter/output/"
cp_dir <- "/mnt/kwiat/well/human/george/globetrotter/input/"
out_dir <- "/mnt/kwiat/well/human/george/admix_sims/globetrotter/lengths/"

for(pop1 in pops)
{
    for(pop2 in pops)
    {
        if(pop1 != pop2 & pop1 != "CEU")
        {
            sim <- paste0(pop1,pop2)
            out_file <- paste0(out_dir,sim,".chunklengths.out")
            print(paste0("generating copying vectors for ",sim))
            xs <- c()
            for(gen in gens)
            {
                for(prop in props)
                {
                    sim_pop <- paste0(sim,gen,"g",prop,"p")
                    file_check <- 0 ## add to this if there's an empty file
                    for(chr in 1:22)
                    {
                        if(chr < 10) chr <- paste0("0",chr)
                        sim_file <- paste0(length_dir,sim_pop,"Chrom",chr,".chunklengths.out")
                        if(file.exists(sim_file)) file_size <- file.info(sim_file)$size
                        if(file_size == 0) file_check <- file_check + 1
                        if(file.exists(sim_file) & file_size > 0)
                        {
                            if(chr == "01")
                            {
                                x1 <- read.table(sim_file,header=T,row.names=1)
                            } else
                            {
                                x1 <- x1+read.table(sim_file,header=T,row.names=1)
                            }
                        }
                    }
                    if(file_check == 0)
                    {
                        x1 <- apply(x1,2,mean)
                        xs <- rbind(xs,x1)
                        rownames(xs)[nrow(xs)] <- sim_pop
                    }
                }
            }
            if(!is.null(xs))
            {
                ## NOW ADD TO THE COPYING VECTORS OF THE ORIGINAL ANALYSIS
                if(pop1 == "MALAWI") original_lengths <- "South_Africa_Niger-Congo"
                if(pop1 == "JOLA") original_lengths <- "Western_Africa_Niger-Congo"
                if(pop1 == "JUHOAN") original_lengths <- "South_Africa_KhoeSan"
                lengths <- paste0(cp_dir,
                                  "MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinal",
                                  original_lengths,".chunklengths.out")
                x <- read.table(lengths,header=T,row.names=1)
                ## REMOVE pop2 COLUMN AS THE SIMS WERENT' PAINTED WITH THIS POP
                x <- x[,colnames(x)!=pop2]
                x <- rbind(x,xs[,colnames(x)])
                x <- x/rowSums(x)
                x <- cbind(rownames(x),x)
                colnames(x)[1] <- "Recipient"
                write.table(x,file=out_file,col.names=T,row.names=F,quote=F)
            }
        }
    }
}

