#!/usr/bin/python

############################################################
## THIS SCRIPT TAKES CHROMOPAINTER OUTPUT AND GENERATES  ##
## A MATRIX OF GENOTYPES THAT CAN BE USED IN THE NATURAL ##
## SELECTION PROGRAMS
############################################################
import h5py # for reading/writing hdf5 files
import gzip # for reading/writing gzipped files
from optparse import OptionParser # for command line options
import numpy as np
import os.path ## for checking if file exists

#file_root = "/mnt/kwiat/well/human/george/chromopainter2/output/FULAInolocalChrom"
#out_root = "/mnt/kwiat/data/bayes/users/george/popgen/analysis3/chromopainter/outputcopyprobs/FULAInolocalAllChroms"
############################################################
# get commandline arguments
#usage = "usage: %prog [options] \n\n **This program takes in ChromoPainter *samples.out file root and generates a matrix of genome-wide paintings with n_haps columns and n_snps rows."
#parser = OptionParser()
#parser = OptionParser(usage=usage, version="%prog 1.0")
#parser.add_option("-i", "--infile", dest="infile",
                  #help="input file root: should be the bit pre *ChromXX.samples.out file from ChromoPainter", metavar="INFILE")
#parser.add_option("-o", "--outfile", dest="outfile", 
                  #help="output file (gzipped)", metavar="OUTROOT")
#parser.add_option("-s", "--samps", dest="numsamps", default = 10,
                  #help="the number of samples generated by ChromoPainter: default is 10", metavar="SAMPS")
#parser.add_option("-q", "--quiet",
                  #action="store_false", dest="verbose", default=True,
                  #help="don't print status messages to stdout")

#(options, args) = parser.parse_args()
#file_root = options.infile
#n_samps = options.numsamps
#out_file = options.outfile

## NON-LOCAL PAINTINGS
file_root1 = '/well/malariagen/malariagen/human/george/chromopainter2/output/AllPopsnolocalChrom'
file_suff1  = "PP.samples.out"

## LOCAL PAINTINGS
file_root2 = '/well/malariagen/malariagen/human/george/chromopainter2/output/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalAllPopsChrom'
file_suff2  = ".samples.out"

## MAPS and RECOMRATES
snp_root = '/data/bayes/users/george/popgen/analysis3/chromopainter/snpfiles/AllPops330KChrom'
snp_suff = 'phased.legend.gz'
map_root = '/well/malariagen/malariagen/human/george/chromopainter2/input/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCPChrom'
map_suff = '.recomrates'

## SAMPLE FILE
sample_file = '/well/malariagen/malariagen/human/george/chromopainter2/analysislists/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCPv2.idfile.txt'
region_file = '/well/malariagen/malariagen/human/george/chromopainter2/analysislists/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCPv2.regions.txt'

## OUTPUT FILE
out_file = '/well/malariagen/malariagen/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
n_samps = 10

   
############################################################
# PROGRAM #

## NB I'VE EDITED THE BELOW FOR THE ARI SIMS 

for chrom in range(1,22):
    if chrom < 10:
        chrom = "0" + str(chrom)
    else:
	chrom = str(chrom)

    in_file1 = file_root1 + chrom + file_suff1
    in_file2 = file_root2 + chrom + file_suff2
    snp_file = snp_root + chrom + snp_suff
    map_file = map_root + chrom + map_suff
    
    hgrp = "paintings/chrom" + chrom
    if os.path.isfile(out_file) == False :
        fout = h5py.File(out_file,'w')
    else:
        fout = h5py.File(out_file,'a')
    
    if hgrp not in fout:
        grp = fout.create_group(hgrp)
    else:
        grp = fout[hgrp]
        
    hdset1 = "nonlocal" 
    hdset2 = "local"
    hdset3 = 'snps'
    hdset4 = 'map'

    print("switching haps to hdf5 for chrom: " + chrom)
############################################################
## FIND NUMBER OF SNPS
    f = open(in_file1)
    for i, line in enumerate(f):
	if i == 2:
	    n_snps = len(line.split())-1
    print('chromosome ' + chrom + ' sample file has: ' + str(n_snps) + ' snps...')
############################################################
## FIND NUMBER OF HAPS
    f = open(in_file1)
    l = [x for x in f.readlines() if x.startswith('HAP')]
    n_haps = len(l)
    print 'sample file has: ' + str(n_haps) + ' haplotypes...'
############################################################
## NOW READ IN ALL NON-LOCALLY COPIED SAMPLES AND MAKE A n_haps * 10 by n_snps matrix
    #print("outputting NON-LOCAL paintings to hdf5 file")
    #mat = np.zeros(shape=(n_snps,n_haps*n_samps),dtype=int)
    #i = 0
    #for line in open(in_file1):
	#l = line.split()
	#if l[0].isdigit():
	    #mat[:,i] =  l[1:len(l)]
	    #i = i + 1
    ##with gzip.open(out_file2,'wb') as f_handle:
        ##np.savetxt(f_handle,mat, fmt='%1s')

    ### OUTPUT TO HDF5 FILE
    #dset = grp.create_dataset(hdset1,data=mat, dtype="i", compression="gzip",compression_opts=9)
    
#############################################################
### NOW READ IN ALL LOCALLY COPIED SAMPLES AND MAKE A n_haps * 10 by n_snps matrix
    #print("outputting LOCAL paintings to hdf5 file")
    #mat = np.zeros(shape=(n_snps,n_haps*n_samps),dtype=int)
    #i = 0
    #for line in open(in_file2):
	#l = line.split()
	#if l[0].isdigit():
	    #mat[:,i] =  l[1:len(l)]
	    #i = i + 1
    ##with gzip.open(out_file2,'wb') as f_handle:
        ##np.savetxt(f_handle,mat, fmt='%1s')

    ### OUTPUT TO HDF5 FILE
    #dset = grp.create_dataset(hdset2,data=mat, dtype="i", compression="gzip",compression_opts=9)
    
############################################################
## NOW READ IN SNPS and MAPS files for each chromosome
    print("outputting MAP and RECOMRATES to hdf5 file")

    snps = np.loadtxt(snp_file, dtype= 'S')
    maps = np.loadtxt(map_file, dtype= 'S')
    maps = maps[1:,]


    ## OUTPUT TO HDF5 FILE
    dset = grp.create_dataset(hdset3,data=snps, compression="gzip",compression_opts=9)
    dset = grp.create_dataset(hdset4,data=maps, compression="gzip",compression_opts=9)
    
    print("\n!!CHROMOSOME " + chrom + " FINISHED!!\n")
    fout.close()
    
############################################################
## FINALLY, ADD SOME INFO ABOUT THE SAMPLES
hgrp = "paintings/samples"
fout = h5py.File(out_file,'a')
if hgrp not in fout:
    grp = fout.create_group(hgrp)
else:
    grp = fout[hgrp]

hdset1 = 'individuals'
hdset2 = 'regions'

inds = np.loadtxt(sample_file, dtype = 'S')
regs = np.loadtxt(region_file, dtype = 'S')
regs = regs[1:,]
dset = grp.create_dataset(hdset1,data=inds, compression="gzip",compression_opts=9)
dset = grp.create_dataset(hdset2,data=regs, compression="gzip",compression_opts=9)    

fout.close()

    
    
    
    
    
    