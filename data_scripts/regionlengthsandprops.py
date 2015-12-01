#!/usr/bin/python

############################################################
## SCRIPT TO GENERATE LENGTH STATISTICS BASED ON CHROMOPAINTER
## PAINTINGS: THE SCRIPT ALSO CONVERTS PAINTING SAMPLES BY IND
## TO PAINTINGS BY POP AND REGION -- ALL OUTPUT IS STORED IN A
## SINGLE HDF5 FILE
## TO DO: PERHAPS DO THE IND --> POP PAINTING CONVERSION ON THE
## FLY??                            GEORGE BUSBY 24/11/2015
############################################################
import h5py # for reading/writing hdf5 files
import gzip # for reading/writing gzipped files
from optparse import OptionParser # for command line options
import numpy as np
import os.path ## for checking if file exists
import string,sys

pop = 'FULAI'
chrom  = 'TEMPCHROM'

in_file = '/well/malariagen/malariagen/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
file = h5py.File(in_file, 'r')
out_file = '/well/malariagen/malariagen/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings' + pop + '.hdf5'


## DEFINE DATASET
#pop = 'TEMPPOP'
#chrom = '02'
nsamps = 10
nolocal = file['/paintings/chrom' + str(chrom) + '/nonlocal/']
local = file['/paintings/chrom' + str(chrom) + '/local/']
snps = file['/paintings/chrom' + str(chrom) + '/snps/']
maps = file['/paintings/chrom' + str(chrom) + '/map/']
inds = file['/paintings/samples/individuals']
regions = file['/paintings/samples/regions']
popcopied = file['/lengths/chrom' + chrom + '/popcopied']
regcopied = file['/lengths/chrom' + chrom + '/regcopied']
copiercopies = file['/lengths/chrom' + chrom + '/copiercopies']
popcopiercopies = file['/lengths/chrom' + chrom + '/popcopiercopies']
regcopiercopies = file['/lengths/chrom' + chrom + '/regcopiercopies']
switches = file['/lengths/chrom' + chrom + '/switches'][:,0]


##########################################
## FUNCTIONS
## CONVERT A HAP INDEX TO A INDEX OF THE PAINTING SAMPLE MATRIX
def hap2sampindex (x,nsamps=10):
    y = ((x+1)*nsamps) - nsamps
    return(y)
  
## CONVERT A PAINTED HAPPLOTYPE INTO LENGTHS COPIED AT EACH SNP  
def findLengths(paintedhap, maps = maps):
    start = 0
    totallength = 0
    haplengths = np.empty(len(paintedhap))
    while start < len(paintedhap) - 1:
        end = start + 1
        while paintedhap[start] == paintedhap[end] and end+1 < len(paintedhap) :
            end = end + 1
        haplength = end - start
        lengthbp = int(maps[end,0]) - int(maps[start,0])
        lengthgen = sum(maps[start:end,1].astype(np.float))
        chunklength = lengthgen*lengthbp
        totallength = totallength + chunklength
        haplengths[start:end] = chunklength 
        start = end
    haplengths[end] = chunklength
    return(haplengths)

## CONVERT A PAINTED HAPPLOTYPE INTO A SEQUENCE OF UNIQUE INDICES
def findUniqueHaps(paintedhap):
    hapnum = 0
    start = 0
    haplengths = np.zeros(len(paintedhap), dtype='i')
    while start < len(paintedhap) - 1:
        end = start + 1
        while paintedhap[start] == paintedhap[end] and end+1 < len(paintedhap) :
            end = end + 1
        haplengths[start:end] = hapnum
        start = end
        hapnum = hapnum + 1
    haplengths[end] = hapnum - 1 
    return(haplengths)


##########################################  
## FIND OUT THE POPULATION ORIGIN OF EACH HAPLOTYPE
popinfo = {}
inds2 = []
d = [r.split() for r in inds[:,1]]
for i in range(len(d)):
    popinfo[2*i] = d[i][0]
    popinfo[2*i+1] = d[i][0]
    inds2.append(d[i][0])
    inds2.append(d[i][0])


reginfo = {x[0]:x[2] for x in regions}
## ANNOYING DIFFERENCE WITH SEMI-BANTU
reginfo["SEMI.BANTU"] = reginfo["SEMI-BANTU"]
del reginfo["SEMI-BANTU"]
    
## GET A LIST OF SNP POSITIONS
pos = [r.split() for r in snps[:,2]] #read the position that we do not know if we really need that
## GET THE RECOMBINATION DISTANCE BETWEEN SNPS
rec = [r.split() for r in maps[:,1]]

## IDENTIFY INDICES FOR POP
sampindices = [i for i in popinfo if popinfo[i] == pop]

## FOR EACH HAP AND SNP IN data WE WANT TO FIND
## WHO IS BEING COPIED
nsnps = snps.shape[0]
ninds = len(sampindices)

## COMPUTE THE SWITCHES PER BP
pos = [float(i) for i in maps[:,0]]
rec = [float(i) for i in maps[:,1]]
posdiff = np.diff(pos)
switchrate = switches[range(len(posdiff))]/posdiff
switchrate = np.append(switchrate,0)


## SET EVERYTHING UP
region = 'Eurasia'
paintings = regcopied[:,sampindices]
isregion = paintings == region
paintedhaps = np.zeros(shape=paintings.shape,dtype='i')
for hap in np.arange(ninds):
    print('getting data for hap ' + str(hap))
    paintedhaps[:,hap] = findUniqueHaps(paintings[:,hap])

## NOW INFER LENGTH AND PROPORTION AROUND EACH SNP
freqs = np.zeros(shape=(nsnps,1))
switchlengths = np.zeros(shape=(nsnps,1))
reclengths = np.zeros(shape=(nsnps,1))
for snp in np.arange(nsnps):
    print('geting props for snp ' +  str(snp))
    ## SUBSET TO COLS THAT COPY FROM THE REGION
    hapsfromreg = paintedhaps[snp,isregion[snp,:]==True]
    snpfreq = float(len(hapsfromreg))/ninds
    test = paintedhaps[:,isregion[snp,:]==True] == hapsfromreg
    switchlength = np.sum(np.sum(test,axis=1)*switchrate)
    reclength = np.sum(np.sum(test,axis=1)*rec)
    freqs[snp] = snpfreq
    switchlengths[snp] = switchlength
    reclengths[snp] = reclength
 
###################################################################
## NOW WRITE THESE OUT TO THE out_file
hgrp = "lengths/chrom" + str(chrom) + '/' + pop + '/' + region
if os.path.isfile(out_file) == False :
    fout = h5py.File(out_file,'w')
else:
    fout = h5py.File(out_file,'a')

if hgrp not in fout:
    grp = fout.create_group(hgrp)
else:
    grp = fout[hgrp]

hdset1 = 'switches'
out_data = np.hstack((freqs,switchlengths,reclengths))
grp.create_dataset(hdset1,data=out_data, compression="gzip",compression_opts=9, dtype = "float")
fout.close()










