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


in_file = '/well/malariagen/malariagen/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintings.hdf5'
#out_file = '/well/malariagen/malariagen/human/george/copy_selection/hdf5files/MalariaGenSelectionPaintingsLengths.hdf5'
out_file = in_file

file = h5py.File(in_file, 'r')

## DEFINE DATASET
chrom  = 'TEMPCHROM'
nsamps = 10
nolocal = file['/paintings/chrom' + str(chrom) + '/nonlocal/']
local = file['/paintings/chrom' + str(chrom) + '/local/']
snps = file['/paintings/chrom' + str(chrom) + '/snps/']
maps = file['/paintings/chrom' + str(chrom) + '/map/']
inds = file['/paintings/samples/individuals']
regions = file['/paintings/samples/regions']

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

## LET'S DO ALL SAMPLES
sampindices = [hap2sampindex(range(len(popinfo))[i]) for i in popinfo]

## FOR EACH HAP AND SNP IN data WE WANT TO FIND
## WHO IS BEING COPIED

nsnps = snps.shape[0]
ninds = len(sampindices)
popwhocopies = np.empty(shape=(nsnps,ninds),dtype="S30")
regwhocopies = np.empty(shape=(nsnps,ninds),dtype="S30")
regwhocopieslengths = np.zeros(shape=(nsnps,ninds),dtype=float)
whocopiedcopies = np.zeros(shape=(nsnps,ninds),dtype=int)
popcopiedcopies = np.empty(shape=(nsnps,ninds),dtype="S30")
regcopiedcopies = np.empty(shape=(nsnps,ninds),dtype="S30")
regcopiedcopieslengths = np.zeros(shape=(nsnps,ninds),dtype=float)

for hap in np.arange(ninds):
    print('getting data for hap ' + str(hap))
    hap1 = sampindices[hap]
    ## THE POPULATION THAT IS COPIED BY THE RECIPIENT
    who = nolocal[:,hap1]-1
    popwho = np.array(inds2)[who].tolist()
    popwhocopies[:,hap] = popwho
    regwho = [reginfo[i] for i in popwho]
    regwhocopies[:,hap] = regwho
    regwhocopieslengths[:,hap] = findLengths(regwho).tolist()
    ## THE POPULATION THAT IS COPIED BY THE COPIER    
    whocopied = [local[j,hap2sampindex(who[j])] for j in range(nsnps)]
    ##
    whocopiedcopies[:,hap] = whocopied
    pop = np.array(inds2)[np.array(whocopied)-1] 
    popcopiedcopies[:,hap] = pop.tolist()
    reg = [reginfo[i] for i in pop]
    regcopiedcopies[:,hap] = reg
    regcopiedcopieslengths[:,hap] = findLengths(reg).tolist()


## NOW COMPUTE SWITCHING AND CHUNKCOUNT INFO
print("computing switches for chromosome: " + chrom)
switch = np.array(nolocal[0:(nolocal.shape[0]-1),:]) != np.array(nolocal[1:(nolocal.shape[0]),:])
## LET'S GET THREE ESTIMATES OF # SWITCHES
firstsample = [(x,x+1,x+2,x+3,x+4) for x in range(nolocal.shape[1]) if x % 10 == 0]
firstsample = [item for sublist in firstsample for item in sublist]
counts = switch.sum(axis = 1)
counts1 = switch[:,firstsample].sum(axis = 1)
counts = np.append(counts,0)
counts1 = np.append(counts1,0)
switches = np.empty(shape=(nolocal.shape[0],2),dtype="int")
switches[:,0] = counts
switches[:,1] = counts1

## NOW INFER CHUNKCOUNTS PER HAPLOTYPE
chunkcounts = switch.sum(axis=0)
## AVERAGE ACROSS PAINTINGS
chunkcounts = [np.mean(chunkcounts[x:(x+10)]) for x in range(len(chunkcounts)) if x % 10 == 0]
## NB WE CAN USE THIS TO COMPARE TO THE TOTAL NUMBER OF CHUNKS
## THAT AN INDIVIDUAL IS INFERRED TO HAVE FROM CHROMOPAINTER

file.close()

###################################################################
## NOW WRITE THESE OUT TO THE out_file
hgrp = "lengths/chrom" + str(chrom) 
if os.path.isfile(out_file) == False :
    fout = h5py.File(out_file,'w')
else:
    fout = h5py.File(out_file,'a')

if hgrp not in fout:
    grp = fout.create_group(hgrp)
else:
    grp = fout[hgrp]
    
## LIST OUTPUT DATASETS, DATA AND DEFINE TYPES FOR THE DATASET
dsets = ['popcopied', 'regcopied', 'copiercopies', 'popcopiercopies', 'regcopiercopies', 'switches', 'chunkcounts']
ddata = [popwhocopies,regwhocopies,whocopiedcopies,popcopiedcopies,regcopiedcopies,switches,chunkcounts]
dsettypes = ['S30', 'S30', 'i', 'S30', 'S30', 'i' , 'float']

for d in range(len(dsets)):
    dset = hgrp + '/' + dsets[d]
    if dset in fout:
        del fout[dset]
    grp.create_dataset(dsets[d],data=ddata[d],compression="gzip",compression_opts=9,dtype=dsettypes[d])

fout.close()










