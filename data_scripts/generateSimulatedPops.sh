#!/bin/bash

admixsimulator=/well/malariagen/malariagen/human/george/admix_sims/FakeAdmixPopGeneratorEXPSimonMultiChromSwitchInd.pl
dataindir=/well/malariagen/malariagen/human/george/admix_sims/chromopainter/input
simoutdir=/well/malariagen/malariagen/human/george/admix_sims/chromopainter/input_sims
recdir=/well/malariagen/malariagen/human/george/chromopainter2/input
recinroot=MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCP

nhaps=200
pops=(MALAWI JOLA JUHOAN CEU)
gens=(100 200 300 400)
props=(95 90 80)

simulate=0
merge=1

## FIRST SUBSET MAIN FILE TO HAVE ONLY HAPS FROM EACH OF THE POPS IN $pops
## I DID THIS ON CLUSTER USING VARIANTS OF THE FOLLOWING COMMAND
##pop=MALAWI
##for chr in $(seq -w 1 22); do
## echo $chr
##  python /home/george/popgen/scripts/mergechromopainter.py \
##    -p /well/malariagen/malariagen/human/george/admix_sims/chromopainter/input/${pop}.pops \
##    -c ${chr} \
##    -d /well/malariagen/malariagen/human/george/chromopainter2/input/ \
##    -o /well/malariagen/malariagen/human/george/admix_sims/chromopainter/input/ \
##     -n ${pop}
##done


#if [[ $simulate == 1 ]]; then
jobs=/well/malariagen/malariagen/human/george/admix_sims/simjobs.txt
for pop1 in ${pops[@]} ; do
    for pop2 in ${pops[@]} ; do
        for gen in ${gens[@]} ; do
            for prop in ${props[@]} ; do
                if [[ $pop1 != $pop2 && $pop1 != 'CEU' ]] ; then
                    echo "generating files for $pop1 v $pop2 at $gen generations with proportion $prop"
                    pop1file=${dataindir}/${pop1}
                    pop2file=${dataindir}/${pop2}
                    recomfile=${recdir}/${recinroot}
                    outfile=${simoutdir}/${pop1}${pop2}${gen}g${prop}p
                    echo "perl $admixsimulator $outfile $nhaps $gen 0.$prop $pop1file $pop2file $recomfile" >> $jobs
                fi
             done
        done
    done
done
#fi

## sed -i '1s/^/200\n/' /well/malariagen/malariagen/human/george/admix_sims/chromopainter/input_sims/*phase

## NOW MERGE ALL SIMS USING THE SAME POP-PAIRS TOGETHER WITH THE MAIN dataindir
## THERE ARE 7 DIFFERENT POP-PAIRS

## NB - THE ADMIX SIMULATOR GIVES A HEADER WITH ONLY 4 LINES, NEED TO ADD A LINE TO THE BEGINNING
## OF THE FILE EG. sed '1s/^/200\n/ ' $FILE
## THIS NEEDS TO BE REMEDIED LATER TO HAVE ONLY 3 LINES AS PER CHROMOSPAINTER v2 STYLE

## POPLISTS
admixpoplist=/well/malariagen/malariagen/human/george/admix_sims/chromopainter/admix_sims.pops
ls /well/malariagen/malariagen/human/george/admix_sims/chromopainter/input_sims/*Chrom22*phase \
    | awk -F '/' '{print $10}' | sed 's/Chrom/ /1' | awk '{print $1}' \
    | sed '/^Mala/ d' > $admixpoplist
    


#if [[ $merge == 1 ]]; then
for pop1 in ${pops[@]} ; do
    for pop2 in ${pops[@]} ; do
        if [[ $pop1 != $pop2 && $pop1 != 'CEU' ]] ; then
            pop=${pop1}${pop2}
            popfile=/well/malariagen/malariagen/human/george/admix_sims/chromopainter/input_sims/${pop}.pops
            if [[ -f $popfile ]]; then
                rm $popfile
            fi
            echo MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCP > $popfile
            grep ^$pop $admixpoplist >> $popfile
            for chr in $(seq -w 1 22); do
                echo $pop $chr
                python /home/george/popgen/scripts/mergechromopainter.py \
                        -p ${popfile} \
                        -c ${chr} \
                        -d /well/malariagen/malariagen/human/george/admix_sims/chromopainter/input_sims/ \
                        -o /well/malariagen/malariagen/human/george/admix_sims/chromopainter/combined_sims/ \
                        -n ${pop}
            done
        fi
    done
done
#fi

## NOW MAKE ID AND POPLIST FILES 
analydir=/well/malariagen/malariagen/human/george/admix_sims/chromopainter/analysislists
while read i ; do
    for hap in $(seq 1 100); do 
        echo IND_${hap} $i 1 
    done > ${analydir}/${i}.inds ;
done < $admixpoplist

for pop1 in ${pops[@]} ; do
    for pop2 in ${pops[@]} ; do
        if [[ $pop1 != $pop2 && $pop1 != 'CEU' ]] ; then
            ## ID FILE
            outfile=${analydir}/${pop1}${pop2}.idfile.txt
            if [[ -f $outfile ]] ; then
                rm $outfile
            fi
            cat ${analydir}/${pop1}nolocal.idfile.txt \
            ${analydir}/${pop1}${pop2}*inds \
            > $outfile
        fi
    done
done

for pop1 in ${pops[@]} ; do
    popslist=$(grep ^${pop1} $admixpoplist)
    for i in $popslist ; do
        echo $i;
        cp ${analydir}/${pop1}nolocal.pops ${analydir}/${i}nolocal.pops
        ## remove pop1 from being a recipient
        sed -i '/^'${pop1}'/ d' ${analydir}/${i}nolocal.pops
        ## remove the minor pop from being a donor
        pop2=$(echo $i | sed 's/'${pop1}'//' | sed 's/[0-9]//g' | sed 's/gp//' )
        sed -i '/^'${pop2}'/ d' ${analydir}/${i}nolocal.pops
        echo "$i R" >> ${analydir}/${i}nolocal.pops
    done
done


## THE BELOW CODE SWITCHES THE THIRD COLUMN IN A FILE TO 0 WHEN THE SECOND COLUMN MATCHES JUHOAN
## awk '{ if ($2 == "JUHOAN") printf("%s %s 0\n", $1, $2); else printf("%s\n", $0);}'

## NEED TO REMOVE SECOND AND FIFTH LINES AND ADD TOTAL NUMBER OF HAPS IN FILE TO
## FIRST LINE (8966) HAPS
## DO THIS ONLY ONCE!
sed -i -e '2d;5d' -e '1s/0/8966/' /well/malariagen/malariagen/human/george/admix_sims/chromopainter/combined_sims/*phase


## cp commanline

#~george/programs/ChromoPainterv2 \
#-g /well/malariagen/malariagen/human/george/admix_sims/chromopainter/combined_sims/TEMPPOP1TEMPPOP2ChromTEMPCHROMphased.phase \
#-r /well/malariagen/malariagen/human/george/admix_sims/chromopainter/input_sims/MalariaGen23EthnicGroups1KGSouthAfricaNoAmericaFinalCPChromTEMPCHROM.recomrates \
#-t /well/malariagen/malariagen/human/george/admix_sims/chromopainter/analysislists/TEMPPOP1TEMPPOP2.idfile \
#-f /well/malariagen/malariagen/human/george/admix_sims/chromopainter/analysislists/TEMPPOP1TEMPPOP2TEMPANALY 0 0 \
#-o /well/malariagen/malariagen/human/george/admix_sims/chromopainter/output/TEMPPOP1TEMPPOP2TEMPANALYChromTEMPCHROM \
#-s 10 -n 196.853660843533 -M 0.000794001328076925

## on rescomp, run the following



pops=(MALAWI JOLA JUHOAN CEU)
gens=(100 200 300 400)
props=(95 90 80)
commandline=/users/kwiatkowski/george/scripts/admixsimcpcommandline.txt
outparams=/users/kwiatkowski/george/scripts/admixsimcp.txt
if [[ -f $outparams ]]; then
    rm $outparams
fi

for pop1 in ${pops[@]} ; do
    for pop2 in ${pops[@]} ; do
        for gen in ${gens[@]} ; do
            for prop in ${props[@]} ; do
                for chrom in $(seq -w 1 22); do
                    if [[ $pop1 != $pop2 && $pop1 != 'CEU' ]] ; then
                        sed -e 's/TEMPPOP1/'${pop1}'/g' \
                            -e 's/TEMPPOP2/'${pop2}'/g' \
                            -e 's/TEMPCHROM/'${chrom}'/g' \
                            -e 's/TEMPANALY/'${gen}g${prop}p'/g' \
                            $commandline >> $outparams
                    fi
                done
            done
       done
   done
done
      
		

## GENERATE MALDER PARAMFILES
pops=(MALAWI JOLA JUHOAN CEU)
gens=(100 200 300 400)
props=(95 90 80)
temp_param=/well/malariagen/malariagen/human/george/admix_sims/malder/paramfiles/TEMP_alderallrefsmalder.par
temp_run=/well/malariagen/malariagen/human/george/admix_sims/malder/scripts/alderTEMPLATEallrefsmalder.sh
for pop1 in ${pops[@]} ; do
    for pop2 in ${pops[@]} ; do
        ## make a pops file
        if [[ $pop1 != $pop2 && $pop1 != 'CEU' ]] ; then
            awk '{print $2}' /well/malariagen/malariagen/human/george/admix_sims/chromopainter/analysislists/${pop1}${pop2}.idfile.txt | \
                sort | uniq | grep -v ${pop1}${pop2} | grep -v ${pop1} | grep -v ${pop2} > \
                /well/malariagen/malariagen/human/george/admix_sims/malder/idfiles/${pop1}${pop2}.malderpops
        fi
        for gen in ${gens[@]} ; do
            for prop in ${props[@]} ; do
                if [[ $pop1 != $pop2 && $pop1 != 'CEU' ]] ; then
                        param_file=/well/malariagen/malariagen/human/george/admix_sims/malder/paramfiles/${pop1}${pop2}${gen}g${prop}p_alderallrefsmalder.par
                        sed -e 's/TEMPPOP1/'${pop1}'/g' \
                            -e 's/TEMPPOP2/'${pop2}'/g' \
                            -e 's/TEMPCHROM/'${chrom}'/g' \
                            -e 's/TEMPANALY/'${gen}g${prop}p'/g' \
                            $temp_param > $param_file
                        run_file=/well/malariagen/malariagen/human/george/admix_sims/malder/scripts/run/${pop1}${pop2}${gen}g${prop}p_alderallrefsmalder.sh
                        sed -e 's/TEMPPOP1/'${pop1}'/g' \
                            -e 's/TEMPPOP2/'${pop2}'/g' \
                            -e 's/TEMPCHROM/'${chrom}'/g' \
                            -e 's/TEMPANALY/'${gen}g${prop}p'/g' \
                            $temp_run > $run_file    
                            
                fi
            done
       done
   done
done



## GENERATE GLOBETROTTER PARAMFILES
pops=(MALAWI JOLA JUHOAN CEU)
gens=(100 200 300 400)
props=(95 90 80)
temp_param=/well/malariagen/malariagen/human/george/admix_sims/globetrotter/paramfiles/TEMP.paramfile.txt
for pop1 in ${pops[@]} ; do
    for pop2 in ${pops[@]} ; do
        for gen in ${gens[@]} ; do
            for prop in ${props[@]} ; do
                if [[ $pop1 != $pop2 && $pop1 != 'CEU' ]] ; then
                    param_file=/well/malariagen/malariagen/human/george/admix_sims/globetrotter/paramfiles/${pop1}${pop2}${gen}g${prop}p.paramfile.txt
                    sample_file=/well/malariagen/malariagen/human/george/admix_sims/globetrotter/paramfiles/${pop1}${pop2}${gen}g${prop}p.samples.txt
                    run_temp=/well/malariagen/malariagen/human/george/admix_sims/globetrotter/scripts/globetrotterTEMPLATE.sh
                    cvpops=$(awk '{if ($2 == "D") printf("%s ", $1)}' /well/malariagen/malariagen/human/george/admix_sims/chromopainter/analysislists/${pop1}${pop2}${gen}g${prop}pnolocal.pops) 
                    surpops=$cvpops
                    sed -e 's/TEMPPOP1/'${pop1}'/g' \
                        -e 's/TEMPPOP2/'${pop2}'/g' \
                        -e 's/TEMPANALY/'${gen}g${prop}p'/g' \
                        -e "s/TEMPCVPOPS/${cvpops}/g" \
                        -e "s/TEMPSURPOPS/${surpops}/g" \
                        $temp_param > $param_file
                        
                    ls /well/malariagen/malariagen/human/george/admix_sims/chromopainter/output/${pop1}${pop2}${gen}g${prop}p*samples.out > $sample_file
                    nsamps=$(wc -l $sample_file | awk '{print $1}')
                    
                    if [[ $nsamps == 22 ]]; then
                        run_file=/well/malariagen/malariagen/human/george/admix_sims/globetrotter/scripts/run/${pop1}${pop2}${gen}g${prop}p.props.sh
                        sed -e 's/TEMPPOP1/'${pop1}'/g' \
                            -e 's/TEMPPOP2/'${pop2}'/g' \
                            -e 's/TEMPANALY/'${gen}g${prop}p'/g' \
                            $run_temp > $run_file
                    fi
                fi
            done
       done
   done
done


## GENERATE GLOBETROTTER RUNSCRIPTS






