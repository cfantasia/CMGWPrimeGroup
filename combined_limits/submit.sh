#!/bin/bash -f

#FullRange="`seq 150 10 190`"
FullRange="`seq 170 10 250` `seq 275 25 350` `seq 400 50 600` `seq 700 100 2000`"
#FullRange="`seq 400 50 600` `seq 700 100 2000`"
#FullRange="`seq 200 100 2000`"
Range="`seq 200 100 2000`"

echo $Range

for Mass in $FullRange; do
    echo "Mass $Mass"
    ./doCombinedLimits.py -M $Mass --doSepCh --onlyMakeCards
done

for Mass in $FullRange; do
    echo "Mass $Mass"
    #./doCombinedLimits.py -M $Mass --ntoys=100 
    ./doCombinedLimits.py -M $Mass --ntoys=1000 --condor --njobs=10 --doSepCh --doDuplicates
    #./doCombinedLimits.py -M $Mass --observed
done

for Mass in $FullRange; do
    echo "Mass $Mass"
    #./doCombinedLimits.py -M $Mass --observed 
    #./doCombinedLimits.py -M $Mass --observed --doDuplicates
    ./doCombinedLimits.py -M $Mass --observed --doSepCh 
    #./doCombinedLimits.py -M $Mass --observed --doSepCh --doDuplicates
done

for Mass in $FullRange; do
    echo "Mass $Mass"
    #./doCombinedLimits.py -M $Mass --observed --doSepCh
done


#for Mass in 1000; do
#   for LtShift in `seq -100 50 100`; do
#       for WindShift in `seq -100 50 100`; do
#           #./doCombinedLimits.py -M $Mass --LtShift=$LtShift --WindShift=$WindShift
#           ./doCombinedLimits.py --LtShift=$LtShift --WindShift=$WindShift -p
#       done
#   done
#done    
