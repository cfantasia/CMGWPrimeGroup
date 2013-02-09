#!/bin/bash -f

#for Mass in `seq 2000 -100 200`; do
#    #./doCombinedLimits.py -M $Mass --ntoys=100 
#    ./doCombinedLimits.py -M $Mass --ntoys=1000 --condor
#    #./doCombinedLimits.py -M $Mass --observed
#done

#for Mass in `seq 200 100 700`; do
#    ./doCombinedLimits.py -M $Mass --observed
#done

for Mass in 1000; do
   for LtShift in `seq -100 50 100`; do
       for WindShift in `seq -100 50 100`; do
           #./doCombinedLimits.py -M $Mass --LtShift=$LtShift --WindShift=$WindShift
           ./doCombinedLimits.py --LtShift=$LtShift --WindShift=$WindShift -p
       done
   done
done    
