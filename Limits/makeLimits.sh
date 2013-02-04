#!/bin/bash

#source ../../../StatisticalTools/RooStatsRoutines/setup/cmslpc_standalone_setup.sh
export ROOT_INCLUDE=${ROOTSYS}/include

root -b -l -q -n 'ExpectedEvts.C+("../../../WprimeWZ.root","cutValues.wz.dat",-1, "scaleMC")' > ExpectedEvts.log 

if [ ! -e nEvents.txt ]; then
    exit 1
fi

echo "Done with counting events"

root -b -l -q -n 'CalcLimit.C(1)' > CalcLimit.log #CLs
#root -b -l -q -n 'CalcLimit.C(0)' > CalcLimit.log #Bay

echo "Done with calculating limits"

root -b -l -q -n 'PlotLimit.C+("WprimeWZ")' #For W'
root -b -l -q -n 'PlotLimit2D.C+' #For TC
#root -b -l -q -n 'PlotLimit.C+("HadVZ")' #For VZ
root -b -l -q plotExclusion.C  

echo "Done with Plotting limits"

Dir=LatestNumbers-`date +%Y-%m-%d`
if [ ! -d "$Dir" ]; then
    mkdir $Dir
fi
mv *.txt *.log *.pdf $Dir/

