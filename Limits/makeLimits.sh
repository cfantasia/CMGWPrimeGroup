#!/bin/bash

source ../../../StatisticalTools/RooStatsRoutines/setup/cmslpc_standalone_setup.sh

for windFracTenths in -1
  do
  Ones=$(($windFracTenths / 10))
  Tenths=$(($windFracTenths - $Ones*10))
  windFrac=$Ones"p"$Tenths
  echo "Using $windFrac as half window size"

#  root -b -q 'ExpectedEvts.C+("../../../old-2011-12-24-MET30/WprimeWZ.root","cutValues.wz.dat",'${windFracTenths}', "scaleMC")' > ExpectedEvts.log 
  root -b -q -n 'ExpectedEvts.C+("../../../WprimeWZ.root","cutValues.wz.dat",'${windFracTenths}', "scaleMC")' > ExpectedEvts.log 
#  root -b -q -n 'ExpectedEvts.C+("../../../old-2011-11-11/HadVZAnalyzer.root","cutValues.VZ.dat",'${windFracTenths}', "")' > ExpectedEvts.log 

  echo "Done with counting events"

  #/afs/hep.wisc.edu/cern/.root/root_v5.30.00.Linux-slc5_amd64-gcc4.3/bin/root -b -q -n 'CalcLimit.C+(1)' > CalcLimit.log
  #root -b -q -n 'CalcLimit.C+(1)' > CalcLimit.log #CLs
  root -b -q -n 'CalcLimit.C+(0)' > CalcLimit.log #Bay

  echo "Done with calculating limits"

  root -b -q -n 'PlotLimit.C+("WprimeWZ")' #For W'
  root -b -q -n 'PlotLimit2D.C+' #For TC
#  root -b -q -n 'PlotLimit.C+("HadVZ")' #For VZ
  
  echo "Done with Plotting limits"

  Dir=windFrac${windFrac}
  if [ ! -d "$Dir" ]; then
    mkdir $Dir
  fi
  mv *.txt *.log *.pdf $Dir/
done
