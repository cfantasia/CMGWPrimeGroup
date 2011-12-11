#!/bin/bash

for windFracTenths in -1
  do
  Ones=$(($windFracTenths / 10))
  Tenths=$(($windFracTenths - $Ones*10))
  windFrac=$Ones"p"$Tenths
  echo "Using $windFrac as half window size"

  root -b -q 'ExpectedEvts.C+("../../../WprimeWZ.root","cutValues.wz.dat",'${windFracTenths}')' > ExpectedEvts.log 
#  root -b -q 'ExpectedEvts.C+("../../../old-2011-11-11/HadVZAnalyzer.root","cutValues.VZ.dat",'${windFracTenths}', "")' > ExpectedEvts.log 

  echo "Done with counting events"


  /afs/hep.wisc.edu/cern/.root/root_v5.30.00.Linux-slc5_amd64-gcc4.3/bin/root -b -q -n 'CalcLimit.C+(1)' > CalcLimit.log

  echo "Done with calculating limits"

  root -b -q 'PlotLimit.C+("WprimeWZ")' #For W'
  root -b -q 'PlotLimit2D.C+' #For TC
#  root -b -q 'PlotLimit.C+("HadVZ")' #For VZ
  
  echo "Done with Plotting limits"

  Dir=windFrac${windFrac}
  if [ ! -d "$Dir" ]; then
    mkdir $Dir
  fi
  mv *.txt *.log *.pdf $Dir/
done
