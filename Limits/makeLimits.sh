#!/bin/bash

for windFracTenths in -1
  do
  Ones=$(($windFracTenths / 10))
  Tenths=$(($windFracTenths - $Ones*10))
  windFrac=$Ones"p"$Tenths
  echo "Using $windFrac as half window size"

#  root -b -q 'ExpectedEvts.C+("../../../WprimeWZ.root",'${windFracTenths}')' > ExpectedEvts.log 
  root -b -q 'ExpectedEvts.C+("../../../HadVZAnalyzer.root",'${windFracTenths}', "useHists")' > ExpectedEvts.log 

  /afs/hep.wisc.edu/cern/.root/root_v5.30.00.Linux-slc5_amd64-gcc4.3/bin/root -b -q -n CalcLimit.C > CalcLimit.log

#  root -b -q 'PlotLimit.C+("WprimeWZ")' #For W'
#  root -b -q 'PlotLimit.C+("TCWZ")' #For TC
  root -b -q 'PlotLimit.C+("HadVZ")' #For VZ
  
  Dir=windFrac${windFrac}
  if [ ! -d "$Dir" ]; then
    mkdir $Dir
  fi
  mv *.txt *.log *.pdf $Dir/
done
