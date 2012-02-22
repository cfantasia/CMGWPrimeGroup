#!/bin/bash

if [ ! $# == 2 ]; then
  echo "Usage: $0 Input Output"
  exit
fi

Input=$1
Output=`pwd`/$2

if [ -f $Output ]
then
    mv -f $Output{,.bck} 
fi

root -b -q 'printNEvts.C+('\"$Input\"', -1)' >> $Output
root -b -q 'printNEvts.C+('\"$Input\"', 0)' >> $Output
root -b -q 'printNEvts.C+('\"$Input\"', 1)' >> $Output
root -b -q 'printNEvts.C+('\"$Input\"', 2)' >> $Output
root -b -q 'printNEvts.C+('\"$Input\"', 3)' >> $Output

cd ../Limits

#Breakdown of Bkg After Ht Cut
#Don't scale MC for this table
root -b -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "tbl noWind")' >> $Output
root -b -q 'printTable.C+(1)' >> $Output

#Breakdown of Bkg After Window Cut
root -b -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "tbl")' >> $Output
root -b -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "scaleMC tbl")' >> junk.txt
root -b -q 'printTable.C+(0)' >> $Output

