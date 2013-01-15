#!/bin/bash

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 Input <Output>"
  exit
fi

Input=$1
if [ "$#" -gt 1 ]; then
    Output=`pwd`/$2
else
    Output=`pwd`/tables-`date +%Y-%m-%d`.dat
fi

if [ -f $Output ]
then
    mv -f $Output{,.bck} 
fi

root -b -l -q 'printNEvts.C+('\"$Input\"', -1)' >> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -2)' >> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -999)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 0)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 1)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 2)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 3)' >> $Output

cd ../Limits

# #Breakdown of Bkg After Lt Cut
# #Don't scale MC for this table
 root -b -l -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "tbl noWind")' >> $Output
 root -b -l -q 'printTable.C+(1)' >> $Output

# #Breakdown of Bkg After Window Cut
root -b -l -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "tbl")' >> $Output
root -b -l -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "scaleMC tbl")' > /dev/null
root -b -l -q 'printTable.C+(0)' >> $Output

