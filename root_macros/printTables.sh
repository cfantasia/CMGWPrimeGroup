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

#######
cd ../Systematics
root -b -l -q calcRecoilSys.C+ >> $Output
root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZMuPtRes.root",   "SysMuPtRes.dat")' 
root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZMuPtScale.root", "SysMuPtScale.dat")' 
root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZElEnScale.root", "SysElEnScale.dat")' 
root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZPUSys.root",     "SysPU.dat")' 

#######

cd ../root_macros
echo "% 4 Channels Combined Yields">> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -1)' >> $Output
echo "% 2 Channels Combined Yields">> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -2)' >> $Output
echo "% Indiv Yields">> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -999)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 0)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 1)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 2)' >> $Output
#root -b -l -q 'printNEvts.C+('\"$Input\"', 3)' >> $Output

cd ../Limits

# #Breakdown of Bkg After Lt Cut
# #Don't scale MC for this table
echo "% Bkg Breakdown After Lt Cut">> $Output
 root -b -l -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "tbl noWind")' >> $Output
 root -b -l -q 'printTable.C+(1)' >> $Output


echo "% Breakdown of Bkg After LT and Window Cut">> $Output
root -b -l -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "tbl")' >> $Output
root -b -l -q 'ExpectedEvts.C+('\"$Input\"',"cutValues.wz.dat",-1, "scaleMC tbl")' > /dev/null
root -b -l -q 'printTable.C+(0)' >> $Output


##############################
cd ../root_macros
echo "% Systematic Tables" >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMETRes.dat",    "../Systematics/SysBkgMETRes.dat",    "../Systematics/SysMETRes.pdf")'    >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMETScale.dat",  "../Systematics/SysBkgMETScale.dat",  "../Systematics/SysMETScale.pdf")'  >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMuPtRes.dat",   "../Systematics/SysBkgMuPtRes.dat",   "../Systematics/SysMuPtRes.pdf")'   >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMuPtScale.dat", "../Systematics/SysBkgMuPtScale.dat", "../Systematics/SysMuPtScale.pdf")' >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigElEnScale.dat", "../Systematics/SysBkgElEnScale.dat", "../Systematics/SysElEnScale.pdf")' >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigPDF.dat",       "../Systematics/SysBkgPDF.dat",       "../Systematics/SysPDF.pdf")'       >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigPU.dat",        "../Systematics/SysBkgPU.dat",        "../Systematics/SysPU.pdf")'        >> $Output
