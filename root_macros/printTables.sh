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

####################
###  Yields ########
####################

cd ../root_macros
echo "% 4 Channels Combined Yields">> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -1)' >> $Output
echo "% 2 Channels Combined Yields">> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -2)' >> $Output
echo "% Indiv Yields">> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -4)' >> $Output
root -b -l -q 'printNEvts.C+('\"$Input\"', -5)' >> $Output

####################
###  Limits ########
####################

cd ../combined_limits
echo "% Final Yield, Limit Tables" >> $Output
root -b -l -q 'printTables.C+('\"$Input\"')' >> $Output

####################
###  Calc Sys ######
####################

#cd ../Systematics
#echo "% Sys Tables" >> $Output
#root -b -l -q calcRecoilSys.C+ >> $Output
#root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZMuPtRes.root",   "SysMuPtRes.dat")' 
#root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZMuPtScale.root", "SysMuPtScale.dat")' 
#root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZElEnScale.root", "SysElEnScale.dat")' 
#root -b -l -q 'compareYields.C+('\"$Input\"', "../../../WprimeWZPUSys.root",     "SysPU.dat")' 

####################
###  Print Sys #####
####################

cd ../root_macros
echo "% Systematic Tables" >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMETRes.dat",    "../Systematics/SysBkgMETRes.dat",    "../Systematics/SysMETRes.pdf")'    >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMETScale.dat",  "../Systematics/SysBkgMETScale.dat",  "../Systematics/SysMETScale.pdf")'  >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMuPtRes.dat",   "../Systematics/SysBkgMuPtRes.dat",   "../Systematics/SysMuPtRes.pdf")'   >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigMuPtScale.dat", "../Systematics/SysBkgMuPtScale.dat", "../Systematics/SysMuPtScale.pdf")' >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigElEnScale.dat", "../Systematics/SysBkgElEnScale.dat", "../Systematics/SysElEnScale.pdf")' >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigPDF.dat",       "../Systematics/SysBkgPDF.dat",       "../Systematics/SysPDF.pdf")'       >> $Output
root -b -l -q 'printSysTable.C+("../Systematics/SysSigPU.dat",        "../Systematics/SysBkgPU.dat",        "../Systematics/SysPU.pdf")'        >> $Output

#####################
#### Plots ##########
#####################

cd ../root_macros
echo "% Plots and Eff Tables" >> $Output
root -b -l -q 'MakePlots.cc+('\"$Input\"', "", "")'

root -b -l -q 'EffVsMass.C+('\"$Input\"', 0)' >> $Output #EWK Cuts
root -b -l -q 'EffVsMass.C+('\"$Input\"', 1)' >> $Output #Lt  Cuts
root -b -l -q 'EffVsMass.C+('\"$Input\"', 2)' >> $Output #Lt + Wind Cuts

root -b -l -q 'PlotResolution.C+('\"$Input\"', 2000)'
root -b -l -q 'ZptVsWpt.C+("../../../WprimeWZ.root", "")'
