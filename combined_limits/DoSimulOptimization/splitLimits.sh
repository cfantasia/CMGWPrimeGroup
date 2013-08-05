#!/bin/tcsh -f

if ( $# < 1 ) then
  echo "Usage: `basename $0` Code <Offset>"
  exit 1
endif

set inCode=${1}

set maxCode=15

@ Code = ( $inCode % $maxCode )
@ Part = ( $inCode / $maxCode )

echo $inCode " " $Code " " $Part

#never use 0
if ( $# > 1 ) then
    @ Code = ( $Code + $2 )
endif


@ Code = ( $Code * 10 )

if ( $Code == 10 ) then
    @ Code = ( 25 )
endif

@ Base = ( $Code / 10 )
@ Frac = ( $Code % 10 )

echo Final code is $Code " and part is " $Part " and base is " $Base " and frac is " $Frac

source /uscmst1/prod/sw/cms/setup/cshrc prod
setenv SCRAM_ARCH slc5_amd64_gcc434

cd /uscms_data/d2/fantasia/CommonWprime/HEAD2/CMSSW_4_2_5/src/UserCode/CMGWPrimeGroup/Limits
eval `scramv1 runtime -csh`

source ../../../StatisticalTools/RooStatsRoutines/setup/cmslpc_standalone_setup.csh

pwd
date


set inFile=nEvents_${Code}0_part${Part}
set outFile=nLimits_${Code}0_part${Part}.txt
rm -f $inFile.tmp

grep "SignalCode"  nEvents.txt > $inFile
grep "^${Base}.${Frac}" nEvents.txt >> $inFile.tmp

echo There are `wc -l nEvents.txt` lines
echo There are `wc -l $inFile.tmp` lines

split --suffix-length=1 --numeric-suffixes -l 19 $inFile.tmp nEvents_${Code}0_tmp
cat nEvents_${Code}0_tmp${Part} >> $inFile


echo $inFile
cat $inFile

echo root -b -q -n 'CalcLimit.C+(1, '\"$inFile\"', '\"$outFile\"')'
root -b -q -n 'CalcLimit.C+(1, '\"$inFile\"', '\"$outFile\"')' >& /dev/null #Logs/CalcLimit_${Code}0_part${Part}.log

rm -f $inFile.tmp
rm -f nEvents_${Code}0_tmp${Part}
rm -f nEvents_${Code}0_tmp?

echo "Done"  
