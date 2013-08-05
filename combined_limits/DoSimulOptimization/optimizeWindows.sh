#!/bin/tcsh -f

source /uscmst1/prod/sw/cms/setup/cshrc prod
setenv SCRAM_ARCH slc5_amd64_gcc434

cd /uscms_data/d2/fantasia/CommonWprime/HEAD2/CMSSW_4_2_5/src/UserCode/CMGWPrimeGroup/Limits
eval `scramv1 runtime -csh`

source ../../../StatisticalTools/RooStatsRoutines/setup/cmslpc_standalone_setup.csh

pwd
date

root -b -q MakeAllWindows.C  

root -b -q 'ExpectedEvts.C+("../../../WprimeWZ.root","cutValues.wzFull.dat",-1, "scaleMC")'

root -b -q -n 'CalcLimit.C+(1)'

#root -b -q OptimizeWindows.C
