#!/bin/bash -f

echo "Arguments are $*"

if [ $# -ne 5 ]; then
    echo "Usage: $0 mass ltshift windshift NToys dir"
    exit 1
fi
echo "YES!!!"

source /uscmst1/prod/sw/cms/setup/shrc prod
#setenv SCRAM_ARCH slc5_amd64_gcc434

cd $5
eval `scramv1 runtime -sh`

#cd -
cd $_CONDOR_SCRATCH_DIR

pwd
date

#python doCombinedLimits.py -M $1 --LtShift=$2 --WindShift=$3 --ntoys=$4 --doDuplicates
python doCombinedLimits.py -M $1 --LtShift=$2 --WindShift=$3 --ntoys=$4 --doDuplicates --doSepCh

