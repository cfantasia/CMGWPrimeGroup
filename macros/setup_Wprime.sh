#!/bin/bash
if [ -f /etc/bashrc ]; then
      . /etc/bashrc
else
echo -e " No /etc/bashrc file! So confused..."
echo -e " Reflect, repent, and try with a better script next time"
echo -e " Exiting...."
exit
fi

# Definitions
RELEASE_VERSION=CMSSW_5_3_7_patch4
WORKING_AREA=V405
# end definitions

export RELEASE_VERSION WORKING_AREA

echo -e  "******************************************"
echo -en " Setting up new working area: "
echo -e $WORKING_AREA
echo -e  "******************************************"

#setup new working area
mkdir $WORKING_AREA
cd $WORKING_AREA
export SCRAM_ARCH=slc5_amd64_gcc462
scramv1 p CMSSW $RELEASE_VERSION
cd $RELEASE_VERSION/src

echo -e  "\n**************************"
echo -e  " Checking out the code..."
echo -e  "**************************"
cvs co -r V00-00-13 RecoMET/METFilters
# For CSC Beam Halo filter
cvs co -r V00-00-08 RecoMET/METAnalyzers
# For the HBHE filter
# For 52x release, please check out the following branch since it fixes a compiling issue
# cvs co -r BV00-03-20_HBHEfor52X CommonTools/RecoAlgos
cvs co -r V00-03-23 CommonTools/RecoAlgos
# Additional packages for the tracking POG filters
cvs co -r V01-00-11-01 DPGAnalysis/Skims
cvs co -r V00-11-17 DPGAnalysis/SiStripTools
cvs co -r V00-00-08 DataFormats/TrackerCommon
cvs co -r V01-09-05 RecoLocalTracker/SubCollectionProducers
#Sam's package for HEEP
cvs -Q co -r V00-08-01 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer 
#Common Wprime Group Package
cvs -Q co -r V00-04-05 UserCode/CMGWPrimeGroup
#cvs -Q co UserCode/CMGWPrimeGroup

echo -e "\n************************"
echo -e " Done. Now will compile"
echo -e "************************"
sleep 2

source /afs/cern.ch/cms/sw/cmsset_default.sh
# Not sure how to make script aware of cmsenv alias...
#cmsenv
scram runtime -sh
scram b -j 4 # (lots of output here, but nothing to worry about)

echo -e "\n*************************************************"
echo -e " Done compiling UserCode/CMGWPrimeGroup"
echo -e " Now making symbolic link to example config file"
echo -e "*************************************************"
echo -e " \n"
