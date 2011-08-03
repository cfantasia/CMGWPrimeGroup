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
RELEASE_VERSION=CMSSW_4_2_3
WORKING_AREA=V290
# end definitions

export RELEASE_VERSION WORKING_AREA

echo -e  "******************************************"
echo -en " Setting up new working area: "
echo -e $WORKING_AREA
echo -e  "******************************************"

#setup new working area
mkdir $WORKING_AREA
cd $WORKING_AREA
export SCRAM_ARCH=slc5_amd64_gcc434
scramv1 p CMSSW $RELEASE_VERSION
cd $RELEASE_VERSION/src

echo -e  "\n**************************"
echo -e  " Checking out the code..."
echo -e  "**************************"
cvs -Q co -r $RELEASE_VERSION DataFormats/PatCandidates
cvs -Q co -r $RELEASE_VERSION PhysicsTools/PatAlgos
cvs -Q co -r V00-05-00 -d SHarper/HEEPAnalyzer UserCode/SHarper/HEEPAnalyzer 
cvs -Q co -r V08-03-09 PhysicsTools/Utilities
cvs -Q co -r V00-02-90 UserCode/CMGWPrimeGroup
# cvs -Q co UserCode/CMGWPrimeGroup

echo -e  "\n************************************************************"
echo -e  " Hack PAT Muon content to include high-pt reconstructors..."
echo -e  " Fix BuildFiles (ie. make 41x-compatible) for HEEP code..."
echo -e  " Add WP80 selection to HEEP code...."
echo -e  "************************************************************"
mv DataFormats/PatCandidates/interface/Muon.h DataFormats/PatCandidates/interface/Muon.h_original
cp UserCode/CMGWPrimeGroup/PAT_hack/Muon.h DataFormats/PatCandidates/interface/
mv DataFormats/PatCandidates/src/Muon.cc DataFormats/PatCandidates/src/Muon.cc_original
cp UserCode/CMGWPrimeGroup/PAT_hack/Muon.cc DataFormats/PatCandidates/src/
mv PhysicsTools/PatAlgos/plugins/PATMuonProducer.h PhysicsTools/PatAlgos/plugins/PATMuonProducer.h_original
cp UserCode/CMGWPrimeGroup/PAT_hack/PATMuonProducer.h PhysicsTools/PatAlgos/plugins/
mv PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc PhysicsTools/PatAlgos/plugins/PATMuonProducer.cc_original
cp UserCode/CMGWPrimeGroup/PAT_hack/PATMuonProducer.cc PhysicsTools/PatAlgos/plugins/
mv PhysicsTools/PatAlgos/python/producersLayer1/muonProducer_cfi.py PhysicsTools/PatAlgos/python/producersLayer1/muonProducer_cfi.py_original
cp UserCode/CMGWPrimeGroup/PAT_hack/muonProducer_cfi.py PhysicsTools/PatAlgos/python/producersLayer1/
cp UserCode/CMGWPrimeGroup/SHarper_hack/BuildFile1.xml SHarper/HEEPAnalyzer/BuildFile.xml
cp UserCode/CMGWPrimeGroup/SHarper_hack/BuildFile2.xml SHarper/HEEPAnalyzer/plugins/BuildFile.xml
cp UserCode/CMGWPrimeGroup/SHarper_hack/WP80SelectionCuts_cfi.py SHarper/HEEPAnalyzer/python/
mv SHarper/HEEPAnalyzer/src/HEEPEleSelector.cc SHarper/HEEPAnalyzer/src/HEEPEleSelector.cc_original
cp UserCode/CMGWPrimeGroup/SHarper_hack/HEEPEleSelector.cc SHarper/HEEPAnalyzer/src/
mv SHarper/HEEPAnalyzer/interface/HEEPEle.h SHarper/HEEPAnalyzer/interface/HEEPEle.h_original
cp UserCode/CMGWPrimeGroup/SHarper_hack/HEEPEle.h SHarper/HEEPAnalyzer/interface/ 
mv SHarper/HEEPAnalyzer/interface/HEEPCutValues.h  SHarper/HEEPAnalyzer/interface/HEEPCutValues.h_original
cp UserCode/CMGWPrimeGroup/SHarper_hack/HEEPCutValues.h SHarper/HEEPAnalyzer/interface/
mv SHarper/HEEPAnalyzer/src/HEEPCutValues.cc  SHarper/HEEPAnalyzer/src/HEEPCutValues.cc_original
cp UserCode/CMGWPrimeGroup/SHarper_hack/HEEPCutValues.cc  SHarper/HEEPAnalyzer/src/
mv SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h SHarper/HEEPAnalyzer/interface/HEEPCutCodes.h_original
cp UserCode/CMGWPrimeGroup/SHarper_hack/HEEPCutCodes.h SHarper/HEEPAnalyzer/interface/
mv SHarper/HEEPAnalyzer/src/HEEPCutCodes.cc SHarper/HEEPAnalyzer/src/HEEPCutCodes.cc_original
cp UserCode/CMGWPrimeGroup/SHarper_hack/HEEPCutCodes.cc SHarper/HEEPAnalyzer/src/
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
ln -s UserCode/CMGWPrimeGroup/test/patTuple_mumet_data.py .
ln -s UserCode/CMGWPrimeGroup/test/patTuple_elmet_data.py .
ln -s UserCode/CMGWPrimeGroup/root_macros/ZMET_data.root .
