#!/bin/bash

#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV
BaseDir=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV
AnalysisJSON=../JSON/json_190456_207898_analysis.txt
#${BaseDir}/Prompt/Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt
LumiJSON=${BaseDir}/PileUp/pileup_JSON_DCSONLY_190389-208686_corr.txt

#JSON file used to filter events (from DCSOnly or Prompt subdir)

pileupCalc.py \
    -i ${AnalysisJSON} \
    --inputLumiJSON ${LumiJSON} \
    --calcMode true \
    --minBiasXsec 69300 \
    --maxPileupBin 60 \
    --numPileupBins 60  \
    DataPileupHistogram-69p3mb.root

pileupCalc.py \
    -i ${AnalysisJSON} \
    --inputLumiJSON ${LumiJSON} \
    --calcMode true \
    --minBiasXsec 73500 \
    --maxPileupBin 60 \
    --numPileupBins 60  \
    DataPileupHistogram-73p5mb.root
