#!/bin/bash

#/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV
BaseDir=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV
AnalysisJSON=${BaseDir}/Prompt/Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt
LumiJSON=${BaseDir}/PileUp/pileup_JSON_DCSONLY_190389-203002_corr.txt

#JSON file used to filter events (from DCSOnly or Prompt subdir)

pileupCalc.py \
    -i ${AnalysisJSON} \
    --inputLumiJSON ${LumiJSON} \
    --calcMode true \
    --minBiasXsec 69400 \
    --maxPileupBin 60 \
    --numPileupBins 60  \
    DataPileupHistogram.root
