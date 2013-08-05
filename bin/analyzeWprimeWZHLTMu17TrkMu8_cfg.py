import sys
sys.path.append("UserCode/CMGWPrimeGroup/bin")
from analyzeWprimeWZ_cfg import *

process.WprimeAnalyzer.outputFile  = "WprimeWZHLTMu17TrkMu8.root"
process.WprimeAnalyzer.logFile     = "WprimeWZHLTMu17TrkMu8.dat"
process.WprimeAnalyzer.candEvtFile = "WprimeWZHLTMu17TrkMu8.evt"

process.WprimeAnalyzer.triggersToUse = cms.vstring(
#    'HLT_DoubleMu7_v*',
#    'HLT_Mu13_Mu8_v*', #1e33 unprescaled
#    'HLT_Mu17_Mu8_v*', #3e33 unprescaled
    'HLT_Mu17_TkMu8_v*', #2011B
#    'HLT_Mu22_TkMu8_v*', #2012 single l1 seeded
#    'HLT_Mu22_TkMu22_v*', #2012

    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*', #2011A
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',#MC
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',#Data
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*',
    'HLT_DoubleEle33_CaloIdT_v*',
    )
