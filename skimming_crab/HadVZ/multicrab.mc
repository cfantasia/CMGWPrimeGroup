[MULTICRAB]
cfg=crab.cfg

[COMMON]
CRAB.scheduler=condor
CMSSW.lumis_per_job=300
CMSSW.pset=../../test/patTuple_jetlep_MC.py
GRID.se_black_list=T2_UK_London_IC
#CRAB.use_server=1
#CRAB.scheduler=glite

[DiLeptonJet-V400B-DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball]
CMSSW.datasetpath=/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v2/AODSIM

[DiLeptonJet-V400B-TTJets_TuneZ2star_8TeV-madgraph-tauola]
CMSSW.datasetpath=/TTJets_TuneZ2star_8TeV-madgraph-tauola/Summer12-PU_S7_START52_V9-v1/AODSIM

