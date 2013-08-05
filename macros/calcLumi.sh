#!/bin/bash

#Run2012A-13Jul2012.json       193621    808.472 pb-1
#Run2012A-06Aug2012ReReco.json 190949     82.136 pb-1 
#Run2012B-13Jul2012.json       196531  4.429 fb-1
#Run2012C-ReReco.json          198913    495.003 pb-1
#Run2012C-PromptReco-v2.json   203746  6.401 fb -1
#Run2012C-EcalRecover.json     201191    134.242 pb-1
#Run2012D-PromptReco-v1.json   207898  7.274 fb-1
####################################################
#Run2012A                      193621    891.6 fb-1
#Run2012B                      196531  4.429 fb-1
#Run2012C                      203746  7.030 fb-1
#Run2012D                      207898  7.274 fb-1
####################################################
#Total Up to run               207898 18.258 inv fb (doesn't include ecal recover)
#Total Up to run               208686 19.624 inv fb

PROMPTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-208686_8TeV_PromptReco_Collisions12_JSON.txt

MINRUN=190456
MAXRUN=193621
INPUTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012A-13Jul2012.json
  
MINRUN=190782
MAXRUN=190949
INPUTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012A-06Aug2012ReReco.json

MINRUN=193833
MAXRUN=196531
INPUTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012B-13Jul2012.json

MINRUN=198022
MAXRUN=198913
INPUTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_198022-198523_8TeV_24Aug2012ReReco_Collisions12_JSON.txt
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012C-ReReco.json
  
MINRUN=198934
MAXRUN=203746
INPUTJSON=${PROMPTJSON}
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012C-PromptReco-v2.json

MINRUN=201191
MAXRUN=201191
INPUTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_201191-201191_8TeV_11Dec2012ReReco-recover_Collisions12_JSON.txt
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012C-EcalRecover.json

MINRUN=203768 
MAXRUN=208686
INPUTJSON=${PROMPTJSON}
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012D-PromptReco-v1.json

###################

rm -f final.json

compareJSON.py --or Run2012A-13Jul2012.json Run2012A-06Aug2012ReReco.json finalA.json
compareJSON.py --or Run2012B-13Jul2012.json      finalA.json finalB.json
compareJSON.py --or Run2012C-ReReco.json         finalB.json finalC.json
compareJSON.py --or Run2012C-PromptReco-v2.json  finalC.json finalD.json
compareJSON.py --or Run2012C-EcalRecover.json    finalD.json finalE.json
compareJSON.py --or Run2012D-PromptReco-v1.json  finalE.json final.json

lumiCalc2.py -i final.json -o final.csv -b stable overview
#pixelLumi doesn't have all the latest numbers so don't use for now
#pixelLumiCalc.py -i final.json  overview #Don't need stable since pixel only on for stable


cat final.csv | awk -F'[:,]' '{print $1 "\t" $NF}' > final.dat

#TGraph* g = TGraph("final.dat")
#for (int i=1;i<g->GetN();++i) g->GetY()[i] += g->GetY()[i-1]; #Make graph cumlative
#g->Draw("al")
