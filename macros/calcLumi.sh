#!/bin/bash

#Up to run 180252: 4.632 inv fb

PROMPTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Prompt/Cert_190456-207898_8TeV_PromptReco_Collisions12_JSON.txt

MINRUN=190456
MAXRUN=193621
INPUTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190456-196531_8TeV_13Jul2012ReReco_Collisions12_JSON_v2.txt
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON JSON_Run2012A-13Jul2012.json
  
MINRUN=190782
MAXRUN=190949
INPUTJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/Reprocessing/Cert_190782-190949_8TeV_06Aug2012ReReco_Collisions12_JSON.txt
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON JSON_Run2012A-06Aug2012ReReco.json

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

MINRUN=203768 
MAXRUN=207898
INPUTJSON=${PROMPTJSON}
./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON Run2012D-PromptReco-v1.json

###################

rm -f final.json

compareJSON.py --or JSON_Run2012A-13Jul2012.json JSON_Run2012A-06Aug2012ReReco.json finalA.json
compareJSON.py --or Run2012B-13Jul2012.json      finalA.json finalB.json
compareJSON.py --or Run2012C-ReReco.json         finalB.json finalC.json
compareJSON.py --or Run2012C-PromptReco-v2.json  finalC.json finalD.json
compareJSON.py --or Run2012D-PromptReco-v1.json  finalD.json final.json

lumiCalc2.py -i final.json -b stable overview

#now make pu dist
#json_190456_207898_analysis.txt