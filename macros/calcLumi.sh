#!/bin/bash

#Up to run 176023: 2.301 inv fb
#Up to run 176309: 2.510 inv fb

JSON_Run2011A_May10=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt
MIN_Run2011A_May10=160404
MAX_Run2011A_May10=163869
./jsonrunsel.py $MIN_Run2011A_May10   $MAX_Run2011A_May10   $JSON_Run2011A_May10   Run2011A_May10.json

JSON_Run2011A_Prompt4=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-176309_7TeV_PromptReco_Collisions11_JSON.txt
MIN_Run2011A_Prompt4=165088
MAX_Run2011A_Prompt4=167913
./jsonrunsel.py $MIN_Run2011A_Prompt4 $MAX_Run2011A_Prompt4 $JSON_Run2011A_Prompt4 Run2011A_Prompt4.json

JSON_Run2011A_Aug05=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.txt
MIN_Run2011A_Aug05=170249
MAX_Run2011A_Aug05=172619
./jsonrunsel.py $MIN_Run2011A_Aug05   $MAX_Run2011A_Aug05   $JSON_Run2011A_Aug05   Run2011A_Aug05.json

JSON_Run2011A_Prompt6=$JSON_Run2011A_Prompt4
MIN_Run2011A_Prompt6=172620
MAX_Run2011A_Prompt6=173692
./jsonrunsel.py $MIN_Run2011A_Prompt6 $MAX_Run2011A_Prompt6 $JSON_Run2011A_Prompt6 Run2011A_Prompt6.json

JSON_Run2011B_Prompt1=$JSON_Run2011A_Prompt4
MIN_Run2011B_Prompt1=175860
MAX_Run2011B_Prompt1=999999
./jsonrunsel.py $MIN_Run2011B_Prompt1 $MAX_Run2011B_Prompt1 $JSON_Run2011B_Prompt1 Run2011B_Prompt1.json

rm -f final.json

compareJSON.py --or Run2011A_May10.json Run2011A_Prompt4.json finalA.json
compareJSON.py --or Run2011A_Aug05.json finalA.json finalB.json
compareJSON.py --or Run2011A_Prompt6.json finalB.json finalC.json
compareJSON.py --or Run2011B_Prompt1.json finalC.json finalD.json

cp finalD.json final.json

rm finalA.json finalB.json finalC.json finalD.json

lumiCalc2.py -i final.json overview
