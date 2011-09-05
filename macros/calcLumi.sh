#!/bin/bash

JSON_May10=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt
MIN_May10=160404
MAX_May10=163869
./jsonrunsel.py $MIN_May10   $MAX_May10   $JSON_May10   May10.json

JSON_Prompt4=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Prompt/Cert_160404-173664_7TeV_PromptReco_Collisions11_JSON.txt
MIN_Prompt4=165088
MAX_Prompt4=167913
./jsonrunsel.py $MIN_Prompt4 $MAX_Prompt4 $JSON_Prompt4 Prompt4.json

JSON_Aug05=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions11/7TeV/Reprocessing/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.txt
MIN_Aug05=170249
MAX_Aug05=172619
./jsonrunsel.py $MIN_Aug05   $MAX_Aug05   $JSON_Aug05   Aug05.json

JSON_Prompt6=$JSON_Prompt4
MIN_Prompt6=172620
MAX_Prompt6=999999
./jsonrunsel.py $MIN_Prompt6 $MAX_Prompt6 $JSON_Prompt6 Prompt6.json

rm -f final.json

compareJSON.py --or May10.json Prompt4.json finalA.json
compareJSON.py --or Aug05.json Prompt6.json finalB.json
compareJSON.py --or finalA.json finalB.json final.json

rm finalA.json finalB.json

lumiCalc2.py -i final.json overview
