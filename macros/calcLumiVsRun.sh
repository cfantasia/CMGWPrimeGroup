#!/bin/bash

STARTRUN=190456
INCRUN=100
#ENDRUN=190956
ENDRUN=207898


rm lumi.log
for Run in `seq $STARTRUN $INCRUN $ENDRUN`
  do
  
  MINRUN=$Run
  MAXRUN=$(($MINRUN + $INCRUN))
  INPUTJSON=../JSON/json_190456_207898_analysis.txt
  OUTPUTJSON=JSON_${MINRUN}-${MAXRUN}.json
  ./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON $OUTPUTJSON >& /dev/null
  
  echo $OUTPUTJSON >> lumi.log
  lumiCalc2.py -i $OUTPUTJSON -b stable overview | tail --lines=4 >> lumi.log

done