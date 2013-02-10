#!/bin/bash

STARTRUN=190456
  INCRUN=000100
  ENDRUN=207898
#ENDRUN=190956


rm lumi.log
for Run in `seq $STARTRUN $INCRUN $ENDRUN`
  do
  echo $Run
  MINRUN=$STARTRUN
  MAXRUN=$Run
  #MINRUN=$Run
  #MAXRUN=$(($MINRUN + $INCRUN))
  INPUTJSON=../JSON/json_190456_207898_analysis.txt
  OUTPUTJSON=JSON_${MINRUN}-${MAXRUN}.json
  ./jsonrunsel.py $MINRUN $MAXRUN $INPUTJSON $OUTPUTJSON >& /dev/null
  
  echo $OUTPUTJSON >> lumi.log
  lumiCalc2.py -i $OUTPUTJSON -b stable overview | tail --lines=4 >> lumi.log

done

grep "|" lumi.log > lumi2.log
awk '0 == NR % 2'  lumi2.log

