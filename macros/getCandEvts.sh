#!/bin/bash

rm -f fileList.txt

#Reads in format run:lumi:evt
while read LINE
do
  RUN=${LINE%%:*}
  LUMI=${LINE%:*}
  LUMI=${LUMI##*:}
#  echo Run:Lumi  ${RUN}:${LUMI}

for ERA in \
      Run2011A-May10ReReco-v1 \
      Run2011A-PromptReco-v4 \
      Run2011A-05Aug2011-v1 \
      Run2011A-PromptReco-v6 \
      Run2011B-PromptReco-v1 
    do
    for PD in \
        DoubleMu DoubleElectron SingleMu SingleElectron
      do
      
      DATASET=/${PD}/${ERA}/AOD
      echo "Looking in $DATASET"

      LIST=`dbs search --query="find file where run=${RUN} and lumi=${LUMI} and dataset=${DATASET}" | grep .root`
      NFILE=`echo $LIST | wc -w`
      if [ "$NFILE" -ge "1" ]; then
#        echo File num $NFILE
          echo $LIST >> fileList.txt
      fi
    done
  done
done < candEvts.txt

cat fileList.txt
cmsRun skimCandEvents.py
