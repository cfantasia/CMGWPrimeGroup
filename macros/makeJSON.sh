#!/bin/bash

#This script makes a new json file by subtracting runs since the last time run.

#e.g. ./makeJSON.sh #Copies all of dcs only file
#e.g. ./makeJSON.sh 155422 #Copies all of dcs only file since run 155422
#e.g. ./makeJSON.sh Cert_SomeJSON_File #Copies all of dcs only file starting at last run in this json

lastRunReturn=0
#function last () {
lastRun () {
    #Find last 6 digit 
    lastRunReturn=`egrep -o "[0-9][0-9][0-9][0-9][0-9][0-9]" $1  | tail -1`
    
    echo "Found last run listed in file $1 to be $lastRunReturn"
    #lastRunReturn=99983
    return 
}



lastSkimmedRun=0
#find last run
if [ "$#" -gt 0 ]; then #if given, use that
    if ! [[ "$1" =~ ^[0-9]+$ ]] ; then #if not a number find in json
        lastRun $1
        lastSkimmedRun=$(($lastRunReturn + 1))
    else 
        lastSkimmedRun=${1}
    fi
else #Default last run number
    lastSkimmedRun=0
fi

#Grab latest dcs only from afs
oldJSON="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/DCSOnly/json_DCSONLY.txt"

#Name of output file
lastRun $oldJSON
lastAvailableRun=$lastRunReturn
newJSON="../JSON/json_${lastSkimmedRun}-${lastAvailableRun}_DCSonly.txt"

#remove older runs from latest dcs
./jsonrunsel.py $lastSkimmedRun 999999 $oldJSON $newJSON

echo "Output file is $newJSON which covers runs $lastSkimmedRun to ${lastAvailableRun}"