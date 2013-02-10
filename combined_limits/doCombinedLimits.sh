#!/bin/bash

function run {
    Mass=$1 
    NToys=10
    Seed=-1

    BaseName=WprimeWZ
    cardFile=card_${BaseName}_M${Mass}.txt
    cardCoreFile=card_Core_${BaseName}_M${Mass}.txt
    
    if [  ! -f $cardFile ]; then
        echo "$cardFile does not exist, creating."
        root -b -l -q 'makeLimitCard.C+("../../../WprimeWZ.root", "'$cardCoreFile'", '${Mass}')'
        cat card_Head.txt $cardCoreFile card_Tail.txt > ${cardFile}
    fi
    
    CombineCMD="combine -H $HintMode -M ${Mode} -s ${Seed} -n ${BaseName} -m ${Mass} ${cardFile}"  #MarkovChainMC Line

    #CombineCMD="combine -H $HintMode -M ${Mode} --tries 200 --iteration 30000 --burnInSteps 100 -s ${Seed} -n ${BaseName} -m ${Mass} ${cardFile}"  #MarkovChainMC Line
    #CombineCMD="combine -H $HintMode -M ${Mode} -s ${Seed} -n ${BaseName} -m ${Mass} ${cardFile}"
    #CombineCMD="combine -H $HintMode --rMax 1 -M ${Mode} -s ${Seed} -n ${BaseName} -m ${Mass} ${cardFile}"
    #CombineCMD="combine -H $HintMode -V -v 2 -M ${Mode} -s ${Seed} -n ${BaseName} -m ${Mass} ${cardFile}"
    echo $CombineCMD
    $CombineCMD --toys ${NToys} # >& /dev/null   #Expected Limit
    $CombineCMD                 # >& /dev/null   #Observed Limit

    name=higgsCombineWprimeWZ.${Mode}.mH${Mass}.root
    hadd -f ${name} higgsCombineWprimeWZ.${Mode}.mH${Mass}.*.root >& /dev/null
}  

function makePlots {
    #Now run code to extract limts and make plots
    echo          'extractLimits.C+("../../../WprimeWZ.root",'\"${Mode}\"')'
    root -b -l -q 'extractLimits.C+("../../../WprimeWZ.root",'\"${Mode}\"')' # >& /dev/null
    cat nLimit_${Mode}.txt

    cd ../Limits/
    root -b -l -q '../Limits/PlotLimit.C+("WprimeWZ", "../combined_limits/nLimit_'${Mode}'.txt")' #For W'
    scp  limitVsMass_WZ.pdf buphy.bu.edu:~/public_html

}
 
    #HintMode=ProfileLikelihood
    HintMode=Asymptotic
    #Mode=BayesianToyMC
    Mode=MarkovChainMC      #Bayesian
    #Mode=ProfileLikelihood  #?
    #Mode=Asymptotic          #Asymptotic CLs

if [ "$#" -eq 0 ]; then
    makePlots
    exit 0
fi

if [ "$#" -gt 1 ]; then

    for M in $*
      do
      echo Running mass $M
      run $M
    done

    makePlots
else
    run $1
fi

