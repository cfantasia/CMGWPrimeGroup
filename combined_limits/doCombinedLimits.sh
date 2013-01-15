#!/bin/bash

function run {
    Mass=$1 
    NToys=100
    Seed=-1

    BaseName=WprimeWZ
    cardFile=card_${BaseName}_M${Mass}.txt
    
    if [  ! -f $cardFile ]; then
        echo "$cardFile does not exist, creating."
        root -b -l -q 'makeLimitCard.C+("../../../WprimeWZ.root", '${Mass}')'
        cat card_Head.txt card_Core_${BaseName}_M${Mass}.txt card_Tail.txt > ${cardFile}
    fi
    
    #CombineCMD="combine -H ProfileLikelihood -M ${Mode} --tries 200 -i 30000 -b 100 -t ${NToys} -s ${Seed} -n ${BaseName} -m ${Mass} ${cardFile}"
    CombineCMD="combine -H ProfileLikelihood  -M ${Mode} -t ${NToys} -s ${Seed} -n ${BaseName} -m ${Mass} ${cardFile}"
    echo $CombineCMD
    $CombineCMD # >& /dev/null

    name=higgsCombineWprimeWZ.${Mode}.mH${Mass}.root
    hadd -f ${name} higgsCombineWprimeWZ.${Mode}.mH${Mass}.*.root >& /dev/null
}  

function makePlots {
    #Now run code to extract limts and make plots
    echo          'extractLimits.C+("../../../WprimeWZ.root",'\"${Mode}\"')'
    root -b -l -q 'extractLimits.C+("../../../WprimeWZ.root",'\"${Mode}\"')' # >& /dev/null
    cat nLimit_${Mode}.txt

    cd ../Limits/
    root -b -l -q '../Limits/PlotLimit.C+("WprimeWZ", "../root_macros/nLimit_'${Mode}'.txt")' #For W'
    scp  limitVsMass_WZ.pdf buphy.bu.edu:~/public_html

}

    #Mode=BayesianToyMC
    #Mode=MarkovChainMC
    Mode=Asymptotic 

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

